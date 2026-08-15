#include "atac_spill_materializer.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "atac_mergeable_spill.h"
#include "barcode_correction.h"
#include "mapping_processor.h"
#include "mapping_writer.h"
#include "sequence_batch.h"

namespace chromap {
namespace {

AtacSpillMaterializationResult Failure(const std::string &message) {
  AtacSpillMaterializationResult result;
  result.message = message;
  return result;
}

bool SameReferences(const std::vector<AtacMergeableSpillReference> &a,
                    const std::vector<AtacMergeableSpillReference> &b) {
  if (a.size() != b.size()) {
    return false;
  }
  for (size_t i = 0; i < a.size(); ++i) {
    if (a[i].name != b[i].name || a[i].length != b[i].length) {
      return false;
    }
  }
  return true;
}

bool SameMaterializationContract(const AtacMergeableSpillMetadata &a,
                                 const AtacMergeableSpillMetadata &b) {
  return a.schema_mask == b.schema_mask && a.flags == b.flags &&
         a.shard_count == b.shard_count &&
         a.barcode_length == b.barcode_length &&
         a.mapq_threshold == b.mapq_threshold &&
         a.tn5_forward_shift == b.tn5_forward_shift &&
         a.tn5_reverse_shift == b.tn5_reverse_shift &&
         a.barcode_whitelist_fingerprint ==
             b.barcode_whitelist_fingerprint &&
         a.barcode_correction_error_threshold ==
             b.barcode_correction_error_threshold &&
         a.barcode_correction_probability_threshold ==
             b.barcode_correction_probability_threshold &&
         a.multi_mapping_allocation_distance ==
             b.multi_mapping_allocation_distance &&
         a.multi_mapping_allocation_seed == b.multi_mapping_allocation_seed &&
         a.max_num_best_mappings == b.max_num_best_mappings &&
         a.cache_size == b.cache_size &&
         a.k_for_minhash == b.k_for_minhash &&
         a.frip_est_coefficients == b.frip_est_coefficients &&
         a.sample_id == b.sample_id && a.input_id == b.input_id &&
         SameReferences(a.references, b.references) &&
         a.barcode_abundance_entries.size() ==
             b.barcode_abundance_entries.size() &&
         std::equal(a.barcode_abundance_entries.begin(),
                    a.barcode_abundance_entries.end(),
                    b.barcode_abundance_entries.begin(),
                    [](const AtacBarcodeAbundanceEntry &x,
                       const AtacBarcodeAbundanceEntry &y) {
                      return x.barcode_key == y.barcode_key;
                    });
}

}  // namespace

AtacSpillMaterializationResult MaterializeAtacSpillRecords(
    const std::vector<std::string> &spill_paths,
    const MappingParameters &output_parameters) {
  if (spill_paths.empty()) {
    return Failure("ATAC spill materializer input set is empty");
  }
  if (output_parameters.mapping_output_file_path.empty()) {
    return Failure("ATAC spill materializer output path is empty");
  }
  if ((output_parameters.mapping_output_format == MAPPINGFORMAT_BAM ||
       output_parameters.mapping_output_format == MAPPINGFORMAT_CRAM) &&
      output_parameters.atac_fragment_output_file_path.empty()) {
    return Failure(
        "ATAC BAM/CRAM materialization requires a fragments output path");
  }
  if (output_parameters.mapping_output_format == MAPPINGFORMAT_CRAM &&
      output_parameters.reference_file_path.empty()) {
    return Failure("ATAC CRAM materialization requires a reference FASTA");
  }
  if ((output_parameters.emit_noY_stream &&
       output_parameters.noY_output_path.empty()) ||
      (output_parameters.emit_Y_stream &&
       output_parameters.Y_output_path.empty())) {
    return Failure("ATAC Y-routing output path is empty");
  }
  if (output_parameters.mapping_output_format != MAPPINGFORMAT_BED &&
      output_parameters.mapping_output_format != MAPPINGFORMAT_BAM &&
      output_parameters.mapping_output_format != MAPPINGFORMAT_CRAM) {
    return Failure("ATAC spill materializer supports BED, BAM, or CRAM output");
  }

  std::vector<std::unique_ptr<AtacMergeableSpillReader>> readers;
  std::vector<std::string> paths_by_ordinal;
  AtacMergeableSpillMetadata contract;
  bool have_contract = false;
  uint64_t total_spill_records = 0;
  std::string error;

  for (const std::string &path : spill_paths) {
    std::unique_ptr<AtacMergeableSpillReader> reader(
        new AtacMergeableSpillReader());
    if (!reader->Open(path, &error)) {
      return Failure(error);
    }
    const AtacMergeableSpillMetadata &metadata = reader->metadata();
    if (!have_contract) {
      contract = metadata;
      readers.resize(contract.shard_count);
      paths_by_ordinal.resize(contract.shard_count);
      have_contract = true;
    } else if (!SameMaterializationContract(contract, metadata)) {
      return Failure("ATAC spill materializer input contracts disagree: " +
                     path);
    }
    if (metadata.shard_count != spill_paths.size()) {
      return Failure(
          "ATAC spill materializer input count does not match shard_count");
    }
    if (readers[metadata.shard_ordinal]) {
      return Failure("ATAC spill materializer has duplicate shard ordinal " +
                     std::to_string(metadata.shard_ordinal));
    }
    total_spill_records += reader->expected_record_count();
    paths_by_ordinal[metadata.shard_ordinal] = path;
    readers[metadata.shard_ordinal] = std::move(reader);
  }
  const bool read_ranges_late_bound =
      (contract.flags & kAtacMergeableReadRangeLateBound) != 0;
  std::vector<uint64_t> global_read_prefixes(readers.size(), 0);
  uint64_t next_global_read = 0;
  for (uint32_t ordinal = 0; ordinal < readers.size(); ++ordinal) {
    if (!readers[ordinal]) {
      return Failure("ATAC spill materializer is missing shard ordinal " +
                     std::to_string(ordinal));
    }
    if (readers[ordinal]->metadata().shard_ordinal != ordinal) {
      return Failure("ATAC spill materializer ordinal ordering failed");
    }
    if (read_ranges_late_bound) {
      if (readers[ordinal]->metadata().first_global_read_ordinal != 0) {
        return Failure(
            "late-bound ATAC spill set contains a precomputed read prefix");
      }
      global_read_prefixes[ordinal] = next_global_read;
      if (readers[ordinal]->metadata().input_record_count >
          std::numeric_limits<uint64_t>::max() - next_global_read) {
        return Failure("ATAC spill materializer read range overflows uint64");
      }
      next_global_read += readers[ordinal]->metadata().input_record_count;
    } else if (ordinal > 0) {
      const auto &previous = readers[ordinal - 1]->metadata();
      const auto &current = readers[ordinal]->metadata();
      if (previous.first_global_read_ordinal >
              std::numeric_limits<uint64_t>::max() -
                  previous.input_record_count ||
          current.first_global_read_ordinal !=
              previous.first_global_read_ordinal +
                  previous.input_record_count) {
        return Failure(
            "ATAC spill materializer shard read ranges are not contiguous");
      }
    }
    const auto &metadata = readers[ordinal]->metadata();
    if (!read_ranges_late_bound) {
      global_read_prefixes[ordinal] = metadata.first_global_read_ordinal;
    }
    if (metadata.input_record_count >
        std::numeric_limits<uint64_t>::max() -
            metadata.first_global_read_ordinal) {
      return Failure("ATAC spill materializer read range overflows uint64");
    }
  }
  const bool is_bulk =
      (contract.schema_mask & kAtacSpillSchemaIsBulk) != 0;
  const bool has_raw_barcode_evidence =
      (contract.schema_mask & kAtacSpillSchemaHasRawBarcodeEvidence) != 0;
  if (!is_bulk && !has_raw_barcode_evidence) {
    return Failure(
        "barcoded ATAC spill set lacks raw barcode correction evidence");
  }
  const bool allocate_multi_mappings =
      (contract.flags & kAtacMergeableAllocateMultiMappings) != 0;
  if (allocate_multi_mappings &&
      (contract.multi_mapping_allocation_distance < 0 ||
       contract.max_num_best_mappings == 0 ||
       contract.max_num_best_mappings >
           static_cast<uint32_t>(std::numeric_limits<int>::max()))) {
    return Failure("ATAC spill materializer has invalid multimapping policy");
  }

  std::unordered_map<uint64_t, uint64_t> global_barcode_abundance;
  uint64_t global_num_sample_barcodes = 0;
  if (!is_bulk) {
    global_barcode_abundance.reserve(
        contract.barcode_abundance_entries.size() * 2u);
    for (const auto &entry : contract.barcode_abundance_entries) {
      global_barcode_abundance.emplace(entry.barcode_key, 0);
    }
    for (const auto &reader : readers) {
      const auto &metadata = reader->metadata();
      if (global_num_sample_barcodes >
          std::numeric_limits<uint64_t>::max() -
              metadata.local_num_sample_barcodes) {
        return Failure("global barcode abundance count overflows");
      }
      global_num_sample_barcodes += metadata.local_num_sample_barcodes;
      for (const auto &entry : metadata.barcode_abundance_entries) {
        uint64_t &count = global_barcode_abundance[entry.barcode_key];
        if (count > std::numeric_limits<uint64_t>::max() - entry.count) {
          return Failure("global barcode abundance entry overflows");
        }
        count += entry.count;
      }
    }
    if (global_num_sample_barcodes == 0) {
      return Failure("global barcode correction model is empty");
    }
  }

  const auto barcode_abundance_lookup = [&](uint64_t key, bool *found) {
    const auto it = global_barcode_abundance.find(key);
    *found = it != global_barcode_abundance.end();
    return *found ? it->second : uint64_t{0};
  };
  auto correct_summary_barcode = [&](uint64_t raw_key, uint32_t n_mask,
                                     const std::string &quality,
                                     uint64_t *corrected_key) {
    *corrected_key = raw_key;
    if (is_bulk) {
      return BarcodeCorrectionStatus::kInWhitelist;
    }
    return CorrectPackedBarcode(
        raw_key, contract.barcode_length, n_mask, quality,
        contract.barcode_correction_error_threshold,
        contract.barcode_correction_probability_threshold,
        global_num_sample_barcodes, barcode_abundance_lookup, corrected_key);
  };

  struct SummaryAggregate {
    uint64_t total = 0;
    uint64_t cache_hit = 0;
    std::set<uint32_t> smallest_cache_slots;
  };
  std::unordered_map<uint64_t, SummaryAggregate> summary_aggregates;
  uint64_t nonwhitelist_total = 0;
  const bool synthesize_summary =
      !output_parameters.summary_metadata_file_path.empty();
  if (synthesize_summary) {
    summary_aggregates.reserve(is_bulk ? 1u
                                       : global_barcode_abundance.size() * 2u);
  }
  for (const auto &reader : readers) {
    if (!synthesize_summary) {
      if (!reader->SkipRemainingSummaryEvidence(&error)) {
        return Failure(error);
      }
      continue;
    }
    AtacSummaryEvidence evidence;
    while (true) {
      bool eof = false;
      if (!reader->ReadNextSummaryEvidence(&evidence, &eof, &error)) {
        return Failure(error);
      }
      if (eof) {
        break;
      }
      uint64_t corrected_key = evidence.raw_barcode_key;
      const BarcodeCorrectionStatus status = correct_summary_barcode(
          evidence.raw_barcode_key, evidence.raw_barcode_n_mask,
          evidence.raw_barcode_qual, &corrected_key);
      if (status == BarcodeCorrectionStatus::kRejected) {
        ++nonwhitelist_total;
        continue;
      }
      SummaryAggregate &aggregate = summary_aggregates[corrected_key];
      ++aggregate.total;
      if (evidence.cache_slot1 >= 0 || evidence.cache_slot2 >= 0) {
        ++aggregate.cache_hit;
      }
      auto add_cache_slot = [&](int32_t slot) {
        if (slot < 0 ||
            (contract.flags & kAtacMergeableSummaryCardinality) == 0) {
          return;
        }
        aggregate.smallest_cache_slots.insert(static_cast<uint32_t>(slot));
        if (aggregate.smallest_cache_slots.size() > contract.k_for_minhash) {
          auto largest = aggregate.smallest_cache_slots.end();
          --largest;
          aggregate.smallest_cache_slots.erase(largest);
        }
      };
      add_cache_slot(evidence.cache_slot1);
      add_cache_slot(evidence.cache_slot2);
    }
  }

  MappingParameters parameters = output_parameters;
  parameters.atac_spill_materialization_mode = true;
  parameters.is_bulk_data = is_bulk;
  parameters.remove_pcr_duplicates =
      (contract.flags & kAtacMergeableRemovePcrDuplicates) != 0;
  parameters.remove_pcr_duplicates_at_bulk_level =
      (contract.flags & kAtacMergeableBulkLevelDedup) != 0;
  parameters.Tn5_shift =
      (contract.flags & kAtacMergeableTn5Shift) != 0;
  parameters.Tn5_forward_shift = contract.tn5_forward_shift;
  parameters.Tn5_reverse_shift = contract.tn5_reverse_shift;
  parameters.mapq_threshold = contract.mapq_threshold;
  parameters.allocate_multi_mappings = allocate_multi_mappings;
  parameters.only_output_unique_mappings = !allocate_multi_mappings;
  parameters.multi_mapping_allocation_distance =
      contract.multi_mapping_allocation_distance;
  parameters.multi_mapping_allocation_seed =
      contract.multi_mapping_allocation_seed;
  parameters.max_num_best_mappings =
      static_cast<int>(contract.max_num_best_mappings);
  parameters.output_mappings_not_in_whitelist =
      (contract.flags & kAtacMergeableOutputMappingsNotInWhitelist) != 0;
  parameters.barcode_whitelist_file_path =
      is_bulk ? std::string() : std::string("<spill-contract>");
  parameters.create_mergeable_spill_record_path.clear();

  SequenceBatch reference(static_cast<uint32_t>(contract.references.size()),
                          SequenceEffectiveRange());
  for (uint32_t rid = 0; rid < contract.references.size(); ++rid) {
    reference.AssignLoadedReferenceMetadata(rid, contract.references[rid].name,
                                            contract.references[rid].length);
  }

  struct HeapRecord {
    uint32_t rid = 0;
    AtacSpillRecord mapping;
    uint32_t shard_ordinal = 0;

    bool operator<(const HeapRecord &other) const {
      if (rid != other.rid) {
        return rid > other.rid;
      }
      const bool a_less_b = mapping < other.mapping;
      const bool b_less_a = other.mapping < mapping;
      if (!a_less_b && !b_less_a) {
        return shard_ordinal > other.shard_ordinal;
      }
      return !a_less_b;
    }
  };

  auto read_next = [&](uint32_t ordinal, HeapRecord *output, bool *eof) {
    uint32_t rid = 0;
    AtacSpillRecord mapping;
    const bool read_ok =
        output_parameters.mapping_output_format == MAPPINGFORMAT_BED
            ? readers[ordinal]->ReadNextBed(&rid, &mapping, eof, &error)
            : readers[ordinal]->ReadNext(&rid, &mapping, eof, &error);
    if (!read_ok) {
      return false;
    }
    if (*eof) {
      return true;
    }
    const auto &metadata = readers[ordinal]->metadata();
    if (!GlobalizeAtacSpillReadId(&mapping,
                                  global_read_prefixes[ordinal],
                                  metadata.input_record_count, &error)) {
      return false;
    }
    if (mapping.num_dups_ != 1) {
      error = "ATAC spill materializer received an already-deduplicated record";
      return false;
    }
    output->rid = rid;
    output->mapping = std::move(mapping);
    output->shard_ordinal = ordinal;
    return true;
  };

  uint64_t corrected_barcode_record_count = 0;
  uint64_t rejected_barcode_record_count = 0;
  const bool output_unresolved_barcodes =
      (contract.flags & kAtacMergeableOutputMappingsNotInWhitelist) != 0;
  auto correct_barcode = [&](AtacSpillRecord *mapping, bool *keep) {
    *keep = true;
    if (is_bulk) {
      return true;
    }
    if (!mapping->HasRawBarcodeEvidence() ||
        mapping->raw_barcode_qual_.size() != contract.barcode_length) {
      error = "ATAC spill record lacks complete raw barcode evidence";
      return false;
    }
    uint64_t corrected_key = mapping->cell_barcode_;
    const BarcodeCorrectionStatus status = CorrectPackedBarcode(
        mapping->cell_barcode_, contract.barcode_length,
        mapping->raw_barcode_n_mask_, mapping->raw_barcode_qual_,
        contract.barcode_correction_error_threshold,
        contract.barcode_correction_probability_threshold,
        global_num_sample_barcodes, barcode_abundance_lookup,
        &corrected_key);
    if (status == BarcodeCorrectionStatus::kRejected) {
      ++rejected_barcode_record_count;
      *keep = output_unresolved_barcodes;
      return true;
    }
    if (status == BarcodeCorrectionStatus::kCorrected) {
      ++corrected_barcode_record_count;
    }
    mapping->cell_barcode_ = corrected_key;
    mapping->sam1.cell_barcode_ = corrected_key;
    mapping->sam2.cell_barcode_ = corrected_key;
    return true;
  };

  uint64_t output_fragment_count = 0;
  {
    MappingWriter<AtacSpillRecord> writer(parameters, contract.barcode_length,
                                          std::vector<int>());
    if (parameters.emit_noY_stream || parameters.emit_Y_stream) {
      writer.OpenYFilterStreams();
    }
    writer.OutputHeader(static_cast<uint32_t>(contract.references.size()),
                        reference);
    auto update_summary_count = [&](uint64_t barcode, int field,
                                    uint64_t value) {
      while (value > 0) {
        const int increment = static_cast<int>(std::min<uint64_t>(
            value, static_cast<uint64_t>(std::numeric_limits<int>::max())));
        writer.UpdateSummaryMetadata(barcode, field, increment);
        value -= static_cast<uint64_t>(increment);
      }
    };
    if (synthesize_summary) {
      for (const auto &entry : summary_aggregates) {
        update_summary_count(entry.first, SUMMARY_METADATA_TOTAL,
                             entry.second.total);
        update_summary_count(entry.first, SUMMARY_METADATA_CACHEHIT,
                             entry.second.cache_hit);
        if ((contract.flags & kAtacMergeableSummaryCardinality) != 0 &&
            entry.second.smallest_cache_slots.size() >=
                contract.k_for_minhash) {
          const uint64_t maximum =
              *entry.second.smallest_cache_slots.rbegin();
          const uint64_t cardinality =
              maximum == 0
                  ? uint64_t{0}
                  : (static_cast<uint64_t>(contract.k_for_minhash) *
                         contract.cache_size) /
                            maximum -
                        1u;
          update_summary_count(entry.first, SUMMARY_METADATA_CARDINALITY,
                               cardinality);
        }
      }
      uint64_t remaining_nonwhitelist = nonwhitelist_total;
      while (remaining_nonwhitelist > 0) {
        const int increment = static_cast<int>(std::min<uint64_t>(
            remaining_nonwhitelist,
            static_cast<uint64_t>(std::numeric_limits<int>::max())));
        writer.UpdateSpeicalCategorySummaryMetadata(
            /*nonwhitelist=*/0, SUMMARY_METADATA_TOTAL, increment);
        remaining_nonwhitelist -= static_cast<uint64_t>(increment);
      }
    }

    std::priority_queue<HeapRecord> heap;
    for (uint32_t ordinal = 0; ordinal < readers.size(); ++ordinal) {
      HeapRecord record;
      bool eof = false;
      if (!read_next(ordinal, &record, &eof)) {
        return Failure(error);
      }
      if (!eof) {
        heap.push(std::move(record));
      }
    }

    bool have_last = false;
    uint32_t last_rid = 0;
    AtacSpillRecord last_mapping;
    uint64_t duplicate_count = 0;
    std::vector<std::vector<AtacSpillRecord>> buffered_mappings(
        contract.references.size());
    std::vector<HeapRecord> group;
    std::vector<AtacSpillRecord> barcode_representatives;

    auto emit_last = [&]() {
      if (!have_last) {
        return;
      }
      last_mapping.num_dups_ = static_cast<uint8_t>(std::min<uint64_t>(
          std::numeric_limits<uint8_t>::max(), duplicate_count));
      if (parameters.Tn5_shift) {
        last_mapping.Tn5Shift(parameters.Tn5_forward_shift,
                              parameters.Tn5_reverse_shift);
      }
      if (synthesize_summary && !allocate_multi_mappings) {
        update_summary_count(last_mapping.GetBarcode(),
                             SUMMARY_METADATA_MAPPED, duplicate_count);
        if (last_mapping.mapq_ >= parameters.mapq_threshold) {
          update_summary_count(last_mapping.GetBarcode(),
                               SUMMARY_METADATA_DUP,
                               duplicate_count == 0 ? 0
                                                    : duplicate_count - 1);
        } else {
          update_summary_count(last_mapping.GetBarcode(),
                               SUMMARY_METADATA_LOWMAPQ, duplicate_count);
        }
      }
      if (allocate_multi_mappings) {
        buffered_mappings[last_rid].push_back(std::move(last_mapping));
      } else if (last_mapping.mapq_ >= parameters.mapq_threshold) {
        writer.AppendMaterializedMapping(last_rid, reference, last_mapping);
        ++output_fragment_count;
      }
    };

    while (!heap.empty()) {
      const uint32_t group_rid = heap.top().rid;
      const uint32_t group_start = heap.top().mapping.fragment_start_position_;
      const uint16_t group_length = heap.top().mapping.fragment_length_;
      group.clear();
      while (!heap.empty() && heap.top().rid == group_rid &&
             heap.top().mapping.fragment_start_position_ == group_start &&
             heap.top().mapping.fragment_length_ == group_length) {
        HeapRecord current =
            std::move(const_cast<HeapRecord &>(heap.top()));
        heap.pop();
        const uint32_t source_ordinal = current.shard_ordinal;
        bool keep = true;
        if (!correct_barcode(&current.mapping, &keep)) {
          return Failure(error);
        }
        if (keep) {
          group.push_back(std::move(current));
        }

        HeapRecord next;
        bool eof = false;
        if (!read_next(source_ordinal, &next, &eof)) {
          return Failure(error);
        }
        if (!eof) {
          heap.push(std::move(next));
        }
      }

      std::sort(group.begin(), group.end(),
                [](const HeapRecord &a, const HeapRecord &b) {
                  const bool a_less_b = a.mapping < b.mapping;
                  const bool b_less_a = b.mapping < a.mapping;
                  if (!a_less_b && !b_less_a) {
                    return a.shard_ordinal < b.shard_ordinal;
                  }
                  return a_less_b;
                });
      const bool bulk_level_cell_dedup =
          parameters.remove_pcr_duplicates && !is_bulk &&
          parameters.remove_pcr_duplicates_at_bulk_level;
      if (bulk_level_cell_dedup && !group.empty()) {
        barcode_representatives.clear();
        uint64_t total_coordinate_duplicates = 0;
        size_t begin = 0;
        while (begin < group.size()) {
          size_t end = begin + 1;
          while (end < group.size() &&
                 group[end].mapping == group[end - 1].mapping) {
            ++end;
          }
          AtacSpillRecord representative =
              std::move(group[end - 1].mapping);
          representative.num_dups_ = static_cast<uint8_t>(
              std::min<size_t>(std::numeric_limits<uint8_t>::max(),
                               end - begin));
          barcode_representatives.push_back(std::move(representative));
          total_coordinate_duplicates += end - begin;
          begin = end;
        }
        size_t best = 0;
        uint64_t best_abundance = 0;
        for (size_t i = 0; i < barcode_representatives.size(); ++i) {
          const auto abundance_it = global_barcode_abundance.find(
              barcode_representatives[i].GetBarcode());
          const uint64_t abundance =
              abundance_it == global_barcode_abundance.end()
                  ? uint64_t{0}
                  : abundance_it->second;
          if (i == 0 ||
              barcode_representatives[i].num_dups_ >
                  barcode_representatives[best].num_dups_ ||
              (barcode_representatives[i].num_dups_ ==
                   barcode_representatives[best].num_dups_ &&
               abundance > best_abundance)) {
            best = i;
            best_abundance = abundance;
          }
        }
        emit_last();
        have_last = true;
        last_rid = group_rid;
        last_mapping = std::move(barcode_representatives[best]);
        duplicate_count = total_coordinate_duplicates;
        continue;
      }
      for (HeapRecord &current : group) {
        const bool duplicate =
            have_last && current.rid == last_rid &&
            (current.mapping == last_mapping ||
             (!is_bulk && parameters.remove_pcr_duplicates_at_bulk_level &&
              current.mapping.IsSamePosition(last_mapping)));
        if (parameters.remove_pcr_duplicates && duplicate) {
          ++duplicate_count;
          // Match Chromap's canonical policy: the last record in sort order is
          // the retained representative (highest MAPQ/read id tie-break).
          last_mapping = std::move(current.mapping);
        } else {
          emit_last();
          have_last = true;
          last_rid = current.rid;
          last_mapping = std::move(current.mapping);
          duplicate_count = 1;
        }
      }
    }
    emit_last();
    if (allocate_multi_mappings) {
      uint64_t num_multi_mappings = 0;
      for (const auto &mappings : buffered_mappings) {
        for (const auto &mapping : mappings) {
          if (mapping.mapq_ < 4) {
            ++num_multi_mappings;
          }
        }
      }
      MappingProcessor<AtacSpillRecord> processor(parameters, 4);
      if (num_multi_mappings > 0) {
        processor.AllocateMultiMappings(
            static_cast<uint32_t>(contract.references.size()),
            num_multi_mappings, parameters.multi_mapping_allocation_distance,
            buffered_mappings);
      }
      processor.SortOutputMappings(
          static_cast<uint32_t>(contract.references.size()),
          buffered_mappings);
      for (const auto &mappings : buffered_mappings) {
        for (const auto &mapping : mappings) {
          if (mapping.mapq_ >= parameters.mapq_threshold) {
            ++output_fragment_count;
          }
        }
      }
      writer.OutputMappings(
          static_cast<uint32_t>(contract.references.size()), reference,
          buffered_mappings);
    }
    if (synthesize_summary) {
      writer.OutputSummaryMetadata(
          contract.frip_est_coefficients,
          (contract.flags & kAtacMergeableSummaryCardinality) != 0);
    }
    writer.FinalizeSortedOutput();
    writer.CloseYFilterStreams();
  }

  uint64_t total_input_records = 0;
  for (const auto &reader : readers) {
    if (reader->metadata().input_record_count >
        std::numeric_limits<uint64_t>::max() - total_input_records) {
      return Failure("ATAC spill materializer input record count overflows");
    }
    total_input_records += reader->metadata().input_record_count;
  }
  AtacSpillMaterializationResult result;
  result.ok = true;
  result.message = "ok";
  result.sample_id = contract.sample_id;
  result.input_id = contract.input_id;
  result.shard_count = contract.shard_count;
  result.input_record_count = total_input_records;
  result.spill_record_count = total_spill_records;
  result.corrected_barcode_record_count = corrected_barcode_record_count;
  result.rejected_barcode_record_count = rejected_barcode_record_count;
  result.output_fragment_count = output_fragment_count;
  return result;
}

}  // namespace chromap
