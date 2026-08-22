#include "atac_spill_materializer.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <unistd.h>

#include "atac_materialized_binary.h"
#include "atac_hot_spill.h"
#include "atac_mergeable_spill.h"
#include "barcode_correction.h"
#include "mapping_processor.h"
#include "mapping_writer.h"
#include "sequence_batch.h"
#include "utils.h"

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

struct HotPartitionOutcome {
  std::string path;
  std::string error;
  uint64_t corrected = 0;
  uint64_t rejected = 0;
  uint64_t output = 0;
};

// Ephemeral per-reference materializer row for barcodes that do not fit in
// ATMBLK1's direct 32-bit barcode_value. The ordered assembler replaces the
// raw key with a deterministic dictionary id before writing durable ATMBLK1.
#pragma pack(push, 1)
struct HotPartitionWideBarcodeRecord {
  uint64_t barcode_key;
  uint32_t start;
  uint16_t rid;
  uint16_t fragment_length;
  uint8_t duplicate_count;
  uint8_t mapq;
  uint8_t flags;
  uint8_t reserved;
};
#pragma pack(pop)

static_assert(sizeof(HotPartitionWideBarcodeRecord) == 20,
              "wide hot partition row must be 20 bytes");

bool MaterializeHotPartitions(
    const std::vector<std::unique_ptr<AtacHotSpillReader>> &hot_readers,
    const std::vector<std::unique_ptr<AtacMergeableSpillReader>> &parents,
    const std::vector<uint64_t> &global_read_prefixes,
    const AtacMergeableSpillMetadata &contract,
    const MappingParameters &parameters,
    const std::unordered_map<uint64_t, uint64_t> &global_barcode_abundance,
    uint64_t global_num_sample_barcodes,
    const std::string &partition_base,
    AtacMaterializedBinaryWriter *binary_writer,
    uint64_t *corrected_barcode_record_count,
    uint64_t *rejected_barcode_record_count,
    uint64_t *output_fragment_count, std::string *error) {
  if (hot_readers.size() != parents.size() || hot_readers.empty() ||
      binary_writer == nullptr || corrected_barcode_record_count == nullptr ||
      rejected_barcode_record_count == nullptr ||
      output_fragment_count == nullptr || contract.references.size() >
          static_cast<uint64_t>(std::numeric_limits<uint16_t>::max()) + 1u) {
    if (error != nullptr) {
      *error = "invalid parallel ATAC hot partition request";
    }
    return false;
  }
  const bool is_bulk =
      (contract.schema_mask & kAtacSpillSchemaIsBulk) != 0;
  const bool use_wide_barcode_partition =
      !is_bulk && contract.barcode_length > 16;
  const bool output_unresolved_barcodes =
      (contract.flags & kAtacMergeableOutputMappingsNotInWhitelist) != 0;
  const bool bulk_level_cell_dedup =
      parameters.remove_pcr_duplicates && !is_bulk &&
      parameters.remove_pcr_duplicates_at_bulk_level;
  const int hot_threads = std::max(1, parameters.num_threads);
  std::vector<HotPartitionOutcome> outcomes(contract.references.size());
  std::atomic<bool> failed(false);

#pragma omp parallel for schedule(dynamic, 1) num_threads(hot_threads)
  for (int64_t rid_index = 0;
       rid_index < static_cast<int64_t>(contract.references.size());
       ++rid_index) {
    if (failed.load(std::memory_order_relaxed)) {
      continue;
    }
    const uint32_t rid = static_cast<uint32_t>(rid_index);
    HotPartitionOutcome &outcome = outcomes[rid];
    uint64_t partition_input_records = 0;
    for (const auto &reader : hot_readers) {
      const uint64_t count = reader->PartitionRecordCount(rid);
      if (count > std::numeric_limits<uint64_t>::max() -
                      partition_input_records) {
        outcome.error = "ATAC hot partition record count overflows";
        failed.store(true, std::memory_order_relaxed);
        break;
      }
      partition_input_records += count;
    }
    if (failed.load(std::memory_order_relaxed)) {
      continue;
    }
    if (partition_input_records == 0) {
      continue;
    }
    outcome.path = partition_base + ".hotpart." +
                   std::to_string(static_cast<uint64_t>(getpid())) + "." +
                   std::to_string(rid);
    FILE *partition = fopen(outcome.path.c_str(), "wb");
    if (partition == nullptr) {
      outcome.error = "cannot create ATAC hot materializer partition " +
                      outcome.path + ": " + std::strerror(errno);
      failed.store(true, std::memory_order_relaxed);
      continue;
    }
    (void)setvbuf(partition, nullptr, _IOFBF, 8u * 1024u * 1024u);

    struct HeapRecord {
      AtacHotSpillDecodedRecord mapping;
      uint32_t shard_ordinal = 0;

      bool operator<(const HeapRecord &other) const {
        const bool a_less_b = mapping < other.mapping;
        const bool b_less_a = other.mapping < mapping;
        if (!a_less_b && !b_less_a) {
          return shard_ordinal > other.shard_ordinal;
        }
        return !a_less_b;
      }
    };

    std::vector<std::unique_ptr<AtacHotSpillPartitionCursor>> cursors(
        hot_readers.size());
    std::priority_queue<HeapRecord> heap;
    bool partition_ok = true;
    auto read_next = [&](uint32_t ordinal, HeapRecord *next, bool *eof) {
      AtacHotSpillDecodedRecord mapping;
      if (!cursors[ordinal]->ReadNext(&mapping, eof, &outcome.error)) {
        return false;
      }
      if (*eof) {
        return true;
      }
      const uint64_t local_record_count =
          parents[ordinal]->metadata().input_record_count;
      if (mapping.read_id_ >= local_record_count ||
          mapping.read_id_ > std::numeric_limits<uint64_t>::max() -
                                 global_read_prefixes[ordinal]) {
        outcome.error = "ATAC hot spill local read id is invalid";
        return false;
      }
      mapping.read_id_ += global_read_prefixes[ordinal];
      next->mapping = std::move(mapping);
      next->shard_ordinal = ordinal;
      return true;
    };
    for (uint32_t ordinal = 0; ordinal < hot_readers.size(); ++ordinal) {
      cursors[ordinal].reset(new AtacHotSpillPartitionCursor());
      if (!hot_readers[ordinal]->OpenPartition(
              rid, cursors[ordinal].get(), &outcome.error)) {
        partition_ok = false;
        break;
      }
      HeapRecord first;
      bool eof = false;
      if (!read_next(ordinal, &first, &eof)) {
        partition_ok = false;
        break;
      }
      if (!eof) {
        heap.push(std::move(first));
      }
    }

    const auto barcode_abundance_lookup = [&](uint64_t key, bool *found) {
      const auto it = global_barcode_abundance.find(key);
      *found = it != global_barcode_abundance.end();
      return *found ? it->second : uint64_t{0};
    };
    auto correct_barcode = [&](AtacHotSpillDecodedRecord *mapping,
                               bool *keep) {
      *keep = true;
      if (is_bulk) {
        return true;
      }
      if (!mapping->has_raw_barcode_evidence ||
          mapping->raw_barcode_qual.size() != contract.barcode_length) {
        outcome.error =
            "ATAC hot spill record lacks complete barcode evidence";
        return false;
      }
      uint64_t corrected_key = mapping->cell_barcode_;
      const BarcodeCorrectionStatus status = CorrectPackedBarcode(
          mapping->cell_barcode_, contract.barcode_length,
          mapping->raw_barcode_n_mask, mapping->raw_barcode_qual,
          contract.barcode_correction_error_threshold,
          contract.barcode_correction_probability_threshold,
          global_num_sample_barcodes, barcode_abundance_lookup,
          &corrected_key);
      if (status == BarcodeCorrectionStatus::kRejected) {
        ++outcome.rejected;
        *keep = output_unresolved_barcodes;
        return true;
      }
      if (status == BarcodeCorrectionStatus::kCorrected) {
        ++outcome.corrected;
      }
      mapping->cell_barcode_ = corrected_key;
      return true;
    };

    bool have_last = false;
    AtacHotSpillDecodedRecord last_mapping;
    uint64_t duplicate_count = 0;
    std::vector<HeapRecord> group;
    std::vector<AtacHotSpillDecodedRecord> barcode_representatives;
    auto emit_last = [&]() {
      if (!have_last) {
        return true;
      }
      last_mapping.num_dups_ = static_cast<uint8_t>(std::min<uint64_t>(
          std::numeric_limits<uint8_t>::max(), duplicate_count));
      if (parameters.Tn5_shift) {
        last_mapping.Tn5Shift(parameters.Tn5_forward_shift,
                              parameters.Tn5_reverse_shift);
      }
      if (last_mapping.mapq_ < parameters.mapq_threshold) {
        return true;
      }
      if ((!is_bulk && !use_wide_barcode_partition &&
           last_mapping.cell_barcode_ >
               std::numeric_limits<uint32_t>::max()) ||
          !IsAtacSpillFragmentLengthRepresentable(
              last_mapping.fragment_length_)) {
        outcome.error = "ATAC hot partition output exceeds compact encoding";
        return false;
      }
      AtacMaterializedBinaryRecordV1 encoded = {};
      encoded.start = last_mapping.fragment_start_position_;
      encoded.barcode_value =
          is_bulk ? 0u : static_cast<uint32_t>(last_mapping.cell_barcode_);
      encoded.rid = static_cast<uint16_t>(rid);
      encoded.fragment_length = last_mapping.fragment_length_;
      encoded.duplicate_count = static_cast<uint8_t>(std::min<uint64_t>(
          std::numeric_limits<uint8_t>::max(), duplicate_count));
      encoded.mapq = last_mapping.mapq_;
      encoded.flags = last_mapping.IsPositiveStrand()
                          ? kAtacMaterializedBinaryPositiveStrand
                          : uint8_t{0};
      bool wrote = false;
      if (use_wide_barcode_partition) {
        HotPartitionWideBarcodeRecord wide = {};
        wide.start = encoded.start;
        wide.barcode_key = last_mapping.cell_barcode_;
        wide.rid = encoded.rid;
        wide.fragment_length = encoded.fragment_length;
        wide.duplicate_count = encoded.duplicate_count;
        wide.mapq = encoded.mapq;
        wide.flags = encoded.flags;
        wrote = fwrite(&wide, sizeof(wide), 1, partition) == 1;
      } else {
        wrote = fwrite(&encoded, sizeof(encoded), 1, partition) == 1;
      }
      if (!wrote) {
        outcome.error = "cannot write ATAC hot materializer partition " +
                        outcome.path;
        return false;
      }
      ++outcome.output;
      return true;
    };

    while (partition_ok && !heap.empty()) {
      const uint32_t group_start =
          heap.top().mapping.fragment_start_position_;
      const uint16_t group_length = heap.top().mapping.fragment_length_;
      group.clear();
      while (!heap.empty() &&
             heap.top().mapping.fragment_start_position_ == group_start &&
             heap.top().mapping.fragment_length_ == group_length) {
        HeapRecord current =
            std::move(const_cast<HeapRecord &>(heap.top()));
        heap.pop();
        const uint32_t source_ordinal = current.shard_ordinal;
        bool keep = true;
        if (!correct_barcode(&current.mapping, &keep)) {
          partition_ok = false;
          break;
        }
        if (keep) {
          group.push_back(std::move(current));
        }
        HeapRecord next;
        bool eof = false;
        if (!read_next(source_ordinal, &next, &eof)) {
          partition_ok = false;
          break;
        }
        if (!eof) {
          heap.push(std::move(next));
        }
      }
      if (!partition_ok) {
        break;
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
          AtacHotSpillDecodedRecord representative =
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
        if (!emit_last()) {
          partition_ok = false;
          break;
        }
        have_last = true;
        last_mapping = std::move(barcode_representatives[best]);
        duplicate_count = total_coordinate_duplicates;
        continue;
      }
      for (HeapRecord &current : group) {
        const bool duplicate =
            have_last &&
            (current.mapping == last_mapping ||
             (!is_bulk && parameters.remove_pcr_duplicates_at_bulk_level &&
              current.mapping.IsSamePosition(last_mapping)));
        if (parameters.remove_pcr_duplicates && duplicate) {
          ++duplicate_count;
          last_mapping = std::move(current.mapping);
        } else {
          if (!emit_last()) {
            partition_ok = false;
            break;
          }
          have_last = true;
          last_mapping = std::move(current.mapping);
          duplicate_count = 1;
        }
      }
    }
    if (partition_ok && !emit_last()) {
      partition_ok = false;
    }
    if (fclose(partition) != 0 && partition_ok) {
      outcome.error = "cannot close ATAC hot materializer partition " +
                      outcome.path;
      partition_ok = false;
    }
    if (!partition_ok) {
      failed.store(true, std::memory_order_relaxed);
    }
  }

  auto cleanup = [&]() {
    for (const auto &outcome : outcomes) {
      if (!outcome.path.empty()) {
        unlink(outcome.path.c_str());
      }
    }
  };
  if (failed.load(std::memory_order_relaxed)) {
    if (error != nullptr) {
      *error = "parallel ATAC hot partition merge failed";
      for (const auto &outcome : outcomes) {
        if (!outcome.error.empty()) {
          *error = outcome.error;
          break;
        }
      }
    }
    cleanup();
    return false;
  }

  std::vector<AtacMaterializedBinaryRecordV1> encoded(
      kAtacMaterializedBinaryTargetRecordsPerBlock);
  std::vector<HotPartitionWideBarcodeRecord> wide_encoded;
  std::vector<uint64_t> wide_barcode_keys;
  if (use_wide_barcode_partition) {
    wide_encoded.resize(kAtacMaterializedBinaryTargetRecordsPerBlock);
    wide_barcode_keys.resize(kAtacMaterializedBinaryTargetRecordsPerBlock);
  }
  for (const auto &outcome : outcomes) {
    if (outcome.path.empty()) {
      continue;
    }
    FILE *partition = fopen(outcome.path.c_str(), "rb");
    if (partition == nullptr) {
      if (error != nullptr) {
        *error = "cannot reopen ATAC hot materializer partition " +
                 outcome.path;
      }
      cleanup();
      return false;
    }
    while (true) {
      const size_t count =
          use_wide_barcode_partition
              ? fread(wide_encoded.data(), sizeof(wide_encoded[0]),
                      wide_encoded.size(), partition)
              : fread(encoded.data(), sizeof(encoded[0]), encoded.size(),
                      partition);
      if (count != 0 && use_wide_barcode_partition) {
        for (size_t i = 0; i < count; ++i) {
          encoded[i].start = wide_encoded[i].start;
          encoded[i].barcode_value = 0;
          encoded[i].rid = wide_encoded[i].rid;
          encoded[i].fragment_length = wide_encoded[i].fragment_length;
          encoded[i].duplicate_count = wide_encoded[i].duplicate_count;
          encoded[i].mapq = wide_encoded[i].mapq;
          encoded[i].flags = wide_encoded[i].flags;
          encoded[i].reserved = wide_encoded[i].reserved;
          wide_barcode_keys[i] = wide_encoded[i].barcode_key;
        }
      }
      const bool appended =
          count == 0 ||
          (use_wide_barcode_partition
               ? binary_writer->AppendEncodedRecordsWithRawBarcodeKeys(
                     encoded.data(), wide_barcode_keys.data(), count, error)
               : binary_writer->AppendEncodedRecordsWithRawBarcodes(
                     encoded.data(), count, error));
      if (!appended) {
        fclose(partition);
        cleanup();
        return false;
      }
      const size_t capacity = use_wide_barcode_partition
                                  ? wide_encoded.size()
                                  : encoded.size();
      if (count != capacity) {
        if (ferror(partition) != 0) {
          if (error != nullptr) {
            *error = "cannot read ATAC hot materializer partition " +
                     outcome.path;
          }
          fclose(partition);
          cleanup();
          return false;
        }
        break;
      }
    }
    fclose(partition);
    *corrected_barcode_record_count += outcome.corrected;
    *rejected_barcode_record_count += outcome.rejected;
    *output_fragment_count += outcome.output;
  }
  cleanup();
  return true;
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
  const bool has_hot_sidecar =
      (contract.flags & kAtacMergeableHasHotSidecar) != 0;
  std::vector<std::unique_ptr<AtacHotSpillReader>> hot_readers;
  if (has_hot_sidecar &&
      output_parameters.mapping_output_format == MAPPINGFORMAT_BED) {
    hot_readers.resize(readers.size());
    for (uint32_t ordinal = 0; ordinal < readers.size(); ++ordinal) {
      hot_readers[ordinal].reset(new AtacHotSpillReader());
      if (!hot_readers[ordinal]->Open(
              paths_by_ordinal[ordinal], readers[ordinal]->metadata(),
              readers[ordinal]->expected_record_count(), &error)) {
        return Failure(error);
      }
    }
  }
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
  bool used_parallel_hot_spill = false;
  double merge_output_seconds = 0.0;
  double terminal_bed_export_seconds = 0.0;
  {
    const double merge_output_start = GetRealTime();
    const bool binary_bed_output =
        parameters.mapping_output_format == MAPPINGFORMAT_BED;
    const bool preserve_materialized_binary =
        !parameters.atac_materialized_binary_output_file_path.empty();
    const std::string materialized_binary_path =
        preserve_materialized_binary
            ? parameters.atac_materialized_binary_output_file_path
            : parameters.mapping_output_file_path + ".materialized." +
                  std::to_string(static_cast<uint64_t>(getpid())) + ".atmb1";
    AtacMaterializedBinaryWriter binary_writer;
    if (binary_bed_output) {
      AtacMaterializedBinaryMetadata binary_metadata;
      binary_metadata.is_bulk = is_bulk;
      binary_metadata.use_barcode_dictionary =
          !parameters.barcode_translate_table_file_path.empty() ||
          (!is_bulk && contract.barcode_length > 16);
      binary_metadata.barcode_length = contract.barcode_length;
      binary_metadata.shard_count = contract.shard_count;
      binary_metadata.sample_id = contract.sample_id;
      binary_metadata.input_id = contract.input_id;
      binary_metadata.references = contract.references;
      if (!binary_writer.Open(materialized_binary_path, binary_metadata,
                              &error)) {
        return Failure(error);
      }
    }
    MappingParameters writer_parameters = parameters;
    if (binary_bed_output) {
      // Summary synthesis still uses the canonical MappingWriter counters,
      // but the BED merge hot path must never construct text rows.
      writer_parameters.mapping_output_file_path = "/dev/null";
    }
    MappingWriter<AtacSpillRecord> writer(writer_parameters,
                                          contract.barcode_length,
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

    const bool use_parallel_hot_spill =
        binary_bed_output && has_hot_sidecar && !allocate_multi_mappings &&
        !synthesize_summary;
    if (use_parallel_hot_spill) {
      used_parallel_hot_spill = true;
      if (!MaterializeHotPartitions(
              hot_readers, readers, global_read_prefixes, contract, parameters,
              global_barcode_abundance, global_num_sample_barcodes,
              materialized_binary_path, &binary_writer,
              &corrected_barcode_record_count,
              &rejected_barcode_record_count, &output_fragment_count,
              &error)) {
        return Failure(error);
      }
    } else {
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
        return true;
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
        if (binary_bed_output) {
          if (!binary_writer.Append(last_rid, last_mapping, duplicate_count,
                                    &error)) {
            return false;
          }
        } else {
          writer.AppendMaterializedMapping(last_rid, reference, last_mapping);
        }
        ++output_fragment_count;
      }
      return true;
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
        if (!emit_last()) {
          return Failure(error);
        }
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
          if (!emit_last()) {
            return Failure(error);
          }
          have_last = true;
          last_rid = current.rid;
          last_mapping = std::move(current.mapping);
          duplicate_count = 1;
        }
      }
    }
    if (!emit_last()) {
      return Failure(error);
    }
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
      for (uint32_t rid = 0; rid < buffered_mappings.size(); ++rid) {
        for (const auto &mapping : buffered_mappings[rid]) {
          if (mapping.mapq_ >= parameters.mapq_threshold) {
            if (binary_bed_output &&
                !binary_writer.Append(rid, mapping, mapping.num_dups_,
                                      &error)) {
              return Failure(error);
            }
            ++output_fragment_count;
          }
        }
      }
      if (!binary_bed_output) {
        writer.OutputMappings(
            static_cast<uint32_t>(contract.references.size()), reference,
            buffered_mappings);
      }
    }
    }
    if (synthesize_summary) {
      writer.OutputSummaryMetadata(
          contract.frip_est_coefficients,
          (contract.flags & kAtacMergeableSummaryCardinality) != 0);
    }
    writer.FinalizeSortedOutput();
    writer.CloseYFilterStreams();
    if (binary_bed_output) {
      if (!binary_writer.Finalize(&error)) {
        return Failure(error);
      }
      merge_output_seconds = GetRealTime() - merge_output_start;
      const double terminal_bed_start = GetRealTime();
      if (!ExportAtacMaterializedBinaryToBed(
              materialized_binary_path, parameters.mapping_output_file_path,
              parameters, &error)) {
        return Failure(error);
      }
      terminal_bed_export_seconds = GetRealTime() - terminal_bed_start;
      if (!preserve_materialized_binary) {
        unlink(materialized_binary_path.c_str());
      }
    } else {
      merge_output_seconds = GetRealTime() - merge_output_start;
    }
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
  result.used_parallel_hot_spill = used_parallel_hot_spill;
  result.merge_output_seconds = merge_output_seconds;
  result.terminal_bed_export_seconds = terminal_bed_export_seconds;
  return result;
}

}  // namespace chromap
