#include "atac_spill_compactor.h"

#include <algorithm>
#include <limits>
#include <memory>
#include <queue>
#include <string>
#include <utility>
#include <vector>
#include <unistd.h>

#include "atac_hot_spill.h"
#include "atac_mergeable_spill.h"

namespace chromap {
namespace {

AtacSpillCompactionResult Failure(const std::string &message) {
  AtacSpillCompactionResult result;
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

bool SameCompactionContract(const AtacMergeableSpillMetadata &a,
                            const AtacMergeableSpillMetadata &b) {
  if (a.schema_mask != b.schema_mask || a.flags != b.flags ||
      a.shard_count != b.shard_count ||
      a.barcode_length != b.barcode_length ||
      a.mapq_threshold != b.mapq_threshold ||
      a.tn5_forward_shift != b.tn5_forward_shift ||
      a.tn5_reverse_shift != b.tn5_reverse_shift ||
      a.barcode_whitelist_fingerprint != b.barcode_whitelist_fingerprint ||
      a.barcode_correction_error_threshold !=
          b.barcode_correction_error_threshold ||
      a.barcode_correction_probability_threshold !=
          b.barcode_correction_probability_threshold ||
      a.multi_mapping_allocation_distance !=
          b.multi_mapping_allocation_distance ||
      a.multi_mapping_allocation_seed != b.multi_mapping_allocation_seed ||
      a.max_num_best_mappings != b.max_num_best_mappings ||
      a.cache_size != b.cache_size || a.k_for_minhash != b.k_for_minhash ||
      a.frip_est_coefficients != b.frip_est_coefficients ||
      a.sample_id != b.sample_id || a.input_id != b.input_id ||
      !SameReferences(a.references, b.references) ||
      a.barcode_abundance_entries.size() !=
          b.barcode_abundance_entries.size()) {
    return false;
  }
  for (size_t i = 0; i < a.barcode_abundance_entries.size(); ++i) {
    if (a.barcode_abundance_entries[i].barcode_key !=
        b.barcode_abundance_entries[i].barcode_key) {
      return false;
    }
  }
  return true;
}

struct InputRun {
  std::string path;
  std::unique_ptr<AtacMergeableSpillReader> reader;
  std::unique_ptr<AtacHotSpillReader> hot_reader;
  uint64_t output_read_prefix = 0;
};

}  // namespace

AtacSpillCompactionResult CompactAtacSpillRecords(
    const std::vector<std::string> &spill_paths,
    const std::string &output_path) {
  if (spill_paths.empty()) {
    return Failure("ATAC raw compactor input set is empty");
  }
  if (output_path.empty() || output_path == "-" ||
      output_path == "/dev/stdout" || output_path == "/dev/stderr") {
    return Failure("ATAC raw compactor requires a regular output path");
  }
  if (access(output_path.c_str(), F_OK) == 0 ||
      access(AtacHotSpillSidecarPath(output_path).c_str(), F_OK) == 0) {
    return Failure("ATAC raw compactor refuses to replace an existing output");
  }
  for (const std::string &path : spill_paths) {
    if (path == output_path || AtacHotSpillSidecarPath(path) == output_path) {
      return Failure("ATAC raw compactor output aliases an input");
    }
  }

  std::vector<InputRun> inputs;
  inputs.reserve(spill_paths.size());
  std::string error;
  for (const std::string &path : spill_paths) {
    InputRun input;
    input.path = path;
    input.reader.reset(new AtacMergeableSpillReader());
    if (!input.reader->Open(path, &error)) {
      return Failure(error);
    }
    inputs.push_back(std::move(input));
  }
  std::sort(inputs.begin(), inputs.end(),
            [](const InputRun &a, const InputRun &b) {
              return a.reader->metadata().shard_ordinal <
                     b.reader->metadata().shard_ordinal;
            });

  const AtacMergeableSpillMetadata &contract =
      inputs.front().reader->metadata();
  AtacMergeableSpillMetadata output_metadata = contract;
  output_metadata.shard_ordinal = contract.shard_ordinal;
  output_metadata.shard_span = 0;
  output_metadata.input_record_count = 0;
  output_metadata.summary_evidence_count = 0;
  output_metadata.local_num_sample_barcodes = 0;
  for (auto &entry : output_metadata.barcode_abundance_entries) {
    entry.count = 0;
  }

  const bool late_bound =
      (contract.flags & kAtacMergeableReadRangeLateBound) != 0;
  const bool has_hot_sidecar =
      (contract.flags & kAtacMergeableHasHotSidecar) != 0;
  uint32_t next_ordinal = contract.shard_ordinal;
  uint64_t next_precomputed_read = contract.first_global_read_ordinal;
  uint64_t total_spill_records = 0;
  for (InputRun &input : inputs) {
    const AtacMergeableSpillMetadata &metadata = input.reader->metadata();
    if (!SameCompactionContract(contract, metadata)) {
      return Failure("ATAC raw compactor input contracts disagree: " +
                     input.path);
    }
    if (metadata.shard_ordinal != next_ordinal) {
      return Failure("ATAC raw compactor inputs are not adjacent ordinal "
                     "ranges");
    }
    if (!late_bound &&
        metadata.first_global_read_ordinal != next_precomputed_read) {
      return Failure("ATAC raw compactor read ranges are not contiguous");
    }
    if (metadata.shard_span >
            std::numeric_limits<uint32_t>::max() -
                output_metadata.shard_span ||
        metadata.input_record_count >
            std::numeric_limits<uint64_t>::max() -
                output_metadata.input_record_count ||
        metadata.summary_evidence_count >
            std::numeric_limits<uint64_t>::max() -
                output_metadata.summary_evidence_count ||
        metadata.local_num_sample_barcodes >
            std::numeric_limits<uint64_t>::max() -
                output_metadata.local_num_sample_barcodes ||
        input.reader->expected_record_count() >
            std::numeric_limits<uint64_t>::max() - total_spill_records) {
      return Failure("ATAC raw compactor aggregate count overflows");
    }
    input.output_read_prefix = output_metadata.input_record_count;
    output_metadata.shard_span += metadata.shard_span;
    output_metadata.input_record_count += metadata.input_record_count;
    output_metadata.summary_evidence_count +=
        metadata.summary_evidence_count;
    output_metadata.local_num_sample_barcodes +=
        metadata.local_num_sample_barcodes;
    total_spill_records += input.reader->expected_record_count();
    next_ordinal += metadata.shard_span;
    if (!late_bound) {
      next_precomputed_read += metadata.input_record_count;
    }
    for (size_t i = 0;
         i < output_metadata.barcode_abundance_entries.size(); ++i) {
      const uint64_t count = metadata.barcode_abundance_entries[i].count;
      uint64_t &aggregate =
          output_metadata.barcode_abundance_entries[i].count;
      if (count > std::numeric_limits<uint64_t>::max() - aggregate) {
        return Failure("ATAC raw compactor barcode abundance overflows");
      }
      aggregate += count;
    }
    if (has_hot_sidecar) {
      input.hot_reader.reset(new AtacHotSpillReader());
      if (!input.hot_reader->Open(
              input.path, metadata, input.reader->expected_record_count(),
              &error)) {
        return Failure(error);
      }
    }
  }
  if (next_ordinal > contract.shard_count ||
      output_metadata.shard_span == 0) {
    return Failure("ATAC raw compactor ordinal coverage is invalid");
  }
  if (late_bound) {
    output_metadata.first_global_read_ordinal = 0;
  }

  AtacMergeableSpillWriter writer;
  if (!writer.Open(output_path, output_metadata, &error)) {
    return Failure(error);
  }
  std::unique_ptr<AtacHotSpillWriter> hot_writer;
  if (has_hot_sidecar) {
    hot_writer.reset(new AtacHotSpillWriter());
    if (!hot_writer->Open(AtacHotSpillSidecarPath(output_path),
                          output_metadata, &error)) {
      return Failure(error);
    }
  }

  for (InputRun &input : inputs) {
    while (true) {
      AtacSummaryEvidence evidence;
      bool eof = false;
      if (!input.reader->ReadNextSummaryEvidence(&evidence, &eof, &error)) {
        return Failure(error);
      }
      if (eof) {
        break;
      }
      if (!writer.AppendSummaryEvidence(evidence, &error)) {
        return Failure(error);
      }
    }
  }

  struct HeapRecord {
    uint32_t rid = 0;
    AtacSpillRecord record;
    size_t input_index = 0;

    bool operator<(const HeapRecord &other) const {
      if (rid != other.rid) {
        return rid > other.rid;
      }
      return other.record < record;
    }
  };
  auto read_next = [&](size_t input_index, HeapRecord *output, bool *eof) {
    InputRun &input = inputs[input_index];
    uint32_t rid = 0;
    AtacSpillRecord record;
    if (!input.reader->ReadNext(&rid, &record, eof, &error)) {
      return false;
    }
    if (*eof) {
      return true;
    }
    const AtacMergeableSpillMetadata &metadata = input.reader->metadata();
    if (record.num_dups_ != 1 ||
        !GlobalizeAtacSpillReadId(&record, input.output_read_prefix,
                                  metadata.input_record_count, &error)) {
      if (record.num_dups_ != 1 && error.empty()) {
        error = "ATAC raw compactor received an already-deduplicated record";
      }
      return false;
    }
    output->rid = rid;
    output->record = std::move(record);
    output->input_index = input_index;
    return true;
  };

  std::priority_queue<HeapRecord> heap;
  for (size_t i = 0; i < inputs.size(); ++i) {
    HeapRecord record;
    bool eof = false;
    if (!read_next(i, &record, &eof)) {
      return Failure(error);
    }
    if (!eof) {
      heap.push(std::move(record));
    }
  }
  uint64_t records_written = 0;
  while (!heap.empty()) {
    HeapRecord current = std::move(const_cast<HeapRecord &>(heap.top()));
    heap.pop();
    const size_t source = current.input_index;
    if (!writer.Append(current.rid, current.record, &error) ||
        (hot_writer &&
         !hot_writer->Append(current.rid, current.record, &error))) {
      return Failure(error);
    }
    ++records_written;

    HeapRecord next;
    bool eof = false;
    if (!read_next(source, &next, &eof)) {
      return Failure(error);
    }
    if (!eof) {
      heap.push(std::move(next));
    }
  }
  if (records_written != total_spill_records) {
    return Failure("ATAC raw compactor record count disagrees with inputs");
  }
  if (hot_writer && !hot_writer->Finalize(&error)) {
    return Failure(error);
  }
  if (!writer.Finalize(&error)) {
    if (hot_writer) {
      (void)unlink(AtacHotSpillSidecarPath(output_path).c_str());
    }
    return Failure(error);
  }

  AtacSpillCompactionResult result;
  result.ok = true;
  result.sample_id = output_metadata.sample_id;
  result.input_id = output_metadata.input_id;
  result.first_shard_ordinal = output_metadata.shard_ordinal;
  result.shard_span = output_metadata.shard_span;
  result.shard_count = output_metadata.shard_count;
  result.input_record_count = output_metadata.input_record_count;
  result.spill_record_count = records_written;
  result.wrote_hot_sidecar = static_cast<bool>(hot_writer);
  return result;
}

}  // namespace chromap
