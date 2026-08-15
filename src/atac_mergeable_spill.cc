#include "atac_mergeable_spill.h"

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <unistd.h>

namespace chromap {
namespace {

constexpr uint32_t kMaximumMetadataBytes = 64u * 1024u * 1024u;
constexpr uint32_t kMaximumRecordBytes = 64u * 1024u * 1024u;
constexpr size_t kBedStreamBufferBytes = 8u * 1024u * 1024u;

bool WriteBytes(FILE *file, const void *data, size_t size) {
  return size == 0 || fwrite(data, 1, size, file) == size;
}

bool ReadBytes(FILE *file, void *data, size_t size) {
  return size == 0 || fread(data, 1, size, file) == size;
}

bool AddMetadataBytes(uint64_t add, uint64_t *total) {
  if (add > kMaximumMetadataBytes ||
      *total > static_cast<uint64_t>(kMaximumMetadataBytes) - add) {
    return false;
  }
  *total += add;
  return true;
}

class MemoryCursor {
 public:
  MemoryCursor(const char *data, size_t size)
      : cursor_(data), end_(data + size) {}

  template <typename T>
  bool Read(T *value) {
    if (value == nullptr || Remaining() < sizeof(T)) {
      return false;
    }
    memcpy(value, cursor_, sizeof(T));
    cursor_ += sizeof(T);
    return true;
  }

  bool Skip(size_t size) {
    if (Remaining() < size) {
      return false;
    }
    cursor_ += size;
    return true;
  }

  const char *current() const { return cursor_; }
  size_t Remaining() const { return static_cast<size_t>(end_ - cursor_); }

 private:
  const char *cursor_;
  const char *end_;
};

bool SkipSerializedSamMapping(MemoryCursor *cursor, uint64_t expected_read_id) {
  static_assert(sizeof(int) == sizeof(uint32_t),
                "ATAC spill SAM codec requires 32-bit int");
  uint64_t read_id = 0;
  uint16_t read_name_bytes = 0;
  if (!cursor->Read(&read_id) || read_id != expected_read_id ||
      !cursor->Read(&read_name_bytes) || !cursor->Skip(read_name_bytes)) {
    return false;
  }

  uint64_t cell_barcode = 0;
  uint8_t num_dups = 0;
  int64_t pos = 0;
  int rid = 0;
  int64_t mate_pos = 0;
  int mate_rid = 0;
  int template_length = 0;
  int flags = 0;
  uint32_t packed_fields = 0;
  int num_cigar = 0;
  if (!cursor->Read(&cell_barcode) || !cursor->Read(&num_dups) ||
      !cursor->Read(&pos) || !cursor->Read(&rid) ||
      !cursor->Read(&mate_pos) || !cursor->Read(&mate_rid) ||
      !cursor->Read(&template_length) || !cursor->Read(&flags) ||
      !cursor->Read(&packed_fields) || !cursor->Read(&num_cigar) ||
      num_cigar < 0 ||
      static_cast<uint64_t>(num_cigar) >
          std::numeric_limits<size_t>::max() / sizeof(uint32_t) ||
      !cursor->Skip(static_cast<size_t>(num_cigar) * sizeof(uint32_t))) {
    return false;
  }

  uint16_t md_bytes = 0;
  uint16_t sequence_bytes = 0;
  if (!cursor->Read(&md_bytes) || !cursor->Skip(md_bytes) ||
      !cursor->Read(&sequence_bytes) ||
      static_cast<size_t>(sequence_bytes) >
          std::numeric_limits<size_t>::max() / 2u ||
      !cursor->Skip(static_cast<size_t>(sequence_bytes) * 2u)) {
    return false;
  }
  return true;
}

bool DecodeBedPayload(const char *payload, size_t payload_bytes,
                      uint16_t file_schema, uint32_t barcode_length,
                      AtacSpillRecord *record) {
  MemoryCursor cursor(payload, payload_bytes);
  uint32_t magic = 0;
  uint16_t version = 0;
  uint16_t payload_mask = 0;
  if (!cursor.Read(&magic) || !cursor.Read(&version) ||
      !cursor.Read(&payload_mask) || magic != kAtacSpillPayloadMagic ||
      version != kAtacSpillRecordCodecVersion) {
    return false;
  }
  const uint16_t optional_mask = static_cast<uint16_t>(
      kAtacSpillSchemaHasBamPair | kAtacSpillSchemaHasRawBarcodeEvidence);
  if ((payload_mask & optional_mask) != (file_schema & optional_mask)) {
    return false;
  }

  uint8_t mapq = 0;
  uint8_t direction = 0;
  uint8_t is_unique = 0;
  if (!cursor.Read(&record->read_id_) ||
      !cursor.Read(&record->cell_barcode_) ||
      !cursor.Read(&record->fragment_start_position_) ||
      !cursor.Read(&record->fragment_length_) || !cursor.Read(&mapq) ||
      !cursor.Read(&direction) || !cursor.Read(&is_unique) ||
      !cursor.Read(&record->num_dups_) ||
      !cursor.Read(&record->positive_alignment_length_) ||
      !cursor.Read(&record->negative_alignment_length_)) {
    return false;
  }
  record->mapq_ = mapq;
  record->direction_ = direction;
  record->is_unique_ = is_unique;

  if ((file_schema & kAtacSpillSchemaHasBamPair) != 0 &&
      (!SkipSerializedSamMapping(&cursor, record->read_id_) ||
       !SkipSerializedSamMapping(&cursor, record->read_id_))) {
    return false;
  }

  record->raw_barcode_n_mask_ = 0;
  record->raw_barcode_qual_.clear();
  if ((file_schema & kAtacSpillSchemaHasRawBarcodeEvidence) != 0) {
    uint32_t quality_bytes = 0;
    if (!cursor.Read(&record->raw_barcode_n_mask_) ||
        !cursor.Read(&quality_bytes) || quality_bytes != barcode_length ||
        quality_bytes > AtacBarcodeQuality::kCapacity ||
        cursor.Remaining() < quality_bytes) {
      return false;
    }
    record->raw_barcode_qual_.assign(cursor.current(), quality_bytes);
    if (!cursor.Skip(quality_bytes)) {
      return false;
    }
  }
  if (cursor.Remaining() != 0) {
    return false;
  }

  // The compact BED record deliberately does not claim to own decoded BAM
  // mates. Row-level Y and raw-barcode flags remain available; bulk comes from
  // the authoritative file schema.
  record->prefix_flags_ = static_cast<uint16_t>(
      (payload_mask & kAtacSpillSchemaHasYHit) |
      (file_schema &
       (kAtacSpillSchemaIsBulk | kAtacSpillSchemaHasRawBarcodeEvidence)));
  return true;
}

}  // namespace

AtacMergeableSpillWriter::~AtacMergeableSpillWriter() {
  if (file_ != nullptr) {
    fclose(file_);
    file_ = nullptr;
  }
  if (!temporary_path_.empty() && !finalized_) {
    unlink(temporary_path_.c_str());
  }
}

bool AtacMergeableSpillWriter::Fail(const std::string &message,
                                    std::string *error) {
  if (error != nullptr) {
    *error = message;
  }
  return false;
}

bool AtacMergeableSpillWriter::Open(
    const std::string &path, const AtacMergeableSpillMetadata &metadata,
    std::string *error) {
  if (file_ != nullptr || finalized_) {
    return Fail("ATAC mergeable spill writer is already open", error);
  }
  if (path.empty() || path == "-" || path == "/dev/stdout" ||
      path == "/dev/stderr") {
    return Fail("ATAC mergeable spill requires a regular output path", error);
  }
  if (metadata.sample_id.empty() || metadata.input_id.empty()) {
    return Fail("ATAC mergeable spill sample and input ids must be non-empty",
                error);
  }
  if (metadata.shard_count == 0 ||
      metadata.shard_ordinal >= metadata.shard_count) {
    return Fail("ATAC mergeable spill shard ordinal/count is invalid", error);
  }
  if ((metadata.schema_mask & kAtacSpillSchemaHasBamPair) == 0) {
    return Fail("ATAC mergeable spill v3 requires BAM-pair payloads", error);
  }
  const bool is_bulk =
      (metadata.schema_mask & kAtacSpillSchemaIsBulk) != 0;
  const bool has_raw_barcode =
      (metadata.schema_mask & kAtacSpillSchemaHasRawBarcodeEvidence) != 0;
  if (!is_bulk &&
      (metadata.barcode_length == 0 || metadata.barcode_length > 32 ||
       metadata.barcode_correction_error_threshold > 2 ||
       !std::isfinite(metadata.barcode_correction_probability_threshold) ||
       metadata.barcode_correction_probability_threshold < 0.0 ||
       metadata.barcode_correction_probability_threshold >= 1.0)) {
    return Fail("invalid ATAC mergeable spill barcode correction policy",
                error);
  }
  if (!is_bulk && (!has_raw_barcode ||
                   metadata.barcode_abundance_entries.empty() ||
                   metadata.barcode_whitelist_fingerprint == 0)) {
    return Fail("barcoded ATAC mergeable spill v3 requires raw barcode "
                "evidence and a local whitelist-abundance table", error);
  }
  if (metadata.references.empty()) {
    return Fail("ATAC mergeable spill reference dictionary is empty", error);
  }
  if ((metadata.flags & kAtacMergeableAllocateMultiMappings) != 0 &&
      (metadata.multi_mapping_allocation_distance < 0 ||
       metadata.max_num_best_mappings == 0 ||
       metadata.max_num_best_mappings >
           static_cast<uint32_t>(std::numeric_limits<int>::max()))) {
    return Fail("invalid ATAC mergeable spill multimapping policy", error);
  }
  const bool has_summary_evidence =
      (metadata.flags & kAtacMergeableHasSummaryEvidence) != 0;
  if ((metadata.flags & kAtacMergeableReadRangeLateBound) != 0 &&
      metadata.first_global_read_ordinal != 0) {
    return Fail("late-bound ATAC mergeable spill has a precomputed prefix",
                error);
  }
  if (!has_summary_evidence ||
      metadata.summary_evidence_count != metadata.input_record_count ||
      metadata.cache_size == 0 ||
      metadata.cache_size >
          static_cast<uint32_t>(std::numeric_limits<int32_t>::max()) ||
      metadata.k_for_minhash == 0 ||
      metadata.frip_est_coefficients.size() != 5) {
    return Fail("incomplete ATAC mergeable spill summary policy", error);
  }
  for (double coefficient : metadata.frip_est_coefficients) {
    if (!std::isfinite(coefficient)) {
      return Fail("invalid ATAC mergeable spill summary coefficient", error);
    }
  }
  uint64_t metadata_bytes = metadata.sample_id.size() + metadata.input_id.size();
  for (const auto &reference : metadata.references) {
    if (reference.name.empty()) {
      return Fail("ATAC mergeable spill reference name is empty", error);
    }
    if (!AddMetadataBytes(sizeof(uint32_t) * 2u + reference.name.size(),
                          &metadata_bytes)) {
      return Fail("ATAC mergeable spill metadata is too large", error);
    }
  }
  uint64_t previous_barcode = 0;
  bool have_previous_barcode = false;
  uint64_t abundance_sum = 0;
  for (const auto &entry : metadata.barcode_abundance_entries) {
    if (have_previous_barcode && entry.barcode_key <= previous_barcode) {
      return Fail("ATAC mergeable spill barcode-abundance keys are not "
                  "strictly sorted", error);
    }
    if (abundance_sum > std::numeric_limits<uint64_t>::max() - entry.count) {
      return Fail("ATAC mergeable spill barcode abundance overflows", error);
    }
    abundance_sum += entry.count;
    previous_barcode = entry.barcode_key;
    have_previous_barcode = true;
    if (!AddMetadataBytes(sizeof(uint64_t) * 2u, &metadata_bytes)) {
      return Fail("ATAC mergeable spill metadata is too large", error);
    }
  }
  if (!is_bulk && abundance_sum != metadata.local_num_sample_barcodes) {
    return Fail("ATAC mergeable spill local barcode abundance total "
                "disagrees with its header", error);
  }
  if (metadata_bytes > kMaximumMetadataBytes ||
      metadata.sample_id.size() > std::numeric_limits<uint32_t>::max() ||
      metadata.input_id.size() > std::numeric_limits<uint32_t>::max() ||
      metadata.references.size() > std::numeric_limits<uint32_t>::max() ||
      metadata.barcode_abundance_entries.size() >
          std::numeric_limits<uint32_t>::max()) {
    return Fail("ATAC mergeable spill metadata exceeds v3 limits", error);
  }

  output_path_ = path;
  temporary_path_ = path + ".tmp." + std::to_string(getpid());
  metadata_ = metadata;
  file_ = fopen(temporary_path_.c_str(), "wb+");
  if (file_ == nullptr) {
    return Fail("cannot open ATAC mergeable spill temporary output " +
                    temporary_path_ + ": " + std::strerror(errno),
                error);
  }

  AtacMergeableSpillHeaderV3 header = {};
  memcpy(header.magic, kAtacMergeableSpillMagicV3, sizeof(header.magic));
  header.format_version = kAtacMergeableSpillFormatVersion;
  header.fixed_header_bytes = sizeof(header);
  header.metadata_bytes = static_cast<uint32_t>(metadata_bytes);
  header.record_codec_version = kAtacSpillRecordCodecVersion;
  header.schema_mask = metadata.schema_mask;
  header.flags = metadata.flags;
  header.shard_ordinal = metadata.shard_ordinal;
  header.shard_count = metadata.shard_count;
  header.first_global_read_ordinal = metadata.first_global_read_ordinal;
  header.input_record_count = metadata.input_record_count;
  header.record_count = 0;
  header.barcode_length = metadata.barcode_length;
  header.num_reference_sequences =
      static_cast<uint32_t>(metadata.references.size());
  header.mapq_threshold = metadata.mapq_threshold;
  header.tn5_forward_shift = metadata.tn5_forward_shift;
  header.tn5_reverse_shift = metadata.tn5_reverse_shift;
  header.barcode_correction_error_threshold =
      metadata.barcode_correction_error_threshold;
  header.sample_id_bytes = static_cast<uint32_t>(metadata.sample_id.size());
  header.input_id_bytes = static_cast<uint32_t>(metadata.input_id.size());
  header.num_barcode_abundance_entries = static_cast<uint32_t>(
      metadata.barcode_abundance_entries.size());
  header.local_num_sample_barcodes = metadata.local_num_sample_barcodes;
  header.barcode_whitelist_fingerprint =
      metadata.barcode_whitelist_fingerprint;
  memcpy(&header.barcode_correction_probability_bits,
         &metadata.barcode_correction_probability_threshold,
         sizeof(header.barcode_correction_probability_bits));
  header.multi_mapping_allocation_distance =
      metadata.multi_mapping_allocation_distance;
  header.multi_mapping_allocation_seed = metadata.multi_mapping_allocation_seed;
  header.max_num_best_mappings = metadata.max_num_best_mappings;
  header.summary_evidence_count = metadata.summary_evidence_count;
  header.cache_size = metadata.cache_size;
  header.k_for_minhash = metadata.k_for_minhash;
  for (size_t i = 0; i < metadata.frip_est_coefficients.size(); ++i) {
    memcpy(&header.frip_est_coefficient_bits[i],
           &metadata.frip_est_coefficients[i],
           sizeof(header.frip_est_coefficient_bits[i]));
  }

  if (!WriteBytes(file_, &header, sizeof(header)) ||
      !WriteBytes(file_, metadata.sample_id.data(), metadata.sample_id.size()) ||
      !WriteBytes(file_, metadata.input_id.data(), metadata.input_id.size())) {
    return Fail("cannot write ATAC mergeable spill header", error);
  }
  for (const auto &reference : metadata.references) {
    const uint32_t name_bytes = static_cast<uint32_t>(reference.name.size());
    if (!WriteBytes(file_, &name_bytes, sizeof(name_bytes)) ||
        !WriteBytes(file_, &reference.length, sizeof(reference.length)) ||
        !WriteBytes(file_, reference.name.data(), reference.name.size())) {
      return Fail("cannot write ATAC mergeable spill reference dictionary",
                  error);
    }
  }
  for (const auto &entry : metadata.barcode_abundance_entries) {
    if (!WriteBytes(file_, &entry.barcode_key, sizeof(entry.barcode_key)) ||
        !WriteBytes(file_, &entry.count, sizeof(entry.count))) {
      return Fail("cannot write ATAC mergeable spill barcode abundance", error);
    }
  }
  return true;
}

bool AtacMergeableSpillWriter::Append(uint32_t rid,
                                      const AtacSpillRecord &record,
                                      std::string *error) {
  if (file_ == nullptr || finalized_) {
    return Fail("ATAC mergeable spill writer is not open", error);
  }
  if (summary_evidence_count_ != metadata_.summary_evidence_count) {
    return Fail("ATAC mergeable spill mapping records precede complete "
                "summary evidence", error);
  }
  if (rid >= metadata_.references.size()) {
    return Fail("ATAC mergeable spill record has an invalid reference id",
                error);
  }
  if (!record.HasBamPairSection()) {
    return Fail("ATAC mergeable spill record is missing BAM-pair payloads",
                error);
  }
  if ((metadata_.schema_mask & kAtacSpillSchemaHasRawBarcodeEvidence) != 0 &&
      (!record.HasRawBarcodeEvidence() ||
       record.raw_barcode_qual_.size() != metadata_.barcode_length)) {
    return Fail("ATAC mergeable spill record is missing raw barcode evidence",
                error);
  }
  const uint32_t valid_n_mask =
      metadata_.barcode_length == 32
          ? std::numeric_limits<uint32_t>::max()
          : ((uint32_t{1} << metadata_.barcode_length) - 1u);
  if (record.HasRawBarcodeEvidence() &&
      (record.raw_barcode_n_mask_ & ~valid_n_mask) != 0) {
    return Fail("ATAC mergeable spill record has an invalid barcode N mask",
                error);
  }
  if (have_previous_ &&
      (rid < previous_rid_ ||
       (rid == previous_rid_ && record < previous_record_))) {
    return Fail("ATAC mergeable spill records are not globally sorted", error);
  }

  char *payload = nullptr;
  size_t payload_size = 0;
  FILE *memory = open_memstream(&payload, &payload_size);
  if (memory == nullptr) {
    return Fail("cannot allocate ATAC mergeable spill record buffer", error);
  }
  const size_t expected = record.SerializedSize();
  const size_t written = record.WriteToFile(memory);
  const bool close_ok = fclose(memory) == 0;
  if (!close_ok || written != expected || payload_size != expected ||
      payload_size > kMaximumRecordBytes ||
      payload_size > std::numeric_limits<uint32_t>::max()) {
    free(payload);
    return Fail("cannot serialize ATAC mergeable spill record", error);
  }
  const uint32_t payload_bytes = static_cast<uint32_t>(payload_size);
  const bool output_ok =
      WriteBytes(file_, &rid, sizeof(rid)) &&
      WriteBytes(file_, &payload_bytes, sizeof(payload_bytes)) &&
      WriteBytes(file_, payload, payload_size);
  free(payload);
  if (!output_ok) {
    return Fail("cannot write ATAC mergeable spill record", error);
  }
  previous_rid_ = rid;
  previous_record_ = record;
  have_previous_ = true;
  ++record_count_;
  return true;
}

bool AtacMergeableSpillWriter::AppendSummaryEvidence(
    const AtacSummaryEvidence &evidence, std::string *error) {
  if (file_ == nullptr || finalized_ || record_count_ != 0) {
    return Fail("ATAC mergeable spill summary evidence is out of order", error);
  }
  if (summary_evidence_count_ >= metadata_.summary_evidence_count) {
    return Fail("ATAC mergeable spill has excess summary evidence", error);
  }
  const bool is_bulk =
      (metadata_.schema_mask & kAtacSpillSchemaIsBulk) != 0;
  if ((is_bulk && (evidence.raw_barcode_key != 0 ||
                   evidence.raw_barcode_n_mask != 0 ||
                   !evidence.raw_barcode_qual.empty())) ||
      (!is_bulk &&
       evidence.raw_barcode_qual.size() != metadata_.barcode_length)) {
    return Fail("ATAC mergeable spill summary barcode evidence is invalid",
                error);
  }
  const uint32_t valid_n_mask =
      metadata_.barcode_length == 32
          ? std::numeric_limits<uint32_t>::max()
          : (metadata_.barcode_length == 0
                 ? uint32_t{0}
                 : ((uint32_t{1} << metadata_.barcode_length) - 1u));
  if ((evidence.raw_barcode_n_mask & ~valid_n_mask) != 0 ||
      evidence.cache_slot1 < -1 || evidence.cache_slot2 < -1 ||
      evidence.cache_slot1 >= static_cast<int32_t>(metadata_.cache_size) ||
      evidence.cache_slot2 >= static_cast<int32_t>(metadata_.cache_size)) {
    return Fail("ATAC mergeable spill summary evidence value is invalid",
                error);
  }
  if (!WriteBytes(file_, &evidence.raw_barcode_key,
                  sizeof(evidence.raw_barcode_key)) ||
      !WriteBytes(file_, &evidence.raw_barcode_n_mask,
                  sizeof(evidence.raw_barcode_n_mask)) ||
      !WriteBytes(file_, &evidence.cache_slot1,
                  sizeof(evidence.cache_slot1)) ||
      !WriteBytes(file_, &evidence.cache_slot2,
                  sizeof(evidence.cache_slot2)) ||
      !WriteBytes(file_, evidence.raw_barcode_qual.data(),
                  evidence.raw_barcode_qual.size())) {
    return Fail("cannot write ATAC mergeable spill summary evidence", error);
  }
  ++summary_evidence_count_;
  return true;
}

bool AtacMergeableSpillWriter::Finalize(std::string *error) {
  if (file_ == nullptr || finalized_) {
    return Fail("ATAC mergeable spill writer cannot be finalized", error);
  }
  if (summary_evidence_count_ != metadata_.summary_evidence_count) {
    return Fail("ATAC mergeable spill summary evidence is incomplete", error);
  }
  if (fflush(file_) != 0) {
    return Fail("cannot flush ATAC mergeable spill output", error);
  }
  const int descriptor = fileno(file_);
  const off_t count_offset =
      static_cast<off_t>(offsetof(AtacMergeableSpillHeaderV3, record_count));
  if (descriptor < 0 ||
      pwrite(descriptor, &record_count_, sizeof(record_count_), count_offset) !=
          static_cast<ssize_t>(sizeof(record_count_)) ||
      fsync(descriptor) != 0) {
    return Fail("cannot commit ATAC mergeable spill header", error);
  }
  if (fclose(file_) != 0) {
    file_ = nullptr;
    return Fail("cannot close ATAC mergeable spill output", error);
  }
  file_ = nullptr;
  if (rename(temporary_path_.c_str(), output_path_.c_str()) != 0) {
    return Fail("cannot publish ATAC mergeable spill output " + output_path_ +
                    ": " + std::strerror(errno),
                error);
  }
  finalized_ = true;
  temporary_path_.clear();
  return true;
}

AtacMergeableSpillReader::~AtacMergeableSpillReader() {
  if (file_ != nullptr) {
    fclose(file_);
  }
}

bool AtacMergeableSpillReader::Fail(const std::string &message,
                                    std::string *error) {
  if (error != nullptr) {
    *error = message;
  }
  return false;
}

bool AtacMergeableSpillReader::Open(const std::string &path,
                                    std::string *error) {
  if (file_ != nullptr) {
    return Fail("ATAC mergeable spill reader is already open", error);
  }
  file_ = fopen(path.c_str(), "rb");
  if (file_ == nullptr) {
    return Fail("cannot open ATAC mergeable spill " + path + ": " +
                    std::strerror(errno),
                error);
  }
  // Eight shard readers are interleaved by the gather heap. A multi-megabyte
  // stdio buffer prevents that access pattern from degenerating into millions
  // of tiny Lustre reads.
  (void)setvbuf(file_, nullptr, _IOFBF, 4u * 1024u * 1024u);
  path_ = path;
  AtacMergeableSpillHeaderV3 header = {};
  if (!ReadBytes(file_, &header, sizeof(header)) ||
      memcmp(header.magic, kAtacMergeableSpillMagicV3,
             sizeof(header.magic)) != 0 ||
      header.format_version != kAtacMergeableSpillFormatVersion ||
      header.fixed_header_bytes != sizeof(header) ||
      header.record_codec_version != kAtacSpillRecordCodecVersion) {
    return Fail("invalid or unsupported ATAC mergeable spill header in " + path,
                error);
  }
  if (header.metadata_bytes > kMaximumMetadataBytes ||
      header.shard_count == 0 || header.shard_ordinal >= header.shard_count ||
      header.sample_id_bytes == 0 || header.input_id_bytes == 0 ||
      header.num_reference_sequences == 0 ||
      (header.schema_mask & kAtacSpillSchemaHasBamPair) == 0) {
    return Fail("invalid ATAC mergeable spill metadata fields in " + path,
                error);
  }
  if ((header.flags & kAtacMergeableAllocateMultiMappings) != 0 &&
      (header.multi_mapping_allocation_distance < 0 ||
       header.max_num_best_mappings == 0 ||
       header.max_num_best_mappings >
           static_cast<uint32_t>(std::numeric_limits<int>::max()))) {
    return Fail("invalid ATAC mergeable spill multimapping metadata in " +
                    path,
                error);
  }
  if ((header.flags & kAtacMergeableHasSummaryEvidence) == 0 ||
      header.summary_evidence_count != header.input_record_count ||
      header.cache_size == 0 ||
      header.cache_size >
          static_cast<uint32_t>(std::numeric_limits<int32_t>::max()) ||
      header.k_for_minhash == 0) {
    return Fail("incomplete ATAC mergeable spill summary metadata in " + path,
                error);
  }
  if ((header.flags & kAtacMergeableReadRangeLateBound) != 0 &&
      header.first_global_read_ordinal != 0) {
    return Fail("late-bound ATAC mergeable spill has a precomputed prefix in " +
                    path,
                error);
  }

  uint64_t consumed = 0;
  auto read_string = [&](uint32_t length, std::string *value) {
    if (!AddMetadataBytes(length, &consumed)) {
      return false;
    }
    value->assign(length, '\0');
    return ReadBytes(file_, length == 0 ? nullptr : &(*value)[0], length);
  };
  if (!read_string(header.sample_id_bytes, &metadata_.sample_id) ||
      !read_string(header.input_id_bytes, &metadata_.input_id)) {
    return Fail("truncated ATAC mergeable spill identity metadata in " + path,
                error);
  }
  metadata_.references.reserve(header.num_reference_sequences);
  for (uint32_t i = 0; i < header.num_reference_sequences; ++i) {
    uint32_t name_bytes = 0;
    uint32_t length = 0;
    if (!AddMetadataBytes(sizeof(name_bytes) + sizeof(length), &consumed) ||
        !ReadBytes(file_, &name_bytes, sizeof(name_bytes)) ||
        !ReadBytes(file_, &length, sizeof(length)) || name_bytes == 0) {
      return Fail("invalid ATAC mergeable spill reference metadata in " + path,
                  error);
    }
    AtacMergeableSpillReference reference;
    reference.length = length;
    if (!read_string(name_bytes, &reference.name)) {
      return Fail("truncated ATAC mergeable spill reference name in " + path,
                  error);
    }
    metadata_.references.push_back(std::move(reference));
  }
  metadata_.barcode_abundance_entries.reserve(
      header.num_barcode_abundance_entries);
  uint64_t abundance_sum = 0;
  uint64_t previous_barcode = 0;
  for (uint32_t i = 0; i < header.num_barcode_abundance_entries; ++i) {
    if (!AddMetadataBytes(sizeof(uint64_t) * 2u, &consumed)) {
      return Fail("ATAC mergeable spill barcode metadata is too large in " +
                      path,
                  error);
    }
    AtacBarcodeAbundanceEntry entry;
    if (!ReadBytes(file_, &entry.barcode_key, sizeof(entry.barcode_key)) ||
        !ReadBytes(file_, &entry.count, sizeof(entry.count)) ||
        (i > 0 && entry.barcode_key <= previous_barcode) ||
        abundance_sum > std::numeric_limits<uint64_t>::max() - entry.count) {
      return Fail("invalid ATAC mergeable spill barcode abundance in " + path,
                  error);
    }
    abundance_sum += entry.count;
    previous_barcode = entry.barcode_key;
    metadata_.barcode_abundance_entries.push_back(entry);
  }
  if (consumed != header.metadata_bytes) {
    return Fail("ATAC mergeable spill metadata length mismatch in " + path,
                error);
  }

  metadata_.schema_mask = header.schema_mask;
  metadata_.flags = header.flags;
  metadata_.shard_ordinal = header.shard_ordinal;
  metadata_.shard_count = header.shard_count;
  metadata_.first_global_read_ordinal = header.first_global_read_ordinal;
  metadata_.input_record_count = header.input_record_count;
  metadata_.barcode_length = header.barcode_length;
  metadata_.mapq_threshold = header.mapq_threshold;
  metadata_.tn5_forward_shift = header.tn5_forward_shift;
  metadata_.tn5_reverse_shift = header.tn5_reverse_shift;
  metadata_.local_num_sample_barcodes = header.local_num_sample_barcodes;
  metadata_.barcode_whitelist_fingerprint =
      header.barcode_whitelist_fingerprint;
  metadata_.barcode_correction_error_threshold =
      header.barcode_correction_error_threshold;
  memcpy(&metadata_.barcode_correction_probability_threshold,
         &header.barcode_correction_probability_bits,
         sizeof(header.barcode_correction_probability_bits));
  metadata_.multi_mapping_allocation_distance =
      header.multi_mapping_allocation_distance;
  metadata_.multi_mapping_allocation_seed =
      header.multi_mapping_allocation_seed;
  metadata_.max_num_best_mappings = header.max_num_best_mappings;
  metadata_.summary_evidence_count = header.summary_evidence_count;
  metadata_.cache_size = header.cache_size;
  metadata_.k_for_minhash = header.k_for_minhash;
  metadata_.frip_est_coefficients.resize(5);
  for (size_t i = 0; i < metadata_.frip_est_coefficients.size(); ++i) {
    memcpy(&metadata_.frip_est_coefficients[i],
           &header.frip_est_coefficient_bits[i],
           sizeof(header.frip_est_coefficient_bits[i]));
    if (!std::isfinite(metadata_.frip_est_coefficients[i])) {
      return Fail("invalid ATAC mergeable spill summary coefficient in " +
                      path,
                  error);
    }
  }
  const bool is_bulk =
      (metadata_.schema_mask & kAtacSpillSchemaIsBulk) != 0;
  if (!is_bulk &&
      ((metadata_.schema_mask & kAtacSpillSchemaHasRawBarcodeEvidence) == 0 ||
       metadata_.barcode_length == 0 || metadata_.barcode_length > 32 ||
       metadata_.barcode_correction_error_threshold > 2 ||
       !std::isfinite(metadata_.barcode_correction_probability_threshold) ||
       metadata_.barcode_correction_probability_threshold < 0.0 ||
       metadata_.barcode_correction_probability_threshold >= 1.0 ||
       metadata_.barcode_abundance_entries.empty() ||
       metadata_.barcode_whitelist_fingerprint == 0 ||
       abundance_sum != metadata_.local_num_sample_barcodes)) {
    return Fail("incomplete ATAC mergeable spill barcode correction metadata "
                "in " + path, error);
  }
  expected_record_count_ = header.record_count;
  expected_summary_evidence_count_ = header.summary_evidence_count;
  return true;
}

bool AtacMergeableSpillReader::ReadNextSummaryEvidence(
    AtacSummaryEvidence *evidence, bool *eof, std::string *error) {
  if (file_ == nullptr || evidence == nullptr || eof == nullptr ||
      records_read_ != 0) {
    return Fail("ATAC mergeable spill summary reader is not initialized",
                error);
  }
  *eof = false;
  if (summary_evidence_read_ == expected_summary_evidence_count_) {
    *eof = true;
    return true;
  }
  const bool is_bulk =
      (metadata_.schema_mask & kAtacSpillSchemaIsBulk) != 0;
  evidence->raw_barcode_qual.assign(is_bulk ? 0 : metadata_.barcode_length,
                                    '\0');
  if (!ReadBytes(file_, &evidence->raw_barcode_key,
                 sizeof(evidence->raw_barcode_key)) ||
      !ReadBytes(file_, &evidence->raw_barcode_n_mask,
                 sizeof(evidence->raw_barcode_n_mask)) ||
      !ReadBytes(file_, &evidence->cache_slot1,
                 sizeof(evidence->cache_slot1)) ||
      !ReadBytes(file_, &evidence->cache_slot2,
                 sizeof(evidence->cache_slot2)) ||
      !ReadBytes(file_, evidence->raw_barcode_qual.empty()
                            ? nullptr
                            : &evidence->raw_barcode_qual[0],
                 evidence->raw_barcode_qual.size())) {
    return Fail("truncated ATAC mergeable spill summary evidence in " + path_,
                error);
  }
  const uint32_t valid_n_mask =
      metadata_.barcode_length == 32
          ? std::numeric_limits<uint32_t>::max()
          : (metadata_.barcode_length == 0
                 ? uint32_t{0}
                 : ((uint32_t{1} << metadata_.barcode_length) - 1u));
  if ((is_bulk && (evidence->raw_barcode_key != 0 ||
                   evidence->raw_barcode_n_mask != 0)) ||
      (evidence->raw_barcode_n_mask & ~valid_n_mask) != 0 ||
      evidence->cache_slot1 < -1 || evidence->cache_slot2 < -1 ||
      evidence->cache_slot1 >= static_cast<int32_t>(metadata_.cache_size) ||
      evidence->cache_slot2 >= static_cast<int32_t>(metadata_.cache_size)) {
    return Fail("invalid ATAC mergeable spill summary evidence in " + path_,
                error);
  }
  ++summary_evidence_read_;
  return true;
}

bool AtacMergeableSpillReader::SkipRemainingSummaryEvidence(
    std::string *error) {
  if (file_ == nullptr || records_read_ != 0 || full_decoder_started_ ||
      bed_decoder_started_) {
    return Fail("ATAC mergeable spill summary skipper is not initialized",
                error);
  }
  if (summary_evidence_read_ > expected_summary_evidence_count_) {
    return Fail("ATAC mergeable spill summary count is invalid in " + path_,
                error);
  }
  const uint64_t remaining =
      expected_summary_evidence_count_ - summary_evidence_read_;
  const bool is_bulk =
      (metadata_.schema_mask & kAtacSpillSchemaIsBulk) != 0;
  const uint64_t record_bytes =
      sizeof(uint64_t) + sizeof(uint32_t) + 2u * sizeof(int32_t) +
      (is_bulk ? uint64_t{0} : metadata_.barcode_length);
  if (record_bytes == 0 ||
      remaining > std::numeric_limits<uint64_t>::max() / record_bytes) {
    return Fail("ATAC mergeable spill summary skip overflows in " + path_,
                error);
  }
  const uint64_t bytes = remaining * record_bytes;
  if (bytes > static_cast<uint64_t>(std::numeric_limits<off_t>::max()) ||
      fseeko(file_, static_cast<off_t>(bytes), SEEK_CUR) != 0) {
    return Fail("cannot skip ATAC mergeable spill summary evidence in " + path_ +
                    ": " + std::strerror(errno),
                error);
  }
  summary_evidence_read_ = expected_summary_evidence_count_;
  return true;
}

bool AtacMergeableSpillReader::EnsureBedStreamBytes(size_t count,
                                                     std::string *error) {
  if (count > static_cast<size_t>(kMaximumRecordBytes)) {
    return Fail("ATAC mergeable spill BED decoder request is too large in " +
                    path_,
                error);
  }
  size_t available = bed_stream_end_ - bed_stream_begin_;
  if (available >= count) {
    return true;
  }
  if (bed_stream_begin_ != 0 && available != 0) {
    memmove(bed_stream_buffer_.data(),
            bed_stream_buffer_.data() + bed_stream_begin_, available);
  }
  bed_stream_begin_ = 0;
  bed_stream_end_ = available;
  const size_t required_capacity = std::max(kBedStreamBufferBytes, count);
  if (bed_stream_buffer_.size() < required_capacity) {
    bed_stream_buffer_.resize(required_capacity);
  }
  while (bed_stream_end_ < count) {
    const size_t capacity = bed_stream_buffer_.size() - bed_stream_end_;
    const size_t got = fread(bed_stream_buffer_.data() + bed_stream_end_, 1,
                             capacity, file_);
    bed_stream_end_ += got;
    if (got == 0) {
      if (ferror(file_) != 0) {
        return Fail("cannot read ATAC mergeable spill BED stream in " + path_ +
                        ": " + std::strerror(errno),
                    error);
      }
      break;
    }
  }
  return bed_stream_end_ - bed_stream_begin_ >= count;
}

bool AtacMergeableSpillReader::ReadNext(uint32_t *rid,
                                        AtacSpillRecord *record, bool *eof,
                                        std::string *error) {
  if (file_ == nullptr || rid == nullptr || record == nullptr || eof == nullptr) {
    return Fail("ATAC mergeable spill reader is not initialized", error);
  }
  if (bed_decoder_started_) {
    return Fail("ATAC mergeable spill decoder mode changed after reading " +
                    path_,
                error);
  }
  full_decoder_started_ = true;
  if (summary_evidence_read_ != expected_summary_evidence_count_) {
    return Fail("ATAC mergeable spill mapping read precedes summary evidence",
                error);
  }
  *eof = false;
  if (records_read_ == expected_record_count_) {
    if (!verified_eof_) {
      const int trailing = fgetc(file_);
      if (trailing != EOF || ferror(file_) != 0) {
        return Fail("ATAC mergeable spill has trailing or unreadable bytes in " +
                        path_,
                    error);
      }
      verified_eof_ = true;
    }
    *eof = true;
    return true;
  }

  uint32_t payload_bytes = 0;
  if (!ReadBytes(file_, rid, sizeof(*rid)) ||
      !ReadBytes(file_, &payload_bytes, sizeof(payload_bytes)) ||
      payload_bytes == 0 || payload_bytes > kMaximumRecordBytes ||
      *rid >= metadata_.references.size()) {
    return Fail("invalid or truncated ATAC mergeable spill record header in " +
                    path_,
                error);
  }
  std::vector<char> payload(payload_bytes);
  if (!ReadBytes(file_, payload.data(), payload.size())) {
    return Fail("truncated ATAC mergeable spill record payload in " + path_,
                error);
  }
  FILE *memory = fmemopen(payload.data(), payload.size(), "rb");
  if (memory == nullptr) {
    return Fail("cannot open ATAC mergeable spill record buffer", error);
  }
  record->LoadFromFile(memory, metadata_.schema_mask);
  const long consumed = ftell(memory);
  fclose(memory);
  if (consumed < 0 || static_cast<size_t>(consumed) != payload.size() ||
      record->SerializedSize() != payload.size() ||
      !record->HasBamPairSection()) {
    return Fail("ATAC mergeable spill record codec length mismatch in " + path_,
                error);
  }
  if ((metadata_.schema_mask & kAtacSpillSchemaHasRawBarcodeEvidence) != 0 &&
      (!record->HasRawBarcodeEvidence() ||
       record->raw_barcode_qual_.size() != metadata_.barcode_length)) {
    return Fail("ATAC mergeable spill raw barcode evidence mismatch in " +
                    path_,
                error);
  }
  const uint32_t valid_n_mask =
      metadata_.barcode_length == 32
          ? std::numeric_limits<uint32_t>::max()
          : ((uint32_t{1} << metadata_.barcode_length) - 1u);
  if (record->HasRawBarcodeEvidence() &&
      (record->raw_barcode_n_mask_ & ~valid_n_mask) != 0) {
    return Fail("ATAC mergeable spill barcode N mask is invalid in " + path_,
                error);
  }
  const PairedEndMappingWithBarcode &mapping_key = *record;
  if (have_previous_ &&
      (*rid < previous_rid_ ||
       (*rid == previous_rid_ && mapping_key < previous_mapping_))) {
    return Fail("ATAC mergeable spill records are not sorted in " + path_,
                error);
  }
  previous_rid_ = *rid;
  previous_mapping_ = mapping_key;
  have_previous_ = true;
  ++records_read_;
  return true;
}

bool AtacMergeableSpillReader::ReadNextBed(uint32_t *rid,
                                           AtacSpillRecord *record, bool *eof,
                                           std::string *error) {
  if (file_ == nullptr || rid == nullptr || record == nullptr || eof == nullptr) {
    return Fail("ATAC mergeable spill BED reader is not initialized", error);
  }
  if (full_decoder_started_) {
    return Fail("ATAC mergeable spill decoder mode changed after reading " +
                    path_,
                error);
  }
  bed_decoder_started_ = true;
  if (summary_evidence_read_ != expected_summary_evidence_count_) {
    return Fail("ATAC mergeable spill mapping read precedes summary evidence",
                error);
  }
  *eof = false;
  if (records_read_ == expected_record_count_) {
    if (!verified_eof_) {
      if (bed_stream_end_ != bed_stream_begin_) {
        return Fail("ATAC mergeable spill has trailing bytes in " + path_,
                    error);
      }
      const int trailing = fgetc(file_);
      if (trailing != EOF || ferror(file_) != 0) {
        return Fail("ATAC mergeable spill has trailing or unreadable bytes in " +
                        path_,
                    error);
      }
      verified_eof_ = true;
    }
    *eof = true;
    return true;
  }

  constexpr size_t kRecordEnvelopeBytes = 2u * sizeof(uint32_t);
  if (!EnsureBedStreamBytes(kRecordEnvelopeBytes, error)) {
    return Fail("truncated ATAC mergeable spill record header in " + path_,
                error);
  }
  uint32_t payload_bytes = 0;
  memcpy(rid, bed_stream_buffer_.data() + bed_stream_begin_, sizeof(*rid));
  memcpy(&payload_bytes,
         bed_stream_buffer_.data() + bed_stream_begin_ + sizeof(*rid),
         sizeof(payload_bytes));
  bed_stream_begin_ += kRecordEnvelopeBytes;
  if (payload_bytes == 0 || payload_bytes > kMaximumRecordBytes ||
      *rid >= metadata_.references.size()) {
    return Fail("invalid ATAC mergeable spill record header in " + path_,
                error);
  }
  if (!EnsureBedStreamBytes(payload_bytes, error)) {
    return Fail("truncated ATAC mergeable spill record payload in " + path_,
                error);
  }
  const char *payload = bed_stream_buffer_.data() + bed_stream_begin_;
  if (!DecodeBedPayload(payload, payload_bytes, metadata_.schema_mask,
                        metadata_.barcode_length, record)) {
    return Fail("invalid ATAC mergeable spill BED payload in " + path_, error);
  }
  bed_stream_begin_ += payload_bytes;

  const uint32_t valid_n_mask =
      metadata_.barcode_length == 32
          ? std::numeric_limits<uint32_t>::max()
          : ((uint32_t{1} << metadata_.barcode_length) - 1u);
  if (record->HasRawBarcodeEvidence() &&
      (record->raw_barcode_n_mask_ & ~valid_n_mask) != 0) {
    return Fail("ATAC mergeable spill barcode N mask is invalid in " + path_,
                error);
  }
  const PairedEndMappingWithBarcode &mapping_key = *record;
  if (have_previous_ &&
      (*rid < previous_rid_ ||
       (*rid == previous_rid_ && mapping_key < previous_mapping_))) {
    return Fail("ATAC mergeable spill records are not sorted in " + path_,
                error);
  }
  previous_rid_ = *rid;
  previous_mapping_ = mapping_key;
  have_previous_ = true;
  ++records_read_;
  return true;
}

bool GlobalizeAtacSpillReadId(AtacSpillRecord *record,
                              uint64_t first_global_read_ordinal,
                              uint64_t input_record_count,
                              std::string *error) {
  if (record == nullptr) {
    if (error != nullptr) {
      *error = "null ATAC spill record";
    }
    return false;
  }
  const uint64_t local_id = record->read_id_;
  if (local_id >= input_record_count ||
      local_id > std::numeric_limits<uint64_t>::max() -
                     first_global_read_ordinal) {
    if (error != nullptr) {
      *error = "ATAC spill read id lies outside its declared shard range";
    }
    return false;
  }
  const uint64_t global_id = first_global_read_ordinal + local_id;
  record->read_id_ = global_id;
  if (record->HasBamPairSection()) {
    if (record->sam1.read_id_ != local_id || record->sam2.read_id_ != local_id) {
      if (error != nullptr) {
        *error = "ATAC spill fragment/BAM read ids disagree";
      }
      return false;
    }
    record->sam1.read_id_ = global_id;
    record->sam2.read_id_ = global_id;
  }
  return true;
}

}  // namespace chromap
