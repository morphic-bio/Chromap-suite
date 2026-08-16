#include "atac_hot_spill.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <sys/stat.h>
#include <unistd.h>

namespace chromap {
namespace {

constexpr size_t kHotReadBufferBytes = 8u * 1024u * 1024u;

bool WriteBytes(FILE *file, const void *data, size_t size) {
  return size == 0 || fwrite(data, 1, size, file) == size;
}

bool PreadAll(int fd, void *data, size_t size, uint64_t offset) {
  char *cursor = static_cast<char *>(data);
  size_t remaining = size;
  while (remaining != 0) {
    const ssize_t got = pread(fd, cursor, remaining,
                              static_cast<off_t>(offset));
    if (got < 0 && errno == EINTR) {
      continue;
    }
    if (got <= 0) {
      return false;
    }
    cursor += got;
    remaining -= static_cast<size_t>(got);
    offset += static_cast<uint64_t>(got);
  }
  return true;
}

uint16_t HotFlagsForMetadata(const AtacMergeableSpillMetadata &metadata) {
  uint16_t flags = 0;
  if ((metadata.schema_mask & kAtacSpillSchemaIsBulk) != 0) {
    flags = static_cast<uint16_t>(flags | kAtacHotSpillIsBulk);
  }
  if ((metadata.schema_mask & kAtacSpillSchemaHasRawBarcodeEvidence) != 0) {
    flags = static_cast<uint16_t>(flags |
                                  kAtacHotSpillHasRawBarcodeEvidence);
  }
  if ((metadata.schema_mask & kAtacSpillSchemaHasYHit) != 0) {
    flags = static_cast<uint16_t>(flags | kAtacHotSpillHasYHit);
  }
  return flags;
}

}  // namespace

std::string AtacHotSpillSidecarPath(const std::string &spill_path) {
  return spill_path + ".hot";
}

AtacHotSpillWriter::~AtacHotSpillWriter() {
  if (file_ != nullptr) {
    fclose(file_);
    file_ = nullptr;
  }
  if (!temporary_path_.empty() && !finalized_) {
    unlink(temporary_path_.c_str());
  }
}

bool AtacHotSpillWriter::Fail(const std::string &message,
                              std::string *error) {
  if (error != nullptr) {
    *error = message;
  }
  return false;
}

bool AtacHotSpillWriter::Open(
    const std::string &path, const AtacMergeableSpillMetadata &metadata,
    std::string *error) {
  if (file_ != nullptr || finalized_) {
    return Fail("ATAC hot spill writer is already open", error);
  }
  const bool is_bulk =
      (metadata.schema_mask & kAtacSpillSchemaIsBulk) != 0;
  const bool has_raw_barcode =
      (metadata.schema_mask & kAtacSpillSchemaHasRawBarcodeEvidence) != 0;
  if (path.empty() || metadata.shard_count == 0 ||
      metadata.shard_ordinal >= metadata.shard_count ||
      metadata.references.empty() || metadata.references.size() >
          static_cast<uint64_t>(std::numeric_limits<uint32_t>::max()) ||
      (!is_bulk && (!has_raw_barcode || metadata.barcode_length == 0 ||
                    metadata.barcode_length > 32)) ||
      (is_bulk && (has_raw_barcode || metadata.barcode_length != 0))) {
    return Fail("invalid ATAC hot spill metadata", error);
  }

  const uint64_t record_bytes = sizeof(AtacHotSpillRecordPrefixV1) +
                                (is_bulk ? 0u : metadata.barcode_length);
  if (record_bytes > std::numeric_limits<uint32_t>::max() ||
      metadata.references.size() >
          (std::numeric_limits<uint64_t>::max() - sizeof(header_)) /
              sizeof(AtacHotSpillDirectoryV1)) {
    return Fail("ATAC hot spill layout overflows", error);
  }

  metadata_ = metadata;
  output_path_ = path;
  temporary_path_ = path + ".tmp." + std::to_string(getpid());
  file_ = fopen(temporary_path_.c_str(), "wb+");
  if (file_ == nullptr) {
    return Fail("cannot open ATAC hot spill temporary output " +
                    temporary_path_ + ": " + std::strerror(errno),
                error);
  }
  (void)setvbuf(file_, nullptr, _IOFBF, 8u * 1024u * 1024u);

  memcpy(header_.magic, kAtacHotSpillMagicV1, sizeof(header_.magic));
  header_.format_version = kAtacHotSpillFormatVersion;
  header_.fixed_header_bytes = sizeof(header_);
  header_.record_prefix_bytes = sizeof(AtacHotSpillRecordPrefixV1);
  header_.flags = HotFlagsForMetadata(metadata);
  header_.endian_marker = kAtacHotSpillEndianMarker;
  header_.shard_ordinal = metadata.shard_ordinal;
  header_.shard_count = metadata.shard_count;
  header_.barcode_length = metadata.barcode_length;
  header_.num_reference_sequences =
      static_cast<uint32_t>(metadata.references.size());
  header_.first_global_read_ordinal = metadata.first_global_read_ordinal;
  header_.input_record_count = metadata.input_record_count;
  header_.barcode_whitelist_fingerprint =
      metadata.barcode_whitelist_fingerprint;
  header_.directory_offset = sizeof(header_);
  header_.data_offset = sizeof(header_) +
                        metadata.references.size() *
                            sizeof(AtacHotSpillDirectoryV1);
  header_.record_bytes = static_cast<uint32_t>(record_bytes);
  directory_.resize(metadata.references.size());
  encoded_record_.resize(header_.record_bytes);
  if (!WriteBytes(file_, &header_, sizeof(header_)) ||
      !WriteBytes(file_, directory_.data(),
                  directory_.size() * sizeof(directory_[0]))) {
    return Fail("cannot initialize ATAC hot spill output", error);
  }
  return true;
}

bool AtacHotSpillWriter::Append(uint32_t rid,
                                const AtacSpillRecord &record,
                                std::string *error) {
  if (file_ == nullptr || finalized_ || rid >= metadata_.references.size() ||
      record.num_dups_ != 1 ||
      record.read_id_ >= metadata_.input_record_count ||
      record.fragment_length_ == 0 ||
      static_cast<uint64_t>(record.fragment_start_position_) +
              record.fragment_length_ >
          metadata_.references[rid].length) {
    return Fail("invalid ATAC hot spill record", error);
  }
  const bool is_bulk =
      (metadata_.schema_mask & kAtacSpillSchemaIsBulk) != 0;
  if ((!is_bulk &&
       (!record.HasRawBarcodeEvidence() ||
        record.raw_barcode_qual_.size() != metadata_.barcode_length)) ||
      (is_bulk && (record.cell_barcode_ != 0 ||
                   record.raw_barcode_n_mask_ != 0 ||
                   !record.raw_barcode_qual_.empty()))) {
    return Fail("ATAC hot spill barcode evidence is invalid", error);
  }
  const uint32_t valid_n_mask =
      metadata_.barcode_length == 32
          ? std::numeric_limits<uint32_t>::max()
          : (metadata_.barcode_length == 0
                 ? uint32_t{0}
                 : (uint32_t{1} << metadata_.barcode_length) - 1u);
  const uint64_t valid_barcode_max =
      metadata_.barcode_length == 32
          ? std::numeric_limits<uint64_t>::max()
          : (metadata_.barcode_length == 0
                 ? uint64_t{0}
                 : (uint64_t{1} << (2u * metadata_.barcode_length)) - 1u);
  if ((record.raw_barcode_n_mask_ & ~valid_n_mask) != 0 ||
      record.cell_barcode_ > valid_barcode_max || record.mapq_ > 63) {
    return Fail("ATAC hot spill barcode or MAPQ value is invalid", error);
  }
  if (have_previous_ &&
      (rid < previous_rid_ ||
       (rid == previous_rid_ &&
        static_cast<const PairedEndMappingWithBarcode &>(record) <
            previous_mapping_))) {
    return Fail("ATAC hot spill records are not globally sorted", error);
  }

  AtacHotSpillRecordPrefixV1 prefix = {};
  prefix.local_read_id = record.read_id_;
  prefix.raw_barcode_key = record.cell_barcode_;
  prefix.start = record.fragment_start_position_;
  prefix.raw_barcode_n_mask = record.raw_barcode_n_mask_;
  prefix.fragment_length = record.fragment_length_;
  prefix.positive_alignment_length = record.positive_alignment_length_;
  prefix.negative_alignment_length = record.negative_alignment_length_;
  prefix.mapq_direction_unique = static_cast<uint8_t>(
      record.mapq_ | (static_cast<uint8_t>(record.direction_) << 6) |
      (static_cast<uint8_t>(record.is_unique_) << 7));
  prefix.row_flags =
      record.HasYHit() ? kAtacHotSpillRowHasYHit : uint8_t{0};
  memcpy(encoded_record_.data(), &prefix, sizeof(prefix));
  if (!is_bulk) {
    memcpy(encoded_record_.data() + sizeof(prefix),
           record.raw_barcode_qual_.data(), record.raw_barcode_qual_.size());
  }
  if (!WriteBytes(file_, encoded_record_.data(), encoded_record_.size())) {
    return Fail("cannot write ATAC hot spill record", error);
  }

  AtacHotSpillDirectoryV1 &entry = directory_[rid];
  if (entry.record_count == 0) {
    entry.first_start = record.fragment_start_position_;
  }
  entry.last_start = record.fragment_start_position_;
  ++entry.record_count;
  ++record_count_;
  previous_rid_ = rid;
  previous_mapping_ =
      static_cast<const PairedEndMappingWithBarcode &>(record);
  have_previous_ = true;
  return true;
}

bool AtacHotSpillWriter::Finalize(std::string *error) {
  if (file_ == nullptr || finalized_) {
    return Fail("ATAC hot spill writer cannot be finalized", error);
  }
  uint64_t offset = header_.data_offset;
  uint64_t total = 0;
  for (auto &entry : directory_) {
    entry.offset = offset;
    if (entry.record_count >
            (std::numeric_limits<uint64_t>::max() - offset) /
                header_.record_bytes ||
        total > std::numeric_limits<uint64_t>::max() - entry.record_count) {
      return Fail("ATAC hot spill directory overflows", error);
    }
    offset += entry.record_count * header_.record_bytes;
    total += entry.record_count;
  }
  if (total != record_count_ || fflush(file_) != 0) {
    return Fail("cannot flush ATAC hot spill output", error);
  }
  header_.record_count = record_count_;
  const int descriptor = fileno(file_);
  if (descriptor < 0 ||
      pwrite(descriptor, &header_, sizeof(header_), 0) !=
          static_cast<ssize_t>(sizeof(header_)) ||
      pwrite(descriptor, directory_.data(),
             directory_.size() * sizeof(directory_[0]),
             static_cast<off_t>(header_.directory_offset)) !=
          static_cast<ssize_t>(directory_.size() * sizeof(directory_[0])) ||
      fsync(descriptor) != 0) {
    return Fail("cannot commit ATAC hot spill directory", error);
  }
  if (fclose(file_) != 0) {
    file_ = nullptr;
    return Fail("cannot close ATAC hot spill output", error);
  }
  file_ = nullptr;
  if (rename(temporary_path_.c_str(), output_path_.c_str()) != 0) {
    return Fail("cannot publish ATAC hot spill output " + output_path_ +
                    ": " + std::strerror(errno),
                error);
  }
  finalized_ = true;
  temporary_path_.clear();
  return true;
}

AtacHotSpillReader::~AtacHotSpillReader() {
  if (fd_ >= 0) {
    close(fd_);
  }
}

bool AtacHotSpillReader::Fail(const std::string &message,
                              std::string *error) {
  if (error != nullptr) {
    *error = message;
  }
  return false;
}

bool AtacHotSpillReader::Open(
    const std::string &parent_spill_path,
    const AtacMergeableSpillMetadata &parent_metadata,
    uint64_t parent_record_count, std::string *error) {
  if (fd_ >= 0) {
    return Fail("ATAC hot spill reader is already open", error);
  }
  path_ = AtacHotSpillSidecarPath(parent_spill_path);
  fd_ = open(path_.c_str(), O_RDONLY);
  if (fd_ < 0) {
    return Fail("cannot open ATAC hot spill " + path_ + ": " +
                    std::strerror(errno),
                error);
  }
  if (!PreadAll(fd_, &header_, sizeof(header_), 0) ||
      memcmp(header_.magic, kAtacHotSpillMagicV1,
             sizeof(header_.magic)) != 0 ||
      header_.format_version != kAtacHotSpillFormatVersion ||
      header_.fixed_header_bytes != sizeof(header_) ||
      header_.record_prefix_bytes != sizeof(AtacHotSpillRecordPrefixV1) ||
      header_.endian_marker != kAtacHotSpillEndianMarker ||
      header_.reserved != 0 ||
      (header_.flags & ~(kAtacHotSpillIsBulk |
                         kAtacHotSpillHasRawBarcodeEvidence |
                         kAtacHotSpillHasYHit)) != 0) {
    return Fail("invalid or unsupported ATAC hot spill header in " + path_,
                error);
  }
  const bool is_bulk =
      (parent_metadata.schema_mask & kAtacSpillSchemaIsBulk) != 0;
  const uint32_t expected_record_bytes =
      sizeof(AtacHotSpillRecordPrefixV1) +
      (is_bulk ? 0u : parent_metadata.barcode_length);
  if (header_.flags != HotFlagsForMetadata(parent_metadata) ||
      header_.shard_ordinal != parent_metadata.shard_ordinal ||
      header_.shard_count != parent_metadata.shard_count ||
      header_.barcode_length != parent_metadata.barcode_length ||
      header_.num_reference_sequences != parent_metadata.references.size() ||
      header_.first_global_read_ordinal !=
          parent_metadata.first_global_read_ordinal ||
      header_.input_record_count != parent_metadata.input_record_count ||
      header_.barcode_whitelist_fingerprint !=
          parent_metadata.barcode_whitelist_fingerprint ||
      header_.record_count != parent_record_count ||
      header_.record_bytes != expected_record_bytes ||
      header_.directory_offset != sizeof(header_) ||
      header_.data_offset !=
          sizeof(header_) + parent_metadata.references.size() *
                                sizeof(AtacHotSpillDirectoryV1)) {
    return Fail("ATAC hot spill disagrees with its parent shard " + path_,
                error);
  }

  directory_.resize(parent_metadata.references.size());
  if (!PreadAll(fd_, directory_.data(),
                directory_.size() * sizeof(directory_[0]),
                header_.directory_offset)) {
    return Fail("truncated ATAC hot spill directory in " + path_, error);
  }
  reference_lengths_.reserve(parent_metadata.references.size());
  uint64_t expected_offset = header_.data_offset;
  uint64_t total = 0;
  for (size_t rid = 0; rid < directory_.size(); ++rid) {
    const auto &entry = directory_[rid];
    reference_lengths_.push_back(parent_metadata.references[rid].length);
    if (entry.offset != expected_offset ||
        (entry.record_count == 0 &&
         (entry.first_start != 0 || entry.last_start != 0)) ||
        (entry.record_count != 0 &&
         (entry.first_start > entry.last_start ||
          entry.last_start >= parent_metadata.references[rid].length)) ||
        entry.record_count >
            (std::numeric_limits<uint64_t>::max() - expected_offset) /
                header_.record_bytes ||
        total > std::numeric_limits<uint64_t>::max() - entry.record_count) {
      return Fail("invalid ATAC hot spill partition directory in " + path_,
                  error);
    }
    expected_offset += entry.record_count * header_.record_bytes;
    total += entry.record_count;
  }
  struct stat status = {};
  if (total != header_.record_count || fstat(fd_, &status) != 0 ||
      status.st_size < 0 ||
      static_cast<uint64_t>(status.st_size) != expected_offset) {
    return Fail("ATAC hot spill size or record count mismatch in " + path_,
                error);
  }
  return true;
}

bool AtacHotSpillReader::OpenPartition(
    uint32_t rid, AtacHotSpillPartitionCursor *cursor,
    std::string *error) const {
  if (fd_ < 0 || cursor == nullptr || rid >= directory_.size()) {
    if (error != nullptr) {
      *error = "invalid ATAC hot spill partition request";
    }
    return false;
  }
  const auto &entry = directory_[rid];
  cursor->fd_ = fd_;
  cursor->path_ = path_;
  cursor->next_offset_ = entry.offset;
  cursor->remaining_records_ = entry.record_count;
  cursor->expected_records_ = entry.record_count;
  cursor->records_read_ = 0;
  cursor->input_record_count_ = header_.input_record_count;
  cursor->record_bytes_ = header_.record_bytes;
  cursor->barcode_length_ = header_.barcode_length;
  cursor->reference_length_ = reference_lengths_[rid];
  cursor->is_bulk_ = (header_.flags & kAtacHotSpillIsBulk) != 0;
  cursor->has_raw_barcode_evidence_ =
      (header_.flags & kAtacHotSpillHasRawBarcodeEvidence) != 0;
  cursor->has_y_hit_ = (header_.flags & kAtacHotSpillHasYHit) != 0;
  cursor->expected_first_start_ = entry.first_start;
  cursor->expected_last_start_ = entry.last_start;
  cursor->buffer_begin_ = 0;
  cursor->buffer_end_ = 0;
  cursor->have_previous_ = false;
  const size_t records_per_buffer = std::max<size_t>(
      1u, kHotReadBufferBytes / cursor->record_bytes_);
  cursor->buffer_.resize(records_per_buffer * cursor->record_bytes_);
  return true;
}

bool AtacHotSpillPartitionCursor::Fail(const std::string &message,
                                       std::string *error) {
  if (error != nullptr) {
    *error = message;
  }
  return false;
}

bool AtacHotSpillPartitionCursor::Refill(std::string *error) {
  if (buffer_begin_ != buffer_end_ || remaining_records_ == 0) {
    return true;
  }
  const uint64_t capacity_records = buffer_.size() / record_bytes_;
  const uint64_t records = std::min(remaining_records_, capacity_records);
  const size_t bytes = static_cast<size_t>(records * record_bytes_);
  if (!PreadAll(fd_, buffer_.data(), bytes, next_offset_)) {
    return Fail("cannot read ATAC hot spill partition in " + path_, error);
  }
  next_offset_ += bytes;
  buffer_begin_ = 0;
  buffer_end_ = bytes;
  return true;
}

bool AtacHotSpillPartitionCursor::ReadNext(
    AtacHotSpillDecodedRecord *record, bool *eof, std::string *error) {
  if (fd_ < 0 || record == nullptr || eof == nullptr) {
    return Fail("ATAC hot spill partition cursor is not initialized", error);
  }
  *eof = false;
  if (remaining_records_ == 0) {
    if (buffer_begin_ != buffer_end_ || records_read_ != expected_records_ ||
        (records_read_ != 0 &&
         (previous_mapping_.fragment_start_position_ !=
          expected_last_start_))) {
      return Fail("ATAC hot spill partition termination mismatch in " + path_,
                  error);
    }
    *eof = true;
    return true;
  }
  if (!Refill(error) || buffer_end_ - buffer_begin_ < record_bytes_) {
    return Fail("truncated ATAC hot spill record in " + path_, error);
  }
  AtacHotSpillRecordPrefixV1 prefix = {};
  memcpy(&prefix, buffer_.data() + buffer_begin_, sizeof(prefix));
  const char *quality =
      buffer_.data() + buffer_begin_ + sizeof(prefix);
  buffer_begin_ += record_bytes_;

  const uint8_t mapq = prefix.mapq_direction_unique & 0x3fu;
  const uint8_t direction =
      static_cast<uint8_t>((prefix.mapq_direction_unique >> 6) & 1u);
  const uint8_t is_unique =
      static_cast<uint8_t>((prefix.mapq_direction_unique >> 7) & 1u);
  const uint32_t valid_n_mask =
      barcode_length_ == 32
          ? std::numeric_limits<uint32_t>::max()
          : (barcode_length_ == 0 ? uint32_t{0}
                                  : (uint32_t{1} << barcode_length_) - 1u);
  const uint64_t valid_barcode_max =
      barcode_length_ == 32
          ? std::numeric_limits<uint64_t>::max()
          : (barcode_length_ == 0
                 ? uint64_t{0}
                 : (uint64_t{1} << (2u * barcode_length_)) - 1u);
  if (prefix.local_read_id >= input_record_count_ ||
      prefix.fragment_length == 0 ||
      static_cast<uint64_t>(prefix.start) + prefix.fragment_length >
          reference_length_ ||
      prefix.raw_barcode_key > valid_barcode_max ||
      (prefix.raw_barcode_n_mask & ~valid_n_mask) != 0 ||
      (prefix.row_flags & ~kAtacHotSpillRowHasYHit) != 0 ||
      ((prefix.row_flags & kAtacHotSpillRowHasYHit) != 0 &&
       !has_y_hit_)) {
    return Fail("invalid ATAC hot spill record in " + path_, error);
  }

  record->read_id_ = prefix.local_read_id;
  record->cell_barcode_ = prefix.raw_barcode_key;
  record->fragment_start_position_ = prefix.start;
  record->fragment_length_ = prefix.fragment_length;
  record->mapq_ = mapq;
  record->direction_ = direction;
  record->is_unique_ = is_unique;
  record->num_dups_ = 1;
  record->positive_alignment_length_ = prefix.positive_alignment_length;
  record->negative_alignment_length_ = prefix.negative_alignment_length;
  record->raw_barcode_n_mask = 0;
  record->raw_barcode_qual.clear();
  record->has_raw_barcode_evidence = has_raw_barcode_evidence_;
  record->has_y_hit =
      (prefix.row_flags & kAtacHotSpillRowHasYHit) != 0;
  if (has_raw_barcode_evidence_) {
    record->raw_barcode_n_mask = prefix.raw_barcode_n_mask;
    record->raw_barcode_qual.assign(quality, barcode_length_);
  }
  if ((records_read_ == 0 && prefix.start != expected_first_start_) ||
      (have_previous_ &&
       static_cast<const PairedEndMappingWithBarcode &>(*record) <
           previous_mapping_)) {
    return Fail("unsorted ATAC hot spill partition in " + path_, error);
  }
  previous_mapping_ =
      static_cast<const PairedEndMappingWithBarcode &>(*record);
  have_previous_ = true;
  --remaining_records_;
  ++records_read_;
  return true;
}

}  // namespace chromap
