#include "atac_materialized_binary.h"

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <sys/stat.h>
#include <unistd.h>

#include "barcode_translator.h"

namespace chromap {
namespace {

constexpr uint64_t kMaximumMetadataBytes = 256u * 1024u * 1024u;

bool WriteBytes(FILE *file, const void *data, size_t bytes) {
  return bytes == 0 || fwrite(data, 1, bytes, file) == bytes;
}

bool PreadAll(int fd, void *data, size_t bytes, uint64_t offset) {
  char *cursor = static_cast<char *>(data);
  size_t remaining = bytes;
  while (remaining != 0) {
    const ssize_t got = pread(fd, cursor, remaining, static_cast<off_t>(offset));
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

bool PwriteAll(int fd, const void *data, size_t bytes, uint64_t offset) {
  const char *cursor = static_cast<const char *>(data);
  size_t remaining = bytes;
  while (remaining != 0) {
    const ssize_t put =
        pwrite(fd, cursor, remaining, static_cast<off_t>(offset));
    if (put < 0 && errno == EINTR) {
      continue;
    }
    if (put <= 0) {
      return false;
    }
    cursor += put;
    remaining -= static_cast<size_t>(put);
    offset += static_cast<uint64_t>(put);
  }
  return true;
}

void AppendUnsigned(std::string *output, uint64_t value) {
  char digits[32];
  char *end = digits + sizeof(digits);
  char *begin = end;
  do {
    *--begin = static_cast<char>('0' + value % 10u);
    value /= 10u;
  } while (value != 0);
  output->append(begin, static_cast<size_t>(end - begin));
}

void AppendPackedBarcode(std::string *output, uint32_t key,
                         uint32_t barcode_length) {
  for (uint32_t i = 0; i < barcode_length; ++i) {
    const uint32_t shift = 2u * (barcode_length - 1u - i);
    output->push_back(Uint8ToChar((key >> shift) & 3u));
  }
}

struct BinaryInput {
  int fd = -1;
  AtacMaterializedBinaryHeaderV1 header{};
  AtacMaterializedBinaryMetadata metadata;
  std::vector<uint64_t> barcode_dictionary;
  std::vector<AtacMaterializedBinaryBlockV1> blocks;

  ~BinaryInput() {
    if (fd >= 0) {
      close(fd);
    }
  }
};

template <typename T>
bool ReadMetadataValue(const std::vector<char> &bytes, size_t *offset,
                       T *value) {
  if (*offset > bytes.size() || bytes.size() - *offset < sizeof(T)) {
    return false;
  }
  memcpy(value, bytes.data() + *offset, sizeof(T));
  *offset += sizeof(T);
  return true;
}

bool ReadMetadataString(const std::vector<char> &bytes, size_t *offset,
                        size_t count, std::string *value) {
  if (*offset > bytes.size() || bytes.size() - *offset < count) {
    return false;
  }
  value->assign(bytes.data() + *offset, count);
  *offset += count;
  return true;
}

bool OpenBinaryInput(const std::string &path, BinaryInput *input,
                     std::string *error) {
  input->fd = open(path.c_str(), O_RDONLY);
  if (input->fd < 0) {
    *error = "cannot open ATAC materialized binary " + path + ": " +
             std::strerror(errno);
    return false;
  }
  struct stat info = {};
  if (fstat(input->fd, &info) != 0 || info.st_size < 0 ||
      !PreadAll(input->fd, &input->header, sizeof(input->header), 0)) {
    *error = "cannot read ATAC materialized binary header from " + path;
    return false;
  }
  const auto &header = input->header;
  if (memcmp(header.magic, kAtacMaterializedBinaryMagicV1,
             sizeof(header.magic)) != 0 ||
      header.format_version != kAtacMaterializedBinaryFormatVersion ||
      header.fixed_header_bytes != sizeof(header) ||
      header.record_bytes != sizeof(AtacMaterializedBinaryRecordV1) ||
      header.metadata_bytes > kMaximumMetadataBytes ||
      header.num_reference_sequences == 0 || header.shard_count == 0 ||
      header.sample_id_bytes == 0 || header.input_id_bytes == 0 ||
      header.target_records_per_block == 0 ||
      header.endian_marker != kAtacMaterializedBinaryEndianMarker ||
      (header.flags & ~(kAtacMaterializedBinaryIsBulk |
                        kAtacMaterializedBinaryHasBarcodeDictionary)) != 0 ||
      header.data_offset != sizeof(header) + header.metadata_bytes ||
      header.barcode_dictionary_offset < header.data_offset ||
      header.directory_offset < header.barcode_dictionary_offset ||
      header.directory_offset > static_cast<uint64_t>(info.st_size)) {
    *error = "invalid ATAC materialized binary header in " + path;
    return false;
  }
  const uint64_t dictionary_bytes =
      header.barcode_dictionary_count * sizeof(uint64_t);
  const uint64_t file_bytes = static_cast<uint64_t>(info.st_size);
  if (header.num_reference_sequences >
          static_cast<uint32_t>(std::numeric_limits<uint16_t>::max()) + 1u ||
      header.barcode_dictionary_count >
          std::numeric_limits<uint64_t>::max() / sizeof(uint64_t) ||
      dictionary_bytes !=
          header.directory_offset - header.barcode_dictionary_offset ||
      header.block_count >
          std::numeric_limits<uint64_t>::max() /
              sizeof(AtacMaterializedBinaryBlockV1) ||
      header.block_count * sizeof(AtacMaterializedBinaryBlockV1) !=
          file_bytes - header.directory_offset) {
    *error = "invalid ATAC materialized binary directory bounds in " + path;
    return false;
  }

  std::vector<char> metadata(static_cast<size_t>(header.metadata_bytes));
  if (!PreadAll(input->fd, metadata.data(), metadata.size(), sizeof(header))) {
    *error = "truncated ATAC materialized binary metadata in " + path;
    return false;
  }
  size_t offset = 0;
  if (!ReadMetadataString(metadata, &offset, header.sample_id_bytes,
                          &input->metadata.sample_id) ||
      !ReadMetadataString(metadata, &offset, header.input_id_bytes,
                          &input->metadata.input_id)) {
    *error = "truncated ATAC materialized binary identity in " + path;
    return false;
  }
  input->metadata.is_bulk =
      (header.flags & kAtacMaterializedBinaryIsBulk) != 0;
  input->metadata.use_barcode_dictionary =
      (header.flags & kAtacMaterializedBinaryHasBarcodeDictionary) != 0;
  input->metadata.barcode_length = header.barcode_length;
  input->metadata.shard_count = header.shard_count;
  input->metadata.references.reserve(header.num_reference_sequences);
  for (uint32_t rid = 0; rid < header.num_reference_sequences; ++rid) {
    uint32_t name_bytes = 0;
    AtacMergeableSpillReference reference;
    if (!ReadMetadataValue(metadata, &offset, &name_bytes) ||
        !ReadMetadataValue(metadata, &offset, &reference.length) ||
        name_bytes == 0 ||
        !ReadMetadataString(metadata, &offset, name_bytes, &reference.name)) {
      *error = "invalid ATAC materialized binary reference metadata in " +
               path;
      return false;
    }
    input->metadata.references.push_back(std::move(reference));
  }
  if (offset != metadata.size()) {
    *error = "ATAC materialized binary metadata length mismatch in " + path;
    return false;
  }

  const bool is_bulk = input->metadata.is_bulk;
  const bool use_barcode_dictionary =
      input->metadata.use_barcode_dictionary;
  if ((!is_bulk &&
       (header.barcode_length == 0 || header.barcode_length > 32)) ||
      (is_bulk && (use_barcode_dictionary ||
                   header.barcode_dictionary_count != 0)) ||
      (!is_bulk && header.barcode_length > 16 &&
       !use_barcode_dictionary) ||
      (!use_barcode_dictionary && header.barcode_dictionary_count != 0) ||
      (use_barcode_dictionary && header.record_count != 0 &&
       header.barcode_dictionary_count == 0) ||
      header.barcode_dictionary_count >
          static_cast<uint64_t>(std::numeric_limits<uint32_t>::max()) + 1u) {
    *error = "invalid ATAC materialized binary barcode dictionary in " + path;
    return false;
  }
  input->barcode_dictionary.resize(
      static_cast<size_t>(header.barcode_dictionary_count));
  if (!PreadAll(input->fd, input->barcode_dictionary.data(),
                input->barcode_dictionary.size() * sizeof(uint64_t),
                header.barcode_dictionary_offset)) {
    *error = "truncated ATAC materialized binary barcode dictionary in " +
             path;
    return false;
  }

  input->blocks.resize(static_cast<size_t>(header.block_count));
  if (!PreadAll(input->fd, input->blocks.data(),
                input->blocks.size() * sizeof(input->blocks[0]),
                header.directory_offset)) {
    *error = "truncated ATAC materialized binary directory in " + path;
    return false;
  }
  uint64_t records = 0;
  uint64_t expected_offset = header.data_offset;
  bool have_previous = false;
  uint32_t previous_rid = 0;
  uint32_t previous_start = 0;
  for (const auto &block : input->blocks) {
    if (block.record_count == 0 || block.offset != expected_offset ||
        block.record_count >
            std::numeric_limits<uint64_t>::max() /
                sizeof(AtacMaterializedBinaryRecordV1)) {
      *error = "invalid ATAC materialized binary block in " + path;
      return false;
    }
    const uint64_t bytes =
        block.record_count * sizeof(AtacMaterializedBinaryRecordV1);
    if (expected_offset > header.barcode_dictionary_offset ||
        bytes > header.barcode_dictionary_offset - expected_offset ||
        block.first_rid >= header.num_reference_sequences ||
        block.last_rid >= header.num_reference_sequences ||
        block.first_rid > block.last_rid ||
        (have_previous &&
         (block.first_rid < previous_rid ||
          (block.first_rid == previous_rid &&
           block.first_start < previous_start)))) {
      *error = "invalid ATAC materialized binary block range in " + path;
      return false;
    }
    expected_offset += bytes;
    records += block.record_count;
    previous_rid = block.last_rid;
    previous_start = block.last_start;
    have_previous = true;
  }
  if (expected_offset != header.barcode_dictionary_offset ||
      records != header.record_count ||
      ((header.record_count == 0) != input->blocks.empty())) {
    *error = "ATAC materialized binary record count mismatch in " + path;
    return false;
  }
  return true;
}

bool FormatBlock(const BinaryInput &input,
                 const AtacMaterializedBinaryBlockV1 &block,
                 const std::vector<std::string> &barcode_text,
                 BarcodeTranslator *translator, bool translate_barcode,
                 std::string *output,
                 std::string *error) {
  std::vector<AtacMaterializedBinaryRecordV1> records(
      static_cast<size_t>(block.record_count));
  if (!PreadAll(input.fd, records.data(),
                records.size() * sizeof(records[0]), block.offset)) {
    *error = "cannot read ATAC materialized binary block";
    return false;
  }
  output->clear();
  output->reserve(records.size() * (input.metadata.is_bulk ? 40u : 42u));
  const uint32_t direct_barcode_max =
      input.metadata.barcode_length >= 16
          ? std::numeric_limits<uint32_t>::max()
          : (input.metadata.barcode_length == 0
                 ? uint32_t{0}
                 : (uint32_t{1} << (2u * input.metadata.barcode_length)) - 1u);
  for (const auto &record : records) {
    const uint64_t end =
        static_cast<uint64_t>(record.start) + record.fragment_length;
    if (record.rid >= input.metadata.references.size() ||
        record.fragment_length == 0 ||
        end > input.metadata.references[record.rid].length ||
        record.duplicate_count == 0 ||
        record.mapq > 63 || record.reserved != 0 ||
        (record.flags & ~kAtacMaterializedBinaryPositiveStrand) != 0 ||
        (input.metadata.is_bulk && record.barcode_value != 0) ||
        (input.metadata.use_barcode_dictionary &&
         record.barcode_value >= input.barcode_dictionary.size()) ||
        (!input.metadata.is_bulk &&
         !input.metadata.use_barcode_dictionary &&
         record.barcode_value > direct_barcode_max)) {
      *error = "invalid ATAC materialized binary record";
      return false;
    }
    output->append(input.metadata.references[record.rid].name);
    output->push_back('\t');
    AppendUnsigned(output, record.start);
    output->push_back('\t');
    AppendUnsigned(output, end);
    output->push_back('\t');
    if (input.metadata.is_bulk) {
      output->append("N\t");
      AppendUnsigned(output, record.mapq);
      output->push_back('\t');
      output->push_back(
          (record.flags & kAtacMaterializedBinaryPositiveStrand) != 0 ? '+'
                                                                      : '-');
      output->push_back('\t');
    } else {
      if (input.metadata.use_barcode_dictionary) {
        output->append(barcode_text[record.barcode_value]);
      } else if (translate_barcode) {
        output->append(translator->Translate(record.barcode_value,
                                             input.metadata.barcode_length));
      } else {
        AppendPackedBarcode(output, record.barcode_value,
                            input.metadata.barcode_length);
      }
      output->push_back('\t');
    }
    AppendUnsigned(output, record.duplicate_count);
    output->push_back('\n');
  }
  return true;
}

}  // namespace

AtacMaterializedBinaryWriter::~AtacMaterializedBinaryWriter() {
  if (file_ != nullptr) {
    fclose(file_);
  }
  if (!temporary_path_.empty() && !finalized_) {
    unlink(temporary_path_.c_str());
  }
}

bool AtacMaterializedBinaryWriter::Fail(const std::string &message,
                                        std::string *error) {
  if (error != nullptr) {
    *error = message;
  }
  return false;
}

bool AtacMaterializedBinaryWriter::Open(
    const std::string &path, const AtacMaterializedBinaryMetadata &metadata,
    std::string *error) {
  if (file_ != nullptr || finalized_ || path.empty()) {
    return Fail("ATAC materialized binary writer is not initializable", error);
  }
  if (metadata.sample_id.empty() || metadata.input_id.empty() ||
      metadata.references.empty() || metadata.shard_count == 0 ||
      (!metadata.is_bulk &&
       (metadata.barcode_length == 0 || metadata.barcode_length > 32))) {
    return Fail("invalid ATAC materialized binary metadata", error);
  }
  uint64_t metadata_bytes = metadata.sample_id.size() + metadata.input_id.size();
  for (const auto &reference : metadata.references) {
    if (reference.name.empty() ||
        metadata_bytes > std::numeric_limits<uint64_t>::max() -
                             2u * sizeof(uint32_t) - reference.name.size()) {
      return Fail("ATAC materialized binary metadata overflows", error);
    }
    metadata_bytes += 2u * sizeof(uint32_t) + reference.name.size();
  }
  if (metadata_bytes > kMaximumMetadataBytes ||
      metadata.sample_id.size() > std::numeric_limits<uint32_t>::max() ||
      metadata.input_id.size() > std::numeric_limits<uint32_t>::max() ||
      metadata.references.size() >
          static_cast<uint64_t>(std::numeric_limits<uint16_t>::max()) + 1u) {
    return Fail("ATAC materialized binary metadata is too large", error);
  }

  output_path_ = path;
  temporary_path_ = path + ".tmp";
  file_ = fopen(temporary_path_.c_str(), "wb+");
  if (file_ == nullptr) {
    return Fail("cannot open ATAC materialized binary temporary output " +
                    temporary_path_ + ": " + std::strerror(errno),
                error);
  }
  (void)setvbuf(file_, nullptr, _IOFBF, 8u * 1024u * 1024u);
  metadata_ = metadata;
  metadata_.use_barcode_dictionary =
      !metadata.is_bulk &&
      (metadata.use_barcode_dictionary || metadata.barcode_length > 16);
  memcpy(header_.magic, kAtacMaterializedBinaryMagicV1,
         sizeof(header_.magic));
  header_.format_version = kAtacMaterializedBinaryFormatVersion;
  header_.fixed_header_bytes = sizeof(header_);
  header_.record_bytes = sizeof(AtacMaterializedBinaryRecordV1);
  header_.flags = metadata.is_bulk ? kAtacMaterializedBinaryIsBulk : 0;
  if (metadata_.use_barcode_dictionary) {
    header_.flags |= kAtacMaterializedBinaryHasBarcodeDictionary;
  }
  header_.barcode_length = metadata.barcode_length;
  header_.num_reference_sequences =
      static_cast<uint32_t>(metadata.references.size());
  header_.shard_count = metadata.shard_count;
  header_.sample_id_bytes = static_cast<uint32_t>(metadata.sample_id.size());
  header_.input_id_bytes = static_cast<uint32_t>(metadata.input_id.size());
  header_.metadata_bytes = metadata_bytes;
  header_.data_offset = sizeof(header_) + metadata_bytes;
  header_.target_records_per_block =
      kAtacMaterializedBinaryTargetRecordsPerBlock;
  header_.endian_marker = kAtacMaterializedBinaryEndianMarker;
  if (!WriteBytes(file_, &header_, sizeof(header_)) ||
      !WriteBytes(file_, metadata.sample_id.data(), metadata.sample_id.size()) ||
      !WriteBytes(file_, metadata.input_id.data(), metadata.input_id.size())) {
    return Fail("cannot write ATAC materialized binary header", error);
  }
  for (const auto &reference : metadata.references) {
    const uint32_t name_bytes = static_cast<uint32_t>(reference.name.size());
    if (!WriteBytes(file_, &name_bytes, sizeof(name_bytes)) ||
        !WriteBytes(file_, &reference.length, sizeof(reference.length)) ||
        !WriteBytes(file_, reference.name.data(), reference.name.size())) {
      return Fail("cannot write ATAC materialized binary reference metadata",
                  error);
    }
  }
  block_records_.reserve(kAtacMaterializedBinaryTargetRecordsPerBlock);
  return true;
}

bool AtacMaterializedBinaryWriter::Append(uint32_t rid,
                                          const AtacSpillRecord &mapping,
                                          uint64_t duplicate_count,
                                          std::string *error) {
  const uint64_t end = static_cast<uint64_t>(mapping.GetStartPosition()) +
                       mapping.fragment_length_;
  if (file_ == nullptr || finalized_ || rid >= metadata_.references.size() ||
      mapping.fragment_length_ == 0 ||
      end > metadata_.references[rid].length || duplicate_count == 0) {
    return Fail("invalid ATAC materialized binary record", error);
  }
  if (have_previous_ &&
      (rid < previous_rid_ ||
       (rid == previous_rid_ &&
        mapping.GetStartPosition() < previous_start_))) {
    return Fail("ATAC materialized binary records are not sorted", error);
  }
  AtacMaterializedBinaryRecordV1 record = {};
  record.start = mapping.GetStartPosition();
  record.rid = static_cast<uint16_t>(rid);
  record.fragment_length = mapping.fragment_length_;
  record.duplicate_count = static_cast<uint8_t>(std::min<uint64_t>(
      std::numeric_limits<uint8_t>::max(), duplicate_count));
  record.mapq = mapping.mapq_;
  record.flags = mapping.IsPositiveStrand()
                     ? kAtacMaterializedBinaryPositiveStrand
                     : uint8_t{0};
  return AppendEncodedRecordWithRawBarcodeKey(
      record, metadata_.is_bulk ? uint64_t{0} : mapping.GetBarcode(), error);
}

bool AtacMaterializedBinaryWriter::AppendEncodedRecords(
    const AtacMaterializedBinaryRecordV1 *records, size_t count,
    std::string *error) {
  if (count != 0 && records == nullptr) {
    return Fail("null ATAC materialized encoded-record batch", error);
  }
  for (size_t i = 0; i < count; ++i) {
    if (!AppendEncodedRecord(records[i], /*allow_dictionary=*/false, error)) {
      return false;
    }
  }
  return true;
}

bool AtacMaterializedBinaryWriter::AppendEncodedRecordsWithRawBarcodes(
    const AtacMaterializedBinaryRecordV1 *records, size_t count,
    std::string *error) {
  if (count != 0 && records == nullptr) {
    return Fail("null ATAC materialized raw-barcode record batch", error);
  }
  if (!metadata_.is_bulk && metadata_.barcode_length > 16) {
    return Fail("raw-barcode record batch exceeds direct 32-bit encoding",
                error);
  }
  for (size_t i = 0; i < count; ++i) {
    if (!AppendEncodedRecordWithRawBarcodeKey(
            records[i], static_cast<uint64_t>(records[i].barcode_value),
            error)) {
      return false;
    }
  }
  return true;
}

bool AtacMaterializedBinaryWriter::AppendEncodedRecordsWithRawBarcodeKeys(
    const AtacMaterializedBinaryRecordV1 *records,
    const uint64_t *barcode_keys, size_t count, std::string *error) {
  if (count != 0 && (records == nullptr || barcode_keys == nullptr)) {
    return Fail("null ATAC materialized wide-barcode record batch", error);
  }
  for (size_t i = 0; i < count; ++i) {
    if (!AppendEncodedRecordWithRawBarcodeKey(records[i], barcode_keys[i],
                                               error)) {
      return false;
    }
  }
  return true;
}

bool AtacMaterializedBinaryWriter::AppendEncodedRecordWithRawBarcodeKey(
    AtacMaterializedBinaryRecordV1 record, uint64_t barcode_key,
    std::string *error) {
  const uint64_t packed_barcode_max =
      metadata_.is_bulk || metadata_.barcode_length == 0
          ? uint64_t{0}
          : (metadata_.barcode_length == 32
                 ? std::numeric_limits<uint64_t>::max()
                 : (uint64_t{1} << (2u * metadata_.barcode_length)) - 1u);
  if ((metadata_.is_bulk && barcode_key != 0) ||
      (!metadata_.is_bulk && barcode_key > packed_barcode_max)) {
    return Fail("invalid ATAC materialized raw barcode", error);
  }
  if (!metadata_.is_bulk && metadata_.use_barcode_dictionary) {
    const auto found = barcode_ids_.find(barcode_key);
    if (found == barcode_ids_.end()) {
      if (barcode_dictionary_.size() >
          std::numeric_limits<uint32_t>::max()) {
        return Fail("ATAC materialized binary has more than 2^32 barcodes",
                    error);
      }
      record.barcode_value =
          static_cast<uint32_t>(barcode_dictionary_.size());
      barcode_ids_.emplace(barcode_key, record.barcode_value);
      barcode_dictionary_.push_back(barcode_key);
    } else {
      record.barcode_value = found->second;
    }
  } else {
    if (barcode_key > std::numeric_limits<uint32_t>::max()) {
      return Fail("ATAC materialized barcode exceeds direct 32-bit encoding",
                  error);
    }
    record.barcode_value = static_cast<uint32_t>(barcode_key);
  }
  return AppendEncodedRecord(record, metadata_.use_barcode_dictionary, error);
}

bool AtacMaterializedBinaryWriter::AppendEncodedRecord(
    const AtacMaterializedBinaryRecordV1 &record, bool allow_dictionary,
    std::string *error) {
  const uint64_t end =
      static_cast<uint64_t>(record.start) + record.fragment_length;
  const uint32_t direct_barcode_max =
      metadata_.barcode_length >= 16
          ? std::numeric_limits<uint32_t>::max()
          : (metadata_.barcode_length == 0
                 ? uint32_t{0}
                 : (uint32_t{1} << (2u * metadata_.barcode_length)) - 1u);
  if (file_ == nullptr || finalized_ ||
      record.rid >= metadata_.references.size() ||
      record.fragment_length == 0 ||
      end > metadata_.references[record.rid].length ||
      record.duplicate_count == 0 || record.mapq > 63 ||
      record.reserved != 0 ||
      (record.flags & ~kAtacMaterializedBinaryPositiveStrand) != 0 ||
      (metadata_.is_bulk && record.barcode_value != 0) ||
      (!metadata_.is_bulk && metadata_.use_barcode_dictionary &&
       (!allow_dictionary ||
        record.barcode_value >= barcode_dictionary_.size())) ||
      (!metadata_.is_bulk && !metadata_.use_barcode_dictionary &&
       record.barcode_value > direct_barcode_max) ||
      (have_previous_ &&
       (record.rid < previous_rid_ ||
        (record.rid == previous_rid_ && record.start < previous_start_)))) {
    return Fail("invalid or unsorted ATAC materialized encoded record", error);
  }
  block_records_.push_back(record);
  ++record_count_;
  previous_rid_ = record.rid;
  previous_start_ = record.start;
  have_previous_ = true;
  if (block_records_.size() ==
      kAtacMaterializedBinaryTargetRecordsPerBlock) {
    return FlushBlock(error);
  }
  return true;
}

bool AtacMaterializedBinaryWriter::FlushBlock(std::string *error) {
  if (block_records_.empty()) {
    return true;
  }
  const off_t offset = ftello(file_);
  if (offset < 0) {
    return Fail("cannot locate ATAC materialized binary block", error);
  }
  AtacMaterializedBinaryBlockV1 block = {};
  block.offset = static_cast<uint64_t>(offset);
  block.record_count = block_records_.size();
  block.first_rid = block_records_.front().rid;
  block.last_rid = block_records_.back().rid;
  block.first_start = block_records_.front().start;
  block.last_start = block_records_.back().start;
  block.first_end = block_records_.front().start +
                    block_records_.front().fragment_length;
  block.last_end =
      block_records_.back().start + block_records_.back().fragment_length;
  if (!WriteBytes(file_, block_records_.data(),
                  block_records_.size() * sizeof(block_records_[0]))) {
    return Fail("cannot write ATAC materialized binary block", error);
  }
  blocks_.push_back(block);
  block_records_.clear();
  return true;
}

bool AtacMaterializedBinaryWriter::Finalize(std::string *error) {
  if (file_ == nullptr || finalized_ || !FlushBlock(error)) {
    return false;
  }
  const off_t barcode_dictionary_offset = ftello(file_);
  if (barcode_dictionary_offset < 0 ||
      !WriteBytes(file_, barcode_dictionary_.data(),
                  barcode_dictionary_.size() * sizeof(uint64_t))) {
    return Fail("cannot write ATAC materialized binary barcode dictionary",
                error);
  }
  const off_t directory_offset = ftello(file_);
  if (directory_offset < 0 ||
      !WriteBytes(file_, blocks_.data(),
                  blocks_.size() * sizeof(blocks_[0])) ||
      fflush(file_) != 0) {
    return Fail("cannot finalize ATAC materialized binary directory", error);
  }
  header_.record_count = record_count_;
  header_.block_count = blocks_.size();
  header_.barcode_dictionary_count = barcode_dictionary_.size();
  header_.barcode_dictionary_offset =
      static_cast<uint64_t>(barcode_dictionary_offset);
  header_.directory_offset = static_cast<uint64_t>(directory_offset);
  if (fseeko(file_, 0, SEEK_SET) != 0 ||
      !WriteBytes(file_, &header_, sizeof(header_)) || fflush(file_) != 0 ||
      fclose(file_) != 0) {
    file_ = nullptr;
    return Fail("cannot commit ATAC materialized binary header", error);
  }
  file_ = nullptr;
  if (rename(temporary_path_.c_str(), output_path_.c_str()) != 0) {
    return Fail("cannot publish ATAC materialized binary " + output_path_ +
                    ": " + std::strerror(errno),
                error);
  }
  finalized_ = true;
  return true;
}

bool ExportAtacMaterializedBinaryToBed(
    const std::string &binary_path, const std::string &bed_path,
    const MappingParameters &parameters, std::string *error) {
  BinaryInput input;
  if (!OpenBinaryInput(binary_path, &input, error)) {
    return false;
  }
  BarcodeTranslator translator;
  if (!parameters.barcode_translate_table_file_path.empty()) {
    translator.SetTranslateTable(
        parameters.barcode_translate_table_file_path,
        parameters.barcode_translate_from_first_column);
  }
  std::vector<std::string> barcode_text(input.barcode_dictionary.size());
  for (size_t i = 0; i < input.barcode_dictionary.size(); ++i) {
    barcode_text[i] = translator.Translate(input.barcode_dictionary[i],
                                           input.metadata.barcode_length);
  }
  const bool translate_barcode =
      !parameters.barcode_translate_table_file_path.empty();
  const std::string temporary_path = bed_path + ".tmp";
  const int output_fd =
      open(temporary_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0666);
  if (output_fd < 0) {
    *error = "cannot open terminal BED output " + temporary_path + ": " +
             std::strerror(errno);
    return false;
  }

  const int effective_threads = std::max(1, parameters.num_threads);
  const size_t requested_threads = static_cast<size_t>(effective_threads);
  const size_t wave_size = std::max<size_t>(
      1, std::min(requested_threads, input.blocks.size()));
  uint64_t output_offset = 0;
  bool ok = true;
  for (size_t begin = 0; begin < input.blocks.size() && ok;
       begin += wave_size) {
    const size_t count = std::min(wave_size, input.blocks.size() - begin);
    std::vector<std::string> buffers(count);
    std::vector<std::string> errors(count);
    std::atomic<bool> format_failed(false);
#pragma omp parallel for schedule(dynamic, 1) num_threads(effective_threads)
    for (size_t i = 0; i < count; ++i) {
      if (!FormatBlock(input, input.blocks[begin + i], barcode_text,
                       &translator, translate_barcode,
                       &buffers[i], &errors[i])) {
        format_failed.store(true, std::memory_order_relaxed);
      }
    }
    if (format_failed.load(std::memory_order_relaxed)) {
      for (const auto &message : errors) {
        if (!message.empty()) {
          *error = message;
          break;
        }
      }
      ok = false;
      break;
    }
    std::vector<uint64_t> offsets(count, output_offset);
    for (size_t i = 0; i < count; ++i) {
      offsets[i] = output_offset;
      output_offset += buffers[i].size();
    }
    std::atomic<bool> write_failed(false);
#pragma omp parallel for schedule(static) num_threads(effective_threads)
    for (size_t i = 0; i < count; ++i) {
      if (!PwriteAll(output_fd, buffers[i].data(), buffers[i].size(),
                     offsets[i])) {
        write_failed.store(true, std::memory_order_relaxed);
      }
    }
    if (write_failed.load(std::memory_order_relaxed)) {
      *error = "cannot write terminal BED output " + temporary_path;
      ok = false;
    }
  }
  if (ok && ftruncate(output_fd, static_cast<off_t>(output_offset)) != 0) {
    *error = "cannot size terminal BED output " + temporary_path;
    ok = false;
  }
  if (close(output_fd) != 0 && ok) {
    *error = "cannot close terminal BED output " + temporary_path;
    ok = false;
  }
  if (ok && rename(temporary_path.c_str(), bed_path.c_str()) != 0) {
    *error = "cannot publish terminal BED output " + bed_path + ": " +
             std::strerror(errno);
    ok = false;
  }
  if (!ok) {
    unlink(temporary_path.c_str());
  }
  return ok;
}

}  // namespace chromap
