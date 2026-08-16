#ifndef ATAC_HOT_SPILL_H_
#define ATAC_HOT_SPILL_H_

#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>

#include "atac_mergeable_spill.h"

namespace chromap {

// Fixed-width, coordinate-partitioned companion to an ATACMS3 worker shard.
// It contains only fields needed by BED correction/deduplication. The full
// ATACMS3 remains the compatibility source for BAM/CRAM materialization.
#pragma pack(push, 1)
struct AtacHotSpillHeaderV1 {
  char magic[8];
  uint16_t format_version;
  uint16_t fixed_header_bytes;
  uint16_t record_prefix_bytes;
  uint16_t flags;
  uint32_t endian_marker;
  uint32_t shard_ordinal;
  uint32_t shard_count;
  uint32_t barcode_length;
  uint32_t num_reference_sequences;
  uint64_t first_global_read_ordinal;
  uint64_t input_record_count;
  uint64_t record_count;
  uint64_t barcode_whitelist_fingerprint;
  uint64_t directory_offset;
  uint64_t data_offset;
  uint32_t record_bytes;
  uint32_t reserved;
};

struct AtacHotSpillDirectoryV1 {
  uint64_t offset;
  uint64_t record_count;
  uint32_t first_start;
  uint32_t last_start;
};

struct AtacHotSpillRecordPrefixV1 {
  uint64_t local_read_id;
  uint64_t raw_barcode_key;
  uint32_t start;
  uint32_t raw_barcode_n_mask;
  uint16_t fragment_length;
  uint16_t positive_alignment_length;
  uint16_t negative_alignment_length;
  // Bits 0..5 MAPQ, bit 6 direction, bit 7 unique.
  uint8_t mapq_direction_unique;
  // Bit 0 records the row-level Y-hit state.
  uint8_t row_flags;
};
#pragma pack(pop)

static_assert(sizeof(AtacHotSpillHeaderV1) == 92,
              "ATAC hot spill header must be 92 bytes");
static_assert(sizeof(AtacHotSpillDirectoryV1) == 24,
              "ATAC hot spill directory entry must be 24 bytes");
static_assert(sizeof(AtacHotSpillRecordPrefixV1) == 32,
              "ATAC hot spill prefix must be 32 bytes");

static constexpr char kAtacHotSpillMagicV1[8] = {
    'A', 'T', 'A', 'C', 'H', 'O', 'T', '1'};
static constexpr uint16_t kAtacHotSpillFormatVersion = 1;
static constexpr uint32_t kAtacHotSpillEndianMarker = 0x01020304u;
static constexpr uint16_t kAtacHotSpillIsBulk = 1u << 0;
static constexpr uint16_t kAtacHotSpillHasRawBarcodeEvidence = 1u << 1;
static constexpr uint16_t kAtacHotSpillHasYHit = 1u << 2;
static constexpr uint8_t kAtacHotSpillRowHasYHit = 1u << 0;

std::string AtacHotSpillSidecarPath(const std::string &spill_path);

// In-memory view used by the BED materializer. Keep this separate from
// AtacSpillRecord: the latter owns two complete SAMMapping payloads that are
// intentionally absent from ATACHOT1 and expensive to construct per record.
struct AtacHotSpillDecodedRecord : public PairedEndMappingWithBarcode {
  uint32_t raw_barcode_n_mask = 0;
  AtacBarcodeQuality raw_barcode_qual;
  bool has_raw_barcode_evidence = false;
  bool has_y_hit = false;

  bool operator<(const AtacHotSpillDecodedRecord &other) const {
    return static_cast<const PairedEndMappingWithBarcode &>(*this) <
           static_cast<const PairedEndMappingWithBarcode &>(other);
  }
  bool operator==(const AtacHotSpillDecodedRecord &other) const {
    return static_cast<const PairedEndMappingWithBarcode &>(*this) ==
           static_cast<const PairedEndMappingWithBarcode &>(other);
  }
  bool IsSamePosition(const AtacHotSpillDecodedRecord &other) const {
    return static_cast<const PairedEndMappingWithBarcode &>(*this)
        .IsSamePosition(
            static_cast<const PairedEndMappingWithBarcode &>(other));
  }
};

class AtacHotSpillWriter {
 public:
  AtacHotSpillWriter() = default;
  ~AtacHotSpillWriter();

  bool Open(const std::string &path,
            const AtacMergeableSpillMetadata &metadata,
            std::string *error);
  bool Append(uint32_t rid, const AtacSpillRecord &record,
              std::string *error);
  bool Finalize(std::string *error);
  uint64_t record_count() const { return record_count_; }

 private:
  bool Fail(const std::string &message, std::string *error);

  FILE *file_ = nullptr;
  std::string output_path_;
  std::string temporary_path_;
  AtacMergeableSpillMetadata metadata_;
  AtacHotSpillHeaderV1 header_{};
  std::vector<AtacHotSpillDirectoryV1> directory_;
  std::vector<char> encoded_record_;
  uint64_t record_count_ = 0;
  bool finalized_ = false;
  bool have_previous_ = false;
  uint32_t previous_rid_ = 0;
  PairedEndMappingWithBarcode previous_mapping_;
};

class AtacHotSpillPartitionCursor;

class AtacHotSpillReader {
 public:
  AtacHotSpillReader() = default;
  ~AtacHotSpillReader();

  bool Open(const std::string &parent_spill_path,
            const AtacMergeableSpillMetadata &parent_metadata,
            uint64_t parent_record_count, std::string *error);
  bool OpenPartition(uint32_t rid, AtacHotSpillPartitionCursor *cursor,
                     std::string *error) const;
  uint64_t record_count() const { return header_.record_count; }
  uint64_t PartitionRecordCount(uint32_t rid) const {
    return rid < directory_.size() ? directory_[rid].record_count : 0;
  }
  const std::string &path() const { return path_; }

 private:
  bool Fail(const std::string &message, std::string *error);

  int fd_ = -1;
  std::string path_;
  AtacHotSpillHeaderV1 header_{};
  std::vector<AtacHotSpillDirectoryV1> directory_;
  std::vector<uint32_t> reference_lengths_;
};

class AtacHotSpillPartitionCursor {
 public:
  AtacHotSpillPartitionCursor() = default;

  bool ReadNext(AtacHotSpillDecodedRecord *record, bool *eof,
                std::string *error);

 private:
  friend class AtacHotSpillReader;
  bool Fail(const std::string &message, std::string *error);
  bool Refill(std::string *error);

  int fd_ = -1;
  std::string path_;
  uint64_t next_offset_ = 0;
  uint64_t remaining_records_ = 0;
  uint64_t expected_records_ = 0;
  uint64_t records_read_ = 0;
  uint64_t input_record_count_ = 0;
  uint32_t record_bytes_ = 0;
  uint32_t barcode_length_ = 0;
  uint32_t reference_length_ = 0;
  bool is_bulk_ = false;
  bool has_raw_barcode_evidence_ = false;
  bool has_y_hit_ = false;
  uint32_t expected_first_start_ = 0;
  uint32_t expected_last_start_ = 0;
  std::vector<char> buffer_;
  size_t buffer_begin_ = 0;
  size_t buffer_end_ = 0;
  bool have_previous_ = false;
  PairedEndMappingWithBarcode previous_mapping_;
};

}  // namespace chromap

#endif  // ATAC_HOT_SPILL_H_
