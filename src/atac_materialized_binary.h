#ifndef ATAC_MATERIALIZED_BINARY_H_
#define ATAC_MATERIALIZED_BINARY_H_

#include <cstdio>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "atac_mergeable_spill.h"
#include "atac_spill_record.h"
#include "mapping_parameters.h"

namespace chromap {

// Canonical post-correction/post-dedup ATAC fragment container.  The gather
// hot path writes only these fixed-width binary records.  Text formats such as
// BED are terminal exports and are never used as intermediate merge state.
#pragma pack(push, 1)
struct AtacMaterializedBinaryHeaderV1 {
  char magic[8];
  uint16_t format_version;
  uint16_t fixed_header_bytes;
  uint32_t record_bytes;
  uint32_t flags;
  uint32_t barcode_length;
  uint32_t num_reference_sequences;
  uint32_t shard_count;
  uint32_t sample_id_bytes;
  uint32_t input_id_bytes;
  uint64_t metadata_bytes;
  uint64_t record_count;
  uint64_t block_count;
  uint64_t barcode_dictionary_count;
  uint64_t barcode_dictionary_offset;
  uint64_t directory_offset;
  uint64_t data_offset;
  uint32_t target_records_per_block;
  uint32_t endian_marker;
};

struct AtacMaterializedBinaryRecordV1 {
  uint32_t start;
  // Packed barcode key when the header has no barcode-dictionary flag;
  // otherwise a dense id into the footer dictionary.
  uint32_t barcode_value;
  uint16_t rid;
  uint16_t fragment_length;
  uint8_t duplicate_count;
  uint8_t mapq;
  uint8_t flags;
  uint8_t reserved;
};

struct AtacMaterializedBinaryBlockV1 {
  uint64_t offset;
  uint64_t record_count;
  uint32_t first_rid;
  uint32_t last_rid;
  uint32_t first_start;
  uint32_t last_start;
  uint32_t first_end;
  uint32_t last_end;
};
#pragma pack(pop)

static_assert(sizeof(AtacMaterializedBinaryHeaderV1) == 104,
              "ATAC materialized binary header must be 104 bytes");
static_assert(sizeof(AtacMaterializedBinaryRecordV1) == 16,
              "ATAC materialized binary record must be 16 bytes");
static_assert(sizeof(AtacMaterializedBinaryBlockV1) == 40,
              "ATAC materialized binary block descriptor must be 40 bytes");

static constexpr char kAtacMaterializedBinaryMagicV1[8] = {
    'A', 'T', 'M', 'B', 'L', 'K', '1', '\0'};
static constexpr uint16_t kAtacMaterializedBinaryFormatVersion = 1;
static constexpr uint32_t kAtacMaterializedBinaryTargetRecordsPerBlock =
    262144;
static constexpr uint32_t kAtacMaterializedBinaryEndianMarker = 0x01020304u;
static constexpr uint32_t kAtacMaterializedBinaryIsBulk = 1u << 0;
static constexpr uint32_t kAtacMaterializedBinaryHasBarcodeDictionary =
    1u << 1;
static constexpr uint8_t kAtacMaterializedBinaryPositiveStrand = 1u << 0;

struct AtacMaterializedBinaryMetadata {
  bool is_bulk = false;
  bool use_barcode_dictionary = false;
  uint32_t barcode_length = 0;
  uint32_t shard_count = 0;
  std::string sample_id;
  std::string input_id;
  std::vector<AtacMergeableSpillReference> references;
};

class AtacMaterializedBinaryWriter {
 public:
  AtacMaterializedBinaryWriter() = default;
  ~AtacMaterializedBinaryWriter();

  bool Open(const std::string &path,
            const AtacMaterializedBinaryMetadata &metadata,
            std::string *error);
  bool Append(uint32_t rid, const AtacSpillRecord &mapping,
              uint64_t duplicate_count, std::string *error);
  // Appends already encoded direct-barcode records produced by parallel hot
  // partitions.
  bool AppendEncodedRecords(const AtacMaterializedBinaryRecordV1 *records,
                            size_t count, std::string *error);
  // Appends hot-partition records whose barcode_value still contains the raw
  // <=16-base packed key. In dictionary mode, dense ids are assigned here in
  // final reference/coordinate order, keeping translation out of worker tasks
  // while producing the same deterministic dictionary as Append().
  bool AppendEncodedRecordsWithRawBarcodes(
      const AtacMaterializedBinaryRecordV1 *records, size_t count,
      std::string *error);
  // Appends compact rows plus their full 64-bit packed barcode keys. This is
  // used when parallel hot partitions carry 17-32-base barcodes that cannot
  // be narrowed before the ordered dictionary is assembled.
  bool AppendEncodedRecordsWithRawBarcodeKeys(
      const AtacMaterializedBinaryRecordV1 *records,
      const uint64_t *barcode_keys, size_t count, std::string *error);
  bool Finalize(std::string *error);
  uint64_t record_count() const { return record_count_; }

 private:
  bool Fail(const std::string &message, std::string *error);
  bool FlushBlock(std::string *error);
  bool AppendEncodedRecordWithRawBarcodeKey(
      AtacMaterializedBinaryRecordV1 record, uint64_t barcode_key,
      std::string *error);
  bool AppendEncodedRecord(const AtacMaterializedBinaryRecordV1 &record,
                           bool allow_dictionary, std::string *error);

  FILE *file_ = nullptr;
  std::string output_path_;
  std::string temporary_path_;
  AtacMaterializedBinaryHeaderV1 header_{};
  AtacMaterializedBinaryMetadata metadata_;
  std::vector<AtacMaterializedBinaryRecordV1> block_records_;
  std::vector<AtacMaterializedBinaryBlockV1> blocks_;
  std::unordered_map<uint64_t, uint32_t> barcode_ids_;
  std::vector<uint64_t> barcode_dictionary_;
  uint64_t record_count_ = 0;
  bool finalized_ = false;
  bool have_previous_ = false;
  uint32_t previous_rid_ = 0;
  uint32_t previous_start_ = 0;
};

// Parallel terminal text export. The input and all materialization
// intermediates remain binary; text exists only in bounded in-memory block
// buffers immediately before disjoint pwrite calls to the final BED.
bool ExportAtacMaterializedBinaryToBed(
    const std::string &binary_path, const std::string &bed_path,
    const MappingParameters &parameters, std::string *error);

}  // namespace chromap

#endif  // ATAC_MATERIALIZED_BINARY_H_
