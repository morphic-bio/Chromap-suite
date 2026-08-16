#ifndef ATAC_KWAY_SPILL_H_
#define ATAC_KWAY_SPILL_H_

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "atac_spill_record.h"

namespace chromap {

// Ephemeral, worker-local ATAC overflow format used by the k-way merge.  It is
// deliberately independent of the durable ATACMS cross-process envelope.
#pragma pack(push, 1)
struct AtacKwaySpillFileHeaderV1 {
  char magic[8];
  uint16_t format_version;
  uint16_t fixed_header_bytes;
  uint16_t schema_mask;
  uint16_t flags;
  uint32_t record_codec_version;
  uint32_t reference_id;
  uint32_t endian_marker;
  uint32_t reserved;
};

struct AtacKwaySpillBlockHeaderV1 {
  uint32_t magic;
  uint32_t record_count;
  uint32_t payload_bytes;
  uint32_t reserved;
};

struct AtacKwaySpillRecordHeaderV1 {
  uint32_t magic;
  uint16_t codec_version;
  uint16_t fixed_header_bytes;
  uint64_t local_read_id;
  uint64_t packed_barcode_key;
  uint32_t fragment_start;
  uint32_t raw_barcode_n_mask;
  uint16_t fragment_length;
  uint16_t positive_alignment_length;
  uint16_t negative_alignment_length;
  uint8_t duplicate_count;
  uint8_t mapq_direction_unique;
  uint8_t row_flags;
  uint8_t barcode_quality_bytes;
  uint32_t bam_pair_bytes;
  uint16_t reserved;
};

struct AtacKwayBamPairHeaderV1 {
  uint16_t flags;
  uint16_t fixed_header_bytes;
  uint32_t template_length_value1;
  uint32_t template_length_value2;
  uint16_t query_name1_bytes;
  uint16_t query_name2_bytes;
  uint32_t mate1_bytes;
  uint32_t mate2_bytes;
};

struct AtacKwayBamMateHeaderV1 {
  uint32_t position;
  uint32_t edit_distance;
  uint16_t reference_id;
  uint16_t sam_flags;
  uint16_t cigar_count;
  uint16_t sequence_length;
  uint16_t md_event_count;
  uint16_t md_event_bytes;
  uint8_t mapq;
  uint8_t attributes;
  uint8_t quality_mode;
  uint8_t reserved;
};
#pragma pack(pop)

static_assert(sizeof(AtacKwaySpillFileHeaderV1) == 32,
              "ATAC k-way file header must be 32 bytes");
static_assert(sizeof(AtacKwaySpillBlockHeaderV1) == 16,
              "ATAC k-way block header must be 16 bytes");
static_assert(sizeof(AtacKwaySpillRecordHeaderV1) == 48,
              "ATAC k-way record header must be 48 bytes");
static_assert(sizeof(AtacKwayBamPairHeaderV1) == 24,
              "ATAC k-way BAM pair header must be 24 bytes");
static_assert(sizeof(AtacKwayBamMateHeaderV1) == 24,
              "ATAC k-way BAM mate header must be 24 bytes");

static constexpr char kAtacKwaySpillMagicV1[8] = {
    'A', 'T', 'K', 'W', 'S', '1', '\0', '\0'};
static constexpr uint16_t kAtacKwaySpillFormatVersion = 1;
static constexpr uint32_t kAtacKwaySpillRecordCodecVersion = 1;
static constexpr uint32_t kAtacKwaySpillEndianMarker = 0x01020304u;
static constexpr uint32_t kAtacKwaySpillBlockMagic = 0x3142574bu;  // KWB1
static constexpr uint32_t kAtacKwaySpillRecordMagic = 0x3152574bu; // KWR1
static constexpr size_t kAtacKwaySpillTargetBlockBytes = 4u * 1024u * 1024u;

static constexpr uint8_t kAtacKwayRowHasYHit = 1u << 0;
static constexpr uint16_t kAtacKwayPairCanonicalTemplateLength = 1u << 0;
static constexpr uint16_t kAtacKwayPairMate1PositiveTemplateLength = 1u << 1;
static constexpr uint16_t kAtacKwayPairExplicitMateCoordinates = 1u << 2;
static constexpr uint8_t kAtacKwayMateReverse = 1u << 0;
static constexpr uint8_t kAtacKwayMateAlt = 1u << 1;
static constexpr uint8_t kAtacKwayMateUnique = 1u << 2;
static constexpr uint8_t kAtacKwayQualityMissing = 0;
static constexpr uint8_t kAtacKwayQualityPhred33 = 1;

bool EncodeAtacKwaySpillRecord(const AtacSpillRecord &record,
                               uint16_t file_schema_mask,
                               std::vector<uint8_t> *encoded,
                               std::string *error);

bool DecodeAtacKwaySpillRecord(const void *bytes, size_t byte_count,
                               uint16_t file_schema_mask,
                               AtacSpillRecord *record,
                               std::string *error);

}  // namespace chromap

#endif  // ATAC_KWAY_SPILL_H_
