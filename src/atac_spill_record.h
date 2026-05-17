#ifndef ATAC_SPILL_RECORD_H_
#define ATAC_SPILL_RECORD_H_

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <string>
#include <utility>

#include "bed_mapping.h"
#include "sam_mapping.h"
#include "utils.h"

namespace chromap {

// Per-run optional sections for ATAC low-memory spill files (file header +
// length-prefixed payloads). Optional BAM pair section is present only when
// HAS_BAM_PAIR is set in the file header schema_mask.
enum AtacSpillSchemaMask : uint16_t {
  kAtacSpillSchemaNone = 0,
  kAtacSpillSchemaHasBamPair = 1u << 0,
  kAtacSpillSchemaHasYHit = 1u << 1,
  kAtacSpillSchemaIsBulk = 1u << 2,
};

// File-level header written once at the start of each overflow temp file when
// ATAC spill mode is enabled (see OverflowWriter::EnableAtacSpillFileHeader).
struct AtacSpillFileHeader {
  uint32_t magic;                 // kAtacSpillFileMagic
  uint16_t format_version;        // 1
  uint16_t schema_mask;           // AtacSpillSchemaMask bits
  uint32_t record_codec_version;  // 1
  uint32_t reserved0;             // 0 (alignment / future use)
};

static_assert(sizeof(AtacSpillFileHeader) == 16, "AtacSpillFileHeader size");

static constexpr uint32_t kAtacSpillFileMagic = 0x63417331u;  // 'cAs1'
static constexpr uint32_t kAtacSpillPayloadMagic = 0x61746131u;  // 'at1'
static constexpr uint32_t kAtacSpillRecordCodecVersion = 1u;
// LoadFromFile: use this sentinel so optional BAM sections are decoded from
// the per-payload mask (legacy temp files). Overflow merge passes the file
// header schema_mask instead so decoding follows the run-level contract.
static constexpr uint16_t kAtacSpillPayloadMaskAuthoritative = 0xFFFFu;

// Single in-memory / spill record for scATAC paired-end fragment path with
// optional BAM pair payload. Sort/dedup use only PairedEndMappingWithBarcode
// fields (via operator< / == / IsSamePosition delegating to the base).
struct AtacSpillRecord : public PairedEndMappingWithBarcode {
  // Row-level flags (e.g. HAS_Y_HIT). HAS_BAM_PAIR / IS_BULK are negotiated from
  // the overflow file header when loading spill payloads; see LoadFromFile.
  uint16_t prefix_flags_ = 0;
  SAMMapping sam1;
  SAMMapping sam2;

  AtacSpillRecord() = default;

  explicit AtacSpillRecord(PairedEndMappingWithBarcode bed_only)
      : PairedEndMappingWithBarcode(std::move(bed_only)), prefix_flags_(0) {}

  AtacSpillRecord(PairedEndMappingWithBarcode bed, SAMMapping a, SAMMapping b)
      : PairedEndMappingWithBarcode(std::move(bed)),
        prefix_flags_(kAtacSpillSchemaHasBamPair),
        sam1(std::move(a)),
        sam2(std::move(b)) {}

  bool operator<(const AtacSpillRecord &m) const {
    return static_cast<const PairedEndMappingWithBarcode &>(*this) <
           static_cast<const PairedEndMappingWithBarcode &>(m);
  }
  bool operator==(const AtacSpillRecord &m) const {
    return static_cast<const PairedEndMappingWithBarcode &>(*this) ==
           static_cast<const PairedEndMappingWithBarcode &>(m);
  }
  bool IsSamePosition(const AtacSpillRecord &m) const {
    return static_cast<const PairedEndMappingWithBarcode &>(*this)
        .IsSamePosition(static_cast<const PairedEndMappingWithBarcode &>(m));
  }

  bool HasBamPairSection() const {
    return (prefix_flags_ & kAtacSpillSchemaHasBamPair) != 0;
  }

  // Conservative in-RAM footprint for low-memory byte accounting.
  size_t HeapMemoryBytes() const {
    size_t n = sizeof(AtacSpillRecord);
    if (HasBamPairSection()) {
      n += sam1.SerializedSize() + sam2.SerializedSize();
    }
    return n;
  }

  size_t SerializedSize() const {
    size_t n = sizeof(uint32_t) + sizeof(uint16_t) + sizeof(uint16_t) +
               BedSerializedPayloadSize();  // magic + ver + mask + bed
    if (HasBamPairSection()) {
      n += sam1.SerializedSize() + sam2.SerializedSize();
    }
    return n;
  }

  size_t WriteToFile(FILE *fp) const {
    const size_t expect = SerializedSize();
    uint32_t magic = kAtacSpillPayloadMagic;
    uint16_t ver = static_cast<uint16_t>(kAtacSpillRecordCodecVersion);
    uint16_t mask = prefix_flags_;
    if (fwrite(&magic, sizeof(magic), 1, fp) != 1 ||
        fwrite(&ver, sizeof(ver), 1, fp) != 1 ||
        fwrite(&mask, sizeof(mask), 1, fp) != 1) {
      return 0;
    }
    if (fwrite(&read_id_, sizeof(read_id_), 1, fp) != 1 ||
        fwrite(&cell_barcode_, sizeof(cell_barcode_), 1, fp) != 1 ||
        fwrite(&fragment_start_position_, sizeof(fragment_start_position_), 1,
               fp) != 1 ||
        fwrite(&fragment_length_, sizeof(fragment_length_), 1, fp) != 1) {
      return 0;
    }
    uint8_t mapq8 = mapq_;
    uint8_t dir8 = direction_;
    uint8_t uniq8 = is_unique_;
    if (fwrite(&mapq8, 1, 1, fp) != 1 || fwrite(&dir8, 1, 1, fp) != 1 ||
        fwrite(&uniq8, 1, 1, fp) != 1 ||
        fwrite(&num_dups_, sizeof(num_dups_), 1, fp) != 1 ||
        fwrite(&positive_alignment_length_, sizeof(positive_alignment_length_),
               1, fp) != 1 ||
        fwrite(&negative_alignment_length_, sizeof(negative_alignment_length_),
               1, fp) != 1) {
      return 0;
    }
    if (HasBamPairSection()) {
      if (sam1.WriteToFile(fp) != sam1.SerializedSize() ||
          sam2.WriteToFile(fp) != sam2.SerializedSize()) {
        return 0;
      }
    }
    return expect;
  }

  // `overflow_file_schema_mask` is the AtacSpillFileHeader.schema_mask for
  // overflow payloads; pass kAtacSpillPayloadMaskAuthoritative for temp files
  // that have no spill header (payload mask alone selects optional sections).
  void LoadFromFile(FILE *fp, uint16_t overflow_file_schema_mask) {
    uint32_t magic = 0;
    uint16_t ver = 0;
    uint16_t payload_mask = 0;
    if (fread(&magic, sizeof(magic), 1, fp) != 1 ||
        magic != kAtacSpillPayloadMagic) {
      ExitWithMessage("Invalid AtacSpillRecord payload magic");
    }
    if (fread(&ver, sizeof(ver), 1, fp) != 1 ||
        fread(&payload_mask, sizeof(payload_mask), 1, fp) != 1 ||
        ver != 1u) {
      ExitWithMessage("Invalid AtacSpillRecord payload header");
    }
    const bool use_file_schema =
        overflow_file_schema_mask != kAtacSpillPayloadMaskAuthoritative;
    const uint16_t file_bam =
        static_cast<uint16_t>(overflow_file_schema_mask &
                              kAtacSpillSchemaHasBamPair);
    const uint16_t payload_bam =
        static_cast<uint16_t>(payload_mask & kAtacSpillSchemaHasBamPair);
    if (use_file_schema) {
      if (payload_bam != file_bam) {
        ExitWithMessage(
            "AtacSpillRecord payload schema_mask disagrees with overflow file "
            "header (HAS_BAM_PAIR mismatch; file may be corrupt)");
      }
    }
    const bool decode_bam =
        use_file_schema ? (file_bam != 0) : (payload_bam != 0);
    if (fread(&read_id_, sizeof(read_id_), 1, fp) != 1 ||
        fread(&cell_barcode_, sizeof(cell_barcode_), 1, fp) != 1 ||
        fread(&fragment_start_position_, sizeof(fragment_start_position_), 1,
              fp) != 1 ||
        fread(&fragment_length_, sizeof(fragment_length_), 1, fp) != 1) {
      ExitWithMessage("Truncated AtacSpillRecord bed payload");
    }
    uint8_t mapq8 = 0, dir8 = 0, uniq8 = 0;
    if (fread(&mapq8, 1, 1, fp) != 1 || fread(&dir8, 1, 1, fp) != 1 ||
        fread(&uniq8, 1, 1, fp) != 1 ||
        fread(&num_dups_, sizeof(num_dups_), 1, fp) != 1 ||
        fread(&positive_alignment_length_, sizeof(positive_alignment_length_), 1,
              fp) != 1 ||
        fread(&negative_alignment_length_, sizeof(negative_alignment_length_), 1,
              fp) != 1) {
      ExitWithMessage("Truncated AtacSpillRecord bed tail");
    }
    mapq_ = mapq8;
    direction_ = dir8;
    is_unique_ = uniq8;
    if (decode_bam) {
      sam1.LoadFromFile(fp);
      sam2.LoadFromFile(fp);
    } else {
      chromap::SAMMapping z1{};
      chromap::SAMMapping z2{};
      sam1 = std::move(z1);
      sam2 = std::move(z2);
    }
    // Canonical flags for consumers: run-level BAM/bulk from file header,
    // row-level Y from payload when present.
    prefix_flags_ = static_cast<uint16_t>(
        (payload_mask & kAtacSpillSchemaHasYHit) |
        (use_file_schema
             ? static_cast<uint16_t>(overflow_file_schema_mask &
                                     (kAtacSpillSchemaHasBamPair |
                                      kAtacSpillSchemaIsBulk))
             : static_cast<uint16_t>(payload_mask &
                                     (kAtacSpillSchemaHasBamPair |
                                      kAtacSpillSchemaIsBulk))));
  }

 private:
  static constexpr size_t BedSerializedPayloadSize() {
    return sizeof(uint32_t) + sizeof(uint64_t) + sizeof(uint32_t) +
           sizeof(uint16_t) + 3 * sizeof(uint8_t) + sizeof(uint8_t) +
           sizeof(uint16_t) + sizeof(uint16_t);
  }
};

inline uint16_t AtacSpillSchemaMaskForParameters(bool dual_bam_fragments,
                                                 bool is_bulk) {
  uint16_t m = 0;
  if (dual_bam_fragments) {
    m = static_cast<uint16_t>(m | kAtacSpillSchemaHasBamPair);
  }
  if (is_bulk) {
    m = static_cast<uint16_t>(m | kAtacSpillSchemaIsBulk);
  }
  return m;
}

}  // namespace chromap

#endif  // ATAC_SPILL_RECORD_H_
