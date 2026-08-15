#ifndef ATAC_SPILL_RECORD_H_
#define ATAC_SPILL_RECORD_H_

#include <algorithm>
#include <array>
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
  kAtacSpillSchemaHasRawBarcodeEvidence = 1u << 3,
};

// File-level header written once at the start of each overflow temp file when
// ATAC spill mode is enabled (see OverflowWriter::EnableAtacSpillFileHeader).
struct AtacSpillFileHeader {
  uint32_t magic;                 // kAtacSpillFileMagic
  uint16_t format_version;        // kAtacSpillFileFormatVersion
  uint16_t schema_mask;           // AtacSpillSchemaMask bits
  uint32_t record_codec_version;  // kAtacSpillRecordCodecVersion
  uint32_t reserved0;             // 0 (alignment / future use)
};

static_assert(sizeof(AtacSpillFileHeader) == 16, "AtacSpillFileHeader size");

static constexpr uint32_t kAtacSpillFileMagic = 0x63417331u;  // 'cAs1'
static constexpr uint32_t kAtacSpillPayloadMagic = 0x61746131u;  // 'at1'
static constexpr uint16_t kAtacSpillFileFormatVersion = 2u;
static constexpr uint32_t kAtacSpillRecordCodecVersion = 2u;
// LoadFromFile: use this sentinel so optional BAM sections are decoded from
// the per-payload mask (legacy temp files). Overflow merge passes the file
// header schema_mask instead so decoding follows the run-level contract.
static constexpr uint16_t kAtacSpillPayloadMaskAuthoritative = 0xFFFFu;

// Barcode qualities in the ATAC spill contract are bounded by the packed
// barcode width (32 bases). Keeping them inline avoids one heap allocation for
// every decoded mapping record during gather. This changes only the in-memory
// representation; the versioned on-disk payload remains unchanged.
class AtacBarcodeQuality {
 public:
  static constexpr size_t kCapacity = 32;

  AtacBarcodeQuality() = default;
  AtacBarcodeQuality(const AtacBarcodeQuality &) = default;
  AtacBarcodeQuality &operator=(const AtacBarcodeQuality &) = default;

  AtacBarcodeQuality &operator=(const std::string &quality) {
    assign(quality.data(), quality.size());
    return *this;
  }

  void assign(size_t count, char value) {
    if (count > kCapacity) {
      ExitWithMessage("ATAC barcode quality exceeds 32 bases");
    }
    size_ = static_cast<uint8_t>(count);
    std::fill(bytes_.begin(), bytes_.begin() + count, value);
  }

  void assign(const char *data, size_t count) {
    if (count > kCapacity || (count != 0 && data == nullptr)) {
      ExitWithMessage("Invalid ATAC barcode quality payload");
    }
    size_ = static_cast<uint8_t>(count);
    if (count != 0) {
      memcpy(bytes_.data(), data, count);
    }
  }

  void clear() { size_ = 0; }
  bool empty() const { return size_ == 0; }
  size_t size() const { return size_; }
  const char *data() const { return bytes_.data(); }
  char *data() { return bytes_.data(); }
  const char &operator[](size_t index) const { return bytes_[index]; }
  char &operator[](size_t index) { return bytes_[index]; }

  bool operator==(const char *text) const {
    return text != nullptr && strlen(text) == size_ &&
           memcmp(bytes_.data(), text, size_) == 0;
  }
  bool operator!=(const char *text) const { return !(*this == text); }

 private:
  std::array<char, kCapacity> bytes_{};
  uint8_t size_ = 0;
};

// Single in-memory / spill record for scATAC paired-end fragment path with
// optional BAM pair payload. Sort/dedup use only PairedEndMappingWithBarcode
// fields (via operator< / == / IsSamePosition delegating to the base).
struct AtacSpillRecord : public PairedEndMappingWithBarcode {
  // Row-level flags (e.g. HAS_Y_HIT). HAS_BAM_PAIR / IS_BULK are negotiated from
  // the overflow file header when loading spill payloads; see LoadFromFile.
  uint16_t prefix_flags_ = 0;
  uint32_t raw_barcode_n_mask_ = 0;
  AtacBarcodeQuality raw_barcode_qual_;
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

  bool HasYHit() const {
    return (prefix_flags_ & kAtacSpillSchemaHasYHit) != 0;
  }

  void SetYHit(bool value) {
    if (value) {
      prefix_flags_ = static_cast<uint16_t>(prefix_flags_ |
                                            kAtacSpillSchemaHasYHit);
    } else {
      prefix_flags_ = static_cast<uint16_t>(prefix_flags_ &
                                            ~kAtacSpillSchemaHasYHit);
    }
  }

  bool HasRawBarcodeEvidence() const {
    return (prefix_flags_ & kAtacSpillSchemaHasRawBarcodeEvidence) != 0;
  }

  void SetRawBarcodeEvidence(uint32_t n_mask, const std::string &quality) {
    prefix_flags_ = static_cast<uint16_t>(
        prefix_flags_ | kAtacSpillSchemaHasRawBarcodeEvidence);
    raw_barcode_n_mask_ = n_mask;
    raw_barcode_qual_ = quality;
  }

  // Conservative in-RAM footprint for low-memory byte accounting.
  size_t HeapMemoryBytes() const {
    size_t n = sizeof(AtacSpillRecord);
    if (HasBamPairSection()) {
      n += sam1.SerializedSize() + sam2.SerializedSize();
    }
    if (HasRawBarcodeEvidence()) {
      n += sizeof(raw_barcode_n_mask_) + sizeof(uint32_t) +
           raw_barcode_qual_.size();
    }
    return n;
  }

  size_t SerializedSize() const {
    size_t n = sizeof(uint32_t) + sizeof(uint16_t) + sizeof(uint16_t) +
               BedSerializedPayloadSize();  // magic + ver + mask + bed
    if (HasBamPairSection()) {
      n += sam1.SerializedSize() + sam2.SerializedSize();
    }
    if (HasRawBarcodeEvidence()) {
      n += sizeof(raw_barcode_n_mask_) + sizeof(uint32_t) +
           raw_barcode_qual_.size();
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
    if (HasRawBarcodeEvidence()) {
      const uint32_t quality_bytes =
          static_cast<uint32_t>(raw_barcode_qual_.size());
      if (fwrite(&raw_barcode_n_mask_, sizeof(raw_barcode_n_mask_), 1, fp) !=
              1 ||
          fwrite(&quality_bytes, sizeof(quality_bytes), 1, fp) != 1 ||
          (quality_bytes > 0 &&
           fwrite(raw_barcode_qual_.data(), 1, quality_bytes, fp) !=
               quality_bytes)) {
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
        ver != kAtacSpillRecordCodecVersion) {
      ExitWithMessage("Invalid AtacSpillRecord payload header");
    }
    const bool use_file_schema =
        overflow_file_schema_mask != kAtacSpillPayloadMaskAuthoritative;
    const uint16_t file_bam =
        static_cast<uint16_t>(overflow_file_schema_mask &
                              kAtacSpillSchemaHasBamPair);
    const uint16_t payload_bam =
        static_cast<uint16_t>(payload_mask & kAtacSpillSchemaHasBamPair);
    const uint16_t file_raw = static_cast<uint16_t>(
        overflow_file_schema_mask & kAtacSpillSchemaHasRawBarcodeEvidence);
    const uint16_t payload_raw = static_cast<uint16_t>(
        payload_mask & kAtacSpillSchemaHasRawBarcodeEvidence);
    if (use_file_schema) {
      if (payload_bam != file_bam || payload_raw != file_raw) {
        ExitWithMessage(
            "AtacSpillRecord payload schema_mask disagrees with overflow file "
            "header (optional-section mismatch; file may be corrupt)");
      }
    }
    const bool decode_bam =
        use_file_schema ? (file_bam != 0) : (payload_bam != 0);
    const bool decode_raw =
        use_file_schema ? (file_raw != 0) : (payload_raw != 0);
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
      auto clear_sam = [](chromap::SAMMapping *sam) {
        sam->read_id_ = 0;
        sam->read_name_.clear();
        sam->cell_barcode_ = 0;
        sam->num_dups_ = 0;
        sam->pos_ = 0;
        sam->rid_ = -1;
        sam->mpos_ = 0;
        sam->mrid_ = -1;
        sam->tlen_ = 0;
        sam->flag_ = 0;
        sam->is_rev_ = 0;
        sam->is_alt_ = 0;
        sam->is_unique_ = 0;
        sam->mapq_ = 0;
        sam->NM_ = 0;
        sam->n_cigar_ = 0;
        sam->cigar_.clear();
        sam->MD_.clear();
        sam->sequence_.clear();
        sam->sequence_qual_.clear();
      };
      clear_sam(&sam1);
      clear_sam(&sam2);
    }
    raw_barcode_n_mask_ = 0;
    raw_barcode_qual_.clear();
    if (decode_raw) {
      uint32_t quality_bytes = 0;
      if (fread(&raw_barcode_n_mask_, sizeof(raw_barcode_n_mask_), 1, fp) !=
              1 ||
          fread(&quality_bytes, sizeof(quality_bytes), 1, fp) != 1 ||
          quality_bytes > 32) {
        ExitWithMessage("Invalid AtacSpillRecord raw barcode evidence");
      }
      raw_barcode_qual_.assign(quality_bytes, '\0');
      if (quality_bytes > 0 &&
          fread(&raw_barcode_qual_[0], 1, quality_bytes, fp) !=
              quality_bytes) {
        ExitWithMessage("Truncated AtacSpillRecord raw barcode quality");
      }
    }
    // Canonical flags for consumers: run-level BAM/bulk from file header,
    // row-level Y from payload when present.
    prefix_flags_ = static_cast<uint16_t>(
        (payload_mask & kAtacSpillSchemaHasYHit) |
        (use_file_schema
             ? static_cast<uint16_t>(overflow_file_schema_mask &
                                     (kAtacSpillSchemaHasBamPair |
                                      kAtacSpillSchemaIsBulk |
                                      kAtacSpillSchemaHasRawBarcodeEvidence))
             : static_cast<uint16_t>(payload_mask &
                                     (kAtacSpillSchemaHasBamPair |
                                      kAtacSpillSchemaIsBulk |
                                      kAtacSpillSchemaHasRawBarcodeEvidence))));
  }

 private:
  static constexpr size_t BedSerializedPayloadSize() {
    return sizeof(uint64_t) + sizeof(uint64_t) + sizeof(uint32_t) +
           sizeof(uint16_t) + 3 * sizeof(uint8_t) + sizeof(uint8_t) +
           sizeof(uint16_t) + sizeof(uint16_t);
  }
};

inline uint16_t AtacSpillSchemaMaskForParameters(bool dual_bam_fragments,
                                                 bool is_bulk,
                                                 bool raw_barcode_evidence = false) {
  uint16_t m = 0;
  if (dual_bam_fragments) {
    m = static_cast<uint16_t>(m | kAtacSpillSchemaHasBamPair);
  }
  if (is_bulk) {
    m = static_cast<uint16_t>(m | kAtacSpillSchemaIsBulk);
  }
  if (raw_barcode_evidence) {
    m = static_cast<uint16_t>(m |
                              kAtacSpillSchemaHasRawBarcodeEvidence);
  }
  return m;
}

}  // namespace chromap

#endif  // ATAC_SPILL_RECORD_H_
