#ifndef ATAC_DUAL_MAPPING_H_
#define ATAC_DUAL_MAPPING_H_

#include <cstdio>
#include <cstdint>
#include <string>
#include <utility>

#include "bed_mapping.h"
#include "sam_mapping.h"
#include "utils.h"

namespace chromap {

// Paired-end scATAC: fragment coordinates (Tn5/dedup/collapse identical to
// PairedEndMappingWithBarcode) plus two SAM rows per pair for BAM in one pass.
struct PairedEndAtacDualMapping : public PairedEndMappingWithBarcode {
  SAMMapping sam1;
  SAMMapping sam2;

  PairedEndAtacDualMapping() = default;
  PairedEndAtacDualMapping(PairedEndMappingWithBarcode bed, SAMMapping a,
                           SAMMapping b)
      : PairedEndMappingWithBarcode(bed),
        sam1(std::move(a)),
        sam2(std::move(b)) {}

  bool operator<(const PairedEndAtacDualMapping &m) const {
    return static_cast<const PairedEndMappingWithBarcode &>(*this) <
           static_cast<const PairedEndMappingWithBarcode &>(m);
  }
  bool operator==(const PairedEndAtacDualMapping &m) const {
    return static_cast<const PairedEndMappingWithBarcode &>(*this) ==
           static_cast<const PairedEndMappingWithBarcode &>(m);
  }
  bool IsSamePosition(const PairedEndAtacDualMapping &m) const {
    return static_cast<const PairedEndMappingWithBarcode &>(*this)
        .IsSamePosition(static_cast<const PairedEndMappingWithBarcode &>(m));
  }

  static constexpr uint32_t kDualSerdeMagic = 0x63706164u;  // 'cpad'

  size_t SerializedSize() const {
    return sizeof(uint32_t) + BedSerializedPayloadSize() +
           sam1.SerializedSize() + sam2.SerializedSize();
  }

  size_t WriteToFile(FILE *fp) const {
    const size_t expect = SerializedSize();
    uint32_t magic = kDualSerdeMagic;
    if (fwrite(&magic, sizeof(magic), 1, fp) != 1) {
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
    if (sam1.WriteToFile(fp) != sam1.SerializedSize() ||
        sam2.WriteToFile(fp) != sam2.SerializedSize()) {
      return 0;
    }
    return expect;
  }

  void LoadFromFile(FILE *fp) {
    uint32_t magic = 0;
    if (fread(&magic, sizeof(magic), 1, fp) != 1 ||
        magic != kDualSerdeMagic) {
      ExitWithMessage("Invalid PairedEndAtacDualMapping serde header");
    }
    if (fread(&read_id_, sizeof(read_id_), 1, fp) != 1 ||
        fread(&cell_barcode_, sizeof(cell_barcode_), 1, fp) != 1 ||
        fread(&fragment_start_position_, sizeof(fragment_start_position_), 1,
              fp) != 1 ||
        fread(&fragment_length_, sizeof(fragment_length_), 1, fp) != 1) {
      ExitWithMessage("Truncated PairedEndAtacDualMapping bed payload");
    }
    uint8_t mapq8 = 0, dir8 = 0, uniq8 = 0;
    if (fread(&mapq8, 1, 1, fp) != 1 || fread(&dir8, 1, 1, fp) != 1 ||
        fread(&uniq8, 1, 1, fp) != 1 ||
        fread(&num_dups_, sizeof(num_dups_), 1, fp) != 1 ||
        fread(&positive_alignment_length_, sizeof(positive_alignment_length_), 1,
              fp) != 1 ||
        fread(&negative_alignment_length_, sizeof(negative_alignment_length_), 1,
              fp) != 1) {
      ExitWithMessage("Truncated PairedEndAtacDualMapping bed payload");
    }
    mapq_ = mapq8;
    direction_ = dir8;
    is_unique_ = uniq8;
    sam1.LoadFromFile(fp);
    sam2.LoadFromFile(fp);
  }

 private:
  static constexpr size_t BedSerializedPayloadSize() {
    return sizeof(uint32_t) + sizeof(uint64_t) + sizeof(uint32_t) +
           sizeof(uint16_t) + 3 * sizeof(uint8_t) + sizeof(uint8_t) +
           sizeof(uint16_t) + sizeof(uint16_t);
  }
};

}  // namespace chromap

#endif  // ATAC_DUAL_MAPPING_H_
