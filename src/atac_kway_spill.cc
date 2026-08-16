#include "atac_kway_spill.h"

#include <algorithm>
#include <cctype>
#include <climits>
#include <cstring>
#include <limits>

namespace chromap {
namespace {

template <typename T>
void AppendValue(std::vector<uint8_t> *output, const T &value) {
  const size_t old_size = output->size();
  output->resize(old_size + sizeof(T));
  memcpy(output->data() + old_size, &value, sizeof(T));
}

void AppendBytes(std::vector<uint8_t> *output, const void *bytes,
                 size_t count) {
  const size_t old_size = output->size();
  output->resize(old_size + count);
  if (count != 0) {
    memcpy(output->data() + old_size, bytes, count);
  }
}

class Cursor {
 public:
  Cursor(const void *bytes, size_t count)
      : begin_(static_cast<const uint8_t *>(bytes)), cursor_(begin_),
        end_(begin_ + count) {}

  template <typename T>
  bool Read(T *value) {
    if (Remaining() < sizeof(T)) {
      return false;
    }
    memcpy(value, cursor_, sizeof(T));
    cursor_ += sizeof(T);
    return true;
  }

  bool ReadBytes(void *output, size_t count) {
    if (Remaining() < count) {
      return false;
    }
    if (count != 0) {
      memcpy(output, cursor_, count);
      cursor_ += count;
    }
    return true;
  }

  bool Skip(size_t count) {
    if (Remaining() < count) {
      return false;
    }
    cursor_ += count;
    return true;
  }

  const uint8_t *data() const { return cursor_; }
  size_t Remaining() const { return static_cast<size_t>(end_ - cursor_); }
  size_t Consumed() const { return static_cast<size_t>(cursor_ - begin_); }

 private:
  const uint8_t *begin_;
  const uint8_t *cursor_;
  const uint8_t *end_;
};

bool Fail(const std::string &message, std::string *error) {
  if (error != nullptr) {
    *error = message;
  }
  return false;
}

uint8_t IupacCode(char base) {
  static const char kIupac[] = "=ACMGRSVTWYHKDBN";
  const char upper = static_cast<char>(
      std::toupper(static_cast<unsigned char>(base)));
  for (uint8_t code = 0; code < 16; ++code) {
    if (kIupac[code] == upper) {
      return code;
    }
  }
  return std::numeric_limits<uint8_t>::max();
}

char IupacBase(uint8_t code) {
  static const char kIupac[] = "=ACMGRSVTWYHKDBN";
  return kIupac[code & 15u];
}

bool AppendPackedBases(const std::string &bases, std::vector<uint8_t> *out) {
  for (size_t i = 0; i < bases.size(); i += 2) {
    const uint8_t high = IupacCode(bases[i]);
    const uint8_t low = i + 1 < bases.size() ? IupacCode(bases[i + 1]) : 0;
    if (high > 15 || low > 15) {
      return false;
    }
    out->push_back(static_cast<uint8_t>((high << 4) | low));
  }
  return true;
}

bool ReadPackedBases(Cursor *cursor, size_t count, std::string *bases) {
  const size_t packed_count = (count + 1u) / 2u;
  if (cursor->Remaining() < packed_count) {
    return false;
  }
  bases->resize(count);
  for (size_t i = 0; i < packed_count; ++i) {
    uint8_t packed = 0;
    if (!cursor->Read(&packed)) {
      return false;
    }
    (*bases)[2u * i] = IupacBase(packed >> 4);
    if (2u * i + 1u < count) {
      (*bases)[2u * i + 1u] = IupacBase(packed);
    }
  }
  return true;
}

enum MdEventType : uint8_t {
  kMdMatch = 0,
  kMdMismatch = 1,
  kMdDeletion = 2,
};

bool EncodeMd(const std::string &md, std::vector<uint8_t> *events,
              uint16_t *event_count, std::string *error) {
  events->clear();
  uint32_t count = 0;
  size_t i = 0;
  while (i < md.size()) {
    if (std::isdigit(static_cast<unsigned char>(md[i]))) {
      uint64_t value = 0;
      do {
        value = value * 10u + static_cast<unsigned>(md[i] - '0');
        if (value > std::numeric_limits<uint32_t>::max()) {
          return Fail("ATAC k-way MD match length overflows uint32", error);
        }
        ++i;
      } while (i < md.size() &&
               std::isdigit(static_cast<unsigned char>(md[i])));
      events->push_back(kMdMatch);
      const uint32_t value32 = static_cast<uint32_t>(value);
      AppendValue(events, value32);
    } else if (md[i] == '^') {
      const size_t start = ++i;
      while (i < md.size() &&
             std::isalpha(static_cast<unsigned char>(md[i]))) {
        ++i;
      }
      const size_t length = i - start;
      if (length == 0 || length > std::numeric_limits<uint16_t>::max()) {
        return Fail("invalid ATAC k-way MD deletion", error);
      }
      events->push_back(kMdDeletion);
      const uint16_t length16 = static_cast<uint16_t>(length);
      AppendValue(events, length16);
      if (!AppendPackedBases(md.substr(start, length), events)) {
        return Fail("unsupported IUPAC base in ATAC k-way MD deletion",
                    error);
      }
    } else if (std::isalpha(static_cast<unsigned char>(md[i]))) {
      const uint8_t base = IupacCode(md[i]);
      if (base > 15) {
        return Fail("unsupported IUPAC base in ATAC k-way MD mismatch",
                    error);
      }
      events->push_back(kMdMismatch);
      events->push_back(base);
      ++i;
    } else {
      return Fail("unsupported character in ATAC k-way MD tag", error);
    }
    if (++count > std::numeric_limits<uint16_t>::max()) {
      return Fail("ATAC k-way MD event count overflows uint16", error);
    }
  }
  if (events->size() > std::numeric_limits<uint16_t>::max()) {
    return Fail("ATAC k-way MD event payload overflows uint16", error);
  }
  *event_count = static_cast<uint16_t>(count);
  return true;
}

void AppendUnsigned(std::string *output, uint32_t value) {
  char buffer[16];
  char *end = buffer + sizeof(buffer);
  char *begin = end;
  do {
    *--begin = static_cast<char>('0' + value % 10u);
    value /= 10u;
  } while (value != 0);
  output->append(begin, static_cast<size_t>(end - begin));
}

bool DecodeMd(Cursor *cursor, uint16_t event_count, uint16_t event_bytes,
              std::string *md) {
  if (cursor->Remaining() < event_bytes) {
    return false;
  }
  Cursor events(cursor->data(), event_bytes);
  md->clear();
  for (uint16_t i = 0; i < event_count; ++i) {
    uint8_t type = 0;
    if (!events.Read(&type)) {
      return false;
    }
    if (type == kMdMatch) {
      uint32_t length = 0;
      if (!events.Read(&length)) {
        return false;
      }
      AppendUnsigned(md, length);
    } else if (type == kMdMismatch) {
      uint8_t base = 0;
      if (!events.Read(&base) || base > 15) {
        return false;
      }
      md->push_back(IupacBase(base));
    } else if (type == kMdDeletion) {
      uint16_t length = 0;
      std::string bases;
      if (!events.Read(&length) || length == 0 ||
          !ReadPackedBases(&events, length, &bases)) {
        return false;
      }
      md->push_back('^');
      md->append(bases);
    } else {
      return false;
    }
  }
  if (events.Remaining() != 0) {
    return false;
  }
  return cursor->Skip(event_bytes);
}

bool EncodePosition(int64_t value, uint32_t *encoded) {
  if (value == -1) {
    *encoded = std::numeric_limits<uint32_t>::max();
    return true;
  }
  if (value < 0 || static_cast<uint64_t>(value) >=
                       std::numeric_limits<uint32_t>::max()) {
    return false;
  }
  *encoded = static_cast<uint32_t>(value);
  return true;
}

int64_t DecodePosition(uint32_t value) {
  return value == std::numeric_limits<uint32_t>::max()
             ? int64_t{-1}
             : static_cast<int64_t>(value);
}

bool EncodeReferenceId(int value, uint16_t *encoded) {
  if (value == -1) {
    *encoded = std::numeric_limits<uint16_t>::max();
    return true;
  }
  if (value < 0 || value >= std::numeric_limits<uint16_t>::max()) {
    return false;
  }
  *encoded = static_cast<uint16_t>(value);
  return true;
}

int DecodeReferenceId(uint16_t value) {
  return value == std::numeric_limits<uint16_t>::max()
             ? -1
             : static_cast<int>(value);
}

bool EncodeMate(const SAMMapping &mate, bool explicit_mate_coordinates,
                std::vector<uint8_t> *encoded, std::string *error) {
  if (mate.flag_ < 0 || mate.flag_ > std::numeric_limits<uint16_t>::max() ||
      mate.n_cigar_ < 0 ||
      mate.cigar_.size() != static_cast<size_t>(mate.n_cigar_) ||
      mate.cigar_.size() > std::numeric_limits<uint16_t>::max() ||
      mate.sequence_.size() > std::numeric_limits<uint16_t>::max() ||
      (!mate.sequence_qual_.empty() &&
       mate.sequence_qual_.size() != mate.sequence_.size())) {
    return Fail("ATAC k-way BAM mate fields exceed codec bounds", error);
  }
  AtacKwayBamMateHeaderV1 header = {};
  if (!EncodePosition(mate.pos_, &header.position) ||
      !EncodeReferenceId(mate.rid_, &header.reference_id)) {
    return Fail("ATAC k-way BAM mate coordinate exceeds codec bounds", error);
  }
  header.edit_distance = mate.NM_;
  header.sam_flags = static_cast<uint16_t>(mate.flag_);
  header.cigar_count = static_cast<uint16_t>(mate.cigar_.size());
  header.sequence_length = static_cast<uint16_t>(mate.sequence_.size());
  header.mapq = static_cast<uint8_t>(mate.mapq_);
  header.attributes = static_cast<uint8_t>(
      (mate.is_rev_ ? kAtacKwayMateReverse : 0) |
      (mate.is_alt_ ? kAtacKwayMateAlt : 0) |
      (mate.is_unique_ ? kAtacKwayMateUnique : 0));
  header.quality_mode = mate.sequence_qual_.empty()
                            ? kAtacKwayQualityMissing
                            : kAtacKwayQualityPhred33;
  std::vector<uint8_t> md_events;
  if (!EncodeMd(mate.MD_, &md_events, &header.md_event_count, error)) {
    return false;
  }
  header.md_event_bytes = static_cast<uint16_t>(md_events.size());
  AppendValue(encoded, header);
  if (explicit_mate_coordinates) {
    uint32_t mate_position = 0;
    uint16_t mate_reference = 0;
    if (!EncodePosition(mate.mpos_, &mate_position) ||
        !EncodeReferenceId(mate.mrid_, &mate_reference)) {
      return Fail("ATAC k-way explicit mate coordinate exceeds codec bounds",
                  error);
    }
    AppendValue(encoded, mate_position);
    AppendValue(encoded, mate_reference);
  }
  AppendBytes(encoded, mate.cigar_.data(),
              mate.cigar_.size() * sizeof(uint32_t));
  AppendBytes(encoded, md_events.data(), md_events.size());
  if (!AppendPackedBases(mate.sequence_, encoded)) {
    return Fail("unsupported IUPAC base in ATAC k-way sequence", error);
  }
  if (!mate.sequence_qual_.empty()) {
    for (char quality : mate.sequence_qual_) {
      const unsigned char value = static_cast<unsigned char>(quality);
      if (value < 33 || value > 126) {
        return Fail("ATAC k-way BAM quality is not Phred+33", error);
      }
      encoded->push_back(static_cast<uint8_t>(value - 33u));
    }
  }
  return true;
}

bool DecodeMate(Cursor *cursor, bool explicit_mate_coordinates,
                uint32_t expected_bytes, SAMMapping *mate) {
  if (cursor->Remaining() < expected_bytes) {
    return false;
  }
  Cursor encoded(cursor->data(), expected_bytes);
  AtacKwayBamMateHeaderV1 header = {};
  if (!encoded.Read(&header) ||
      header.quality_mode > kAtacKwayQualityPhred33 ||
      (header.attributes & ~(kAtacKwayMateReverse | kAtacKwayMateAlt |
                             kAtacKwayMateUnique)) != 0) {
    return false;
  }
  mate->pos_ = DecodePosition(header.position);
  mate->rid_ = DecodeReferenceId(header.reference_id);
  mate->flag_ = header.sam_flags;
  mate->NM_ = header.edit_distance;
  mate->mapq_ = header.mapq;
  mate->is_rev_ = (header.attributes & kAtacKwayMateReverse) != 0;
  mate->is_alt_ = (header.attributes & kAtacKwayMateAlt) != 0;
  mate->is_unique_ = (header.attributes & kAtacKwayMateUnique) != 0;
  if (explicit_mate_coordinates) {
    uint32_t mate_position = 0;
    uint16_t mate_reference = 0;
    if (!encoded.Read(&mate_position) || !encoded.Read(&mate_reference)) {
      return false;
    }
    mate->mpos_ = DecodePosition(mate_position);
    mate->mrid_ = DecodeReferenceId(mate_reference);
  }
  mate->n_cigar_ = header.cigar_count;
  mate->cigar_.resize(header.cigar_count);
  if (!encoded.ReadBytes(mate->cigar_.data(),
                         mate->cigar_.size() * sizeof(uint32_t)) ||
      !DecodeMd(&encoded, header.md_event_count, header.md_event_bytes,
                &mate->MD_) ||
      !ReadPackedBases(&encoded, header.sequence_length, &mate->sequence_)) {
    return false;
  }
  mate->sequence_qual_.clear();
  if (header.quality_mode == kAtacKwayQualityPhred33) {
    mate->sequence_qual_.resize(header.sequence_length);
    for (uint16_t i = 0; i < header.sequence_length; ++i) {
      uint8_t quality = 0;
      if (!encoded.Read(&quality) || quality > 93) {
        return false;
      }
      mate->sequence_qual_[i] = static_cast<char>(quality + 33u);
    }
  }
  return encoded.Remaining() == 0 && cursor->Skip(expected_bytes);
}

bool EncodeBamPair(const AtacSpillRecord &record,
                   std::vector<uint8_t> *encoded, std::string *error) {
  const SAMMapping &one = record.sam1;
  const SAMMapping &two = record.sam2;
  if (one.read_id_ != record.read_id_ || two.read_id_ != record.read_id_ ||
      one.cell_barcode_ != record.cell_barcode_ ||
      two.cell_barcode_ != record.cell_barcode_ ||
      one.num_dups_ != record.num_dups_ || two.num_dups_ != record.num_dups_ ||
      one.read_name_.size() > std::numeric_limits<uint16_t>::max() ||
      two.read_name_.size() > std::numeric_limits<uint16_t>::max()) {
    return Fail("ATAC k-way BAM pair disagrees with decision row", error);
  }
  AtacKwayBamPairHeaderV1 pair = {};
  pair.fixed_header_bytes = sizeof(pair);
  const int64_t tlen1 = one.tlen_;
  const int64_t tlen2 = two.tlen_;
  const uint64_t magnitude = tlen1 < 0 ? static_cast<uint64_t>(-tlen1)
                                       : static_cast<uint64_t>(tlen1);
  if (tlen1 == -tlen2 && magnitude <= std::numeric_limits<uint32_t>::max()) {
    pair.flags = static_cast<uint16_t>(
        pair.flags | kAtacKwayPairCanonicalTemplateLength |
        (tlen1 >= 0 ? kAtacKwayPairMate1PositiveTemplateLength : 0));
    pair.template_length_value1 = static_cast<uint32_t>(magnitude);
  } else {
    static_assert(sizeof(int) == sizeof(int32_t),
                  "ATAC k-way codec requires 32-bit int SAM fields");
    memcpy(&pair.template_length_value1, &one.tlen_, sizeof(uint32_t));
    memcpy(&pair.template_length_value2, &two.tlen_, sizeof(uint32_t));
  }
  const bool explicit_mates =
      one.mrid_ != two.rid_ || two.mrid_ != one.rid_ ||
      one.mpos_ != two.pos_ || two.mpos_ != one.pos_;
  if (explicit_mates) {
    pair.flags = static_cast<uint16_t>(
        pair.flags | kAtacKwayPairExplicitMateCoordinates);
  }
  pair.query_name1_bytes = static_cast<uint16_t>(one.read_name_.size());
  pair.query_name2_bytes = one.read_name_ == two.read_name_
                               ? 0
                               : static_cast<uint16_t>(two.read_name_.size());

  std::vector<uint8_t> mate1;
  std::vector<uint8_t> mate2;
  if (!EncodeMate(one, explicit_mates, &mate1, error) ||
      !EncodeMate(two, explicit_mates, &mate2, error) ||
      mate1.size() > std::numeric_limits<uint32_t>::max() ||
      mate2.size() > std::numeric_limits<uint32_t>::max()) {
    return false;
  }
  pair.mate1_bytes = static_cast<uint32_t>(mate1.size());
  pair.mate2_bytes = static_cast<uint32_t>(mate2.size());
  AppendValue(encoded, pair);
  AppendBytes(encoded, one.read_name_.data(), one.read_name_.size());
  if (pair.query_name2_bytes != 0) {
    AppendBytes(encoded, two.read_name_.data(), two.read_name_.size());
  }
  AppendBytes(encoded, mate1.data(), mate1.size());
  AppendBytes(encoded, mate2.data(), mate2.size());
  return true;
}

bool DecodeBamPair(Cursor *cursor, uint32_t byte_count,
                   AtacSpillRecord *record) {
  if (cursor->Remaining() < byte_count) {
    return false;
  }
  Cursor pair_bytes(cursor->data(), byte_count);
  AtacKwayBamPairHeaderV1 pair = {};
  if (!pair_bytes.Read(&pair) || pair.fixed_header_bytes != sizeof(pair) ||
      (pair.flags & ~(kAtacKwayPairCanonicalTemplateLength |
                      kAtacKwayPairMate1PositiveTemplateLength |
                      kAtacKwayPairExplicitMateCoordinates)) != 0 ||
      pair.query_name1_bytes == 0) {
    return false;
  }
  record->sam1.read_name_.resize(pair.query_name1_bytes);
  if (!pair_bytes.ReadBytes(&record->sam1.read_name_[0],
                            pair.query_name1_bytes)) {
    return false;
  }
  if (pair.query_name2_bytes == 0) {
    record->sam2.read_name_ = record->sam1.read_name_;
  } else {
    record->sam2.read_name_.resize(pair.query_name2_bytes);
    if (!pair_bytes.ReadBytes(&record->sam2.read_name_[0],
                              pair.query_name2_bytes)) {
      return false;
    }
  }
  const bool explicit_mates =
      (pair.flags & kAtacKwayPairExplicitMateCoordinates) != 0;
  if (!DecodeMate(&pair_bytes, explicit_mates, pair.mate1_bytes,
                  &record->sam1) ||
      !DecodeMate(&pair_bytes, explicit_mates, pair.mate2_bytes,
                  &record->sam2) ||
      pair_bytes.Remaining() != 0) {
    return false;
  }
  if (!explicit_mates) {
    record->sam1.mrid_ = record->sam2.rid_;
    record->sam1.mpos_ = record->sam2.pos_;
    record->sam2.mrid_ = record->sam1.rid_;
    record->sam2.mpos_ = record->sam1.pos_;
  }
  if ((pair.flags & kAtacKwayPairCanonicalTemplateLength) != 0) {
    if (pair.template_length_value1 >
        static_cast<uint32_t>(std::numeric_limits<int32_t>::max())) {
      return false;
    }
    const int32_t magnitude =
        static_cast<int32_t>(pair.template_length_value1);
    record->sam1.tlen_ =
        (pair.flags & kAtacKwayPairMate1PositiveTemplateLength) != 0
            ? magnitude
            : -magnitude;
    record->sam2.tlen_ = -record->sam1.tlen_;
  } else {
    memcpy(&record->sam1.tlen_, &pair.template_length_value1,
           sizeof(uint32_t));
    memcpy(&record->sam2.tlen_, &pair.template_length_value2,
           sizeof(uint32_t));
  }
  record->sam1.read_id_ = record->read_id_;
  record->sam2.read_id_ = record->read_id_;
  record->sam1.cell_barcode_ = record->cell_barcode_;
  record->sam2.cell_barcode_ = record->cell_barcode_;
  record->sam1.num_dups_ = record->num_dups_;
  record->sam2.num_dups_ = record->num_dups_;
  return cursor->Skip(byte_count);
}

void ClearSam(SAMMapping *sam) {
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
}

}  // namespace

bool EncodeAtacKwaySpillRecord(const AtacSpillRecord &record,
                               uint16_t file_schema_mask,
                               std::vector<uint8_t> *encoded,
                               std::string *error) {
  if (encoded == nullptr || record.fragment_length_ == 0 || record.mapq_ > 63) {
    return Fail("invalid ATAC k-way decision row", error);
  }
  const bool file_bam =
      (file_schema_mask & kAtacSpillSchemaHasBamPair) != 0;
  const bool file_raw =
      (file_schema_mask & kAtacSpillSchemaHasRawBarcodeEvidence) != 0;
  if (file_bam != record.HasBamPairSection() ||
      file_raw != record.HasRawBarcodeEvidence() ||
      record.raw_barcode_qual_.size() > 32) {
    return Fail("ATAC k-way record disagrees with file schema", error);
  }

  std::vector<uint8_t> pair;
  if (file_bam && !EncodeBamPair(record, &pair, error)) {
    return false;
  }
  if (pair.size() > std::numeric_limits<uint32_t>::max()) {
    return Fail("ATAC k-way BAM pair payload overflows uint32", error);
  }

  AtacKwaySpillRecordHeaderV1 header = {};
  header.magic = kAtacKwaySpillRecordMagic;
  header.codec_version = kAtacKwaySpillRecordCodecVersion;
  header.fixed_header_bytes = sizeof(header);
  header.local_read_id = record.read_id_;
  header.packed_barcode_key = record.cell_barcode_;
  header.fragment_start = record.fragment_start_position_;
  header.raw_barcode_n_mask = record.raw_barcode_n_mask_;
  header.fragment_length = record.fragment_length_;
  header.positive_alignment_length = record.positive_alignment_length_;
  header.negative_alignment_length = record.negative_alignment_length_;
  header.duplicate_count = record.num_dups_;
  header.mapq_direction_unique = static_cast<uint8_t>(
      record.mapq_ | (static_cast<uint8_t>(record.direction_) << 6) |
      (static_cast<uint8_t>(record.is_unique_) << 7));
  header.row_flags = record.HasYHit() ? kAtacKwayRowHasYHit : 0;
  header.barcode_quality_bytes =
      static_cast<uint8_t>(record.raw_barcode_qual_.size());
  header.bam_pair_bytes = static_cast<uint32_t>(pair.size());

  encoded->clear();
  encoded->reserve(sizeof(header) + header.barcode_quality_bytes + pair.size());
  AppendValue(encoded, header);
  for (size_t i = 0; i < record.raw_barcode_qual_.size(); ++i) {
    const unsigned char value =
        static_cast<unsigned char>(record.raw_barcode_qual_[i]);
    if (value < 33 || value > 126) {
      return Fail("ATAC k-way barcode quality is not Phred+33", error);
    }
    encoded->push_back(static_cast<uint8_t>(value - 33u));
  }
  AppendBytes(encoded, pair.data(), pair.size());
  return true;
}

bool DecodeAtacKwaySpillRecord(const void *bytes, size_t byte_count,
                               uint16_t file_schema_mask,
                               AtacSpillRecord *record,
                               std::string *error) {
  if (bytes == nullptr || record == nullptr) {
    return Fail("invalid ATAC k-way decode target", error);
  }
  Cursor cursor(bytes, byte_count);
  AtacKwaySpillRecordHeaderV1 header = {};
  if (!cursor.Read(&header) || header.magic != kAtacKwaySpillRecordMagic ||
      header.codec_version != kAtacKwaySpillRecordCodecVersion ||
      header.fixed_header_bytes != sizeof(header) ||
      header.fragment_length == 0 ||
      (header.row_flags & ~kAtacKwayRowHasYHit) != 0 ||
      header.barcode_quality_bytes > 32) {
    return Fail("invalid ATAC k-way record header", error);
  }
  const bool file_bam =
      (file_schema_mask & kAtacSpillSchemaHasBamPair) != 0;
  const bool file_raw =
      (file_schema_mask & kAtacSpillSchemaHasRawBarcodeEvidence) != 0;
  if ((file_bam && header.bam_pair_bytes == 0) ||
      (!file_bam && header.bam_pair_bytes != 0) ||
      (file_raw && header.barcode_quality_bytes == 0) ||
      (!file_raw && (header.barcode_quality_bytes != 0 ||
                     header.raw_barcode_n_mask != 0))) {
    return Fail("ATAC k-way record optional sections disagree with schema",
                error);
  }

  record->read_id_ = header.local_read_id;
  record->cell_barcode_ = header.packed_barcode_key;
  record->fragment_start_position_ = header.fragment_start;
  record->raw_barcode_n_mask_ = header.raw_barcode_n_mask;
  record->fragment_length_ = header.fragment_length;
  record->positive_alignment_length_ = header.positive_alignment_length;
  record->negative_alignment_length_ = header.negative_alignment_length;
  record->num_dups_ = header.duplicate_count;
  record->mapq_ = header.mapq_direction_unique & 0x3fu;
  record->direction_ = (header.mapq_direction_unique >> 6) & 1u;
  record->is_unique_ = (header.mapq_direction_unique >> 7) & 1u;
  record->raw_barcode_qual_.clear();
  if (file_raw) {
    record->raw_barcode_qual_.assign(header.barcode_quality_bytes, '\0');
    for (uint8_t i = 0; i < header.barcode_quality_bytes; ++i) {
      uint8_t quality = 0;
      if (!cursor.Read(&quality) || quality > 93) {
        return Fail("invalid ATAC k-way barcode quality", error);
      }
      record->raw_barcode_qual_[i] = static_cast<char>(quality + 33u);
    }
  }
  ClearSam(&record->sam1);
  ClearSam(&record->sam2);
  if (file_bam && !DecodeBamPair(&cursor, header.bam_pair_bytes, record)) {
    return Fail("invalid ATAC k-way BAM pair payload", error);
  }
  if (cursor.Remaining() != 0) {
    return Fail("trailing bytes in ATAC k-way record", error);
  }
  record->prefix_flags_ = static_cast<uint16_t>(
      file_schema_mask & (kAtacSpillSchemaHasBamPair |
                          kAtacSpillSchemaIsBulk |
                          kAtacSpillSchemaHasRawBarcodeEvidence));
  if ((header.row_flags & kAtacKwayRowHasYHit) != 0) {
    record->prefix_flags_ = static_cast<uint16_t>(
        record->prefix_flags_ | kAtacSpillSchemaHasYHit);
  }
  return true;
}

}  // namespace chromap
