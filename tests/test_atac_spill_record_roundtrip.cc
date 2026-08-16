// Serializer roundtrip for AtacSpillRecord (prefix-only and prefix + SAM pair).
#include "atac_spill_record.h"
#include "atac_kway_spill.h"
#include "overflow_reader.h"
#include "overflow_writer.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <unistd.h>

namespace {

void Fail(const char *msg) {
  std::cerr << "FAIL: " << msg << std::endl;
  std::exit(1);
}

void RoundtripPrefixOnly() {
  const uint64_t read_id =
      static_cast<uint64_t>(std::numeric_limits<uint32_t>::max()) + 42u;
  chromap::PairedEndMappingWithBarcode bed(
      read_id, 0xdeadbeefcafeull, 10001u, 250u, 37u, 1u, 1u, 3u, 120u,
      130u);
  chromap::AtacSpillRecord a(bed);
  FILE *fp = std::tmpfile();
  if (!fp) {
    Fail("tmpfile");
  }
  if (a.WriteToFile(fp) != a.SerializedSize()) {
    Fail("WriteToFile size");
  }
  std::rewind(fp);
  chromap::AtacSpillRecord b;
  b.LoadFromFile(fp, chromap::kAtacSpillPayloadMaskAuthoritative);
  std::fclose(fp);
  if (!(static_cast<const chromap::PairedEndMappingWithBarcode &>(a) ==
          static_cast<const chromap::PairedEndMappingWithBarcode &>(b))) {
    Fail("prefix-only bed mismatch");
  }
  if (b.read_id_ != read_id) {
    Fail("64-bit prefix read id mismatch");
  }
  if (b.HasBamPairSection()) {
    Fail("unexpected BAM section");
  }
}

void RoundtripWithBamPair() {
  const uint64_t read_id =
      static_cast<uint64_t>(std::numeric_limits<uint32_t>::max()) + 7u;
  chromap::PairedEndMappingWithBarcode bed(read_id, 12345ull, 500u, 180u,
                                           60u, 0u, 1u, 1u, 90u, 90u);
  chromap::SAMMapping s1;
  s1.read_id_ = read_id;
  s1.read_name_ = "r1";
  s1.cell_barcode_ = 12345ull;
  s1.num_dups_ = 1;
  s1.pos_ = 100;
  s1.rid_ = 0;
  s1.mpos_ = 200;
  s1.mrid_ = 0;
  s1.tlen_ = 180;
  s1.flag_ = 99;
  s1.is_rev_ = 0;
  s1.is_alt_ = 0;
  s1.is_unique_ = 1;
  s1.mapq_ = 60;
  s1.NM_ = 0;
  s1.n_cigar_ = 0;
  s1.sequence_ = "ACGT";
  s1.sequence_qual_ = "IIII";
  chromap::SAMMapping s2 = s1;
  s2.read_name_ = "r2";
  chromap::AtacSpillRecord a(bed, std::move(s1), std::move(s2));
  a.SetRawBarcodeEvidence(/*n_mask=*/5, "I!II");
  FILE *fp = std::tmpfile();
  if (!fp) {
    Fail("tmpfile");
  }
  if (a.WriteToFile(fp) != a.SerializedSize()) {
    Fail("WriteToFile size dual");
  }
  std::rewind(fp);
  chromap::AtacSpillRecord b;
  b.LoadFromFile(fp, chromap::kAtacSpillPayloadMaskAuthoritative);
  std::fclose(fp);
  if (!b.HasBamPairSection()) {
    Fail("expected BAM section");
  }
  if (!b.HasRawBarcodeEvidence() || b.raw_barcode_n_mask_ != 5 ||
      b.raw_barcode_qual_ != "I!II") {
    Fail("raw barcode evidence roundtrip");
  }
  if (b.sam1.read_name_ != "r1" || b.sam2.read_name_ != "r2") {
    Fail("SAM roundtrip names");
  }
  if (b.read_id_ != read_id || b.sam1.read_id_ != read_id ||
      b.sam2.read_id_ != read_id) {
    Fail("64-bit BAM-pair read id mismatch");
  }
}

void RoundtripBulkFileSchema() {
  chromap::PairedEndMappingWithBarcode bed(17u, 0u, 700u, 160u, 60u, 1u,
                                           1u, 1u, 80u, 80u);
  chromap::SAMMapping s1;
  s1.read_id_ = 17u;
  s1.read_name_ = "bulk-r1";
  chromap::SAMMapping s2 = s1;
  s2.read_name_ = "bulk-r2";
  chromap::AtacSpillRecord a(bed, std::move(s1), std::move(s2));
  if (a.HasRawBarcodeEvidence()) {
    Fail("bulk record unexpectedly has raw barcode evidence");
  }
  FILE *fp = std::tmpfile();
  if (!fp) {
    Fail("tmpfile bulk");
  }
  if (a.WriteToFile(fp) != a.SerializedSize()) {
    Fail("WriteToFile size bulk");
  }
  std::rewind(fp);
  chromap::AtacSpillRecord b;
  b.LoadFromFile(fp,
                 static_cast<uint16_t>(chromap::kAtacSpillSchemaHasBamPair |
                                       chromap::kAtacSpillSchemaIsBulk));
  std::fclose(fp);
  if (!b.HasBamPairSection() || b.HasRawBarcodeEvidence() ||
      (b.prefix_flags_ & chromap::kAtacSpillSchemaIsBulk) == 0) {
    Fail("bulk file-schema roundtrip");
  }
}

void FileHeaderLayout() {
  chromap::AtacSpillFileHeader h = {};
  h.magic = chromap::kAtacSpillFileMagic;
  h.format_version = chromap::kAtacSpillFileFormatVersion;
  h.schema_mask = chromap::kAtacSpillSchemaHasBamPair;
  h.record_codec_version = chromap::kAtacSpillRecordCodecVersion;
  h.reserved0 = 0;
  if (sizeof(h) != 16u) {
    Fail("AtacSpillFileHeader size");
  }
  char buf[16];
  std::memcpy(buf, &h, sizeof(h));
  chromap::AtacSpillFileHeader h2;
  std::memcpy(&h2, buf, sizeof(h2));
  if (h2.magic != chromap::kAtacSpillFileMagic ||
      h2.schema_mask != chromap::kAtacSpillSchemaHasBamPair) {
    Fail("header memcpy");
  }
}

void KwayNormalizedCodecRoundtrip() {
  const uint64_t read_id =
      static_cast<uint64_t>(std::numeric_limits<uint32_t>::max()) + 91u;
  chromap::PairedEndMappingWithBarcode bed(
      read_id, 0x12345u, 101u, 179u, 55u, 1u, 1u, 1u, 75u, 76u);
  chromap::SAMMapping s1;
  s1.read_id_ = read_id;
  s1.read_name_ = "pair-name";
  s1.cell_barcode_ = bed.cell_barcode_;
  s1.num_dups_ = 1;
  s1.pos_ = 100;
  s1.rid_ = 2;
  s1.mpos_ = 204;
  s1.mrid_ = 2;
  s1.tlen_ = 179;
  s1.flag_ = 99;
  s1.is_rev_ = 0;
  s1.is_alt_ = 0;
  s1.is_unique_ = 1;
  s1.mapq_ = 55;
  s1.NM_ = 3;
  s1.cigar_.push_back(bam_cigar_gen(4, BAM_CMATCH));
  s1.cigar_.push_back(bam_cigar_gen(1, BAM_CINS));
  s1.cigar_.push_back(bam_cigar_gen(3, BAM_CMATCH));
  s1.n_cigar_ = static_cast<int>(s1.cigar_.size());
  s1.MD_ = "2A1^CG3";
  s1.sequence_ = "ACGTRYSW";
  s1.sequence_qual_ = "IHFEDCBA";
  chromap::SAMMapping s2 = s1;
  s2.flag_ = 147;
  s2.is_rev_ = 1;
  s2.pos_ = 204;
  s2.mpos_ = 100;
  s2.tlen_ = -179;
  s2.sequence_ = "TGCANNNN";
  s2.sequence_qual_ = "ABCDEFGH";
  chromap::AtacSpillRecord input(bed, std::move(s1), std::move(s2));
  input.SetRawBarcodeEvidence(3u, "I!HG");
  input.SetYHit(true);
  const uint16_t schema = static_cast<uint16_t>(
      chromap::kAtacSpillSchemaHasBamPair |
      chromap::kAtacSpillSchemaHasRawBarcodeEvidence);
  std::vector<uint8_t> encoded;
  std::string error;
  if (!chromap::EncodeAtacKwaySpillRecord(input, schema, &encoded, &error)) {
    Fail(error.c_str());
  }
  if (encoded.size() >= input.SerializedSize()) {
    Fail("normalized k-way payload did not shrink legacy SAM-pair payload");
  }
  chromap::AtacSpillRecord output;
  if (!chromap::DecodeAtacKwaySpillRecord(
          encoded.data(), encoded.size(), schema, &output, &error)) {
    Fail(error.c_str());
  }
  if (output.read_id_ != input.read_id_ ||
      output.cell_barcode_ != input.cell_barcode_ ||
      output.fragment_start_position_ != input.fragment_start_position_ ||
      output.fragment_length_ != input.fragment_length_ ||
      output.raw_barcode_n_mask_ != input.raw_barcode_n_mask_ ||
      output.raw_barcode_qual_ != "I!HG" || !output.HasYHit()) {
    Fail("normalized k-way decision row mismatch");
  }
  const chromap::SAMMapping *a[2] = {&input.sam1, &input.sam2};
  const chromap::SAMMapping *b[2] = {&output.sam1, &output.sam2};
  for (int i = 0; i < 2; ++i) {
    if (a[i]->read_id_ != b[i]->read_id_ ||
        a[i]->read_name_ != b[i]->read_name_ ||
        a[i]->cell_barcode_ != b[i]->cell_barcode_ ||
        a[i]->num_dups_ != b[i]->num_dups_ || a[i]->pos_ != b[i]->pos_ ||
        a[i]->rid_ != b[i]->rid_ || a[i]->mpos_ != b[i]->mpos_ ||
        a[i]->mrid_ != b[i]->mrid_ || a[i]->tlen_ != b[i]->tlen_ ||
        a[i]->flag_ != b[i]->flag_ || a[i]->is_rev_ != b[i]->is_rev_ ||
        a[i]->is_alt_ != b[i]->is_alt_ ||
        a[i]->is_unique_ != b[i]->is_unique_ ||
        a[i]->mapq_ != b[i]->mapq_ || a[i]->NM_ != b[i]->NM_ ||
        a[i]->cigar_ != b[i]->cigar_ || a[i]->MD_ != b[i]->MD_ ||
        a[i]->sequence_ != b[i]->sequence_ ||
        a[i]->sequence_qual_ != b[i]->sequence_qual_) {
      Fail("normalized k-way BAM mate mismatch");
    }
  }

  OverflowWriter writer("/tmp", "chromap_kway_unit");
  writer.EnableAtacSpillFileHeader(schema);
  writer.WriteAtac(2, input);
  const std::vector<std::string> paths = writer.Close();
  if (paths.size() != 1) {
    Fail("normalized k-way writer path count");
  }
  OverflowReader reader(paths[0]);
  if (!reader.IsValid() || !reader.FileHasAtacKwayHeader() ||
      reader.AtacKwayReferenceIdFromFileHeader() != 2 ||
      reader.AtacSpillSchemaFromFileHeader() != schema) {
    Fail("normalized k-way file header roundtrip");
  }
  uint32_t rid = 0;
  std::string payload;
  if (!reader.ReadNext(rid, payload) || rid != 2 ||
      payload.size() != encoded.size() ||
      memcmp(payload.data(), encoded.data(), encoded.size()) != 0 ||
      reader.ReadNext(rid, payload)) {
    Fail("normalized k-way block roundtrip");
  }
  unlink(paths[0].c_str());
}

}  // namespace

int main() {
  RoundtripPrefixOnly();
  RoundtripWithBamPair();
  RoundtripBulkFileSchema();
  FileHeaderLayout();
  KwayNormalizedCodecRoundtrip();
  std::cerr << "PASS: AtacSpillRecord roundtrip\n";
  return 0;
}
