// Serializer roundtrip for AtacSpillRecord (prefix-only and prefix + SAM pair).
#include "atac_spill_record.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

namespace {

void Fail(const char *msg) {
  std::cerr << "FAIL: " << msg << std::endl;
  std::exit(1);
}

void RoundtripPrefixOnly() {
  chromap::PairedEndMappingWithBarcode bed(
      42u, 0xdeadbeefcafeull, 10001u, 250u, 37u, 1u, 1u, 3u, 120u, 130u);
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
  if (b.HasBamPairSection()) {
    Fail("unexpected BAM section");
  }
}

void RoundtripWithBamPair() {
  chromap::PairedEndMappingWithBarcode bed(7u, 12345ull, 500u, 180u, 60u, 0u, 1u,
                                           1u, 90u, 90u);
  chromap::SAMMapping s1;
  s1.read_id_ = 7;
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
  if (b.sam1.read_name_ != "r1" || b.sam2.read_name_ != "r2") {
    Fail("SAM roundtrip names");
  }
}

void FileHeaderLayout() {
  chromap::AtacSpillFileHeader h = {};
  h.magic = chromap::kAtacSpillFileMagic;
  h.format_version = 1;
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

}  // namespace

int main() {
  RoundtripPrefixOnly();
  RoundtripWithBamPair();
  FileHeaderLayout();
  std::cerr << "PASS: AtacSpillRecord roundtrip\n";
  return 0;
}
