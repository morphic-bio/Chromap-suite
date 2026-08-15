#include <cstdio>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "atac_mergeable_spill.h"
#include "atac_spill_materializer.h"
#include "summary_metadata.h"

namespace {

void Check(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL: " << message << "\n";
    std::exit(1);
  }
}

chromap::SAMMapping MakeSam(uint64_t read_id, bool read1) {
  chromap::SAMMapping sam;
  sam.read_id_ = read_id;
  sam.read_name_ = "read" + std::to_string(read_id);
  sam.cell_barcode_ = 1;
  sam.num_dups_ = 1;
  sam.pos_ = 100;
  sam.rid_ = 0;
  sam.mpos_ = 140;
  sam.mrid_ = 0;
  sam.tlen_ = 80;
  sam.flag_ = BAM_FPAIRED | BAM_FPROPER_PAIR |
              (read1 ? BAM_FREAD1 : BAM_FREAD2);
  sam.is_rev_ = read1 ? 0 : 1;
  sam.is_alt_ = 0;
  sam.is_unique_ = 1;
  sam.mapq_ = 60;
  sam.NM_ = 0;
  sam.n_cigar_ = 1;
  sam.cigar_.push_back(bam_cigar_gen(20, BAM_CMATCH));
  sam.MD_ = "20";
  sam.sequence_ = "ACGTACGTACGTACGTACGT";
  sam.sequence_qual_ = "IIIIIIIIIIIIIIIIIIII";
  return sam;
}

chromap::AtacSpillRecord MakeRecord(uint64_t local_read_id, uint32_t start,
                                    uint8_t mapq, uint64_t barcode,
                                    bool y_hit = false) {
  chromap::PairedEndMappingWithBarcode fragment(
      local_read_id, barcode, start, 80, mapq, 1, 1, 1, 20, 20);
  chromap::SAMMapping sam1 = MakeSam(local_read_id, true);
  chromap::SAMMapping sam2 = MakeSam(local_read_id, false);
  sam1.cell_barcode_ = barcode;
  sam2.cell_barcode_ = barcode;
  sam1.pos_ = start;
  sam1.mpos_ = start + 40;
  sam2.pos_ = start + 40;
  sam2.mpos_ = start;
  chromap::AtacSpillRecord record(fragment, std::move(sam1),
                                  std::move(sam2));
  record.SetRawBarcodeEvidence(/*n_mask=*/0, "IIII");
  record.SetYHit(y_hit);
  return record;
}

chromap::AtacMergeableSpillMetadata Metadata(uint32_t ordinal,
                                             uint64_t first) {
  chromap::AtacMergeableSpillMetadata metadata;
  metadata.schema_mask = static_cast<uint16_t>(
      chromap::kAtacSpillSchemaHasBamPair |
      chromap::kAtacSpillSchemaHasRawBarcodeEvidence);
  metadata.flags = static_cast<uint16_t>(
      chromap::kAtacMergeableRemovePcrDuplicates |
      chromap::kAtacMergeableHasSummaryEvidence);
  metadata.shard_ordinal = ordinal;
  metadata.shard_count = 2;
  metadata.first_global_read_ordinal = first;
  metadata.input_record_count = 2;
  metadata.summary_evidence_count = 2;
  metadata.barcode_length = 4;
  metadata.mapq_threshold = 30;
  metadata.barcode_whitelist_fingerprint = 0x123456789abcdef0ULL;
  metadata.barcode_correction_error_threshold = 1;
  metadata.barcode_correction_probability_threshold = 0.90;
  metadata.frip_est_coefficients = {0.0, 0.0, 0.0, 0.0, 0.0};
  metadata.barcode_abundance_entries.push_back(
      chromap::AtacBarcodeAbundanceEntry{0, ordinal == 0 ? 1u : 100u});
  metadata.barcode_abundance_entries.push_back(
      chromap::AtacBarcodeAbundanceEntry{5, ordinal == 0 ? 10u : 0u});
  metadata.local_num_sample_barcodes = ordinal == 0 ? 11 : 100;
  metadata.sample_id = "sample-A";
  metadata.input_id = "atac-triplet-A";
  chromap::AtacMergeableSpillReference reference;
  reference.name = "chr1";
  reference.length = 1000000;
  metadata.references.push_back(reference);
  return metadata;
}

void WriteShard(const std::string &path,
                const chromap::AtacMergeableSpillMetadata &metadata,
                const std::vector<chromap::AtacSpillRecord> &records) {
  chromap::AtacMergeableSpillWriter writer;
  std::string error;
  Check(writer.Open(path, metadata, &error), error);
  for (uint64_t i = 0; i < metadata.input_record_count; ++i) {
    chromap::AtacSummaryEvidence evidence;
    evidence.raw_barcode_key =
        i < records.size() ? records[i].cell_barcode_ : uint64_t{0};
    evidence.raw_barcode_qual = "IIII";
    Check(writer.AppendSummaryEvidence(evidence, &error), error);
  }
  for (const auto &record : records) {
    Check(writer.Append(0, record, &error), error);
  }
  Check(writer.Finalize(&error), error);
}

std::vector<std::vector<std::string>> ReadBed(const std::string &path) {
  std::ifstream input(path.c_str());
  Check(input.is_open(), "cannot read materialized BED");
  std::vector<std::vector<std::string>> rows;
  std::string line;
  while (std::getline(input, line)) {
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, '\t')) {
      fields.push_back(field);
    }
    rows.push_back(std::move(fields));
  }
  return rows;
}

}  // namespace

int main(int argc, char **argv) {
  Check(argc == 2, "expected HDD artifact directory argument");
  const std::string root = argv[1];
  const std::string shard0 = root + "/shard0.atacms";
  const std::string shard1 = root + "/shard1.atacms";
  const std::string duplicate0 = root + "/duplicate0.atacms";
  const std::string mismatched_model = root + "/mismatched_model.atacms";
  const std::string output = root + "/materialized.bed";
  const std::string bam_output = root + "/materialized.bam";
  const std::string fragments_output = root + "/materialized.fragments.tsv";
  const std::string evidence_output = root + "/materialized.aev1";
  const std::string noy_bam_output = root + "/materialized.noY.bam";
  const std::string y_bam_output = root + "/materialized.Y.bam";
  const std::string cram_output = root + "/materialized.cram";
  const std::string cram_fragments_output =
      root + "/materialized.cram.fragments.tsv";
  const std::string reference_fasta = root + "/reference.fa";
  const std::string summary_output = root + "/materialized.summary.csv";
  const std::string bulk_dedup_shard0 = root + "/bulk_dedup_shard0.atacms";
  const std::string bulk_dedup_shard1 = root + "/bulk_dedup_shard1.atacms";
  const std::string bulk_dedup_output = root + "/bulk_dedup.bed";
  const std::string allocation_shard0 = root + "/allocation_shard0.atacms";
  const std::string allocation_shard1 = root + "/allocation_shard1.atacms";
  const std::string allocation_output = root + "/allocation.bed";
  const std::string late_shard0 = root + "/late_shard0.atacms";
  const std::string late_shard1 = root + "/late_shard1.atacms";
  const std::string late_output = root + "/late.bed";
  const std::string wide_summary_output = root + "/wide.summary.csv";
  const uint64_t wide_prefix =
      static_cast<uint64_t>(std::numeric_limits<uint32_t>::max()) + 100u;

  WriteShard(shard0, Metadata(0, wide_prefix),
             {MakeRecord(0, 100, 20, 0), MakeRecord(1, 200, 50, 1)});
  WriteShard(shard1, Metadata(1, wide_prefix + 2u),
             {MakeRecord(0, 100, 60, 0),
              MakeRecord(1, 300, 50, 5, true)});

  chromap::AtacMergeableSpillReader reader;
  std::string error;
  Check(reader.Open(shard1, &error), error);
  Check(reader.metadata().shard_ordinal == 1 &&
            reader.metadata().first_global_read_ordinal == wide_prefix + 2u &&
            reader.metadata().barcode_whitelist_fingerprint ==
                0x123456789abcdef0ULL &&
            reader.metadata().barcode_abundance_entries.size() == 2 &&
            reader.metadata().local_num_sample_barcodes == 100 &&
            reader.expected_record_count() == 2,
        "shard header did not round-trip");
  uint32_t rid = 0;
  chromap::AtacSpillRecord first;
  bool eof = false;
  for (uint64_t i = 0; i < reader.metadata().summary_evidence_count; ++i) {
    chromap::AtacSummaryEvidence evidence;
    Check(reader.ReadNextSummaryEvidence(&evidence, &eof, &error) && !eof,
          error);
  }
  Check(reader.ReadNext(&rid, &first, &eof, &error) && !eof, error);
  Check(chromap::GlobalizeAtacSpillReadId(
            &first, reader.metadata().first_global_read_ordinal,
            reader.metadata().input_record_count, &error),
        error);
  Check(first.read_id_ == wide_prefix + 2u &&
            first.sam1.read_id_ == wide_prefix + 2u &&
            first.sam2.read_id_ == wide_prefix + 2u,
        "global read ordinal did not propagate through BAM pair");

  chromap::AtacMergeableSpillReader bed_reader;
  Check(bed_reader.Open(shard1, &error), error);
  Check(bed_reader.SkipRemainingSummaryEvidence(&error), error);
  chromap::AtacSpillRecord compact_first;
  Check(bed_reader.ReadNextBed(&rid, &compact_first, &eof, &error) && !eof,
        error);
  Check(!compact_first.HasBamPairSection() &&
            compact_first.read_id_ == 0 &&
            compact_first.fragment_start_position_ == 100 &&
            compact_first.fragment_length_ == 80 &&
            compact_first.mapq_ == 60 &&
            compact_first.raw_barcode_n_mask_ == 0 &&
            compact_first.raw_barcode_qual_ == "IIII",
        "compact BED spill decoder did not preserve required fields");
  Check(chromap::GlobalizeAtacSpillReadId(
            &compact_first, bed_reader.metadata().first_global_read_ordinal,
            bed_reader.metadata().input_record_count, &error),
        error);
  Check(compact_first.read_id_ == wide_prefix + 2u,
        "compact BED spill decoder did not support 64-bit globalization");

  chromap::MappingParameters parameters;
  parameters.mapping_output_format = chromap::MAPPINGFORMAT_BED;
  parameters.mapping_output_file_path = output;
  parameters.summary_metadata_file_path = summary_output;
  const auto result = chromap::MaterializeAtacSpillRecords(
      {shard1, shard0}, parameters);
  Check(result.ok, result.message);
  Check(result.shard_count == 2 && result.input_record_count == 4 &&
            result.spill_record_count == 4 &&
            result.corrected_barcode_record_count == 1 &&
            result.rejected_barcode_record_count == 0 &&
            result.output_fragment_count == 3,
        "materialization counters are incorrect");

  const auto rows = ReadBed(output);
  Check(rows.size() == 3, "materialized BED row count is incorrect");
  Check(rows[0].size() == 5 && rows[0][0] == "chr1" &&
            rows[0][1] == "100" && rows[0][2] == "180" &&
            rows[0][4] == "2",
        "cross-shard duplicate was not globally retained/counted");
  Check(rows[1][1] == "200" && rows[1][4] == "1" &&
            rows[1][3] == "AAAA" &&
            rows[2][1] == "300" && rows[2][4] == "1",
        "post-spill correction or nonduplicate records are incorrect");
  std::ifstream summary_stream(summary_output.c_str());
  Check(summary_stream.is_open(), "materialized summary was not written");
  std::string summary_text((std::istreambuf_iterator<char>(summary_stream)),
                           std::istreambuf_iterator<char>());
  Check(summary_text.find("AAAA,3,1,0,0,0") != std::string::npos &&
            summary_text.find("AACC,1,0,0,0,0") != std::string::npos &&
            summary_text.find("non-whitelist,0,0,0,0,0") !=
                std::string::npos,
        "materialized summary counts are incorrect");

  chromap::MappingParameters bam_parameters;
  bam_parameters.mapping_output_format = chromap::MAPPINGFORMAT_BAM;
  bam_parameters.mapping_output_file_path = bam_output;
  bam_parameters.atac_fragment_output_file_path = fragments_output;
  bam_parameters.atac_fragment_binary_output_file_path = evidence_output;
  bam_parameters.sort_bam = true;
  bam_parameters.temp_directory_path = root;
  bam_parameters.emit_noY_stream = true;
  bam_parameters.noY_output_path = noy_bam_output;
  bam_parameters.emit_Y_stream = true;
  bam_parameters.Y_output_path = y_bam_output;
  const auto bam_result = chromap::MaterializeAtacSpillRecords(
      {shard1, shard0}, bam_parameters);
  Check(bam_result.ok && bam_result.output_fragment_count == 3,
        bam_result.message);
  std::ifstream bam_stream(bam_output.c_str(), std::ios::binary | std::ios::ate);
  Check(bam_stream.is_open() && bam_stream.tellg() > 0,
        "64-bit BAM-sort materialization produced no BAM output");
  std::ifstream evidence_stream(evidence_output.c_str(),
                                std::ios::binary | std::ios::ate);
  std::ifstream evidence_chroms_stream((evidence_output + ".chroms.tsv").c_str(),
                                       std::ios::binary | std::ios::ate);
  Check(evidence_stream.is_open() && evidence_stream.tellg() > 32 &&
            evidence_chroms_stream.is_open() &&
            evidence_chroms_stream.tellg() > 0,
        "AEV1 evidence materialization produced an empty output");
  std::ifstream noy_bam_stream(noy_bam_output.c_str(),
                               std::ios::binary | std::ios::ate);
  std::ifstream y_bam_stream(y_bam_output.c_str(),
                             std::ios::binary | std::ios::ate);
  Check(noy_bam_stream.is_open() && noy_bam_stream.tellg() > 0 &&
            y_bam_stream.is_open() && y_bam_stream.tellg() > 0,
        "Y-routed BAM materialization produced an empty output");

  {
    std::ofstream fasta(reference_fasta.c_str());
    Check(fasta.is_open(), "cannot create bounded CRAM reference");
    fasta << ">chr1\n";
    const std::string bases(1000, 'A');
    for (int i = 0; i < 1000; ++i) {
      fasta << bases << '\n';
    }
  }
  chromap::MappingParameters cram_parameters;
  cram_parameters.mapping_output_format = chromap::MAPPINGFORMAT_CRAM;
  cram_parameters.mapping_output_file_path = cram_output;
  cram_parameters.atac_fragment_output_file_path = cram_fragments_output;
  cram_parameters.reference_file_path = reference_fasta;
  cram_parameters.temp_directory_path = root;
  const auto cram_result = chromap::MaterializeAtacSpillRecords(
      {shard0, shard1}, cram_parameters);
  Check(cram_result.ok && cram_result.output_fragment_count == 3,
        cram_result.message);
  std::ifstream cram_stream(cram_output.c_str(),
                            std::ios::binary | std::ios::ate);
  Check(cram_stream.is_open() && cram_stream.tellg() > 0,
        "CRAM materialization produced no output");

  auto bulk_metadata0 = Metadata(0, wide_prefix + 10u);
  auto bulk_metadata1 = Metadata(1, wide_prefix + 12u);
  bulk_metadata0.flags = static_cast<uint16_t>(
      bulk_metadata0.flags | chromap::kAtacMergeableBulkLevelDedup);
  bulk_metadata1.flags = static_cast<uint16_t>(
      bulk_metadata1.flags | chromap::kAtacMergeableBulkLevelDedup);
  WriteShard(bulk_dedup_shard0, bulk_metadata0,
             {MakeRecord(0, 400, 60, 5)});
  WriteShard(bulk_dedup_shard1, bulk_metadata1,
             {MakeRecord(0, 400, 60, 0)});
  chromap::MappingParameters bulk_dedup_parameters;
  bulk_dedup_parameters.mapping_output_format = chromap::MAPPINGFORMAT_BED;
  bulk_dedup_parameters.mapping_output_file_path = bulk_dedup_output;
  const auto bulk_dedup_result = chromap::MaterializeAtacSpillRecords(
      {bulk_dedup_shard1, bulk_dedup_shard0}, bulk_dedup_parameters);
  Check(bulk_dedup_result.ok &&
            bulk_dedup_result.output_fragment_count == 1,
        bulk_dedup_result.message);
  const auto bulk_dedup_rows = ReadBed(bulk_dedup_output);
  Check(bulk_dedup_rows.size() == 1 &&
            bulk_dedup_rows[0][1] == "400" &&
            bulk_dedup_rows[0][3] == "AAAA" &&
            bulk_dedup_rows[0][4] == "2",
        "global abundance did not resolve barcoded bulk-level dedup");

  auto allocation_metadata0 = Metadata(0, wide_prefix + 20u);
  auto allocation_metadata1 = Metadata(1, wide_prefix + 22u);
  for (auto *metadata : {&allocation_metadata0, &allocation_metadata1}) {
    metadata->flags = static_cast<uint16_t>(
        (metadata->flags & ~chromap::kAtacMergeableRemovePcrDuplicates) |
        chromap::kAtacMergeableAllocateMultiMappings);
    metadata->mapq_threshold = 0;
    metadata->max_num_best_mappings = 2;
    metadata->multi_mapping_allocation_distance = 0;
    metadata->multi_mapping_allocation_seed = 11;
  }
  WriteShard(allocation_shard0, allocation_metadata0,
             {MakeRecord(0, 1000, 60, 0), MakeRecord(1, 1010, 0, 0),
              MakeRecord(1, 5000, 0, 0)});
  WriteShard(allocation_shard1, allocation_metadata1,
             {MakeRecord(0, 1020, 60, 0), MakeRecord(1, 7000, 60, 0)});
  chromap::MappingParameters allocation_parameters;
  allocation_parameters.mapping_output_format = chromap::MAPPINGFORMAT_BED;
  allocation_parameters.mapping_output_file_path = allocation_output;
  const auto allocation_result = chromap::MaterializeAtacSpillRecords(
      {allocation_shard0, allocation_shard1}, allocation_parameters);
  Check(allocation_result.ok && allocation_result.output_fragment_count == 4,
        allocation_result.message);
  const auto allocation_rows = ReadBed(allocation_output);
  bool saw_allocated = false;
  bool saw_unweighted = false;
  for (const auto &row : allocation_rows) {
    saw_allocated = saw_allocated || row[1] == "1010";
    saw_unweighted = saw_unweighted || row[1] == "5000";
  }
  Check(allocation_rows.size() == 4 && saw_allocated && !saw_unweighted,
        "global multimapping allocation did not use gathered unique support");

  auto late_metadata0 = Metadata(0, 0);
  auto late_metadata1 = Metadata(1, 0);
  late_metadata0.flags = static_cast<uint16_t>(
      late_metadata0.flags | chromap::kAtacMergeableReadRangeLateBound);
  late_metadata1.flags = static_cast<uint16_t>(
      late_metadata1.flags | chromap::kAtacMergeableReadRangeLateBound);
  WriteShard(late_shard0, late_metadata0,
             {MakeRecord(0, 100, 20, 0), MakeRecord(1, 200, 50, 1)});
  WriteShard(late_shard1, late_metadata1,
             {MakeRecord(0, 100, 60, 0), MakeRecord(1, 300, 50, 5)});
  chromap::MappingParameters late_parameters;
  late_parameters.mapping_output_format = chromap::MAPPINGFORMAT_BED;
  late_parameters.mapping_output_file_path = late_output;
  const auto late_result = chromap::MaterializeAtacSpillRecords(
      {late_shard1, late_shard0}, late_parameters);
  Check(late_result.ok && late_result.input_record_count == 4 &&
            late_result.output_fragment_count == 3,
        late_result.message);

  {
    chromap::SummaryMetadata wide_summary;
    wide_summary.SetBarcodeLength(4);
    const uint64_t wide_count =
        static_cast<uint64_t>(std::numeric_limits<uint32_t>::max()) + 123u;
    wide_summary.UpdateCount(/*AAAA=*/0, chromap::SUMMARY_METADATA_TOTAL,
                             wide_count);
    wide_summary.Output(wide_summary_output.c_str(), /*has_white_list=*/false,
                        {0.0, 0.0, 0.0, 0.0, 0.0},
                        /*output_num_cache_slots_info=*/false);
    std::ifstream wide_summary_stream(wide_summary_output.c_str());
    Check(wide_summary_stream.is_open(), "64-bit summary was not written");
    const std::string wide_summary_text(
        (std::istreambuf_iterator<char>(wide_summary_stream)),
        std::istreambuf_iterator<char>());
    Check(wide_summary_text.find(
              "AAAA,4294967418,0,4294967418,0,0,0.00000,0.00000") !=
              std::string::npos,
          "summary count above UINT32_MAX did not round-trip");
  }

  const auto incomplete =
      chromap::MaterializeAtacSpillRecords({shard0}, parameters);
  Check(!incomplete.ok, "incomplete ordinal set was accepted");
  WriteShard(duplicate0, Metadata(0, wide_prefix),
             {MakeRecord(0, 400, 50, 0)});
  const auto duplicate = chromap::MaterializeAtacSpillRecords(
      {shard0, duplicate0}, parameters);
  Check(!duplicate.ok, "duplicate shard ordinal was accepted");
  auto mismatched_metadata = Metadata(1, wide_prefix + 2u);
  mismatched_metadata.barcode_whitelist_fingerprint++;
  WriteShard(mismatched_model, mismatched_metadata,
             {MakeRecord(0, 500, 50, 5)});
  const auto mismatch = chromap::MaterializeAtacSpillRecords(
      {shard0, mismatched_model}, parameters);
  Check(!mismatch.ok, "mismatched barcode correction models were accepted");
  const auto empty = chromap::MaterializeAtacSpillRecords({}, parameters);
  Check(!empty.ok, "empty spill input set was accepted");

  std::cout << "PASS: mergeable ATAC spill writer/materializer\n";
  return 0;
}
