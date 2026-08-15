#ifndef ATAC_MERGEABLE_SPILL_H_
#define ATAC_MERGEABLE_SPILL_H_

#include <cstdio>
#include <cstdint>
#include <string>
#include <vector>

#include "atac_spill_record.h"

namespace chromap {

// Durable, cross-process envelope for sorted, pre-deduplication ATAC mapping
// records. The AtacSpillRecord payload codec remains independently versioned;
// this header versions shard identity and materialization policy.
enum AtacMergeableSpillFlags : uint16_t {
  kAtacMergeableRemovePcrDuplicates = 1u << 0,
  kAtacMergeableBulkLevelDedup = 1u << 1,
  kAtacMergeableTn5Shift = 1u << 2,
  kAtacMergeableOutputMappingsNotInWhitelist = 1u << 3,
  kAtacMergeableAllocateMultiMappings = 1u << 4,
  kAtacMergeableHasSummaryEvidence = 1u << 5,
  kAtacMergeableSummaryCardinality = 1u << 6,
  // Worker did not know its global prefix or final synchronized-record count
  // at process start. input_record_count is the count observed at EOF and the
  // materializer derives first_global_read_ordinal from ordered shard counts.
  kAtacMergeableReadRangeLateBound = 1u << 7,
};

#pragma pack(push, 1)
struct AtacMergeableSpillHeaderV3 {
  char magic[8];
  uint16_t format_version;
  uint16_t fixed_header_bytes;
  uint32_t metadata_bytes;
  uint32_t record_codec_version;
  uint16_t schema_mask;
  uint16_t flags;
  uint32_t shard_ordinal;
  uint32_t shard_count;
  uint64_t first_global_read_ordinal;
  uint64_t input_record_count;
  uint64_t record_count;
  uint32_t barcode_length;
  uint32_t num_reference_sequences;
  uint8_t mapq_threshold;
  int8_t tn5_forward_shift;
  int8_t tn5_reverse_shift;
  uint8_t barcode_correction_error_threshold;
  uint32_t sample_id_bytes;
  uint32_t input_id_bytes;
  uint32_t num_barcode_abundance_entries;
  uint64_t local_num_sample_barcodes;
  uint64_t barcode_whitelist_fingerprint;
  uint64_t barcode_correction_probability_bits;
  int32_t multi_mapping_allocation_distance;
  int32_t multi_mapping_allocation_seed;
  uint32_t max_num_best_mappings;
  uint64_t summary_evidence_count;
  uint32_t cache_size;
  uint32_t k_for_minhash;
  uint64_t frip_est_coefficient_bits[5];
};
#pragma pack(pop)

static_assert(sizeof(AtacMergeableSpillHeaderV3) == 172,
              "ATAC mergeable spill header must be 172 bytes");

static constexpr char kAtacMergeableSpillMagicV3[8] = {
    'A', 'T', 'A', 'C', 'M', 'S', '3', '\0'};
static constexpr uint16_t kAtacMergeableSpillFormatVersion = 3;

struct AtacMergeableSpillReference {
  std::string name;
  uint32_t length = 0;
};

struct AtacBarcodeAbundanceEntry {
  uint64_t barcode_key = 0;
  uint64_t count = 0;

  AtacBarcodeAbundanceEntry() = default;
  AtacBarcodeAbundanceEntry(uint64_t key, uint64_t abundance)
      : barcode_key(key), count(abundance) {}
};

struct AtacSummaryEvidence {
  uint64_t raw_barcode_key = 0;
  uint32_t raw_barcode_n_mask = 0;
  int32_t cache_slot1 = -1;
  int32_t cache_slot2 = -1;
  std::string raw_barcode_qual;
};

struct AtacMergeableSpillMetadata {
  uint16_t schema_mask = 0;
  uint16_t flags = 0;
  uint32_t shard_ordinal = 0;
  uint32_t shard_count = 0;
  uint64_t first_global_read_ordinal = 0;
  uint64_t input_record_count = 0;
  uint32_t barcode_length = 0;
  uint8_t mapq_threshold = 0;
  int8_t tn5_forward_shift = 0;
  int8_t tn5_reverse_shift = 0;
  uint64_t local_num_sample_barcodes = 0;
  uint64_t barcode_whitelist_fingerprint = 0;
  uint8_t barcode_correction_error_threshold = 0;
  double barcode_correction_probability_threshold = 0.0;
  int32_t multi_mapping_allocation_distance = 0;
  int32_t multi_mapping_allocation_seed = 11;
  uint32_t max_num_best_mappings = 1;
  uint64_t summary_evidence_count = 0;
  uint32_t cache_size = 4000003;
  uint32_t k_for_minhash = 250;
  std::vector<double> frip_est_coefficients;
  std::string sample_id;
  std::string input_id;
  std::vector<AtacMergeableSpillReference> references;
  std::vector<AtacBarcodeAbundanceEntry> barcode_abundance_entries;
};

class AtacMergeableSpillWriter {
 public:
  AtacMergeableSpillWriter() = default;
  ~AtacMergeableSpillWriter();

  bool Open(const std::string &path,
            const AtacMergeableSpillMetadata &metadata,
            std::string *error);
  bool Append(uint32_t rid, const AtacSpillRecord &record,
              std::string *error);
  bool AppendSummaryEvidence(const AtacSummaryEvidence &evidence,
                             std::string *error);
  bool Finalize(std::string *error);
  uint64_t record_count() const { return record_count_; }
  uint64_t summary_evidence_count() const { return summary_evidence_count_; }

 private:
  bool Fail(const std::string &message, std::string *error);

  FILE *file_ = nullptr;
  std::string output_path_;
  std::string temporary_path_;
  AtacMergeableSpillMetadata metadata_;
  uint64_t record_count_ = 0;
  uint64_t summary_evidence_count_ = 0;
  bool finalized_ = false;
  bool have_previous_ = false;
  uint32_t previous_rid_ = 0;
  AtacSpillRecord previous_record_;
};

class AtacMergeableSpillReader {
 public:
  AtacMergeableSpillReader() = default;
  ~AtacMergeableSpillReader();

  bool Open(const std::string &path, std::string *error);
  // On success, `eof` distinguishes a record from the verified end of stream.
  bool ReadNext(uint32_t *rid, AtacSpillRecord *record, bool *eof,
                std::string *error);
  // BED gather needs only the fixed fragment prefix and raw barcode evidence.
  // This decoder streams large blocks, validates and skips both serialized SAM
  // mates in place, and never constructs their strings/vectors.
  bool ReadNextBed(uint32_t *rid, AtacSpillRecord *record, bool *eof,
                   std::string *error);
  bool ReadNextSummaryEvidence(AtacSummaryEvidence *evidence, bool *eof,
                               std::string *error);
  // Summary evidence precedes mappings. A BED-only materializer that was not
  // asked to emit summary metadata can seek over this fixed-width section.
  bool SkipRemainingSummaryEvidence(std::string *error);

  const AtacMergeableSpillMetadata &metadata() const { return metadata_; }
  uint64_t expected_record_count() const { return expected_record_count_; }
  uint64_t records_read() const { return records_read_; }
  const std::string &path() const { return path_; }

 private:
  bool Fail(const std::string &message, std::string *error);
  bool EnsureBedStreamBytes(size_t count, std::string *error);

  FILE *file_ = nullptr;
  std::string path_;
  AtacMergeableSpillMetadata metadata_;
  uint64_t expected_record_count_ = 0;
  uint64_t records_read_ = 0;
  uint64_t expected_summary_evidence_count_ = 0;
  uint64_t summary_evidence_read_ = 0;
  bool verified_eof_ = false;
  bool have_previous_ = false;
  uint32_t previous_rid_ = 0;
  PairedEndMappingWithBarcode previous_mapping_;
  std::vector<char> bed_stream_buffer_;
  size_t bed_stream_begin_ = 0;
  size_t bed_stream_end_ = 0;
  bool full_decoder_started_ = false;
  bool bed_decoder_started_ = false;
};

// Adds the shard range's global prefix to the fragment and BAM-pair read ids.
// Returns false if the uint64 read-id contract would overflow or if a
// record's local id lies outside the declared input range.
bool GlobalizeAtacSpillReadId(AtacSpillRecord *record,
                              uint64_t first_global_read_ordinal,
                              uint64_t input_record_count,
                              std::string *error);

}  // namespace chromap

#endif  // ATAC_MERGEABLE_SPILL_H_
