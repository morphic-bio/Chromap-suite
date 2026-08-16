#ifndef ATAC_SPILL_MATERIALIZER_H_
#define ATAC_SPILL_MATERIALIZER_H_

#include <cstdint>
#include <string>
#include <vector>

#include "mapping_parameters.h"

namespace chromap {

struct AtacSpillMaterializationResult {
  bool ok = false;
  std::string message;
  std::string sample_id;
  std::string input_id;
  uint32_t shard_count = 0;
  uint64_t input_record_count = 0;
  uint64_t spill_record_count = 0;
  uint64_t corrected_barcode_record_count = 0;
  uint64_t rejected_barcode_record_count = 0;
  uint64_t output_fragment_count = 0;
  bool used_parallel_hot_spill = false;
  // Timed application phases only. Up-front contract validation is not added
  // to either value.
  double merge_output_seconds = 0.0;
  double terminal_bed_export_seconds = 0.0;
};

// Validates an ordinal-complete mergeable-spill set, globally merges and
// deduplicates it, then emits canonical Chromap ATAC outputs. The genome index
// and genome sequence are not loaded; reference names and lengths come from
// the validated spill headers. CRAM output may still name a FASTA for htslib.
AtacSpillMaterializationResult MaterializeAtacSpillRecords(
    const std::vector<std::string> &spill_paths,
    const MappingParameters &output_parameters);

}  // namespace chromap

#endif  // ATAC_SPILL_MATERIALIZER_H_
