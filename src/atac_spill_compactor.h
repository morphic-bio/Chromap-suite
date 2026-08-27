#ifndef ATAC_SPILL_COMPACTOR_H_
#define ATAC_SPILL_COMPACTOR_H_

#include <cstdint>
#include <string>
#include <vector>

namespace chromap {

struct AtacSpillCompactionResult {
  bool ok = false;
  std::string message;
  std::string sample_id;
  std::string input_id;
  uint32_t first_shard_ordinal = 0;
  uint32_t shard_span = 0;
  uint32_t shard_count = 0;
  uint64_t input_record_count = 0;
  uint64_t spill_record_count = 0;
  bool wrote_hot_sidecar = false;
};

// Merge adjacent, sorted ATACMS3/ATACMS4 raw inputs into one immutable
// ATACMS4 run. This operation only rebases run-local read ids and combines raw
// evidence. It deliberately performs no barcode correction, mapping filter,
// Tn5 shift, multimapping allocation, or duplicate collapse.
AtacSpillCompactionResult CompactAtacSpillRecords(
    const std::vector<std::string> &spill_paths,
    const std::string &output_path);

}  // namespace chromap

#endif  // ATAC_SPILL_COMPACTOR_H_
