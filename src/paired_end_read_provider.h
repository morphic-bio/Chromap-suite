#ifndef PAIRED_END_READ_PROVIDER_H_
#define PAIRED_END_READ_PROVIDER_H_

#include <cstddef>
#include <cstdint>
#include <string>

#include "sequence_batch.h"

namespace chromap {

// Executor-neutral paired-read input boundary. Implementations may read a
// regular file, a shared-filesystem range, or another bounded source, but must
// return synchronized records. MappingParameters retains shared ownership for
// the complete run. Chromap remains responsible for effective-range
// processing, mapping, and spill/output semantics.
class PairedEndReadProvider {
 public:
  virtual ~PairedEndReadProvider() = default;

  virtual bool HasBarcode() const = 0;

  // Preserve Chromap's ordinary comma-separated FASTQ contract: each input
  // lane is mapped by one pass through the lane-scoped mapping loop. Providers
  // with more than one lane must select each lane before LoadBatch is called.
  // The defaults retain the original single-lane provider behavior.
  virtual size_t NumInputLanes() const { return 1; }
  virtual bool SelectInputLane(size_t lane_index, std::string &error) {
    if (lane_index != 0) {
      error = "paired-end read provider lane index is out of range";
      return false;
    }
    error.clear();
    return true;
  }

  // Fill at most max_pairs records into the supplied batches. Return true on
  // success, including clean EOF (num_loaded_pairs == 0). Implementations must
  // reset all supplied batches before filling them and preserve record order.
  virtual bool LoadBatch(uint32_t max_pairs, SequenceBatch &read_batch1,
                         SequenceBatch &read_batch2,
                         SequenceBatch &barcode_batch,
                         uint32_t &num_loaded_pairs,
                         std::string &error) = 0;
};

}  // namespace chromap

#endif  // PAIRED_END_READ_PROVIDER_H_
