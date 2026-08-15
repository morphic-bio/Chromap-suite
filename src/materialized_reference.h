#ifndef MATERIALIZED_REFERENCE_H_
#define MATERIALIZED_REFERENCE_H_

#include <cstdint>

namespace chromap {

// Identity recorded in both the optional materialized-reference sidecar and
// the Chromap index footer that binds the two files together.
struct MaterializedReferenceInfo {
  uint64_t fingerprint = 0;
  uint64_t num_bases = 0;
  uint32_t num_sequences = 0;
};

}  // namespace chromap

#endif  // MATERIALIZED_REFERENCE_H_
