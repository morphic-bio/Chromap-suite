#ifndef MACS3_FRAGMENT_BUCKETS_H_
#define MACS3_FRAGMENT_BUCKETS_H_

#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "rapidmacs/fragments.h"

namespace chromap {
namespace peaks {

// Mapping uses reference-sequence ids for its in-memory buckets, so sparse
// inputs leave empty buckets for sequences that were never observed. File FRAG
// input interns only observed sequences. Compact the memory representation to
// that same contract before peak calling.
inline bool CompactNonEmptyMacs3FragmentBuckets(
    std::vector<std::vector<macs3::FragmentRecord>> *buckets,
    std::vector<std::string> *chrom_names,
    std::string *error) {
  if (buckets == nullptr || chrom_names == nullptr) {
    if (error != nullptr) {
      *error = "missing fragment buckets or chromosome names";
    }
    return false;
  }
  if (chrom_names->size() < buckets->size()) {
    if (error != nullptr) {
      *error = "chromosome name table does not cover fragment buckets";
    }
    return false;
  }

  // Validate before mutating so a malformed bucket table fails atomically.
  for (size_t bucket_id = 0; bucket_id < buckets->size(); ++bucket_id) {
    for (const auto &record : (*buckets)[bucket_id]) {
      if (record.chrom_id != static_cast<int32_t>(bucket_id)) {
        if (error != nullptr) {
          *error = "fragment chromosome id does not match its bucket";
        }
        return false;
      }
    }
  }

  size_t write_id = 0;
  for (size_t read_id = 0; read_id < buckets->size(); ++read_id) {
    if ((*buckets)[read_id].empty()) {
      continue;
    }
    if (write_id > static_cast<size_t>(std::numeric_limits<int32_t>::max())) {
      if (error != nullptr) {
        *error = "too many non-empty chromosome buckets";
      }
      return false;
    }
    if (write_id != read_id) {
      (*buckets)[write_id] = std::move((*buckets)[read_id]);
      (*chrom_names)[write_id] = std::move((*chrom_names)[read_id]);
    }
    for (auto &record : (*buckets)[write_id]) {
      record.chrom_id = static_cast<int32_t>(write_id);
    }
    ++write_id;
  }
  buckets->resize(write_id);
  chrom_names->resize(write_id);
  return true;
}

}  // namespace peaks
}  // namespace chromap

#endif  // MACS3_FRAGMENT_BUCKETS_H_
