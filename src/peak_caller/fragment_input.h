#ifndef CHROMAP_PEAK_CALLER_FRAGMENT_INPUT_H_
#define CHROMAP_PEAK_CALLER_FRAGMENT_INPUT_H_

#include <cstdint>
#include <string>
#include <vector>

namespace chromap {
namespace peaks {

// Chromap/ARC style: chrom, start, end, barcode, duplicate_count (0-based,
// [start, end) BED style).
struct Fragment {
  int32_t start = 0;
  int32_t end = 0;
  int32_t count = 1;
};

struct ChromFragments {
  std::string name;
  std::vector<Fragment> frags;
  int32_t max_end = 0;
};

// Reads a tab-separated fragments file (plain or .gz). Lines starting with
// '#', empty, or malformed are skipped. Returns true on full success (file
// opened, no read errors; empty input is an error for callers to decide).
// Fragments are grouped by `chrom` in strict lexicographic order and sorted
// by (start, end) within each group.
bool LoadFragmentsFromTsv(
    const std::string& path, std::vector<ChromFragments>* by_chrom);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_FRAGMENT_INPUT_H_
