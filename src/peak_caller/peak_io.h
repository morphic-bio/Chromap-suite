#ifndef CHROMAP_PEAK_CALLER_PEAK_IO_H_
#define CHROMAP_PEAK_CALLER_PEAK_IO_H_

#include <cstdio>
#include <string>
#include <vector>

#include "call_peaks.h"

namespace chromap {
namespace peaks {

// Writes 10-col narrowPeak (tab); sorts peaks by (chrom, start, end) first.
// Returns false on I/O or invalid intervals.
bool WriteNarrowPeak(const std::string& path,
                     const std::vector<Peak>& peaks);

// Simple tab-separated metrics: peak_count, total_bp, median_width_bp, ...
bool WriteRunSummaryTsv(
    const std::string& path, const std::vector<Peak>& all_peaks,
    const std::string& argv_line, const std::string& program_params);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_PEAK_IO_H_
