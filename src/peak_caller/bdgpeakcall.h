#ifndef CHROMAP_PEAK_CALLER_BDGPEAKCALL_H_
#define CHROMAP_PEAK_CALLER_BDGPEAKCALL_H_

#include <cstdint>
#include <string>
#include <vector>

#include "macs3_bdg_track.h"

namespace chromap {
namespace peaks {

class StageProfileCollector;

// MACS3 `bdgpeakcall` region calling on a score bedGraph (diagnostic).
// Reads the file like MACS3 BedGraphIO.read_bedGraph / bedGraphTrackI.add_loc
// (sequential merge; no implicit gap fill). Runs BedGraph.call_peaks with
// call_summits=False (interval endpoints from merged above-cutoff segments).
//
// Writes BED3 (chrom, start, end), one interval per line, sorted lexicographically
// by chrom then start, end.
bool RunMacs3StyleBdgPeakCallFromBedGraph(const std::string& path_bdg,
                                          const std::string& path_bed_out,
                                          float cutoff, int32_t min_length,
                                          int32_t max_gap);

// MACS3-style narrowPeak + summits: regions from ppois bdg (same as `bdgpeakcall` on
// that track); summit = midpoint of the treat-pileup-maximal above-cutoff subsegment
// (matches MACS3 `__close_peak_wo_subpeaks`). p/q (BH over summit p) and signal from
// treat+lambda at that subsegment. Summits BED5: same peak order; score = -log10 q.
bool RunMacs3FragPpoisNarrowPeaks(
    const std::string& path_ppois_bdg, const std::string& path_treat_bdg,
    const std::string& path_lambda_bdg, float cutoff, int32_t min_length,
    int32_t max_gap, double pseudocount, const std::string& path_narrowpeak,
    const std::string& path_summits_bed, StageProfileCollector* stage_profile = nullptr,
    int64_t profile_input_frag_rows = -1);

bool RunMacs3FragPpoisNarrowPeaksFromTracks(
    const std::vector<std::string>& chrom_names,
    const std::vector<Macs3BdgTrack>& ppois_tracks,
    const std::vector<Macs3BdgTrack>& treat_tracks,
    const std::vector<Macs3BdgTrack>& lambda_tracks,
    float cutoff, int32_t min_length, int32_t max_gap, double pseudocount,
    const std::string& path_narrowpeak, const std::string& path_summits_bed,
    StageProfileCollector* stage_profile = nullptr,
    int64_t profile_input_frag_rows = -1);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_BDGPEAKCALL_H_
