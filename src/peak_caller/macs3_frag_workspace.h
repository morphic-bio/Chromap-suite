#ifndef CHROMAP_PEAK_CALLER_MACS3_FRAG_WORKSPACE_H_
#define CHROMAP_PEAK_CALLER_MACS3_FRAG_WORKSPACE_H_

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "fragment_input.h"
#include "macs3_bdg_track.h"

namespace chromap {
namespace peaks {

struct Macs3FragWorkspaceParams {
  int32_t frag_pileup_max_count = 0;
  bool macs3_uint8_counts = true;
  int64_t effective_genome_size = 2913022398LL;
  int32_t llocal_bp = 10000;
  double score_pseudocount = 0.0;
};

struct Macs3FragEvent {
  Macs3FragEvent() : pos(0), delta(0) {}
  Macs3FragEvent(int32_t p, int32_t d) : pos(p), delta(d) {}

  int32_t pos;
  int32_t delta;
};

struct Macs3FragEventBuffer {
  std::vector<Macs3FragEvent> treat_span_events;
  std::vector<Macs3FragEvent> lambda_end_window_events;
};

struct Macs3FragPeakWorkspace {
  Macs3FragWorkspaceParams params;
  std::vector<std::string> chrom_names;
  std::vector<Macs3FragEventBuffer> events_by_chrom;
  std::vector<Macs3BdgTrack> treat_tracks;
  std::vector<Macs3BdgTrack> lambda_tracks;
  std::vector<Macs3BdgTrack> ppois_tracks;
  int64_t total_fragment_bp = 0;
  int64_t total_fragment_count = 0;
};

bool InitMacs3FragWorkspace(const std::vector<std::string>& chrom_names,
                            const Macs3FragWorkspaceParams& params,
                            Macs3FragPeakWorkspace* workspace);

bool AddMacs3FragWorkspaceFragment(Macs3FragPeakWorkspace* workspace,
                                   uint32_t chrom_id, int32_t start,
                                   int32_t end, int32_t count);

bool BuildMacs3FragWorkspaceFromFragments(
    const std::vector<ChromFragments>& by_chrom,
    const Macs3FragWorkspaceParams& params,
    Macs3FragPeakWorkspace* workspace);

bool FinalizeMacs3FragTreatTracks(Macs3FragPeakWorkspace* workspace);
bool FinalizeMacs3FragLambdaTracks(Macs3FragPeakWorkspace* workspace);
bool FinalizeMacs3FragPpoisTracks(Macs3FragPeakWorkspace* workspace);

bool WriteMacs3BdgTracksBedGraph(const std::string& path,
                                 const std::vector<std::string>& chrom_names,
                                 const std::vector<Macs3BdgTrack>& tracks);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_MACS3_FRAG_WORKSPACE_H_
