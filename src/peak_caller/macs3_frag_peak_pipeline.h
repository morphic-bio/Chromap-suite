#ifndef CHROMAP_PEAK_CALLER_MACS3_FRAG_PEAK_PIPELINE_H_
#define CHROMAP_PEAK_CALLER_MACS3_FRAG_PEAK_PIPELINE_H_

#include <cstdint>
#include <string>
#include <vector>

#include "fragment_input.h"

namespace chromap {
namespace peaks {

class StageProfileCollector;
struct Macs3FragPeakWorkspace;

// Parameters for the validated MACS3-compatible FRAG pileup → lambda → ppois
// → bdgpeakcall → narrowPeak/summits path (see parity harnesses).
struct Macs3FragPeakPipelineParams {
  float bdgpeakcall_cutoff = 5.f;
  int32_t min_length = 200;
  int32_t max_gap = 30;
  int32_t frag_pileup_max_count = 0;
  bool macs3_uint8_counts = true;
  int64_t effective_genome_size = 2913022398LL;
  int32_t llocal_bp = 10000;
  double score_pseudocount = 0.0;
};

// Optional explicit intermediate bedGraph paths (empty = files under work_dir).
struct Macs3FragPeakPipelinePaths {
  std::string treat_bdg;
  std::string lambda_bdg;
  std::string ppois_bdg;
};

// -log10(p) cutoff for bdgpeakcall on the ppois score track (MACS3 -p 1e-5 → 5).
float BdgPeakCallCutoffFromPValue(double p_value);

// Fragment-span pileup, control_lambda, ppois score, then narrowPeak + summits.
// keep_intermediates_dir: if non-empty, writes treat/lambda/ppois bedGraphs there
// (directory is created) and does not delete it. If empty, uses a unique
// directory under work_dir_parent (or /tmp) and removes it after success.
bool RunMacs3FragPeakPipelineFromFragments(
    const std::vector<ChromFragments>& by_chrom,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message,
    StageProfileCollector* stage_profile = nullptr);

bool RunMacs3FragPeakPipelineFromWorkspace(
    Macs3FragPeakWorkspace* workspace,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message);

}  // namespace peaks
}  // namespace chromap

#endif  // CHROMAP_PEAK_CALLER_MACS3_FRAG_PEAK_PIPELINE_H_
