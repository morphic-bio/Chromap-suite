// Standalone fragment-based peak caller (MACS3-inspired baseline; not a full
// MACS3 reimplementation). Next phase: optional pipeline integration in Chromap/STAR.
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "peak_caller/binned_signal.h"
#include "peak_caller/bdgpeakcall.h"
#include "peak_caller/call_peaks.h"
#include "peak_caller/fragment_input.h"
#include "peak_caller/macs3_frag_peak_pipeline.h"
#include "peak_caller/peak_io.h"
#include "peak_caller/stage_profile.h"

namespace {

void Usage() {
  std::cerr
      << "Usage: chromap_callpeaks --input <fragments.tsv|fragments.tsv.gz> "
         "[--out-prefix <path/prefix>] [options]\n"
         "  --input|-i     Chromap/ARC 5-col TSV(/.gz): chrom start end barcode count\n"
         "  --out-prefix   Prefix for <prefix>.narrowPeak and <prefix>.summary.tsv "
         "(omit if only emitting pileup outputs)\n"
         "  --pileup-bdg   Binned Tn5-extended *cut* pileup bedGraph (diagnostic; "
         "differs from MACS3 `pileup -f FRAG`)\n"
         "  --frag-span-pileup-bdg  Fragment-span pileup bedGraph ([start,end) weighted by "
         "col5; matches MACS3 `pileup -f FRAG` when --frag-pileup-macs3-uint8-counts is set)\n"
         "  --frag-pileup-max-count N  Cap col5 like MACS3 --max-count (0 = no cap; applies "
         "to --frag-span-pileup-bdg only)\n"
         "  --frag-pileup-macs3-uint8-counts  Coerce fragment weights through uint8 (MACS3 "
         "PETrackII); recommended for byte-for-byte parity vs macs3 pileup\n"
         "  --frag-lambda-bdg  MACS3 no-control `control_lambda` (diagnostic; FRAG + llocal, "
         "see MACS3 PeakDetect __call_peaks_wo_control)\n"
         "  --frag-lambda-effective-genome N  -g size (default 2913022398 = hs GRCh38)\n"
         "  --frag-lambda-llocal  Large local window in bp (default: 10000, MACS3 llocal)\n"
         "  --frag-lambda-slocal N  (ignored for no-control FRAG; slocal is not applied to\n"
         "                lambda in MACS3; reserved for future -c / parity work)\n"
         "  --frag-score-ppois-bdg PATH  MACS3 `bdgcmp -m ppois` score bedGraph (diagnostic)\n"
         "  --frag-score-fe-bdg PATH     MACS3 `bdgcmp -m FE` score bedGraph (diagnostic)\n"
         "  --frag-score-pseudocount X   Same as `macs3 bdgcmp -p` (default: 0)\n"
         "  --bdgpeakcall-input PATH     Score bedGraph; run MACS3-style bdgpeakcall (diagnostic)\n"
         "  --bdgpeakcall-output PATH    BED3 regions (requires --bdgpeakcall-input)\n"
         "  --bdgpeakcall-cutoff X       -c cutoff (default: 5)\n"
         "  --bdgpeakcall-min-len N      -l min length (default: 200)\n"
         "  --bdgpeakcall-max-gap N      -g max gap (default: 30)\n"
         "  --macs3-frag-narrowpeak PATH  Diagnostic narrowPeak (FRAG: ppois+lambda+span;\n"
         "                pileup/lambda/ppois bedGraphs are written automatically if omitted)\n"
         "  --macs3-frag-summits PATH    Diagnostic one-line summits BED5 (q -log10 as score)\n"
         "  --macs3-frag-ppois-bdg PATH  Override ppois bedGraph for --macs3-frag-* (else --frag-score-ppois-bdg)\n"
         "  --bin-width    Bin width in bp (default: 50)\n"
         "  --ext-size     Tn5 extension in bp, symmetric; default 150 (2x75; shiftTAG-like)\n"
         "  --local-window  Half-window in bins for local lambda (default: 200; ~20kb for 50bp)\n"
         "  --fdr          FDR (BH) threshold (default: 0.05)\n"
         "  --p-value      Raw p upper tail threshold (default: 0.01)\n"
         "  --merge-gap     Merge significant bins separated by up to this many bins (0)\n"
         "  --min-peak-bp   Drop peaks narrower than this width in bp (0 = off)\n"
         "  --min-summit-cuts  Drop peaks whose summit bin count is below this (0 = off)\n"
         "  --profile-stages PATH  Optional stage timing TSV (wall/user/sys/RSS/rows/bytes)\n"
         "  -h|--help      This message\n";
}

std::string JoinArgv(int argc, char** argv) {
  std::ostringstream o;
  for (int i = 0; i < argc; ++i) {
    if (i) {
      o << ' ';
    }
    o << argv[i];
  }
  return o.str();
}

}  // namespace

int main(int argc, char** argv) {
  std::string input;
  std::string out_prefix;
  std::string pileup_bdg;
  std::string frag_span_pileup_bdg;
  std::string frag_lambda_bdg;
  std::string frag_score_ppois_bdg;
  std::string frag_score_fe_bdg;
  int32_t frag_pileup_max_count = 0;
  bool frag_pileup_macs3_uint8 = false;
  int64_t frag_lambda_effective_genome = 2913022398LL;  // MACS3 Constants.EFFECTIVEGS["hs"]
  int32_t frag_lambda_llocal = 10000;
  bool frag_lambda_slocal_set = false;
  double frag_score_pseudocount = 0.0;
  std::string bdgpeakcall_input;
  std::string bdgpeakcall_output;
  float bdgpeakcall_cutoff = 5.f;
  int32_t bdgpeakcall_min_len = 200;
  int32_t bdgpeakcall_max_gap = 30;
  std::string macs3_frag_narrowpeak_out;
  std::string macs3_frag_summits_out;
  std::string macs3_frag_ppois_override;
  int32_t bin_width = 50;
  int32_t ext_size = 150;
  int32_t local_window = 200;
  double fdr = 0.05;
  double pval = 0.01;
  int32_t merge_gap = 0;
  int32_t min_peak_bp = 0;
  int64_t min_summit_cuts = 0;
  std::string profile_stages_path;
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "-h" || a == "--help") {
      Usage();
      return 0;
    }
    if (a == "--input" || a == "-i") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      input = argv[++i];
    } else if (a == "--out-prefix") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      out_prefix = argv[++i];
    } else if (a == "--pileup-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      pileup_bdg = argv[++i];
    } else if (a == "--frag-span-pileup-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_span_pileup_bdg = argv[++i];
    } else if (a == "--frag-pileup-max-count") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_pileup_max_count = std::atoi(argv[++i]);
    } else if (a == "--frag-pileup-macs3-uint8-counts") {
      frag_pileup_macs3_uint8 = true;
    } else if (a == "--frag-lambda-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_lambda_bdg = argv[++i];
    } else if (a == "--frag-lambda-effective-genome") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_lambda_effective_genome = std::atoll(argv[++i]);
    } else if (a == "--frag-lambda-llocal") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_lambda_llocal = std::atoi(argv[++i]);
    } else if (a == "--frag-lambda-slocal") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      (void)std::strtol(argv[++i], nullptr, 10);
      frag_lambda_slocal_set = true;
    } else if (a == "--frag-score-ppois-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_score_ppois_bdg = argv[++i];
    } else if (a == "--frag-score-fe-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_score_fe_bdg = argv[++i];
    } else if (a == "--frag-score-pseudocount") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      frag_score_pseudocount = std::atof(argv[++i]);
    } else if (a == "--bdgpeakcall-input") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bdgpeakcall_input = argv[++i];
    } else if (a == "--bdgpeakcall-output") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bdgpeakcall_output = argv[++i];
    } else if (a == "--bdgpeakcall-cutoff") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bdgpeakcall_cutoff = static_cast<float>(std::atof(argv[++i]));
    } else if (a == "--bdgpeakcall-min-len") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bdgpeakcall_min_len = std::atoi(argv[++i]);
    } else if (a == "--bdgpeakcall-max-gap") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bdgpeakcall_max_gap = std::atoi(argv[++i]);
    } else if (a == "--macs3-frag-narrowpeak") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs3_frag_narrowpeak_out = argv[++i];
    } else if (a == "--macs3-frag-summits") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs3_frag_summits_out = argv[++i];
    } else if (a == "--macs3-frag-ppois-bdg") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      macs3_frag_ppois_override = argv[++i];
    } else if (a == "--bin-width") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      bin_width = std::atoi(argv[++i]);
    } else if (a == "--ext-size") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      ext_size = std::atoi(argv[++i]);
    } else if (a == "--local-window") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      local_window = std::atoi(argv[++i]);
    } else if (a == "--fdr") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      fdr = std::atof(argv[++i]);
    } else if (a == "--p-value") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      pval = std::atof(argv[++i]);
    } else if (a == "--merge-gap") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      merge_gap = std::atoi(argv[++i]);
    } else if (a == "--min-peak-bp") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      min_peak_bp = std::atoi(argv[++i]);
    } else if (a == "--min-summit-cuts") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      min_summit_cuts = std::atoll(argv[++i]);
    } else if (a == "--profile-stages") {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << a << "\n";
        return 2;
      }
      profile_stages_path = argv[++i];
    } else {
      std::cerr << "Unknown argument: " << a << "\n";
      Usage();
      return 2;
    }
  }
  const bool bdgpeakcall_pair =
      !bdgpeakcall_input.empty() || !bdgpeakcall_output.empty();
  if (bdgpeakcall_pair &&
      (bdgpeakcall_input.empty() || bdgpeakcall_output.empty())) {
    std::cerr << "Both --bdgpeakcall-input and --bdgpeakcall-output are required "
                 "together.\n";
    return 2;
  }
  const bool need_macs3_frag_np =
      !macs3_frag_narrowpeak_out.empty() || !macs3_frag_summits_out.empty();
  const bool need_fragment_input =
      !out_prefix.empty() || !pileup_bdg.empty() || !frag_span_pileup_bdg.empty() ||
      !frag_lambda_bdg.empty() || !frag_score_ppois_bdg.empty() ||
      !frag_score_fe_bdg.empty() || need_macs3_frag_np;
  if (need_fragment_input && input.empty()) {
    Usage();
    return 2;
  }
  if (input.empty() && !bdgpeakcall_pair) {
    Usage();
    return 2;
  }
  if (!bdgpeakcall_pair && out_prefix.empty() && pileup_bdg.empty() &&
      frag_span_pileup_bdg.empty() && frag_lambda_bdg.empty() &&
      frag_score_ppois_bdg.empty() && frag_score_fe_bdg.empty() &&
      !need_macs3_frag_np) {
    Usage();
    return 2;
  }
  if (frag_lambda_slocal_set) {
    std::cerr
        << "Note: --frag-lambda-slocal is not used for MACS3 no-control FRAG lambda; "
           "only --frag-lambda-llocal matches macs3 `callpeak` (PeakDetect __call_peaks_wo_control).\n";
  }
  const bool need_bins = !out_prefix.empty() || !pileup_bdg.empty();
  if (need_bins && (bin_width < 1 || ext_size < 2 || (ext_size % 2) != 0)) {
    std::cerr
        << "Invalid --bin-width or --ext-size (ext-size must be even, >=2)\n";
    return 2;
  }
  if (frag_pileup_max_count < 0) {
    std::cerr << "Invalid --frag-pileup-max-count\n";
    return 2;
  }
  if (min_peak_bp < 0 || min_summit_cuts < 0) {
    std::cerr << "Invalid --min-peak-bp or --min-summit-cuts\n";
    return 2;
  }
  if (frag_lambda_llocal < 1 || frag_lambda_effective_genome < 1) {
    std::cerr << "Invalid --frag-lambda-llocal or --frag-lambda-effective-genome\n";
    return 2;
  }
  if (bdgpeakcall_min_len < 1 || bdgpeakcall_max_gap < 0) {
    std::cerr << "Invalid --bdgpeakcall-min-len or --bdgpeakcall-max-gap\n";
    return 2;
  }
  if (fdr < 0.0 || fdr > 1.0 || pval < 0.0 || pval > 1.0) {
    return 2;
  }
  std::unique_ptr<chromap::peaks::StageProfileCollector> profile_coll;
  std::unique_ptr<chromap::peaks::StageProfileTotalGuard> profile_total_guard;
  int64_t profile_frags = 0;
  int64_t profile_out_peaks = 0;
  if (!profile_stages_path.empty()) {
    profile_coll.reset(new chromap::peaks::StageProfileCollector(profile_stages_path));
    if (!profile_coll->IsOpen()) {
      std::cerr << "Failed to open --profile-stages output: " << profile_stages_path << "\n";
      return 1;
    }
    profile_total_guard.reset(new chromap::peaks::StageProfileTotalGuard(
        profile_coll.get(), chromap::peaks::StageProfileCollector::Now(), &profile_frags,
        &profile_out_peaks));
  }
  if (!bdgpeakcall_input.empty()) {
    if (!chromap::peaks::RunMacs3StyleBdgPeakCallFromBedGraph(
            bdgpeakcall_input, bdgpeakcall_output, bdgpeakcall_cutoff,
            bdgpeakcall_min_len, bdgpeakcall_max_gap)) {
      std::cerr << "Failed MACS3-style bdgpeakcall diagnostic\n";
      return 1;
    }
    if (!need_fragment_input) {
      return 0;
    }
  }
  std::vector<chromap::peaks::ChromFragments> chs;
  if (profile_coll) {
    const chromap::peaks::StageProfileCollector::Tick t0 =
        chromap::peaks::StageProfileCollector::Now();
    if (!chromap::peaks::LoadFragmentsFromTsv(input, &chs)) {
      std::cerr << "Failed to read input: " << input << "\n";
      return 1;
    }
    const chromap::peaks::StageProfileCollector::Tick t_load_end =
        chromap::peaks::StageProfileCollector::Now();
    for (const chromap::peaks::ChromFragments& ch : chs) {
      profile_frags += static_cast<int64_t>(ch.frags.size());
    }
    profile_coll->Record("load_fragments", t0, t_load_end, profile_frags,
                         static_cast<int64_t>(chs.size()), 0, "");
    const chromap::peaks::StageProfileCollector::Tick ts =
        chromap::peaks::StageProfileCollector::Now();
    profile_coll->Record("sort_fragments_by_chrom", ts, ts, profile_frags,
                         static_cast<int64_t>(chs.size()), 0,
                         "included_in_load_fragments;not_timed_separately");
  } else {
    if (!chromap::peaks::LoadFragmentsFromTsv(input, &chs)) {
      std::cerr << "Failed to read input: " << input << "\n";
      return 1;
    }
  }
  if (!pileup_bdg.empty()) {
    if (!chromap::peaks::WritePileupBedGraph(pileup_bdg, chs, bin_width, ext_size)) {
      std::cerr << "Failed to write pileup bedGraph: " << pileup_bdg << "\n";
      return 1;
    }
  }
  if (need_macs3_frag_np) {
    chromap::peaks::Macs3FragPeakPipelineParams pr;
    pr.bdgpeakcall_cutoff = bdgpeakcall_cutoff;
    pr.min_length = bdgpeakcall_min_len;
    pr.max_gap = bdgpeakcall_max_gap;
    pr.frag_pileup_max_count = frag_pileup_max_count;
    pr.macs3_uint8_counts = frag_pileup_macs3_uint8;
    pr.effective_genome_size = frag_lambda_effective_genome;
    pr.llocal_bp = frag_lambda_llocal;
    pr.score_pseudocount = frag_score_pseudocount;
    chromap::peaks::Macs3FragPeakPipelinePaths io;
    io.treat_bdg = frag_span_pileup_bdg;
    io.lambda_bdg = frag_lambda_bdg;
    io.ppois_bdg = !macs3_frag_ppois_override.empty() ? macs3_frag_ppois_override
                                                         : frag_score_ppois_bdg;
    std::string err;
    std::string work_dir;
    if (!chromap::peaks::RunMacs3FragPeakPipelineFromFragments(
            chs, pr, io, macs3_frag_narrowpeak_out, macs3_frag_summits_out, "", "",
            &work_dir, &err, profile_coll.get())) {
      std::cerr << err << "\n";
      return 1;
    }
    if (profile_coll && !macs3_frag_narrowpeak_out.empty()) {
      int64_t nb = 0, nl = 0;
      (void)chromap::peaks::StageProfileCollector::FileMetrics(macs3_frag_narrowpeak_out, &nb,
                                                               &nl);
      profile_out_peaks = nl;
    }
    if (!frag_score_fe_bdg.empty()) {
      if (!chromap::peaks::WriteFragMacs3BdgcmpScoreBedGraphs(
              "", frag_score_fe_bdg, chs, frag_pileup_max_count, frag_pileup_macs3_uint8,
              frag_lambda_effective_genome, frag_lambda_llocal, frag_score_pseudocount)) {
        std::cerr << "Failed to write MACS3-style FE score bedGraph\n";
        return 1;
      }
    }
  } else {
    if (!frag_span_pileup_bdg.empty()) {
      if (!chromap::peaks::WriteFragSpanPileupBedGraph(
              frag_span_pileup_bdg, chs, frag_pileup_max_count, frag_pileup_macs3_uint8)) {
        std::cerr << "Failed to write fragment-span pileup bedGraph: "
                  << frag_span_pileup_bdg << "\n";
        return 1;
      }
    }
    if (!frag_lambda_bdg.empty()) {
      if (!chromap::peaks::WriteFragMacs3NoControlLambdaBedGraph(
              frag_lambda_bdg, chs, frag_pileup_max_count, frag_pileup_macs3_uint8,
              frag_lambda_effective_genome, frag_lambda_llocal)) {
        std::cerr << "Failed to write frag lambda bedGraph: " << frag_lambda_bdg << "\n";
        return 1;
      }
    }
    if (!frag_score_ppois_bdg.empty() || !frag_score_fe_bdg.empty()) {
      if (!chromap::peaks::WriteFragMacs3BdgcmpScoreBedGraphs(
              frag_score_ppois_bdg, frag_score_fe_bdg, chs, frag_pileup_max_count,
              frag_pileup_macs3_uint8, frag_lambda_effective_genome, frag_lambda_llocal,
              frag_score_pseudocount)) {
        std::cerr << "Failed to write MACS3-style score bedGraph(s)\n";
        return 1;
      }
    }
  }
  if (out_prefix.empty()) {
    return 0;
  }
  std::vector<chromap::peaks::Peak> all;
  for (const chromap::peaks::ChromFragments& ch : chs) {
    std::vector<int64_t> obs;
    int32_t max_c = 0;
    if (!chromap::peaks::BuildBinnedCutSignal(ch, bin_width, ext_size, &obs, &max_c)) {
      continue;
    }
    std::vector<chromap::peaks::Peak> out;
    chromap::peaks::CallPeaksOnBins(ch.name, bin_width, obs, local_window, fdr, pval,
                                merge_gap, &out);
    for (const auto& p : out) {
      all.push_back(p);
    }
  }
  if (min_peak_bp > 0 || min_summit_cuts > 0) {
    std::vector<chromap::peaks::Peak> keep;
    keep.reserve(all.size());
    for (const auto& p : all) {
      if (p.end <= p.start) {
        continue;
      }
      if (min_peak_bp > 0 && (p.end - p.start) < min_peak_bp) {
        continue;
      }
      if (min_summit_cuts > 0 && p.max_signal < min_summit_cuts) {
        continue;
      }
      keep.push_back(p);
    }
    all.swap(keep);
  }
  const std::string np = out_prefix + ".narrowPeak";
  const std::string sm = out_prefix + ".summary.tsv";
  if (!chromap::peaks::WriteNarrowPeak(np, all)) {
    std::cerr << "Failed to write " << np << "\n";
    return 1;
  }
  std::ostringstream pstr;
  pstr.setf(std::ios::fixed);
  pstr << "bin_width=" << bin_width << " ext_size=" << ext_size
       << " local_window_bins=" << local_window << " fdr=" << std::setprecision(4) << fdr
       << " p_value=" << pval << " merge_bin_gap=" << merge_gap
       << " min_peak_bp=" << min_peak_bp << " min_summit_cuts=" << min_summit_cuts
       << " input=" << input;
  if (!chromap::peaks::WriteRunSummaryTsv(sm, all, JoinArgv(argc, argv), pstr.str())) {
    std::cerr << "Failed to write " << sm << "\n";
    return 1;
  }
  return 0;
}
