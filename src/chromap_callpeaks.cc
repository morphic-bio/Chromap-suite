// Standalone fragment-based peak caller (MACS3-inspired baseline; not a full
// MACS3 reimplementation). Next phase: optional pipeline integration in Chromap/STAR.
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "peak_caller/binned_signal.h"
#include "peak_caller/call_peaks.h"
#include "peak_caller/fragment_input.h"
#include "peak_caller/peak_io.h"

namespace {

void Usage() {
  std::cerr
      << "Usage: chromap_callpeaks --input <fragments.tsv|fragments.tsv.gz> "
         "--out-prefix <path/prefix> [options]\n"
         "  --input|-i     Chromap/ARC 5-col TSV(/.gz): chrom start end barcode count\n"
         "  --out-prefix   Prefix for <prefix>.narrowPeak and <prefix>.summary.tsv\n"
         "  --bin-width    Bin width in bp (default: 50)\n"
         "  --ext-size     Tn5 extension in bp, symmetric; default 150 (2x75; shiftTAG-like)\n"
         "  --local-window  Half-window in bins for local lambda (default: 200; ~20kb for 50bp)\n"
         "  --fdr          FDR (BH) threshold (default: 0.05)\n"
         "  --p-value      Raw p upper tail threshold (default: 0.01)\n"
         "  --merge-gap     Merge significant bins separated by up to this many bins (0)\n"
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
  int32_t bin_width = 50;
  int32_t ext_size = 150;
  int32_t local_window = 200;
  double fdr = 0.05;
  double pval = 0.01;
  int32_t merge_gap = 0;
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
    } else {
      std::cerr << "Unknown argument: " << a << "\n";
      Usage();
      return 2;
    }
  }
  if (input.empty() || out_prefix.empty()) {
    Usage();
    return 2;
  }
  if (bin_width < 1 || ext_size < 2 || (ext_size % 2) != 0) {
    std::cerr
        << "Invalid --bin-width or --ext-size (ext-size must be even, >=2)\n";
    return 2;
  }
  if (fdr < 0.0 || fdr > 1.0 || pval < 0.0 || pval > 1.0) {
    return 2;
  }
  std::vector<chromap::peaks::ChromFragments> chs;
  if (!chromap::peaks::LoadFragmentsFromTsv(input, &chs)) {
    std::cerr << "Failed to read input: " << input << "\n";
    return 1;
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
       << " input=" << input;
  if (!chromap::peaks::WriteRunSummaryTsv(sm, all, JoinArgv(argc, argv), pstr.str())) {
    std::cerr << "Failed to write " << sm << "\n";
    return 1;
  }
  return 0;
}
