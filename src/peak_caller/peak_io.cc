#include "peak_io.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace chromap {
namespace peaks {

namespace {

int ScoreFromQ(double q) {
  if (q < 0.0) {
    q = 0.0;
  }
  if (q > 1.0) {
    q = 1.0;
  }
  const double s = -10.0 * std::log10(std::max(q, 1e-300));
  int out = static_cast<int>(s + 0.5);
  if (out < 0) {
    out = 0;
  }
  if (out > 1000) {
    out = 1000;
  }
  return out;
}

}  // namespace

bool WriteNarrowPeak(const std::string& path, const std::vector<Peak>& peaks) {
  std::vector<Peak> sorted = peaks;
  std::sort(sorted.begin(), sorted.end(), [](const Peak& a, const Peak& b) {
    if (a.chrom != b.chrom) {
      return a.chrom < b.chrom;
    }
    if (a.start != b.start) {
      return a.start < b.start;
    }
    if (a.end != b.end) {
      return a.end < b.end;
    }
    return a.peak_offset < b.peak_offset;
  });
  for (const Peak& p : sorted) {
    if (p.end <= p.start) {
      return false;
    }
  }
  for (size_t i = 1; i < sorted.size(); ++i) {
    if (sorted[i - 1].chrom == sorted[i].chrom) {
      if (sorted[i - 1].start > sorted[i].start) {
        return false;
      }
    } else if (sorted[i - 1].chrom > sorted[i].chrom) {
      return false;
    }
  }
  FILE* fp = std::fopen(path.c_str(), "w");
  if (fp == nullptr) {
    return false;
  }
  for (size_t i = 0; i < sorted.size(); ++i) {
    const Peak& p = sorted[i];
    const std::string name = "peak_" + std::to_string(i + 1);
    const int sc = ScoreFromQ(p.q_value);
    const double sfg = static_cast<double>(p.max_signal);
    const double pv = std::min(1.0, std::max(0.0, p.p_value));
    const double qv = std::min(1.0, std::max(0.0, p.q_value));
    const double logp = -std::log10(std::max(pv, 1e-300));
    const double logq = -std::log10(std::max(qv, 1e-300));
    if (std::fprintf(fp, "%s\t%d\t%d\t%s\t%d\t.\t%g\t%g\t%g\t%d\n", p.chrom.c_str(),
                     p.start, p.end, name.c_str(), sc, sfg, logp, logq,
                     p.peak_offset) < 0) {
      std::fclose(fp);
      return false;
    }
  }
  if (std::fclose(fp) != 0) {
    return false;
  }
  return true;
}

bool WriteRunSummaryTsv(const std::string& path,
                        const std::vector<Peak>& all_peaks,
                        const std::string& argv_line,
                        const std::string& program_params) {
  std::ofstream o(path);
  if (!o) {
    return false;
  }
  o << "metric\tvalue\n";
  o << "peak_count\t" << all_peaks.size() << "\n";
  int64_t tot = 0;
  std::vector<int32_t> widths;
  for (const Peak& p : all_peaks) {
    if (p.end > p.start) {
      const int w = p.end - p.start;
      tot += w;
      widths.push_back(w);
    }
  }
  o << "total_peak_bp\t" << tot << "\n";
  if (widths.empty()) {
    o << "median_peak_width_bp\t0\n";
  } else {
    std::sort(widths.begin(), widths.end());
    const size_t m = widths.size() / 2;
    o << "median_peak_width_bp\t" << widths[m] << "\n";
  }
  o << "argv\t" << argv_line << "\n";
  o << "params\t" << program_params << "\n";
  o.close();
  if (!o) {
    return false;
  }
  return true;
}

}  // namespace peaks
}  // namespace chromap
