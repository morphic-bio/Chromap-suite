#include "call_peaks.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <vector>

#include "htslib/kfunc.h"

namespace chromap {
namespace peaks {

namespace {

const double kEps = 1e-12;

}  // namespace

double PoissonSurvivalUpper(int k, double lambda) {
  if (k <= 0) {
    return 1.0;
  }
  if (lambda < kEps) {
    return k == 0 ? 1.0 : 0.0;
  }
  // P(Poisson(λ) >= k) = P(k, λ) regularized lower incomplete = kf_gammap(k, λ)
  // Verified against R ppois(k-1, lambda, lower.tail=FALSE) for k>=1.
  double p = kf_gammap(static_cast<double>(k), lambda);
  if (p < 0.0) {
    p = 0.0;
  }
  if (p > 1.0) {
    p = 1.0;
  }
  return p;
}

void BenjaminiHochbergFdr(const std::vector<double>& p_values,
                          std::vector<double>* q_values) {
  const size_t m = p_values.size();
  q_values->assign(m, 0.0);
  if (m == 0) {
    return;
  }
  std::vector<size_t> order(m);
  std::iota(order.begin(), order.end(), 0);
  std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
    if (p_values[a] == p_values[b]) {
      return a < b;
    }
    return p_values[a] < p_values[b];
  });
  std::vector<double> q(m);
  const double dmin = m + 0.0;
  for (int64_t k = static_cast<int64_t>(m) - 1; k >= 0; --k) {
    const size_t i = order[static_cast<size_t>(k)];
    const double raw = p_values[i] * dmin / static_cast<double>(k + 1);
    if (k + 1 < static_cast<int64_t>(m)) {
      const size_t nxt = order[static_cast<size_t>(k + 1)];
      q[static_cast<size_t>(i)] = std::min(std::min(raw, 1.0), q[nxt]);
    } else {
      q[static_cast<size_t>(i)] = std::min(raw, 1.0);
    }
  }
  for (size_t i = 0; i < m; ++i) {
    (*q_values)[i] = std::min(std::max(q[i], 0.0), 1.0);
  }
}

void ComputeLocalBackground(const std::vector<int64_t>& obs,
                            int32_t local_window_bins, double* global_per_bin,
                            std::vector<double>* lambda) {
  const int64_t n = static_cast<int64_t>(obs.size());
  lambda->assign(static_cast<size_t>(n), 0.0);
  *global_per_bin = 0.0;
  if (n < 1) {
    return;
  }
  int64_t s_all = 0;
  for (int64_t v : obs) {
    s_all += v;
  }
  *global_per_bin = static_cast<double>(s_all) / static_cast<double>(n);
  if (local_window_bins < 1) {
    for (int64_t i = 0; i < n; ++i) {
      (*lambda)[static_cast<size_t>(i)] =
          std::max(*global_per_bin, kEps);
    }
    return;
  }
  std::vector<int64_t> pref(static_cast<size_t>(n) + 1, 0);
  for (int64_t i = 0; i < n; ++i) {
    pref[static_cast<size_t>(i) + 1] =
        pref[static_cast<size_t>(i)] + obs[static_cast<size_t>(i)];
  }
  const int64_t w = static_cast<int64_t>(local_window_bins);
  for (int64_t i = 0; i < n; ++i) {
    int64_t lo = i - w;
    int64_t hi = i + w;
    if (lo < 0) {
      lo = 0;
    }
    if (hi >= n) {
      hi = n - 1;
    }
    int64_t sumw = pref[static_cast<size_t>(hi) + 1] - pref[static_cast<size_t>(lo)];
    int64_t sum_wo = sumw - obs[static_cast<size_t>(i)];
    int64_t cnt = (hi - lo);  // number of j != i in [lo,hi]
    double local = *global_per_bin;
    if (cnt > 0) {
      local = std::max(local, static_cast<double>(sum_wo) / static_cast<double>(cnt));
    }
    if (local < kEps) {
      local = kEps;
    }
    (*lambda)[static_cast<size_t>(i)] =
        std::max(*global_per_bin, local);
  }
}

int CallPeaksOnBins(const std::string& chrom, int32_t bin_width,
                    const std::vector<int64_t>& obs, int32_t local_window_bins,
                    double fdr, double p_cutoff, int32_t merge_bin_gap,
                    std::vector<Peak>* out_peaks) {
  out_peaks->clear();
  const int64_t n = static_cast<int64_t>(obs.size());
  if (n < 1) {
    return 0;
  }
  double g = 0.0;
  std::vector<double> lambda;
  ComputeLocalBackground(obs, local_window_bins, &g, &lambda);
  std::vector<double> p_values(static_cast<size_t>(n), 1.0);
  for (int64_t i = 0; i < n; ++i) {
    const int o = static_cast<int>(std::min<int64_t>(
        obs[static_cast<size_t>(i)],
        static_cast<int64_t>(std::numeric_limits<int>::max() / 2)));
    if (o < 1) {
      p_values[static_cast<size_t>(i)] = 1.0;
    } else {
      const double lam = lambda[static_cast<size_t>(i)];
      p_values[static_cast<size_t>(i)] = PoissonSurvivalUpper(o, lam);
    }
  }
  // BH only for bins with positive observed count to avoid m inflation.
  std::vector<double> p_sub;
  std::vector<int64_t> sub_ix;
  for (int64_t i = 0; i < n; ++i) {
    if (obs[static_cast<size_t>(i)] >= 1) {
      sub_ix.push_back(i);
      p_sub.push_back(p_values[static_cast<size_t>(i)]);
    }
  }
  std::vector<double> q_sub;
  BenjaminiHochbergFdr(p_sub, &q_sub);
  std::vector<double> q_values(static_cast<size_t>(n), 1.0);
  for (size_t j = 0; j < sub_ix.size(); ++j) {
    q_values[static_cast<size_t>(sub_ix[j])] = q_sub[j];
  }
  const int cap = 1 + std::max(0, merge_bin_gap);
  std::vector<int32_t> sig;
  for (int64_t i = 0; i < n; ++i) {
    if (q_values[static_cast<size_t>(i)] > fdr) {
      continue;
    }
    if (p_values[static_cast<size_t>(i)] > p_cutoff) {
      continue;
    }
    if (obs[static_cast<size_t>(i)] < 1) {
      continue;
    }
    sig.push_back(static_cast<int32_t>(i));
  }
  if (sig.empty()) {
    return 0;
  }
  // Merge: sorted sig; cluster when next - last <= cap.
  size_t a = 0;
  while (a < sig.size()) {
    size_t b = a;
    while (b + 1 < sig.size() &&
           sig[b + 1] - sig[b] <= static_cast<int32_t>(cap)) {
      ++b;
    }
    const int32_t b0 = sig[a];
    const int32_t b1 = sig[b];
    int32_t best = b0;
    int64_t best_v = -1;
    for (int32_t t = b0; t <= b1; ++t) {
      const int64_t v = obs[static_cast<size_t>(t)];
      if (v > best_v || (v == best_v && t < best)) {
        best_v = v;
        best = t;
      }
    }
    const int32_t start = b0 * bin_width;
    const int32_t end = (b1 + 1) * bin_width;
    const int32_t summit_0 = best * bin_width + std::max(0, bin_width / 2);
    Peak p;
    p.chrom = chrom;
    p.start = start;
    p.end = end;
    p.thick_start = start;
    p.thick_end = end;
    p.peak_offset = summit_0 - start;
    p.max_signal_bin = best;
    p.max_signal = best_v;
    p.p_value = p_values[static_cast<size_t>(best)];
    p.q_value = q_values[static_cast<size_t>(best)];
    out_peaks->push_back(p);
    a = b + 1;
  }
  return static_cast<int>(out_peaks->size());
}

}  // namespace peaks
}  // namespace chromap
