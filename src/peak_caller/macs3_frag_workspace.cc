#include "macs3_frag_workspace.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <tuple>
#include <unordered_map>

#include "peak_io.h"

namespace chromap {
namespace peaks {

namespace {

float Quantize5Float(double x) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.5f", x);
  return std::strtof(buf, nullptr);
}

double Round5Score(double x) {
  return std::round(x * 1e5) / 1e5;
}

int32_t EffectiveFragWeight(int32_t raw_count, int32_t max_count,
                            bool macs3_uint8_counts) {
  int64_t w = raw_count;
  if (max_count > 0) {
    w = std::min<int64_t>(w, max_count);
  }
  if (w <= 0) {
    return 0;
  }
  if (!macs3_uint8_counts) {
    if (w > std::numeric_limits<int32_t>::max()) {
      return std::numeric_limits<int32_t>::max();
    }
    return static_cast<int32_t>(w);
  }
  const auto u = static_cast<uint8_t>(static_cast<unsigned int>(w));
  return static_cast<int32_t>(u);
}

int32_t RowCountForMacsTotals(int32_t raw_count, int32_t max_count,
                              bool macs3_uint8_counts) {
  int64_t t = raw_count;
  if (max_count > 0) {
    t = std::min<int64_t>(t, max_count);
  }
  if (t <= 0) {
    return 0;
  }
  if (macs3_uint8_counts) {
    return static_cast<int32_t>(static_cast<uint8_t>(static_cast<unsigned int>(t)));
  }
  if (t > std::numeric_limits<int32_t>::max()) {
    return std::numeric_limits<int32_t>::max();
  }
  return static_cast<int32_t>(t);
}

void ReleaseEvents(std::vector<Macs3FragEvent>* ev) {
  std::vector<Macs3FragEvent>().swap(*ev);
}

int NormalizeThreads(int requested, size_t n_tasks) {
  if (requested < 1 || n_tasks <= 1) {
    return 1;
  }
  const size_t rt = static_cast<size_t>(requested);
  return static_cast<int>(std::min(rt, n_tasks));
}

bool CoalesceEventPositions(std::vector<Macs3FragEvent>* ev) {
  if (ev->empty()) {
    return true;
  }
  std::sort(ev->begin(), ev->end(),
            [](const Macs3FragEvent& a, const Macs3FragEvent& b) {
              if (a.pos != b.pos) {
                return a.pos < b.pos;
              }
              return a.delta < b.delta;
            });
  size_t wri = 0;
  for (size_t i = 0; i < ev->size(); ++i) {
    if (wri > 0 && ev->at(i).pos == ev->at(wri - 1).pos) {
      const int64_t delta = static_cast<int64_t>(ev->at(wri - 1).delta) +
                            static_cast<int64_t>(ev->at(i).delta);
      if (delta < std::numeric_limits<int32_t>::min() ||
          delta > std::numeric_limits<int32_t>::max()) {
        return false;
      }
      ev->at(wri - 1).delta = static_cast<int32_t>(delta);
    } else {
      ev->at(wri++) = ev->at(i);
    }
  }
  ev->resize(wri);
  return true;
}

void TrackFromEvents(const std::vector<Macs3FragEvent>& ev, Macs3BdgTrack* out) {
  out->p.clear();
  out->v.clear();
  if (ev.empty()) {
    return;
  }
  double z = 0.0;
  double pre_z = -1.0e4;
  int32_t s = 0;
  out->p.reserve(ev.size());
  out->v.reserve(ev.size());
  for (const auto& e : ev) {
    const int32_t pos = e.pos;
    const double v = static_cast<double>(e.delta);
    if (pos != s) {
      if (z == pre_z && !out->p.empty()) {
        out->p.back() = pos;
      } else {
        out->p.push_back(pos);
        out->v.push_back(Quantize5Float(z));
        pre_z = z;
      }
    }
    z += v;
    s = pos;
  }
}

void BuildLambdaTrackFromTreatAndControl(const std::vector<int32_t>& t_p,
                                         const Macs3BdgTrack& ctrl,
                                         Macs3BdgTrack* out) {
  out->p.clear();
  out->v.clear();
  if (t_p.empty() || ctrl.p.empty() || ctrl.p.size() != ctrl.v.size()) {
    return;
  }
  size_t it = 0;
  size_t ic = 0;
  bool have_prev = false;
  int32_t prev_p = 0;
  float pre_v = 0.f;
  out->p.reserve(t_p.size() + ctrl.p.size());
  out->v.reserve(t_p.size() + ctrl.p.size());
  while (it < t_p.size() && ic < ctrl.p.size()) {
    int32_t cur_p = 0;
    const float cur_v = ctrl.v[ic];
    if (t_p[it] < ctrl.p[ic]) {
      cur_p = t_p[it];
      ++it;
    } else if (t_p[it] > ctrl.p[ic]) {
      cur_p = ctrl.p[ic];
      ++ic;
    } else {
      cur_p = t_p[it];
      ++it;
      ++ic;
    }
    if (!have_prev) {
      pre_v = cur_v;
      prev_p = cur_p;
      have_prev = true;
      continue;
    }
    if (std::fabs(static_cast<double>(pre_v) - static_cast<double>(cur_v)) > 1e-5) {
      out->p.push_back(prev_p);
      out->v.push_back(Quantize5Float(pre_v));
      pre_v = cur_v;
    }
    prev_p = cur_p;
  }
  if (have_prev) {
    out->p.push_back(prev_p);
    out->v.push_back(Quantize5Float(pre_v));
  }
}

void ScoreRunPush(Macs3BdgTrack* out, int32_t seg_lo, int32_t seg_hi,
                  double val, bool* has, int32_t* lo, int32_t* hi, double* last) {
  if (seg_hi <= seg_lo) {
    return;
  }
  if (!*has) {
    *has = true;
    *lo = seg_lo;
    *hi = seg_hi;
    *last = val;
    return;
  }
  if (std::fabs(val - *last) <= 1e-5 && seg_lo == *hi) {
    *hi = seg_hi;
    return;
  }
  out->p.push_back(*hi);
  out->v.push_back(Quantize5Float(*last));
  *lo = seg_lo;
  *hi = seg_hi;
  *last = val;
}

void ScoreRunFlush(Macs3BdgTrack* out, bool* has, int32_t* lo, int32_t* hi,
                   double* last) {
  (void)lo;
  if (!*has) {
    return;
  }
  out->p.push_back(*hi);
  out->v.push_back(Quantize5Float(*last));
  *has = false;
}

uint32_t FloatBits(float v) {
  uint32_t out = 0;
  std::memcpy(&out, &v, sizeof(out));
  return out;
}

struct PpoisScoreKey {
  uint32_t treat = 0;
  uint32_t lambda = 0;

  bool operator==(const PpoisScoreKey& other) const {
    return treat == other.treat && lambda == other.lambda;
  }
};

struct PpoisScoreKeyHash {
  size_t operator()(const PpoisScoreKey& k) const {
    return (static_cast<size_t>(k.treat) * 1315423911u) ^
           static_cast<size_t>(k.lambda);
  }
};

}  // namespace

bool BuildMacs3FragWorkspaceFromFragments(
    const std::vector<ChromFragments>& by_chrom,
    const Macs3FragWorkspaceParams& params,
    Macs3FragPeakWorkspace* workspace) {
  if (workspace == nullptr || params.frag_pileup_max_count < 0 ||
      params.effective_genome_size < 1 || params.llocal_bp < 1) {
    return false;
  }
  std::vector<std::string> chrom_names;
  chrom_names.reserve(by_chrom.size());
  for (const ChromFragments& ch : by_chrom) {
    chrom_names.push_back(ch.name);
  }
  if (!InitMacs3FragWorkspace(chrom_names, params, workspace)) {
    return false;
  }
  for (size_t i = 0; i < by_chrom.size(); ++i) {
    const ChromFragments& ch = by_chrom[i];
    workspace->events_by_chrom[i].treat_span_events.reserve(ch.frags.size() * 2);
    workspace->events_by_chrom[i].lambda_end_window_events.reserve(ch.frags.size() * 4);
    for (const Fragment& f : ch.frags) {
      if (!AddMacs3FragWorkspaceFragment(
              workspace, static_cast<uint32_t>(i), f.start, f.end, f.count)) {
        return false;
      }
    }
  }
  return true;
}

bool InitMacs3FragWorkspace(const std::vector<std::string>& chrom_names,
                            const Macs3FragWorkspaceParams& params,
                            Macs3FragPeakWorkspace* workspace) {
  if (workspace == nullptr || params.frag_pileup_max_count < 0 ||
      params.effective_genome_size < 1 || params.llocal_bp < 1) {
    return false;
  }
  workspace->params = params;
  workspace->chrom_names = chrom_names;
  workspace->events_by_chrom.clear();
  workspace->treat_tracks.clear();
  workspace->lambda_tracks.clear();
  workspace->ppois_tracks.clear();
  workspace->total_fragment_bp = 0;
  workspace->total_fragment_count = 0;
  workspace->events_by_chrom.resize(chrom_names.size());
  return true;
}

bool AddMacs3FragWorkspaceFragment(Macs3FragPeakWorkspace* workspace,
                                   uint32_t chrom_id, int32_t start,
                                   int32_t end, int32_t count) {
  if (workspace == nullptr || chrom_id >= workspace->events_by_chrom.size() ||
      end <= start || count <= 0) {
    return false;
  }
  const Macs3FragWorkspaceParams& params = workspace->params;
  const int32_t w = EffectiveFragWeight(count, params.frag_pileup_max_count,
                                        params.macs3_uint8_counts);
  if (w == 0) {
    return true;
  }
  Macs3FragEventBuffer& ev = workspace->events_by_chrom[chrom_id];
  ev.treat_span_events.push_back(Macs3FragEvent{start, w});
  ev.treat_span_events.push_back(Macs3FragEvent{end, -w});

  const int32_t half = params.llocal_bp / 2;
  const int32_t ll = start - half;
  const int32_t rl = ll + params.llocal_bp;
  const int32_t lr = end - half;
  const int32_t rr = lr + params.llocal_bp;
  ev.lambda_end_window_events.push_back(Macs3FragEvent{ll, w});
  ev.lambda_end_window_events.push_back(Macs3FragEvent{rl, -w});
  ev.lambda_end_window_events.push_back(Macs3FragEvent{lr, w});
  ev.lambda_end_window_events.push_back(Macs3FragEvent{rr, -w});

  const int32_t c = RowCountForMacsTotals(count, params.frag_pileup_max_count,
                                          params.macs3_uint8_counts);
  if (c > 0) {
    workspace->total_fragment_count += c;
    workspace->total_fragment_bp += static_cast<int64_t>(end - start) * c;
  }
  return true;
}

bool FinalizeMacs3FragTreatTracks(Macs3FragPeakWorkspace* workspace,
                                   int num_threads) {
  if (workspace == nullptr) {
    return false;
  }
  workspace->treat_tracks.clear();
  workspace->treat_tracks.resize(workspace->events_by_chrom.size());
  const size_t n = workspace->events_by_chrom.size();
  std::vector<char> ok(n, 1);
  const int threads = NormalizeThreads(num_threads, n);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
    const size_t i = static_cast<size_t>(ii);
    std::vector<Macs3FragEvent>& ev = workspace->events_by_chrom[i].treat_span_events;
    if (!CoalesceEventPositions(&ev)) {
      ok[i] = 0;
      continue;
    }
    TrackFromEvents(ev, &workspace->treat_tracks[i]);
    ReleaseEvents(&ev);
  }
  return std::find(ok.begin(), ok.end(), 0) == ok.end();
}

bool FinalizeMacs3FragLambdaTracks(Macs3FragPeakWorkspace* workspace,
                                    int num_threads) {
  if (workspace == nullptr || workspace->total_fragment_count <= 0 ||
      workspace->total_fragment_bp <= 0 ||
      workspace->treat_tracks.size() != workspace->events_by_chrom.size()) {
    return false;
  }
  const double lambda_bg =
      static_cast<double>(workspace->total_fragment_bp) /
      static_cast<double>(workspace->params.effective_genome_size);
  const double scale =
      static_cast<double>(workspace->total_fragment_bp) /
      (static_cast<double>(workspace->params.llocal_bp) *
       static_cast<double>(workspace->total_fragment_count) * 2.0);

  workspace->lambda_tracks.clear();
  workspace->lambda_tracks.resize(workspace->events_by_chrom.size());
  const size_t n = workspace->events_by_chrom.size();
  std::vector<char> ok(n, 1);
  const int threads = NormalizeThreads(num_threads, n);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t ii = 0; ii < static_cast<int64_t>(n); ++ii) {
    const size_t i = static_cast<size_t>(ii);
    std::vector<Macs3FragEvent>& ev = workspace->events_by_chrom[i].lambda_end_window_events;
    if (!CoalesceEventPositions(&ev)) {
      ok[i] = 0;
      continue;
    }
    Macs3BdgTrack ctrl;
    TrackFromEvents(ev, &ctrl);
    ReleaseEvents(&ev);
    for (float& v : ctrl.v) {
      double z = static_cast<double>(v) * scale;
      if (z < lambda_bg) {
        z = lambda_bg;
      }
      v = static_cast<float>(z);
    }
    BuildLambdaTrackFromTreatAndControl(workspace->treat_tracks[i].p, ctrl,
                                        &workspace->lambda_tracks[i]);
  }
  return std::find(ok.begin(), ok.end(), 0) == ok.end();
}

bool FinalizeMacs3FragPpoisTracks(Macs3FragPeakWorkspace* workspace,
                                   int num_threads) {
  if (workspace == nullptr ||
      workspace->treat_tracks.size() != workspace->lambda_tracks.size()) {
    return false;
  }
  workspace->ppois_tracks.clear();
  workspace->ppois_tracks.resize(workspace->treat_tracks.size());
  const size_t n = workspace->treat_tracks.size();
  std::vector<char> ok(n, 1);
  const int threads = NormalizeThreads(num_threads, n);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
  for (int64_t cc = 0; cc < static_cast<int64_t>(n); ++cc) {
    const size_t c = static_cast<size_t>(cc);
    if (!BuildMacs3FragPpoisTrack(workspace->treat_tracks[c],
                                  workspace->lambda_tracks[c],
                                  workspace->params.score_pseudocount,
                                  &workspace->ppois_tracks[c])) {
      ok[c] = 0;
    }
  }
  return std::find(ok.begin(), ok.end(), 0) == ok.end();
}

bool BuildMacs3FragPpoisTrack(const Macs3BdgTrack& treat,
                              const Macs3BdgTrack& lambda,
                              double pseudocount, Macs3BdgTrack* out) {
  if (out == nullptr || treat.p.size() != treat.v.size() ||
      lambda.p.size() != lambda.v.size()) {
    return false;
  }
  out->p.clear();
  out->v.clear();
  if (treat.p.empty() || lambda.p.empty()) {
    return true;
  }
  const float pseudocount_f = static_cast<float>(pseudocount);
  std::unordered_map<PpoisScoreKey, double, PpoisScoreKeyHash> score_cache;
  size_t i1 = 0;
  size_t i2 = 0;
  int32_t prev = 0;
  bool has = false;
  int32_t lo = 0;
  int32_t hi = 0;
  double last = 0.0;
  while (i1 < treat.p.size() && i2 < lambda.p.size()) {
    const int32_t p1 = treat.p[i1];
    const int32_t p2 = lambda.p[i2];
    const int32_t nx = std::min(p1, p2);
    if (nx > prev) {
      const int32_t seg_lo = std::max(prev, static_cast<int32_t>(0));
      const int32_t seg_hi = std::max(nx, static_cast<int32_t>(0));
      if (seg_hi > seg_lo) {
        PpoisScoreKey key;
        key.treat = FloatBits(treat.v[i1]);
        key.lambda = FloatBits(lambda.v[i2]);
        auto it = score_cache.find(key);
        double sp = 0.0;
        if (it == score_cache.end()) {
          sp = Round5Score(Macs3BdgcmpPpoisNegLog10P(
              treat.v[i1], lambda.v[i2], static_cast<double>(pseudocount_f)));
          score_cache.emplace(key, sp);
        } else {
          sp = it->second;
        }
        ScoreRunPush(out, seg_lo, seg_hi, sp, &has, &lo, &hi, &last);
      }
    }
    if (p1 <= p2) {
      ++i1;
    }
    if (p2 <= p1) {
      ++i2;
    }
    prev = nx;
  }
  ScoreRunFlush(out, &has, &lo, &hi, &last);
  return true;
}

bool WriteMacs3BdgTracksBedGraph(const std::string& path,
                                 const std::vector<std::string>& chrom_names,
                                 const std::vector<Macs3BdgTrack>& tracks) {
  if (path.empty() || chrom_names.size() != tracks.size()) {
    return false;
  }
  FILE* fp = std::fopen(path.c_str(), "w");
  if (fp == nullptr) {
    return false;
  }
  for (size_t i = 0; i < tracks.size(); ++i) {
    const Macs3BdgTrack& tr = tracks[i];
    int32_t start = 0;
    for (size_t j = 0; j < tr.p.size(); ++j) {
      const int32_t end = tr.p[j];
      if (std::fprintf(fp, "%s\t%d\t%d\t%.5f\n", chrom_names[i].c_str(), start,
                       end, static_cast<double>(tr.v[j])) < 0) {
        std::fclose(fp);
        return false;
      }
      start = end;
    }
  }
  if (std::fclose(fp) != 0) {
    return false;
  }
  return true;
}

}  // namespace peaks
}  // namespace chromap
