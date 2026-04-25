#include "bdgpeakcall.h"

#include "peak_caller/call_peaks.h"
#include "peak_caller/peak_io.h"
#include "peak_caller/stage_profile.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace chromap {
namespace peaks {

namespace {

bool StartsWith(const std::string& s, const char* prefix) {
  const size_t n = std::strlen(prefix);
  return s.size() >= n && s.compare(0, n, prefix) == 0;
}

void AddLocMacs3(std::map<std::string, Macs3BdgTrack>* by_chrom,
                 const std::string& chromosome, int32_t startpos, int32_t endpos,
                 float value, float baseline_value) {
  if (endpos <= 0) {
    return;
  }
  if (startpos < 0) {
    startpos = 0;
  }
  Macs3BdgTrack& ch = (*by_chrom)[chromosome];
  if (ch.p.empty()) {
    if (startpos > 0) {
      ch.p.push_back(startpos);
      ch.v.push_back(baseline_value);
    }
    ch.p.push_back(endpos);
    ch.v.push_back(value);
    return;
  }
  const float pre_v = ch.v.back();
  if (pre_v == value) {
    ch.p.back() = endpos;
  } else {
    ch.p.push_back(endpos);
    ch.v.push_back(value);
  }
}

bool ReadBedGraphMacs3Style(const std::string& path, float baseline_value,
                            std::map<std::string, Macs3BdgTrack>* out) {
  std::ifstream in(path.c_str());
  if (!in) {
    return false;
  }
  out->clear();
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty()) {
      continue;
    }
    if (StartsWith(line, "track") || StartsWith(line, "#") ||
        StartsWith(line, "browse")) {
      continue;
    }
    std::istringstream iss(line);
    std::string chrom;
    std::string s1, s2, s3;
    if (!(iss >> chrom >> s1 >> s2 >> s3)) {
      continue;
    }
    const int32_t start = static_cast<int32_t>(std::strtol(s1.c_str(), nullptr, 10));
    const int32_t end = static_cast<int32_t>(std::strtol(s2.c_str(), nullptr, 10));
    const float val = static_cast<float>(std::strtod(s3.c_str(), nullptr));
    if (end <= start) {
      continue;
    }
    AddLocMacs3(out, chrom, start, end, val, baseline_value);
  }
  return in.eof() || !in.bad();
}

struct PeakPiece {
  int32_t start = 0;
  int32_t end = 0;
  float value = 0.f;
};

void ClosePeakMacs3(const std::vector<PeakPiece>& peak_content, int32_t min_length,
                    std::vector<std::pair<int32_t, int32_t>>* out_ivs) {
  if (peak_content.empty()) {
    return;
  }
  const int32_t peak_length = peak_content.back().end - peak_content.front().start;
  if (peak_length < min_length) {
    return;
  }
  out_ivs->emplace_back(peak_content.front().start, peak_content.back().end);
}

// Mirrors MACS3 Signal/BedGraph.py bedGraphTrackI.call_peaks (call_summits ignored
// for interval bounds; summits are not written in diagnostic BED3 mode).
void CallPeaksOnTrack(const Macs3BdgTrack& tr, float cutoff, int32_t min_length,
                      int32_t max_gap,
                      std::vector<std::pair<int32_t, int32_t>>* out_ivs) {
  out_ivs->clear();
  const size_t n = tr.v.size();
  if (n == 0 || tr.p.size() != n) {
    return;
  }
  std::vector<char> above(n);
  size_t n_above = 0;
  for (size_t i = 0; i < n; ++i) {
    above[i] = (tr.v[i] >= cutoff) ? 1 : 0;
    if (above[i]) {
      ++n_above;
    }
  }
  if (n_above == 0) {
    return;
  }
  std::vector<int32_t> endpos;
  std::vector<int32_t> startpos;
  std::vector<float> vals;
  endpos.reserve(n_above);
  startpos.reserve(n_above);
  vals.reserve(n_above);
  for (size_t i = 0; i < n; ++i) {
    if (!above[i]) {
      continue;
    }
    endpos.push_back(tr.p[i]);
    vals.push_back(tr.v[i]);
  }
  for (size_t i = 0; i < n; ++i) {
    const size_t j = (i + 1 < n) ? (i + 1) : 0;
    if (!above[j]) {
      continue;
    }
    startpos.push_back(tr.p[i]);
  }
  if (startpos.size() != endpos.size()) {
    return;
  }
  if (above[0]) {
    startpos[0] = 0;
  }
  std::vector<PeakPiece> peak_content;
  peak_content.reserve(endpos.size());
  {
    PeakPiece pc;
    pc.start = startpos[0];
    pc.end = endpos[0];
    pc.value = vals[0];
    peak_content.push_back(pc);
  }
  for (size_t i = 1; i < startpos.size(); ++i) {
    if (startpos[i] - peak_content.back().end <= max_gap) {
      PeakPiece pc;
      pc.start = startpos[i];
      pc.end = endpos[i];
      pc.value = vals[i];
      peak_content.push_back(pc);
    } else {
      ClosePeakMacs3(peak_content, min_length, out_ivs);
      peak_content.clear();
      PeakPiece pc;
      pc.start = startpos[i];
      pc.end = endpos[i];
      pc.value = vals[i];
      peak_content.push_back(pc);
    }
  }
  if (!peak_content.empty()) {
    ClosePeakMacs3(peak_content, min_length, out_ivs);
  }
}

}  // namespace

bool RunMacs3StyleBdgPeakCallFromBedGraph(const std::string& path_bdg,
                                          const std::string& path_bed_out,
                                          float cutoff, int32_t min_length,
                                          int32_t max_gap) {
  if (path_bdg.empty() || path_bed_out.empty() || min_length < 1 || max_gap < 0) {
    return false;
  }
  std::map<std::string, Macs3BdgTrack> tracks;
  if (!ReadBedGraphMacs3Style(path_bdg, 0.f, &tracks)) {
    return false;
  }
  std::vector<std::string> chroms;
  chroms.reserve(tracks.size());
  for (const auto& kv : tracks) {
    chroms.push_back(kv.first);
  }
  std::sort(chroms.begin(), chroms.end());

  std::vector<std::pair<std::string, std::pair<int32_t, int32_t>>> all;
  for (const std::string& c : chroms) {
    const Macs3BdgTrack& tr = tracks[c];
    std::vector<std::pair<int32_t, int32_t>> ivs;
    CallPeaksOnTrack(tr, cutoff, min_length, max_gap, &ivs);
    for (const auto& iv : ivs) {
      all.emplace_back(c, iv);
    }
  }
  std::sort(all.begin(), all.end(),
            [](const std::pair<std::string, std::pair<int32_t, int32_t>>& a,
               const std::pair<std::string, std::pair<int32_t, int32_t>>& b) {
              if (a.first != b.first) {
                return a.first < b.first;
              }
              if (a.second.first != b.second.first) {
                return a.second.first < b.second.first;
              }
              return a.second.second < b.second.second;
            });

  FILE* fp = std::fopen(path_bed_out.c_str(), "w");
  if (fp == nullptr) {
    return false;
  }
  for (const auto& e : all) {
    if (e.second.second <= e.second.first) {
      std::fclose(fp);
      return false;
    }
    if (std::fprintf(fp, "%s\t%d\t%d\n", e.first.c_str(), static_cast<int>(e.second.first),
                     static_cast<int>(e.second.second)) < 0) {
      std::fclose(fp);
      return false;
    }
  }
  if (std::fclose(fp) != 0) {
    return false;
  }
  return true;
}

// --- MACS3-style narrowPeak: bdgpeakcall regions on ppois, summit = max(treat) chunk ---

struct PpoisMergePiece {
  int32_t start = 0;
  int32_t end = 0;
  float ppois = 0.f;
};

struct PpoisMergedRegion {
  int32_t start = 0;
  int32_t end = 0;
  std::vector<PpoisMergePiece> pieces;
};

float BdgValueAtOrZero(const Macs3BdgTrack& tr, int32_t pos) {
  if (tr.p.empty() || tr.v.size() != tr.p.size()) {
    return 0.f;
  }
  const auto it = std::upper_bound(tr.p.begin(), tr.p.end(), pos);
  if (it == tr.p.end()) {
    return 0.f;
  }
  const size_t j = static_cast<size_t>(it - tr.p.begin());
  const int32_t left = (j == 0) ? 0 : tr.p[j - 1];
  return (pos >= left) ? tr.v[j] : 0.f;
}

// Same as CallPeaksOnTrack, but also records each above-cutoff segment as a PpoisMergePiece
// (genomic [start, end) from bdg, value = ppois for that run).
void CallPpoisTrackMergedRegions(
    const Macs3BdgTrack& tr, float cutoff, int32_t min_length, int32_t max_gap,
    std::vector<PpoisMergedRegion>* out) {
  out->clear();
  const size_t n = tr.v.size();
  if (n == 0 || tr.p.size() != n) {
    return;
  }
  std::vector<char> above(n);
  size_t n_above = 0;
  for (size_t i = 0; i < n; ++i) {
    above[i] = (tr.v[i] >= cutoff) ? 1 : 0;
    if (above[i]) {
      ++n_above;
    }
  }
  if (n_above == 0) {
    return;
  }
  std::vector<int32_t> endpos;
  std::vector<int32_t> startpos;
  std::vector<float> vals;
  endpos.reserve(n_above);
  startpos.reserve(n_above);
  vals.reserve(n_above);
  for (size_t i = 0; i < n; ++i) {
    if (!above[i]) {
      continue;
    }
    endpos.push_back(tr.p[i]);
    vals.push_back(tr.v[i]);
  }
  for (size_t i = 0; i < n; ++i) {
    const size_t j = (i + 1 < n) ? (i + 1) : 0;
    if (!above[j]) {
      continue;
    }
    startpos.push_back(tr.p[i]);
  }
  if (startpos.size() != endpos.size()) {
    return;
  }
  if (above[0]) {
    startpos[0] = 0;
  }
  struct PeakPiece {
    int32_t start = 0;
    int32_t end = 0;
    float value = 0.f;
  };
  auto close_peak = [&](const std::vector<PeakPiece>& pc,
                          std::vector<PpoisMergedRegion>* dst) {
    if (pc.empty()) {
      return;
    }
    const int32_t peak_length = pc.back().end - pc.front().start;
    if (peak_length < min_length) {
      return;
    }
    PpoisMergedRegion r;
    r.start = pc.front().start;
    r.end = pc.back().end;
    r.pieces.reserve(pc.size());
    for (const auto& a : pc) {
      PpoisMergePiece e;
      e.start = a.start;
      e.end = a.end;
      e.ppois = a.value;
      r.pieces.push_back(e);
    }
    dst->push_back(std::move(r));
  };
  std::vector<PeakPiece> peak_content;
  peak_content.reserve(endpos.size());
  {
    PeakPiece pc;
    pc.start = startpos[0];
    pc.end = endpos[0];
    pc.value = vals[0];
    peak_content.push_back(pc);
  }
  for (size_t i = 1; i < startpos.size(); ++i) {
    if (startpos[i] - peak_content.back().end <= max_gap) {
      PeakPiece pc;
      pc.start = startpos[i];
      pc.end = endpos[i];
      pc.value = vals[i];
      peak_content.push_back(pc);
    } else {
      close_peak(peak_content, out);
      peak_content.clear();
      PeakPiece pc;
      pc.start = startpos[i];
      pc.end = endpos[i];
      pc.value = vals[i];
      peak_content.push_back(pc);
    }
  }
  if (!peak_content.empty()) {
    close_peak(peak_content, out);
  }
}

static int NarrowPeakIntScoreQ(double neg_log10_q) {
  if (neg_log10_q < 0.0) {
    return 0;
  }
  int s = static_cast<int>(10.0 * neg_log10_q + 0.5);
  if (s < 0) {
    s = 0;
  }
  if (s > 1000) {
    s = 1000;
  }
  return s;
}

bool RunMacs3FragPpoisNarrowPeaks(
    const std::string& path_ppois_bdg, const std::string& path_treat_bdg,
    const std::string& path_lambda_bdg, float cutoff, int32_t min_length,
    int32_t max_gap, double pseudocount, const std::string& path_narrowpeak,
    const std::string& path_summits_bed, StageProfileCollector* stage_profile,
    int64_t profile_input_frag_rows) {
  if (path_ppois_bdg.empty() || path_treat_bdg.empty() || path_lambda_bdg.empty() ||
      (path_narrowpeak.empty() && path_summits_bed.empty()) || min_length < 1 ||
      max_gap < 0) {
    return false;
  }
  std::map<std::string, Macs3BdgTrack> m_p;
  std::map<std::string, Macs3BdgTrack> m_t;
  std::map<std::string, Macs3BdgTrack> m_l;
  if (!ReadBedGraphMacs3Style(path_ppois_bdg, 0.f, &m_p) ||
      !ReadBedGraphMacs3Style(path_treat_bdg, 0.f, &m_t) ||
      !ReadBedGraphMacs3Style(path_lambda_bdg, 0.f, &m_l)) {
    return false;
  }
  std::vector<std::string> chrom_names;
  for (const auto& kv : m_p) {
    chrom_names.push_back(kv.first);
  }
  std::sort(chrom_names.begin(), chrom_names.end());
  std::vector<Macs3BdgTrack> ppois_tracks;
  std::vector<Macs3BdgTrack> treat_tracks;
  std::vector<Macs3BdgTrack> lambda_tracks;
  ppois_tracks.reserve(chrom_names.size());
  treat_tracks.reserve(chrom_names.size());
  lambda_tracks.reserve(chrom_names.size());
  for (const std::string& chrom : chrom_names) {
    ppois_tracks.push_back(m_p[chrom]);
    const auto itt = m_t.find(chrom);
    const auto itl = m_l.find(chrom);
    treat_tracks.push_back(itt == m_t.end() ? Macs3BdgTrack() : itt->second);
    lambda_tracks.push_back(itl == m_l.end() ? Macs3BdgTrack() : itl->second);
  }
  return RunMacs3FragPpoisNarrowPeaksFromTracks(
      chrom_names, ppois_tracks, treat_tracks, lambda_tracks, cutoff, min_length,
      max_gap, pseudocount, path_narrowpeak, path_summits_bed, stage_profile,
      profile_input_frag_rows);
}

bool RunMacs3FragPpoisNarrowPeaksFromTracks(
    const std::vector<std::string>& chrom_names,
    const std::vector<Macs3BdgTrack>& ppois_tracks,
    const std::vector<Macs3BdgTrack>& treat_tracks,
    const std::vector<Macs3BdgTrack>& lambda_tracks,
    float cutoff, int32_t min_length, int32_t max_gap, double pseudocount,
    const std::string& path_narrowpeak, const std::string& path_summits_bed,
    StageProfileCollector* stage_profile, int64_t profile_input_frag_rows) {
  if ((path_narrowpeak.empty() && path_summits_bed.empty()) || min_length < 1 ||
      max_gap < 0 || chrom_names.size() != ppois_tracks.size() ||
      chrom_names.size() != treat_tracks.size() ||
      chrom_names.size() != lambda_tracks.size()) {
    return false;
  }
  StageProfileCollector::Tick t_regions_start{};
  if (stage_profile) {
    t_regions_start = StageProfileCollector::Now();
  }
  struct W {
    std::string chrom;
    int32_t start;
    int32_t end;
    int32_t summit;
    int32_t peak_off;
    double fc;
    double p_mlog;
    double q_mlog;
  };
  std::vector<W> w;
  w.reserve(256);
  for (size_t ci = 0; ci < chrom_names.size(); ++ci) {
    const Macs3BdgTrack& ppois = ppois_tracks[ci];
    const Macs3BdgTrack& treat = treat_tracks[ci];
    const Macs3BdgTrack& lambda = lambda_tracks[ci];
    if (ppois.p.empty() || treat.p.empty() || lambda.p.empty()) {
      continue;
    }
    std::vector<PpoisMergedRegion> regs;
    CallPpoisTrackMergedRegions(ppois, cutoff, min_length, max_gap, &regs);
    for (const PpoisMergedRegion& r : regs) {
      if (r.pieces.empty()) {
        continue;
      }
      float best_tv = 0.f;
      bool have = false;
      std::vector<int32_t> midpts;
      std::vector<size_t> pidx;
      for (size_t i = 0; i < r.pieces.size(); ++i) {
        const PpoisMergePiece& pc = r.pieces[i];
        const float tv = BdgValueAtOrZero(treat, pc.start);
        if (!have || tv > best_tv) {
          have = true;
          best_tv = tv;
          midpts.clear();
          pidx.clear();
          midpts.push_back((pc.start + pc.end) / 2);
          pidx.push_back(i);
        } else if (tv == best_tv) {
          midpts.push_back((pc.start + pc.end) / 2);
          pidx.push_back(i);
        }
      }
      if (!have || midpts.empty()) {
        continue;
      }
      const size_t mxi = (midpts.size() + 1) / 2 - 1;
      const int32_t summit = midpts[mxi];
      const PpoisMergePiece& winp = r.pieces[pidx[mxi]];
      const float t_w = BdgValueAtOrZero(treat, winp.start);
      const float c_w = BdgValueAtOrZero(lambda, winp.start);
      const double pml = Macs3BdgcmpPpoisNegLog10P(t_w, c_w, pseudocount);
      const double fc =
          (static_cast<double>(t_w) + pseudocount) /
          (static_cast<double>(c_w) + pseudocount);
      W row;
      row.chrom = chrom_names[ci];
      row.start = r.start;
      row.end = r.end;
      row.summit = summit;
      row.peak_off = static_cast<int32_t>(row.summit - row.start);
      row.fc = fc;
      row.p_mlog = pml;
      row.q_mlog = 0.0;
      w.push_back(std::move(row));
    }
  }
  StageProfileCollector::Tick t_np_start{};
  if (stage_profile) {
    const int64_t in_fr =
        profile_input_frag_rows >= 0 ? profile_input_frag_rows : 0;
    stage_profile->Record("bdgpeakcall_regions", t_regions_start,
                          StageProfileCollector::Now(), in_fr,
                          static_cast<int64_t>(w.size()), 0,
                          "in_memory_tracks");
    t_np_start = StageProfileCollector::Now();
  }
  std::vector<double> p_lin;
  p_lin.reserve(w.size());
  for (const W& e : w) {
    const double p = std::min(1.0, std::max(0.0, std::pow(10.0, -e.p_mlog)));
    p_lin.push_back(p);
  }
  std::vector<double> qv;
  BenjaminiHochbergFdr(p_lin, &qv);
  for (size_t i = 0; i < w.size(); ++i) {
    const double q = std::min(1.0, std::max(0.0, qv[i]));
    w[i].q_mlog = -std::log10(std::max(q, 1e-300));
  }
  std::sort(w.begin(), w.end(), [](const W& a, const W& b) {
    if (a.chrom != b.chrom) {
      return a.chrom < b.chrom;
    }
    if (a.start != b.start) {
      return a.start < b.start;
    }
    if (a.end != b.end) {
      return a.end < b.end;
    }
    return a.peak_off < b.peak_off;
  });
  if (!path_narrowpeak.empty()) {
    FILE* fpn = std::fopen(path_narrowpeak.c_str(), "w");
    if (fpn == nullptr) {
      return false;
    }
    for (size_t i = 0; i < w.size(); ++i) {
      const W& p = w[i];
      const std::string name = "peak_" + std::to_string(i + 1);
      const int sc = NarrowPeakIntScoreQ(p.q_mlog);
      if (std::fprintf(fpn, "%s\t%d\t%d\t%s\t%d\t.\t%.6g\t%.6g\t%.6g\t%d\n",
                       p.chrom.c_str(), static_cast<int>(p.start),
                       static_cast<int>(p.end), name.c_str(), sc, p.fc,
                       p.p_mlog, p.q_mlog, static_cast<int>(p.peak_off)) < 0) {
        std::fclose(fpn);
        return false;
      }
    }
    if (std::fclose(fpn) != 0) {
      return false;
    }
  }
  if (!path_summits_bed.empty()) {
    FILE* fps = std::fopen(path_summits_bed.c_str(), "w");
    if (fps == nullptr) {
      return false;
    }
    for (size_t i = 0; i < w.size(); ++i) {
      const W& p = w[i];
      const std::string name = "peak_" + std::to_string(i + 1);
      if (std::fprintf(fps, "%s\t%d\t%d\t%s\t%.6g\n", p.chrom.c_str(),
                       static_cast<int>(p.summit),
                       static_cast<int>(p.summit + 1), name.c_str(),
                       p.q_mlog) < 0) {
        std::fclose(fps);
        return false;
      }
    }
    if (std::fclose(fps) != 0) {
      return false;
    }
  }
  if (stage_profile) {
    const StageProfileCollector::Tick t_summits_end = StageProfileCollector::Now();
    const StageProfileCollector::Tick tm0 = StageProfileCollector::Now();
    int64_t bbytes = 0, sbytes = 0;
    if (!path_narrowpeak.empty()) {
      (void)StageProfileCollector::FileMetrics(path_narrowpeak, &bbytes, nullptr);
    }
    if (!path_summits_bed.empty()) {
      (void)StageProfileCollector::FileMetrics(path_summits_bed, &sbytes, nullptr);
    }
    const StageProfileCollector::Tick tm1 = StageProfileCollector::Now();
    const int64_t in_fr = profile_input_frag_rows >= 0 ? profile_input_frag_rows : 0;
    stage_profile->Record("narrowpeak_summits", t_np_start, t_summits_end, in_fr,
                          static_cast<int64_t>(w.size()), bbytes + sbytes,
                          "in_memory_tracks");
    stage_profile->Record("profile_file_metrics_narrowpeak_summits", tm0, tm1, 0, 0,
                          bbytes + sbytes, "stat_size_only");
  }
  return true;
}

}  // namespace peaks
}  // namespace chromap
