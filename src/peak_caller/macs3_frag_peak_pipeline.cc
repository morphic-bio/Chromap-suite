#include "macs3_frag_peak_pipeline.h"

#include <algorithm>
#include <cmath>
#include <cerrno>
#include <cstdio>
#include <cstdlib>

#include <ftw.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "bdgpeakcall.h"
#include "macs3_frag_workspace.h"
#include "peak_io.h"
#include "stage_profile.h"

namespace chromap {
namespace peaks {

namespace {

bool DirExists(const std::string& path) {
  struct stat st;
  return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

bool MakeDir(const std::string& path) {
  if (path.empty()) {
    return false;
  }
  if (mkdir(path.c_str(), 0755) == 0) {
    return true;
  }
  return errno == EEXIST && DirExists(path);
}

int UnlinkPath(const char* fpath, const struct stat* sb, int typeflag,
               struct FTW* ftwbuf) {
  (void)sb;
  (void)typeflag;
  (void)ftwbuf;
  return remove(fpath);
}

bool RemoveTree(const std::string& root) {
  if (root.empty()) {
    return true;
  }
  return nftw(root.c_str(), UnlinkPath, 64, FTW_DEPTH | FTW_PHYS) == 0;
}

std::string DefaultWorkParent(const std::string& work_dir_parent) {
  if (!work_dir_parent.empty()) {
    return work_dir_parent;
  }
  const char* tmp = getenv("TMPDIR");
  if (tmp && tmp[0] != '\0') {
    return std::string(tmp);
  }
  return "/tmp";
}

int NormalizeThreads(int requested, size_t n_tasks) {
  if (requested < 1 || n_tasks <= 1) {
    return 1;
  }
  const size_t rt = static_cast<size_t>(requested);
  return static_cast<int>(std::min(rt, n_tasks));
}

}  // namespace

float BdgPeakCallCutoffFromPValue(double p_value) {
  if (p_value <= 0.0 || p_value > 1.0 || std::isnan(p_value)) {
    return -1.f;
  }
  return static_cast<float>(-std::log10(p_value));
}

bool RunMacs3FragPeakPipelineFromFragments(
    std::vector<ChromFragments>* by_chrom,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message, StageProfileCollector* stage_profile,
    bool release_fragments_after_workspace) {
  auto fail = [&](const char* msg) {
    if (error_message) {
      *error_message = msg;
    }
    return false;
  };
  if (by_chrom == nullptr || narrowpeak_out.empty() || summits_out.empty()) {
    return fail("narrowPeak and summits output paths are required");
  }
  if (params.min_length < 1 || params.max_gap < 0 || params.llocal_bp < 1 ||
      params.effective_genome_size < 1 || params.peak_caller_threads < 1) {
    return fail("invalid MACS3 FRAG peak pipeline numeric parameters");
  }
  if (params.bdgpeakcall_cutoff <= 0.f) {
    return fail("invalid bdgpeakcall cutoff (expected -log10 p, e.g. 5 for 1e-5)");
  }
  const int peak_threads =
      NormalizeThreads(params.peak_caller_threads, by_chrom->size());

  int64_t n_frag_rows = 0;
  if (stage_profile) {
    for (const ChromFragments& ch : *by_chrom) {
      n_frag_rows += static_cast<int64_t>(ch.frags.size());
    }
  }

  const bool export_treat = !keep_intermediates_dir.empty() || !user_paths.treat_bdg.empty();
  const bool export_lambda = !keep_intermediates_dir.empty() || !user_paths.lambda_bdg.empty();
  const bool export_ppois = !keep_intermediates_dir.empty() || !user_paths.ppois_bdg.empty();
  const bool any_intermediate_export = export_treat || export_lambda || export_ppois;

  std::string work;
  bool own_workdir = false;
  if (!keep_intermediates_dir.empty()) {
    work = keep_intermediates_dir;
    if (!MakeDir(work)) {
      return fail("failed to create --macs3-frag-keep-intermediates directory");
    }
  } else if (any_intermediate_export) {
    const std::string parent = DefaultWorkParent(work_dir_parent);
    if (!DirExists(parent) && !MakeDir(parent)) {
      return fail("failed to create temp parent directory for MACS3 FRAG workdir");
    }
    work = parent + "/chromap_macs3_frag_XXXXXX";
    std::vector<char> tmpl(work.begin(), work.end() + 1);
    if (mkdtemp(tmpl.data()) == nullptr) {
      return fail("mkdtemp failed for MACS3 FRAG peak pipeline");
    }
    work.assign(tmpl.data());
    own_workdir = true;
  }
  if (work_dir_used) {
    *work_dir_used = work;
  }

  const std::string treat =
      user_paths.treat_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_treat_pileup.bdg"))
                                   : user_paths.treat_bdg;
  const std::string lambda =
      user_paths.lambda_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_control_lambda.bdg"))
                                    : user_paths.lambda_bdg;
  const std::string ppois =
      user_paths.ppois_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_score_ppois.bdg"))
                                   : user_paths.ppois_bdg;

  Macs3FragWorkspaceParams workspace_params;
  workspace_params.frag_pileup_max_count = params.frag_pileup_max_count;
  workspace_params.macs3_uint8_counts = params.macs3_uint8_counts;
  workspace_params.effective_genome_size = params.effective_genome_size;
  workspace_params.llocal_bp = params.llocal_bp;
  workspace_params.score_pseudocount = params.score_pseudocount;
  Macs3FragPeakWorkspace workspace;
  {
    auto build_workspace = [&]() {
      const bool ok =
          BuildMacs3FragWorkspaceFromFragments(*by_chrom, workspace_params, &workspace);
      if (ok && release_fragments_after_workspace) {
        std::vector<ChromFragments>().swap(*by_chrom);
      }
      return ok;
    };
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      if (!build_workspace()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("BuildMacs3FragWorkspaceFromFragments failed");
      }
      stage_profile->Record("build_workspace", t0, StageProfileCollector::Now(),
                            n_frag_rows, static_cast<int64_t>(workspace.chrom_names.size()),
                            0, "in_memory_events");
    } else {
      if (!build_workspace()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("BuildMacs3FragWorkspaceFromFragments failed");
      }
    }
  }

  {
    auto finalize_treat = [&]() {
      return FinalizeMacs3FragTreatTracks(&workspace, peak_threads);
    };
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      if (!finalize_treat()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragTreatTracks failed");
      }
      const StageProfileCollector::Tick t1 = StageProfileCollector::Now();
      stage_profile->Record("finalize_treat_tracks", t0, t1, n_frag_rows,
                            static_cast<int64_t>(workspace.treat_tracks.size()), 0,
                            "in_memory_tracks;threads=" + std::to_string(peak_threads));
    } else {
      if (!finalize_treat()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragTreatTracks failed");
      }
    }
  }
  {
    auto finalize_lam = [&]() {
      return FinalizeMacs3FragLambdaTracks(&workspace, peak_threads);
    };
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      if (!finalize_lam()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragLambdaTracks failed");
      }
      const StageProfileCollector::Tick t1 = StageProfileCollector::Now();
      stage_profile->Record("finalize_lambda_tracks", t0, t1, n_frag_rows,
                            static_cast<int64_t>(workspace.lambda_tracks.size()), 0,
                            "in_memory_tracks;threads=" + std::to_string(peak_threads));
    } else {
      if (!finalize_lam()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragLambdaTracks failed");
      }
    }
  }
  if (export_ppois) {
    auto finalize_score = [&]() {
      return FinalizeMacs3FragPpoisTracks(&workspace, peak_threads);
    };
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      if (!finalize_score()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragPpoisTracks failed");
      }
      const StageProfileCollector::Tick t1 = StageProfileCollector::Now();
      stage_profile->Record("finalize_ppois_tracks", t0, t1, n_frag_rows,
                            static_cast<int64_t>(workspace.ppois_tracks.size()), 0,
                            "in_memory_tracks;threads=" + std::to_string(peak_threads));
    } else {
      if (!finalize_score()) {
        if (own_workdir) {
          RemoveTree(work);
        }
        return fail("FinalizeMacs3FragPpoisTracks failed");
      }
    }
  }
  auto export_track = [&](const char* stage, const std::string& path,
                          const std::vector<Macs3BdgTrack>& tracks) {
    if (path.empty()) {
      return true;
    }
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      if (!WriteMacs3BdgTracksBedGraph(path, workspace.chrom_names, tracks)) {
        return false;
      }
      const StageProfileCollector::Tick t1 = StageProfileCollector::Now();
      const StageProfileCollector::Tick tm0 = StageProfileCollector::Now();
      int64_t bytes = 0, lines = 0;
      (void)StageProfileCollector::FileMetrics(path, &bytes, &lines);
      const StageProfileCollector::Tick tm1 = StageProfileCollector::Now();
      stage_profile->Record(stage, t0, t1, n_frag_rows, lines, bytes, "");
      stage_profile->Record(std::string("profile_file_metrics_") + stage, tm0, tm1,
                            n_frag_rows, lines, bytes, "line_count_and_size_scan");
      return true;
    }
    return WriteMacs3BdgTracksBedGraph(path, workspace.chrom_names, tracks);
  };
  if (!export_track("export_treat_pileup", treat, workspace.treat_tracks)) {
    if (own_workdir) {
      RemoveTree(work);
    }
    return fail("WriteFragSpanPileupBedGraph failed");
  }
  if (!export_track("export_control_lambda", lambda, workspace.lambda_tracks)) {
    if (own_workdir) {
      RemoveTree(work);
    }
    return fail("WriteFragMacs3NoControlLambdaBedGraph failed");
  }
  if (export_ppois) {
    if (!export_track("export_score_ppois", ppois, workspace.ppois_tracks)) {
      if (own_workdir) {
        RemoveTree(work);
      }
      return fail("WriteMacs3BdgTracksBedGraph ppois failed");
    }
    if (!RunMacs3FragPpoisNarrowPeaksFromTracks(
            workspace.chrom_names, workspace.ppois_tracks, workspace.treat_tracks,
            workspace.lambda_tracks, params.bdgpeakcall_cutoff, params.min_length,
            params.max_gap, params.score_pseudocount, narrowpeak_out, summits_out,
            stage_profile, n_frag_rows)) {
      if (own_workdir) {
        RemoveTree(work);
      }
      return fail("RunMacs3FragPpoisNarrowPeaks failed");
    }
  } else {
    StageProfileCollector::Tick t_ppois_start{};
    StageProfileCollector::Tick t_regions_start{};
    if (stage_profile) {
      t_ppois_start = StageProfileCollector::Now();
    }
    const size_t n_chrom = workspace.chrom_names.size();
    std::vector<std::vector<Macs3FragNarrowPeakRow>> rows_by_chrom(n_chrom);
    std::vector<char> ok(n_chrom, 1);
    const int threads = NormalizeThreads(peak_threads, n_chrom);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
    for (int64_t ii = 0; ii < static_cast<int64_t>(n_chrom); ++ii) {
      const size_t i = static_cast<size_t>(ii);
      Macs3BdgTrack ppois_one;
      if (!BuildMacs3FragPpoisTrack(workspace.treat_tracks[i], workspace.lambda_tracks[i],
                                    params.score_pseudocount, &ppois_one)) {
        ok[i] = 0;
        continue;
      }
      if (!CollectMacs3FragNarrowPeakRowsFromTracks(
              workspace.chrom_names[i], ppois_one, workspace.treat_tracks[i],
              workspace.lambda_tracks[i], params.bdgpeakcall_cutoff, params.min_length,
              params.max_gap, params.score_pseudocount, &rows_by_chrom[i])) {
        ok[i] = 0;
        continue;
      }
      if (!export_treat) {
        std::vector<int32_t>().swap(workspace.treat_tracks[i].p);
        std::vector<float>().swap(workspace.treat_tracks[i].v);
      }
      if (!export_lambda) {
        std::vector<int32_t>().swap(workspace.lambda_tracks[i].p);
        std::vector<float>().swap(workspace.lambda_tracks[i].v);
      }
    }
    if (std::find(ok.begin(), ok.end(), 0) != ok.end()) {
      if (own_workdir) {
        RemoveTree(work);
      }
      return fail("chromosome-at-a-time ppois peak collection failed");
    }
    std::vector<Macs3FragNarrowPeakRow> rows;
    size_t row_count = 0;
    for (const auto& chrom_rows : rows_by_chrom) {
      row_count += chrom_rows.size();
    }
    rows.reserve(row_count);
    for (auto& chrom_rows : rows_by_chrom) {
      rows.insert(rows.end(), chrom_rows.begin(), chrom_rows.end());
    }
    if (stage_profile) {
      t_regions_start = StageProfileCollector::Now();
      stage_profile->Record("finalize_ppois_tracks", t_ppois_start, t_regions_start,
                            n_frag_rows, static_cast<int64_t>(workspace.chrom_names.size()),
                            0, "chromosome_at_a_time;includes_region_collection;threads=" +
                                   std::to_string(threads));
      stage_profile->Record("bdgpeakcall_regions", t_regions_start, t_regions_start,
                            n_frag_rows, static_cast<int64_t>(rows.size()), 0,
                            "included_in_finalize_ppois_tracks");
    }
    StageProfileCollector::Tick t_np_start{};
    if (stage_profile) {
      t_np_start = StageProfileCollector::Now();
    }
    if (!WriteMacs3FragNarrowPeakRows(&rows, narrowpeak_out, summits_out)) {
      if (own_workdir) {
        RemoveTree(work);
      }
      return fail("WriteMacs3FragNarrowPeakRows failed");
    }
    if (stage_profile) {
      const StageProfileCollector::Tick t_summits_end = StageProfileCollector::Now();
      const StageProfileCollector::Tick tm0 = StageProfileCollector::Now();
      int64_t bbytes = 0, sbytes = 0;
      if (!narrowpeak_out.empty()) {
        (void)StageProfileCollector::FileMetrics(narrowpeak_out, &bbytes, nullptr);
      }
      if (!summits_out.empty()) {
        (void)StageProfileCollector::FileMetrics(summits_out, &sbytes, nullptr);
      }
      const StageProfileCollector::Tick tm1 = StageProfileCollector::Now();
      stage_profile->Record("narrowpeak_summits", t_np_start, t_summits_end, n_frag_rows,
                            static_cast<int64_t>(rows.size()), bbytes + sbytes,
                            "chromosome_at_a_time");
      stage_profile->Record("profile_file_metrics_narrowpeak_summits", tm0, tm1, 0, 0,
                            bbytes + sbytes, "stat_size_only");
    }
  }
  {
    auto do_cleanup = [&]() {
      if (own_workdir) {
        RemoveTree(work);
      }
    };
    if (stage_profile) {
      const StageProfileCollector::Tick t0 = StageProfileCollector::Now();
      do_cleanup();
      const std::string notes = own_workdir ? "removed_temp_workdir" : "no_temp_workdir";
      stage_profile->Record("cleanup", t0, StageProfileCollector::Now(), 0, 0, 0, notes);
    } else {
      do_cleanup();
    }
  }
  return true;
}

bool RunMacs3FragPeakPipelineFromWorkspace(
    Macs3FragPeakWorkspace* workspace,
    const Macs3FragPeakPipelineParams& params,
    const Macs3FragPeakPipelinePaths& user_paths,
    const std::string& narrowpeak_out, const std::string& summits_out,
    const std::string& keep_intermediates_dir,
    const std::string& work_dir_parent, std::string* work_dir_used,
    std::string* error_message) {
  auto fail = [&](const char* msg) {
    if (error_message) {
      *error_message = msg;
    }
    return false;
  };
  if (workspace == nullptr) {
    return fail("MACS3 FRAG workspace is null");
  }
  if (narrowpeak_out.empty() || summits_out.empty()) {
    return fail("narrowPeak and summits output paths are required");
  }
  if (params.min_length < 1 || params.max_gap < 0 || params.llocal_bp < 1 ||
      params.effective_genome_size < 1 || params.peak_caller_threads < 1) {
    return fail("invalid MACS3 FRAG peak pipeline numeric parameters");
  }
  if (params.bdgpeakcall_cutoff <= 0.f) {
    return fail("invalid bdgpeakcall cutoff (expected -log10 p, e.g. 5 for 1e-5)");
  }
  const int peak_threads =
      NormalizeThreads(params.peak_caller_threads, workspace->chrom_names.size());

  const bool export_treat = !keep_intermediates_dir.empty() || !user_paths.treat_bdg.empty();
  const bool export_lambda = !keep_intermediates_dir.empty() || !user_paths.lambda_bdg.empty();
  const bool export_ppois = !keep_intermediates_dir.empty() || !user_paths.ppois_bdg.empty();
  const bool any_intermediate_export = export_treat || export_lambda || export_ppois;

  std::string work;
  bool own_workdir = false;
  if (!keep_intermediates_dir.empty()) {
    work = keep_intermediates_dir;
    if (!MakeDir(work)) {
      return fail("failed to create --macs3-frag-keep-intermediates directory");
    }
  } else if (any_intermediate_export) {
    const std::string parent = DefaultWorkParent(work_dir_parent);
    if (!DirExists(parent) && !MakeDir(parent)) {
      return fail("failed to create temp parent directory for MACS3 FRAG workdir");
    }
    work = parent + "/chromap_macs3_frag_XXXXXX";
    std::vector<char> tmpl(work.begin(), work.end() + 1);
    if (mkdtemp(tmpl.data()) == nullptr) {
      return fail("mkdtemp failed for MACS3 FRAG peak pipeline");
    }
    work.assign(tmpl.data());
    own_workdir = true;
  }
  if (work_dir_used) {
    *work_dir_used = work;
  }

  const std::string treat =
      user_paths.treat_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_treat_pileup.bdg"))
                                   : user_paths.treat_bdg;
  const std::string lambda =
      user_paths.lambda_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_control_lambda.bdg"))
                                    : user_paths.lambda_bdg;
  const std::string ppois =
      user_paths.ppois_bdg.empty() ? (work.empty() ? "" : (work + "/macs3_frag_score_ppois.bdg"))
                                   : user_paths.ppois_bdg;

  auto cleanup_on_fail = [&]() {
    if (own_workdir) {
      RemoveTree(work);
    }
  };
  if (!FinalizeMacs3FragTreatTracks(workspace, peak_threads)) {
    cleanup_on_fail();
    return fail("FinalizeMacs3FragTreatTracks failed");
  }
  if (!FinalizeMacs3FragLambdaTracks(workspace, peak_threads)) {
    cleanup_on_fail();
    return fail("FinalizeMacs3FragLambdaTracks failed");
  }
  if (!treat.empty() &&
      !WriteMacs3BdgTracksBedGraph(treat, workspace->chrom_names, workspace->treat_tracks)) {
    cleanup_on_fail();
    return fail("WriteFragSpanPileupBedGraph failed");
  }
  if (!lambda.empty() &&
      !WriteMacs3BdgTracksBedGraph(lambda, workspace->chrom_names, workspace->lambda_tracks)) {
    cleanup_on_fail();
    return fail("WriteFragMacs3NoControlLambdaBedGraph failed");
  }
  if (export_ppois) {
    if (!FinalizeMacs3FragPpoisTracks(workspace, peak_threads)) {
      cleanup_on_fail();
      return fail("FinalizeMacs3FragPpoisTracks failed");
    }
    if (!WriteMacs3BdgTracksBedGraph(ppois, workspace->chrom_names,
                                     workspace->ppois_tracks)) {
      cleanup_on_fail();
      return fail("WriteMacs3BdgTracksBedGraph ppois failed");
    }
    if (!RunMacs3FragPpoisNarrowPeaksFromTracks(
            workspace->chrom_names, workspace->ppois_tracks, workspace->treat_tracks,
            workspace->lambda_tracks, params.bdgpeakcall_cutoff, params.min_length,
            params.max_gap, params.score_pseudocount, narrowpeak_out, summits_out,
            nullptr, -1)) {
      cleanup_on_fail();
      return fail("RunMacs3FragPpoisNarrowPeaks failed");
    }
  } else {
    const size_t n_chrom = workspace->chrom_names.size();
    std::vector<std::vector<Macs3FragNarrowPeakRow>> rows_by_chrom(n_chrom);
    std::vector<char> ok(n_chrom, 1);
    const int threads = NormalizeThreads(peak_threads, n_chrom);
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads) if(threads > 1)
    for (int64_t ii = 0; ii < static_cast<int64_t>(n_chrom); ++ii) {
      const size_t i = static_cast<size_t>(ii);
      Macs3BdgTrack ppois_one;
      if (!BuildMacs3FragPpoisTrack(workspace->treat_tracks[i],
                                    workspace->lambda_tracks[i],
                                    params.score_pseudocount, &ppois_one) ||
          !CollectMacs3FragNarrowPeakRowsFromTracks(
              workspace->chrom_names[i], ppois_one, workspace->treat_tracks[i],
              workspace->lambda_tracks[i], params.bdgpeakcall_cutoff, params.min_length,
              params.max_gap, params.score_pseudocount, &rows_by_chrom[i])) {
        ok[i] = 0;
        continue;
      }
      if (!export_treat) {
        std::vector<int32_t>().swap(workspace->treat_tracks[i].p);
        std::vector<float>().swap(workspace->treat_tracks[i].v);
      }
      if (!export_lambda) {
        std::vector<int32_t>().swap(workspace->lambda_tracks[i].p);
        std::vector<float>().swap(workspace->lambda_tracks[i].v);
      }
    }
    if (std::find(ok.begin(), ok.end(), 0) != ok.end()) {
      cleanup_on_fail();
      return fail("chromosome-at-a-time ppois peak collection failed");
    }
    std::vector<Macs3FragNarrowPeakRow> rows;
    size_t row_count = 0;
    for (const auto& chrom_rows : rows_by_chrom) {
      row_count += chrom_rows.size();
    }
    rows.reserve(row_count);
    for (auto& chrom_rows : rows_by_chrom) {
      rows.insert(rows.end(), chrom_rows.begin(), chrom_rows.end());
    }
    if (!WriteMacs3FragNarrowPeakRows(&rows, narrowpeak_out, summits_out)) {
      cleanup_on_fail();
      return fail("WriteMacs3FragNarrowPeakRows failed");
    }
  }
  if (own_workdir) {
    RemoveTree(work);
  }
  return true;
}

}  // namespace peaks
}  // namespace chromap
