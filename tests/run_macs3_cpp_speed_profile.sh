#!/usr/bin/env bash
# Benchmark harness: external macs3 callpeak vs chromap_callpeaks with optional
# --profile-stages. Emits a summary TSV (wall/RSS, parity, artifact paths).
# RUN_MACS3=0 skips external MACS3 even when macs3 is on PATH (C++ only).
# EXPORT_BDGS=0 skips C++ diagnostic treat/lambda/ppois bedGraph exports and
# measures the workspace production shape (narrowPeak/summits only).
# CPP_THREADS controls chromap_callpeaks --peak-caller-threads (default: 1).
#
# 5M example (after tests/build_peak_speed_fragments_5m.sh):
#   FRAGMENTS_TSV_GZ=/path/to/Chromap-suite/out/peak_speed_parity_20260425/fragments_5m.tsv.gz \
#   OUTDIR=/path/to/Chromap-suite/out/peak_speed_parity_20260425/5m RUN_MACS3=1 EXPORT_BDGS=0 \
#   bash tests/run_macs3_cpp_speed_profile.sh
# shellcheck shell=bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

BENCH_ROOT="${CHROMAP_100K_BENCH:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k}"
CHROMAP_PEAK_RUN_ROOT="${CHROMAP_PEAK_RUN_ROOT:-}"
FRAGMENTS_TSV_GZ="${FRAGMENTS_TSV_GZ:-}"
OUTDIR="${OUTDIR:-}"
MACS3_BIN="${MACS3_BIN:-macs3}"
RUN_FULL="${RUN_FULL:-0}"
RUN_MACS3="${RUN_MACS3:-0}"
EXPORT_BDGS="${EXPORT_BDGS:-1}"
CPP_THREADS="${CPP_THREADS:-1}"

CUTOFF="${CUTOFF:-5}"
MINLEN="${MINLEN:-200}"
MAXGAP="${MAXGAP:-30}"
GENOME_SIZE="${GENOME_SIZE:-hs}"

# shellcheck source=peak_caller_100k_common.sh
source "${SCRIPT_DIR}/peak_caller_100k_common.sh"

gnu_time_field() {
  local log=$1
  local key=$2
  grep "^[[:space:]]*${key}" "${log}" 2>/dev/null | head -1 | awk '{print $NF}' || true
}

# GNU /usr/bin/time -v prints m:ss or h:mm:ss in the last field.
time_to_sec() {
  local t=$1
  [[ -n "${t}" ]] || {
    echo ""
    return
  }
  local n
  n="${t//[^:]/}"
  if [[ "${#n}" -eq 1 ]]; then
    local m s
    m="${t%%:*}"
    s="${t#*:}"
    awk -v mm="${m}" -v ss="${s}" 'BEGIN { printf "%.4f", mm * 60 + ss }'
  elif [[ "${#n}" -eq 2 ]]; then
    local h rest m s
    h="${t%%:*}"
    rest="${t#*:}"
    m="${rest%%:*}"
    s="${rest#*:}"
    awk -v hh="${h}" -v mm="${m}" -v ss="${s}" 'BEGIN { printf "%.4f", hh * 3600 + mm * 60 + ss }'
  else
    printf '%s' "${t}"
  fi
}

main() {
  if [[ "${RUN_FULL}" == "1" ]]; then
    if [[ -z "${FRAGMENTS_TSV_GZ}" || ! -f "${FRAGMENTS_TSV_GZ}" ]]; then
      echo "RUN_FULL=1 requires an existing FRAGMENTS_TSV_GZ=/path/to/fragments.tsv.gz" >&2
      return 2
    fi
  fi

  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/macs3_cpp_speed_profile.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}/macs3_callpeak" "${OUTDIR}/cpp"

  (cd "${REPO_ROOT}" && make chromap_callpeaks)

  if ! peak_100k_resolve_inputs; then
    return 1
  fi
  local fr="${PEAK_100K_FRAG}"
  local caller="${REPO_ROOT}/chromap_callpeaks"

  local macs_wall="" macs_rss="" macs_np="" macs_log="${OUTDIR}/macs3_time.log"
  if [[ "${RUN_MACS3}" != "0" ]] && command -v "${MACS3_BIN}" &>/dev/null; then
    echo "[macs3] callpeak -f FRAG (GNU time -v) -> ${macs_log}" >&2
    set +e
    /usr/bin/time -v "${MACS3_BIN}" callpeak -t "${fr}" -f FRAG -g "${GENOME_SIZE}" \
      -n macs3_frag -p 1e-5 --min-length "${MINLEN}" --max-gap "${MAXGAP}" \
      --outdir "${OUTDIR}/macs3_callpeak" >"${OUTDIR}/macs3_stdout.log" 2>"${macs_log}"
    local macs_ec=$?
    set -e
    if [[ "${macs_ec}" -ne 0 ]]; then
      echo "WARN: macs3 callpeak exited ${macs_ec} (see ${macs_log})" >&2
    else
      macs_np="${OUTDIR}/macs3_callpeak/macs3_frag_peaks.narrowPeak"
      macs_wall="$(gnu_time_field "${macs_log}" "Elapsed (wall clock)")"
      macs_rss="$(gnu_time_field "${macs_log}" "Maximum resident set size")"
    fi
  elif [[ "${RUN_MACS3}" == "0" ]]; then
    echo "[macs3] skipped (RUN_MACS3=0)" >&2
  else
    echo "WARN: ${MACS3_BIN} not in PATH; skipping external MACS3 timing" >&2
  fi

  local treat="${OUTDIR}/cpp/cpp_frag_treat_pileup.bdg"
  local lam="${OUTDIR}/cpp/cpp_frag_control_lambda.bdg"
  local cpp_p="${OUTDIR}/cpp/cpp_frag_score_ppois.bdg"
  local cpp_np_off="${OUTDIR}/cpp/cpp_macs3_frag_off.narrowPeak"
  local cpp_smt_off="${OUTDIR}/cpp/cpp_macs3_frag_summits_off.bed"
  local cpp_np_on="${OUTDIR}/cpp/cpp_macs3_frag_on.narrowPeak"
  local cpp_smt_on="${OUTDIR}/cpp/cpp_macs3_frag_summits_on.bed"
  local stage_tsv="${OUTDIR}/cpp/stage_profile.tsv"

  local cpp_common_args=(
    -i "${fr}"
    --frag-pileup-macs3-uint8-counts
    --frag-lambda-effective-genome 2913022398
    --frag-score-pseudocount 0
    --bdgpeakcall-cutoff "${CUTOFF}"
    --bdgpeakcall-min-len "${MINLEN}"
    --bdgpeakcall-max-gap "${MAXGAP}"
    --peak-caller-threads "${CPP_THREADS}"
  )
  if [[ "${EXPORT_BDGS}" != "0" ]]; then
    cpp_common_args+=(
      --frag-span-pileup-bdg "${treat}"
      --frag-lambda-bdg "${lam}"
      --frag-score-ppois-bdg "${cpp_p}"
    )
  fi

  echo "[cpp] chromap_callpeaks (profiling off) -> ${cpp_np_off}" >&2
  local cpp_off_log="${OUTDIR}/cpp/chromap_off_time.log"
  set +e
  /usr/bin/time -v "${caller}" "${cpp_common_args[@]}" \
    --macs3-frag-narrowpeak "${cpp_np_off}" \
    --macs3-frag-summits "${cpp_smt_off}" >"${OUTDIR}/cpp/chromap_off_stdout.log" 2>"${cpp_off_log}"
  local cpp_off_ec=$?
  set -e
  if [[ "${cpp_off_ec}" -ne 0 ]]; then
    echo "ERROR: chromap_callpeaks (no profile) failed (${cpp_off_ec}); see ${cpp_off_log}" >&2
    return 1
  fi
  local cpp_wall_off cpp_rss_off
  cpp_wall_off="$(gnu_time_field "${cpp_off_log}" "Elapsed (wall clock)")"
  cpp_rss_off="$(gnu_time_field "${cpp_off_log}" "Maximum resident set size")"

  echo "[cpp] chromap_callpeaks --profile-stages ${stage_tsv}" >&2
  local cpp_on_log="${OUTDIR}/cpp/chromap_on_time.log"
  set +e
  /usr/bin/time -v "${caller}" "${cpp_common_args[@]}" \
    --profile-stages "${stage_tsv}" \
    --macs3-frag-narrowpeak "${cpp_np_on}" \
    --macs3-frag-summits "${cpp_smt_on}" >"${OUTDIR}/cpp/chromap_on_stdout.log" 2>"${cpp_on_log}"
  local cpp_on_ec=$?
  set -e
  if [[ "${cpp_on_ec}" -ne 0 ]]; then
    echo "ERROR: chromap_callpeaks (with profile) failed (${cpp_on_ec}); see ${cpp_on_log}" >&2
    return 1
  fi
  local cpp_wall_on cpp_rss_on
  cpp_wall_on="$(gnu_time_field "${cpp_on_log}" "Elapsed (wall clock)")"
  cpp_rss_on="$(gnu_time_field "${cpp_on_log}" "Maximum resident set size")"

  local parity_np="NA" parity_smt="NA"
  if cmp -s "${cpp_np_off}" "${cpp_np_on}"; then
    parity_np="identical_bytes"
  else
    parity_np="DIFFER"
  fi
  if cmp -s "${cpp_smt_off}" "${cpp_smt_on}"; then
    parity_smt="identical_bytes"
  else
    parity_smt="DIFFER"
  fi

  local summary="${OUTDIR}/speed_profile_summary.tsv"
  {
    echo -e "metric\tvalue"
    echo -e "fragments_tsv_gz\t${fr}"
    echo -e "outdir\t${OUTDIR}"
    echo -e "export_bdgs\t${EXPORT_BDGS}"
    echo -e "cpp_threads\t${CPP_THREADS}"
    echo -e "macs3_wall_sec\t${macs_wall}"
    echo -e "macs3_max_rss_kb\t${macs_rss}"
    echo -e "macs3_narrowpeak\t${macs_np}"
    echo -e "macs3_time_log\t${macs_log}"
    echo -e "cpp_wall_sec_prof_off\t${cpp_wall_off}"
    echo -e "cpp_max_rss_kb_prof_off\t${cpp_rss_off}"
    echo -e "cpp_wall_sec_prof_on\t${cpp_wall_on}"
    echo -e "cpp_max_rss_kb_prof_on\t${cpp_rss_on}"
    echo -e "cpp_stage_profile_tsv\t${stage_tsv}"
    echo -e "cpp_off_log\t${cpp_off_log}"
    echo -e "cpp_on_log\t${cpp_on_log}"
    echo -e "parity_narrowpeak_prof_on_vs_off\t${parity_np}"
    echo -e "parity_summits_prof_on_vs_off\t${parity_smt}"
    local macs_s cpp_s ratio
    macs_s="$(time_to_sec "${macs_wall}")"
    cpp_s="$(time_to_sec "${cpp_wall_on}")"
    ratio="NA"
    if [[ -n "${macs_s}" && -n "${cpp_s}" ]]; then
      ratio="$(awk -v m="${macs_s}" -v c="${cpp_s}" 'BEGIN { if (m > 0) printf "%.4f", c / m; else print "NA" }')"
    fi
    echo -e "macs3_cpp_wall_ratio_cpp_over_macs3\t${ratio}"
  } | tee "${summary}"

  echo "Summary: ${summary}" >&2
  echo "Stage profile: ${stage_tsv}" >&2
  if [[ "${parity_np}" != "identical_bytes" || "${parity_smt}" != "identical_bytes" ]]; then
    echo "FAIL: profiling changed narrowPeak or summits bytes" >&2
    return 1
  fi
  return 0
}

main "$@"
