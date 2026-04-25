#!/usr/bin/env bash
# Phase 6: MACS3 callpeak FRAG (narrowPeak + summits) vs C++ MACS-pipeline narrowPeak.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

BENCH_ROOT="${CHROMAP_100K_BENCH:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k}"
CHROMAP_PEAK_RUN_ROOT="${CHROMAP_PEAK_RUN_ROOT:-}"
FRAGMENTS_TSV_GZ="${FRAGMENTS_TSV_GZ:-}"
ATAC_BAM="${ATAC_BAM:-}"
OUTDIR="${OUTDIR:-}"
MACS3_BIN="${MACS3_BIN:-macs3}"
export RUN_MACS3="${RUN_MACS3:-0}"

CUTOFF="${CUTOFF:-5}"
MINLEN="${MINLEN:-200}"
MAXGAP="${MAXGAP:-30}"
GENOME_SIZE="${GENOME_SIZE:-hs}"

# shellcheck source=peak_caller_100k_common.sh
source "${SCRIPT_DIR}/peak_caller_100k_common.sh"

bed3_from_narrow() {
  local in_np=$1
  local out_bed=$2
  awk 'BEGIN{FS=OFS="\t"} !/^track/ && !/^#/ && NF>=3 { print $1, $2+0, $3+0 }' "${in_np}" |
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n >"${out_bed}"
}

main() {
  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/chromap_narrowpeak_parity.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}"
  local log="${OUTDIR}/phase6_harness.log"
  log_msg() { echo "$@" | tee -a "${log}"; }

  if ! command -v "${MACS3_BIN}" &>/dev/null; then
    echo "Missing ${MACS3_BIN}" >&2
    exit 2
  fi

  (cd "${REPO_ROOT}" && make chromap_callpeaks) >>"${log}" 2>&1

  if ! peak_100k_resolve_inputs; then
    exit 1
  fi
  local fr
  fr="${PEAK_100K_FRAG}"

  local d="${OUTDIR}/phase6_narrowpeak"
  mkdir -p "${d}/macs3_callpeak" "${d}/macs3_summits"
  # -p 1e-5 matches the ~5 -log10(p) ppois track used with bdgpeakcall -c 5; explicit
  # -l / --max-gap align internal peak caller with the Phase 5 bdg defaults.
  log_msg "[1] MACS3 callpeak -f FRAG -p 1e-5 (reference narrowPeak; p-based like ppois c=${CUTOFF})"
  ${MACS3_BIN} callpeak -t "${fr}" -f FRAG -g "${GENOME_SIZE}" -n macs3_frag \
    -p 1e-5 --min-length "${MINLEN}" --max-gap "${MAXGAP}" \
    --outdir "${d}/macs3_callpeak" >>"${log}" 2>&1
  log_msg "[2] MACS3 callpeak --call-summits (reference summits.bed, same gating as [1])"
  ${MACS3_BIN} callpeak -t "${fr}" -f FRAG -g "${GENOME_SIZE}" -p 1e-5 \
    --min-length "${MINLEN}" --max-gap "${MAXGAP}" --call-summits \
    -n macs3_frag_summits --outdir "${d}/macs3_summits" >>"${log}" 2>&1

  local macs_np="${d}/macs3_callpeak/macs3_frag_peaks.narrowPeak"
  local macs_smt="${d}/macs3_summits/macs3_frag_summits_summits.bed"
  local macs_xls="${d}/macs3_callpeak/macs3_frag_peaks.xls"
  if [[ ! -f "${macs_np}" ]]; then
    echo "Missing ${macs_np}" >&2
    exit 1
  fi
  if [[ ! -f "${macs_smt}" ]]; then
    # alternate naming if -n differed
    macs_smt="$(ls "${d}/macs3_summits/"*_summits.bed 2>/dev/null | head -1 || true)"
  fi

  local mac_bed3="${d}/macs3_frag.bed3"
  local cpp_bed3="${d}/cpp_macs3_frag.bed3"
  bed3_from_narrow "${macs_np}" "${mac_bed3}"

  local treat="${d}/cpp_frag_treat_pileup.bdg"
  local lam="${d}/cpp_frag_control_lambda.bdg"
  local cpp_p="${d}/cpp_frag_score_ppois.bdg"
  local cpp_np="${d}/cpp_macs3_frag_peaks.narrowPeak"
  local cpp_smt="${d}/cpp_macs3_frag_summits.bed"
  local caller="${REPO_ROOT}/chromap_callpeaks"

  log_msg "[3] C++: pileup, lambda, ppois, diagnostic narrowPeak + summits"
  "${caller}" -i "${fr}" \
    --frag-pileup-macs3-uint8-counts \
    --frag-span-pileup-bdg "${treat}" \
    --frag-lambda-bdg "${lam}" \
    --frag-lambda-effective-genome 2913022398 \
    --frag-score-ppois-bdg "${cpp_p}" \
    --frag-score-pseudocount 0 \
    --bdgpeakcall-cutoff "${CUTOFF}" --bdgpeakcall-min-len "${MINLEN}" --bdgpeakcall-max-gap "${MAXGAP}" \
    --macs3-frag-narrowpeak "${cpp_np}" \
    --macs3-frag-summits "${cpp_smt}" >>"${log}" 2>&1

  bed3_from_narrow "${cpp_np}" "${cpp_bed3}"
  local cmp_tsv="${d}/phase6_narrowpeak_compare.tsv"
  log_msg "[4] compare BED3 + columns"
  local cmp_py="${SCRIPT_DIR}/compare_narrowpeak_parity.py"
  chmod +x "${cmp_py}" 2>/dev/null || true
  sm_args=()
  if [[ -f "${macs_smt}" && -f "${cpp_smt}" ]]; then
    sm_args=(--macs3-summits "${macs_smt}" --cpp-summits "${cpp_smt}")
  fi
  python3 "${cmp_py}" \
    --macs3-np "${macs_np}" --cpp-np "${cpp_np}" \
    "${sm_args[@]}" \
    --out-tsv "${cmp_tsv}" 2>>"${log}" | tee -a "${log}"

  local b3_tsv="${d}/bed3_compare.tsv"
  log_msg "[5] BED3 interval metrics (compare_bdgpeakcall_regions.py)"
  python3 "${SCRIPT_DIR}/compare_bdgpeakcall_regions.py" "${mac_bed3}" "${cpp_bed3}" \
    --label-a macs3_narrow_bed3 --label-b cpp_macs3_frag_bed3 | tee "${b3_tsv}" | tee -a "${log}"

  local summary="${d}/phase6_narrowpeak_summary.tsv"
  {
    echo "phase	6_narrowpeak_parity"
    echo "outdir	${d}"
    echo "bed3_compare	${b3_tsv}"
    echo "narrowpeak_metrics	${cmp_tsv}"
    echo "log	${log}"
    echo "macs3_narrowpeak	${macs_np}"
    echo "macs3_summits	${macs_smt}"
    echo "macs3_xls	${macs_xls}"
    echo "cpp_narrowpeak	${cpp_np}"
    echo "cpp_summits	${cpp_smt}"
    if [[ -f "${macs_xls}" ]]; then
      echo "xls_header_macs3	$(head -1 "${macs_xls}" | tr '\t' ',')"
    fi
  } | tee "${summary}"

  if [[ -f "${macs_smt}" && -f "${cpp_smt}" ]]; then
    log_msg "Summits: MACS3 ${macs_smt} vs C++ ${cpp_smt}"
  else
    log_msg "Note: some summits files missing (skipping in compare)"
  fi

  if ! cmp -s "${mac_bed3}" "${cpp_bed3}"; then
    echo "FAIL: sorted BED3 from narrowPeak differs (see ${b3_tsv} and ${cmp_tsv}; fix intervals before summits/scores)." >&2
    return 1
  fi
  log_msg "OK: sorted BED3 identical (callpeak vs C++ diagnostic narrowPeak BED3)."
  return 0
}

main "$@"
