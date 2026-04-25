#!/usr/bin/env bash
# Phase 2B: MACS3 `pileup -f FRAG` vs chromap_callpeaks `--frag-span-pileup-bdg`
# (fragment-span pileup). Uses `--frag-pileup-macs3-uint8-counts` so weights match
# MACS3 PETrackII (uint8 storage). MACS3 `pileup` supports `--max-count` for FRAG
# (same as callpeak). Diagnostic only.
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
RELAX="${CHROMAP_FRAG_PILEUP_RELAX:-0}"

# shellcheck source=peak_caller_100k_common.sh
source "${SCRIPT_DIR}/peak_caller_100k_common.sh"

main() {
  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/chromap_frag_span_pileup.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}"
  local log="${OUTDIR}/phase2b_harness.log"
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

  local d="${OUTDIR}/phase2b_frag_span_pileup"
  mkdir -p "${d}"

  local macs_def cpp_def macs_mc cpp_mc
  macs_def="${d}/macs3_pileup_FRAG_default.bdg"
  cpp_def="${d}/chromap_frag_span_default.bdg"
  macs_mc="${d}/macs3_pileup_FRAG_maxcount1.bdg"
  cpp_mc="${d}/chromap_frag_span_maxcount1.bdg"

  local cmd_macs_def cmd_macs_mc cmd_cpp_def cmd_cpp_mc
  cmd_macs_def="${MACS3_BIN} pileup -i ${fr} -f FRAG -o ${macs_def}"
  log_msg "[1] ${cmd_macs_def}"
  ${MACS3_BIN} pileup -i "${fr}" -f FRAG -o "${macs_def}" >>"${log}" 2>&1

  local caller="${REPO_ROOT}/chromap_callpeaks"
  cmd_cpp_def="${caller} -i ${fr} --frag-span-pileup-bdg ${cpp_def} --frag-pileup-macs3-uint8-counts"
  log_msg "[2] ${cmd_cpp_def}"
  "${caller}" -i "${fr}" --frag-span-pileup-bdg "${cpp_def}" \
    --frag-pileup-macs3-uint8-counts >>"${log}" 2>&1

  cmd_macs_mc="${MACS3_BIN} pileup -i ${fr} -f FRAG --max-count 1 -o ${macs_mc}"
  log_msg "[3] ${cmd_macs_mc}"
  ${MACS3_BIN} pileup -i "${fr}" -f FRAG --max-count 1 -o "${macs_mc}" >>"${log}" 2>&1

  cmd_cpp_mc="${caller} -i ${fr} --frag-span-pileup-bdg ${cpp_mc} --frag-pileup-macs3-uint8-counts --frag-pileup-max-count 1"
  log_msg "[4] ${cmd_cpp_mc}"
  "${caller}" -i "${fr}" --frag-span-pileup-bdg "${cpp_mc}" \
    --frag-pileup-macs3-uint8-counts --frag-pileup-max-count 1 >>"${log}" 2>&1

  local py="${SCRIPT_DIR}/compare_pileup_bdg.py"
  local sum_def sum_mc
  sum_def="${d}/frag_span_pileup_compare_default.tsv"
  sum_mc="${d}/frag_span_pileup_compare_maxcount1.tsv"
  log_msg "[5] compare default (MACS3 vs C++)"
  python3 "${py}" "${macs_def}" "${cpp_def}" --high-prefix "${d}/pileup_high_default" | tee "${sum_def}"
  log_msg "[6] compare --max-count 1"
  python3 "${py}" "${macs_mc}" "${cpp_mc}" --high-prefix "${d}/pileup_high_maxcount1" | tee "${sum_mc}"

  local hj_def hj_mc
  hj_def="NA"
  hj_mc="NA"
  if command -v bedtools &>/dev/null; then
    local ha hb
    if [[ -f "${d}/pileup_high_default_a.bed" && -f "${d}/pileup_high_default_b.bed" ]]; then
      bedtools sort -i "${d}/pileup_high_default_a.bed" >"${d}/high_def_a.s.bed"
      bedtools sort -i "${d}/pileup_high_default_b.bed" >"${d}/high_def_b.s.bed"
      hj_def="$(bedtools jaccard -a "${d}/high_def_a.s.bed" -b "${d}/high_def_b.s.bed" 2>/dev/null | awk 'NR==2 {print $3}')"
    fi
    if [[ -f "${d}/pileup_high_maxcount1_a.bed" && -f "${d}/pileup_high_maxcount1_b.bed" ]]; then
      bedtools sort -i "${d}/pileup_high_maxcount1_a.bed" >"${d}/high_mc_a.s.bed"
      bedtools sort -i "${d}/pileup_high_maxcount1_b.bed" >"${d}/high_mc_b.s.bed"
      hj_mc="$(bedtools jaccard -a "${d}/high_mc_a.s.bed" -b "${d}/high_mc_b.s.bed" 2>/dev/null | awk 'NR==2 {print $3}')"
    fi
  fi
  echo -e "high_signal_jaccard_default\t${hj_def}" | tee -a "${sum_def}"
  echo -e "high_signal_jaccard_maxcount1\t${hj_mc}" | tee -a "${sum_mc}"

  local summary="${d}/phase2b_summary.tsv"
  {
    echo "metric	value"
    echo "phase	2B_frag_span_pileup"
    echo "compare_default_tsv	${sum_def}"
    echo "compare_maxcount1_tsv	${sum_mc}"
    echo "note	MACS3 pileup supports --max-count for FRAG (see macs3 pileup --help)"
    echo "note2	C++ uses MACS3 pileup_PV segment rules; pass --frag-pileup-macs3-uint8-counts for PETrackII count storage parity"
    grep -E '^(bdg_rows_identical_sorted|pearson_bp_weighted_either_nonempty|total_signal_rel_diff|spearman_bp_weighted_either_nonempty)[[:space:]]' "${sum_def}" | sed 's/^/default_/' || true
    echo "high_signal_jaccard_default	${hj_def}"
    grep -E '^(bdg_rows_identical_sorted|pearson_bp_weighted_either_nonempty|total_signal_rel_diff|spearman_bp_weighted_either_nonempty)[[:space:]]' "${sum_mc}" | sed 's/^/maxcount1_/' || true
    echo "high_signal_jaccard_maxcount1	${hj_mc}"
  } | tee "${summary}"

  {
    echo "PEAK_100K_RUN_ROOT=${PEAK_100K_RUN_ROOT:-}"
    echo "FRAGMENTS_TSV_GZ=${fr}"
    echo "CMD_MACS3_PILEUP_DEFAULT=${cmd_macs_def}"
    echo "CMD_CHROMAP_FRAG_SPAN_DEFAULT=${cmd_cpp_def}"
    echo "CMD_MACS3_PILEUP_MAXCOUNT1=${cmd_macs_mc}"
    echo "CMD_CHROMAP_FRAG_SPAN_MAXCOUNT1=${cmd_cpp_mc}"
  } >"${OUTDIR}/phase2b_commands.txt"

  log_msg "Done. Summary TSV: ${summary}"
  log_msg "Output directory: ${OUTDIR}"

  local ok_def ok_mc
  ok_def="$(awk -F'\t' '$1=="bdg_rows_identical_sorted"{print $2}' "${sum_def}")"
  ok_mc="$(awk -F'\t' '$1=="bdg_rows_identical_sorted"{print $2}' "${sum_mc}")"
  if [[ "${RELAX}" != "1" ]]; then
    if [[ "${ok_def}" != "True" ]]; then
      echo "FAIL: default fragment-span pileup rows differ from MACS3 (see ${sum_def})" >&2
      exit 1
    fi
    if [[ "${ok_mc}" != "True" ]]; then
      echo "FAIL: max-count 1 fragment-span pileup rows differ from MACS3 (see ${sum_mc})" >&2
      exit 1
    fi
  fi
}

main "$@"
