#!/usr/bin/env bash
# Phase 2: Compare MACS3 `pileup -f FRAG` bedGraph vs chromap_callpeaks `--pileup-bdg`
# (binned Tn5-extended cut pileup). Diagnostic only; semantics differ by design.
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

# shellcheck source=peak_caller_100k_common.sh
source "${SCRIPT_DIR}/peak_caller_100k_common.sh"

RUN_MACS3=0

main() {
  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/chromap_pileup_parity.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}"
  local log="${OUTDIR}/phase2_harness.log"
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

  local d="${OUTDIR}/phase2_pileup"
  mkdir -p "${d}"

  local macs_bdg="${d}/macs3_pileup_FRAG.bdg"
  local cpp_bdg="${d}/chromap_binned_cut_pileup.bdg"

  local cmd_macs
  cmd_macs="${MACS3_BIN} pileup -i ${fr} -f FRAG -o ${macs_bdg}"
  log_msg "[1] ${cmd_macs}"
  ${MACS3_BIN} pileup -i "${fr}" -f FRAG -o "${macs_bdg}" >>"${log}" 2>&1

  local caller="${REPO_ROOT}/chromap_callpeaks"
  local cmd_cpp
  cmd_cpp="${caller} -i ${fr} --pileup-bdg ${cpp_bdg} --bin-width 50 --ext-size 150"
  log_msg "[2] ${cmd_cpp}"
  "${caller}" -i "${fr}" --pileup-bdg "${cpp_bdg}" --bin-width 50 --ext-size 150 >>"${log}" 2>&1

  local py="${SCRIPT_DIR}/compare_pileup_bdg.py"
  local cmp="${d}/pileup_compare_summary.tsv"
  log_msg "[3] python3 ${py} ..."
  python3 "${py}" "${macs_bdg}" "${cpp_bdg}" --high-prefix "${d}/pileup_high_signal" | tee "${cmp}"

  local hj="NA"
  if command -v bedtools &>/dev/null; then
    local ha hb
    ha="${d}/pileup_high_signal_a.bed"
    hb="${d}/pileup_high_signal_b.bed"
    if [[ -f "${ha}" && -f "${hb}" ]]; then
      bedtools sort -i "${ha}" >"${d}/high_a.s.bed"
      bedtools sort -i "${hb}" >"${d}/high_b.s.bed"
      hj="$(bedtools jaccard -a "${d}/high_a.s.bed" -b "${d}/high_b.s.bed" 2>/dev/null | awk 'NR==2 {print $3}')"
    fi
  fi
  echo -e "high_signal_jaccard_1bp\t${hj}" | tee -a "${cmp}"

  {
    echo "PEAK_100K_RUN_ROOT=${PEAK_100K_RUN_ROOT:-}"
    echo "FRAGMENTS_TSV_GZ=${fr}"
    echo "CMD_MACS3_PILEUP=${cmd_macs}"
    echo "CMD_CHROMAP_PILEUP_BDG=${cmd_cpp}"
  } >"${OUTDIR}/phase2_commands.txt"

  log_msg "Done. Summary: ${cmp}"
}

main "$@"
