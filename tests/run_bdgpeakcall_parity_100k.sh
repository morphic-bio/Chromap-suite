#!/usr/bin/env bash
# Phase 5: MACS3 `bdgpeakcall` vs C++ diagnostic region caller on ppois score bedGraph.
# Primary: same MACS3 ppois input -> MACS3 bdgpeakcall vs C++ (isolates region logic).
# Secondary: C++ ppois bedGraph + C++ bdgpeakcall vs MACS3 ppois + MACS3 bdgpeakcall.
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

# shellcheck source=peak_caller_100k_common.sh
source "${SCRIPT_DIR}/peak_caller_100k_common.sh"

bed3_from_macs_narrowpeak() {
  local in_np=$1
  local out_bed=$2
  awk 'BEGIN{FS=OFS="\t"} !/^track/ && !/^#/ && NF>=3 { print $1, $2+0, $3+0 }' "${in_np}" |
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n >"${out_bed}"
}

main() {
  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/chromap_bdgpeakcall_parity.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}"
  local log="${OUTDIR}/phase5_harness.log"
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

  local d="${OUTDIR}/phase5_bdgpeakcall"
  mkdir -p "${d}/macs3_callpeak_bdg"
  local macs_dir="${d}/macs3_callpeak_bdg"
  local n="macs3_frag_bdg"

  local cmd_macs
  cmd_macs="${MACS3_BIN} callpeak -t ${fr} -f FRAG -g hs -B -n ${n} --outdir ${macs_dir}"
  log_msg "[1] ${cmd_macs}"
  ${MACS3_BIN} callpeak -t "${fr}" -f FRAG -g hs -B -n "${n}" --outdir "${macs_dir}" >>"${log}" 2>&1

  local macs_treat="${macs_dir}/${n}_treat_pileup.bdg"
  local macs_lambda="${macs_dir}/${n}_control_lambda.bdg"
  for p in "${macs_treat}" "${macs_lambda}"; do
    if [[ ! -f "$p" ]]; then
      echo "Missing expected MACS3 output: $p" >&2
      exit 1
    fi
  done

  local macs_score_p="${d}/macs3_scores_ppois.bdg"
  local cmd_bdgcmp
  cmd_bdgcmp="${MACS3_BIN} bdgcmp -t ${macs_treat} -c ${macs_lambda} -m ppois FE --outdir ${d} --o-prefix macs3_scores -p 0"
  log_msg "[2] ${cmd_bdgcmp}"
  ${MACS3_BIN} bdgcmp -t "${macs_treat}" -c "${macs_lambda}" -m ppois FE --outdir "${d}" \
    --o-prefix macs3_scores -p 0 >>"${log}" 2>&1

  if [[ ! -f "${macs_score_p}" ]]; then
    echo "Missing MACS3 ppois bedGraph: ${macs_score_p}" >&2
    exit 1
  fi

  local macs_np="${d}/macs3_bdgpeakcall_on_macs3_ppois.narrowPeak"
  local cmd_bdgpc
  cmd_bdgpc="${MACS3_BIN} bdgpeakcall -i ${macs_score_p} -c ${CUTOFF} -l ${MINLEN} -g ${MAXGAP} -o $(basename "${macs_np}") --outdir ${d}"
  log_msg "[3] ${cmd_bdgpc}"
  ${MACS3_BIN} bdgpeakcall -i "${macs_score_p}" -c "${CUTOFF}" -l "${MINLEN}" -g "${MAXGAP}" \
    -o "$(basename "${macs_np}")" --outdir "${d}" >>"${log}" 2>&1

  if [[ ! -f "${macs_np}" ]]; then
    echo "Missing MACS3 bdgpeakcall output ${macs_np}" >&2
    exit 1
  fi

  local macs_bed3="${d}/macs3_regions_on_macs3_ppois.bed3"
  bed3_from_macs_narrowpeak "${macs_np}" "${macs_bed3}"

  local caller="${REPO_ROOT}/chromap_callpeaks"
  local cpp_on_macs="${d}/cpp_bdgpeakcall_on_macs3_ppois.bed"
  local cmd_cpp_macs
  cmd_cpp_macs="${caller} --bdgpeakcall-input ${macs_score_p} --bdgpeakcall-output ${cpp_on_macs} --bdgpeakcall-cutoff ${CUTOFF} --bdgpeakcall-min-len ${MINLEN} --bdgpeakcall-max-gap ${MAXGAP}"
  log_msg "[4] ${cmd_cpp_macs}"
  "${caller}" --bdgpeakcall-input "${macs_score_p}" --bdgpeakcall-output "${cpp_on_macs}" \
    --bdgpeakcall-cutoff "${CUTOFF}" --bdgpeakcall-min-len "${MINLEN}" \
    --bdgpeakcall-max-gap "${MAXGAP}" >>"${log}" 2>&1

  local cpp_on_macs_s="${d}/cpp_on_macs3_ppois.s.bed3"
  LC_ALL=C sort -k1,1 -k2,2n -k3,3n "${cpp_on_macs}" >"${cpp_on_macs_s}"

  local cmp_py="${SCRIPT_DIR}/compare_bdgpeakcall_regions.py"
  local primary_tsv="${d}/bdgpeakcall_primary_same_input_compare.tsv"
  log_msg "[5] compare primary (MACS3 vs C++ on same MACS3 ppois)"
  python3 "${cmp_py}" "${macs_bed3}" "${cpp_on_macs_s}" \
    --label-a macs3_bdgpeakcall --label-b cpp_bdgpeakcall | tee "${primary_tsv}"

  local ident
  ident="$(awk -F'\t' '$1=="bed3_sorted_identical"{print $2}' "${primary_tsv}")"

  local cpp_p="${d}/cpp_frag_score_ppois.bdg"
  local cmd_cpp_scores
  cmd_cpp_scores="${caller} -i ${fr} --frag-pileup-macs3-uint8-counts --frag-score-ppois-bdg ${cpp_p} --frag-lambda-effective-genome 2913022398 --frag-score-pseudocount 0"
  log_msg "[6] ${cmd_cpp_scores}"
  "${caller}" -i "${fr}" --frag-pileup-macs3-uint8-counts \
    --frag-score-ppois-bdg "${cpp_p}" \
    --frag-lambda-effective-genome 2913022398 --frag-score-pseudocount 0 >>"${log}" 2>&1

  local secondary_tsv="${d}/bdgpeakcall_secondary_e2e_compare.tsv"
  local sec_note="ran"
  if [[ ! -f "${cpp_p}" ]]; then
    sec_note="skipped_missing_cpp_ppois"
    echo "metric	value" >"${secondary_tsv}"
    echo "note	${sec_note}" >>"${secondary_tsv}"
  else
    local cpp_e2e="${d}/cpp_bdgpeakcall_on_cpp_ppois.bed"
    log_msg "[7] C++ bdgpeakcall on C++ ppois -> ${cpp_e2e}"
    "${caller}" --bdgpeakcall-input "${cpp_p}" --bdgpeakcall-output "${cpp_e2e}" \
      --bdgpeakcall-cutoff "${CUTOFF}" --bdgpeakcall-min-len "${MINLEN}" \
      --bdgpeakcall-max-gap "${MAXGAP}" >>"${log}" 2>&1
    local cpp_e2e_s="${d}/cpp_on_cpp_ppois.s.bed3"
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n "${cpp_e2e}" >"${cpp_e2e_s}"
    log_msg "[8] compare secondary (MACS3 regions vs C++ on C++ ppois)"
    python3 "${cmp_py}" "${macs_bed3}" "${cpp_e2e_s}" \
      --label-a macs3_on_macs3_ppois --label-b cpp_on_cpp_ppois | tee "${secondary_tsv}"
  fi

  local summary="${d}/phase5_bdgpeakcall_summary.tsv"
  {
    echo "metric	value"
    echo "phase	5_bdgpeakcall_region_parity"
    echo "primary_compare_tsv	${primary_tsv}"
    echo "secondary_compare_tsv	${secondary_tsv}"
    echo "same_input_bed3_sorted_identical	${ident}"
    echo "secondary_note	${sec_note}"
  } | tee "${summary}"

  {
    echo "PEAK_100K_RUN_ROOT=${PEAK_100K_RUN_ROOT:-}"
    echo "FRAGMENTS_TSV_GZ=${fr}"
    echo "CMD_MACS3_CALLPEAK_B=${cmd_macs}"
    echo "CMD_MACS3_BDGCM=${cmd_bdgcmp}"
    echo "CMD_MACS3_BDGPEAKCALL=${cmd_bdgpc}"
    echo "CMD_CPP_BDGPEAKCALL_ON_MACS3_PPOIS=${cmd_cpp_macs}"
    echo "CMD_CPP_SCORE_PPOIS=${cmd_cpp_scores}"
  } >"${OUTDIR}/phase5_commands.txt"

  log_msg "Done. Summary TSV: ${summary}"
  log_msg "Output directory: ${OUTDIR}"

  if [[ "${ident}" != "True" ]]; then
    echo "FAIL: same-input bdgpeakcall BED3 mismatch (see ${primary_tsv}); diagnose region caller." >&2
    exit 1
  fi
  exit 0
}

main "$@"
