#!/usr/bin/env bash
# Phase 3: MACS3 `callpeak -f FRAG -B` treat_pileup + control_lambda vs chromap_callpeaks
# diagnostic outputs (fragment pileup + no-control lambda). Shared fragment discovery
# with other 100K harnesses; RUN_MACS3=0 allows frags without paired BAM.
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

# shellcheck source=peak_caller_100k_common.sh
source "${SCRIPT_DIR}/peak_caller_100k_common.sh"

main() {
  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/chromap_lambda_parity.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}"
  local log="${OUTDIR}/phase3_harness.log"
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

  local d="${OUTDIR}/phase3_lambda_parity"
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

  local caller="${REPO_ROOT}/chromap_callpeaks"
  local cpp_treat="${d}/cpp_frag_treat_pileup.bdg"
  local cpp_lambda="${d}/cpp_frag_control_lambda.bdg"
  # Match MACS3 `callpeak -g hs` (GRCh38 effective 2913022398) and default llocal 10000.
  local cmd_cpp
  cmd_cpp="${caller} -i ${fr} --frag-span-pileup-bdg ${cpp_treat} --frag-pileup-macs3-uint8-counts --frag-lambda-bdg ${cpp_lambda} --frag-lambda-effective-genome 2913022398"
  log_msg "[2] ${cmd_cpp}"
  "${caller}" -i "${fr}" --frag-span-pileup-bdg "${cpp_treat}" --frag-pileup-macs3-uint8-counts \
    --frag-lambda-bdg "${cpp_lambda}" --frag-lambda-effective-genome 2913022398 >>"${log}" 2>&1

  local py="${SCRIPT_DIR}/compare_pileup_bdg.py"
  local sum_treat sum_lambda
  sum_treat="${d}/lambda_parity_treat_compare.tsv"
  sum_lambda="${d}/lambda_parity_control_lambda_compare.tsv"
  log_msg "[3] compare C++ treat vs MACS3 treat_pileup"
  python3 "${py}" "${macs_treat}" "${cpp_treat}" \
    --high-prefix "${d}/treat_high" --largest-n-diff 5 --largest-out "${d}/largest_treat_diff.tsv" | tee "${sum_treat}"
  log_msg "[4] compare C++ lambda vs MACS3 control_lambda"
  python3 "${py}" "${macs_lambda}" "${cpp_lambda}" \
    --high-prefix "${d}/lambda_high" --largest-n-diff 5 --largest-out "${d}/largest_lambda_diff.tsv" | tee "${sum_lambda}"

  local ht_hj
  ht_hj="NA"
  if command -v bedtools &>/dev/null; then
    if [[ -f "${d}/treat_high_a.bed" && -f "${d}/treat_high_b.bed" ]]; then
      bedtools sort -i "${d}/treat_high_a.bed" >"${d}/treat_ha.s.bed"
      bedtools sort -i "${d}/treat_high_b.bed" >"${d}/treat_hb.s.bed"
      ht_hj="$(bedtools jaccard -a "${d}/treat_ha.s.bed" -b "${d}/treat_hb.s.bed" 2>/dev/null | awk 'NR==2 {print $3}')"
    fi
  fi
  echo -e "high_signal_treat_jaccard\t${ht_hj}" | tee -a "${sum_treat}"

  local hl_hj
  hl_hj="NA"
  if command -v bedtools &>/dev/null; then
    if [[ -f "${d}/lambda_high_a.bed" && -f "${d}/lambda_high_b.bed" ]]; then
      bedtools sort -i "${d}/lambda_high_a.bed" >"${d}/lambda_ha.s.bed"
      bedtools sort -i "${d}/lambda_high_b.bed" >"${d}/lambda_hb.s.bed"
      hl_hj="$(bedtools jaccard -a "${d}/lambda_ha.s.bed" -b "${d}/lambda_hb.s.bed" 2>/dev/null | awk 'NR==2 {print $3}')"
    fi
  fi
  echo -e "high_mass_lambda_jaccard\t${hl_hj}" | tee -a "${sum_lambda}"

  local treat_ok lambda_ok
  treat_ok="$(awk -F'\t' '$1=="bdg_rows_identical_sorted"{print $2}' "${sum_treat}")"
  lambda_ok="$(awk -F'\t' '$1=="bdg_rows_identical_sorted"{print $2}' "${sum_lambda}")"
  local pear_l
  pear_l="$(awk -F'\t' '$1=="pearson_bp_weighted_either_nonempty"{print $2}' "${sum_lambda}")"

  local summary="${d}/phase3_lambda_summary.tsv"
  {
    echo "metric	value"
    echo "phase	3_background_lambda_parity"
    echo "treat_compare_tsv	${sum_treat}"
    echo "lambda_compare_tsv	${sum_lambda}"
    echo "bdg_treat_rows_match	${treat_ok}"
    echo "bdg_lambda_rows_match	${lambda_ok}"
    echo "lambda_pearson_either_nonempty	${pear_l}"
    echo "note	High-lambda BEDs use the same high-frac mass heuristic as high-signal (|value|*bp)"
    echo "conclusion	see harness stdout"
  } | tee "${summary}"

  {
    echo "PEAK_100K_RUN_ROOT=${PEAK_100K_RUN_ROOT:-}"
    echo "FRAGMENTS_TSV_GZ=${fr}"
    echo "CMD_MACS3_CALLPEAK_B=${cmd_macs}"
    echo "CMD_CHROMAP_PILEUP_LAMBDA=${cmd_cpp}"
  } >"${OUTDIR}/phase3_commands.txt"

  log_msg "Done. Summary TSV: ${summary}"
  log_msg "Output directory: ${OUTDIR}"
  if [[ "${treat_ok}" != "True" ]]; then
    echo "FAIL: C++ treat pileup != MACS3 treat_pileup (see ${sum_treat})" >&2
    exit 1
  fi
  if [[ "${lambda_ok}" == "True" ]]; then
    log_msg "Lambda bedGraph: exact row match to MACS3 control_lambda."
  else
    log_msg "Lambda bedGraph: differs from MACS3; see ${sum_lambda} and ${d}/largest_lambda_diff.tsv"
    log_msg "Likely follow-up: port MACS3 __chrom_pair_treat_ctrl (min-merge of treat/ctrl p,v arrays, with tail) then __write 1e-5; a sorted union of breakpoints is not equivalent to MACS3's paired stream."
  fi
  exit 0
}

main "$@"
