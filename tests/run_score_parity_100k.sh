#!/usr/bin/env bash
# Phase 4: MACS3 `bdgcmp -m ppois FE` vs chromap_callpeaks diagnostic score bedGraphs.
# Uses the same 100K fragment discovery as other peak_caller harnesses.
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
    OUTDIR="$(mktemp -d /tmp/chromap_score_parity.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}"
  local log="${OUTDIR}/phase4_harness.log"
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

  local d="${OUTDIR}/phase4_score_parity"
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
  local macs_score_f="${d}/macs3_scores_FE.bdg"
  local cmd_bdgcmp
  # MACS3 argparse: repeat `-m` overwrites; one `-m` with multiple args runs all methods.
  cmd_bdgcmp="${MACS3_BIN} bdgcmp -t ${macs_treat} -c ${macs_lambda} -m ppois FE --outdir ${d} --o-prefix macs3_scores -p 0"
  log_msg "[2] ${cmd_bdgcmp}"
  ${MACS3_BIN} bdgcmp -t "${macs_treat}" -c "${macs_lambda}" -m ppois FE --outdir "${d}" \
    --o-prefix macs3_scores -p 0 >>"${log}" 2>&1

  if [[ ! -f "${macs_score_p}" || ! -f "${macs_score_f}" ]]; then
    echo "Missing MACS3 bdgcmp outputs (expected ${macs_score_p} and ${macs_score_f})" >&2
    exit 1
  fi

  local caller="${REPO_ROOT}/chromap_callpeaks"
  local cpp_p="${d}/cpp_frag_score_ppois.bdg"
  local cpp_f="${d}/cpp_frag_score_FE.bdg"
  local cmd_cpp
  cmd_cpp="${caller} -i ${fr} --frag-pileup-macs3-uint8-counts --frag-score-ppois-bdg ${cpp_p} --frag-score-fe-bdg ${cpp_f} --frag-lambda-effective-genome 2913022398 --frag-score-pseudocount 0"
  log_msg "[3] ${cmd_cpp}"
  "${caller}" -i "${fr}" --frag-pileup-macs3-uint8-counts \
    --frag-score-ppois-bdg "${cpp_p}" --frag-score-fe-bdg "${cpp_f}" \
    --frag-lambda-effective-genome 2913022398 --frag-score-pseudocount 0 >>"${log}" 2>&1

  local py="${SCRIPT_DIR}/compare_pileup_bdg.py"
  local sum_p sum_f
  sum_p="${d}/score_parity_ppois_compare.tsv"
  sum_f="${d}/score_parity_fe_compare.tsv"
  log_msg "[4] compare C++ ppois vs MACS3 bdgcmp"
  python3 "${py}" "${macs_score_p}" "${cpp_p}" \
    --high-prefix "${d}/ppois_high" --largest-n-diff 5 --largest-out "${d}/largest_ppois_diff.tsv" \
    --exact-row-tol 1e-4 | tee "${sum_p}"
  log_msg "[5] compare C++ FE vs MACS3 bdgcmp"
  python3 "${py}" "${macs_score_f}" "${cpp_f}" \
    --high-prefix "${d}/fe_high" --largest-n-diff 5 --largest-out "${d}/largest_fe_diff.tsv" \
    --exact-row-tol 1e-4 | tee "${sum_f}"

  local hj_p hj_f
  hj_p="NA"
  hj_f="NA"
  if command -v bedtools &>/dev/null; then
    if [[ -f "${d}/ppois_high_a.bed" && -f "${d}/ppois_high_b.bed" ]]; then
      bedtools sort -i "${d}/ppois_high_a.bed" >"${d}/ppois_ha.s.bed"
      bedtools sort -i "${d}/ppois_high_b.bed" >"${d}/ppois_hb.s.bed"
      hj_p="$(bedtools jaccard -a "${d}/ppois_ha.s.bed" -b "${d}/ppois_hb.s.bed" 2>/dev/null | awk 'NR==2 {print $3}')"
    fi
    if [[ -f "${d}/fe_high_a.bed" && -f "${d}/fe_high_b.bed" ]]; then
      bedtools sort -i "${d}/fe_high_a.bed" >"${d}/fe_ha.s.bed"
      bedtools sort -i "${d}/fe_high_b.bed" >"${d}/fe_hb.s.bed"
      hj_f="$(bedtools jaccard -a "${d}/fe_ha.s.bed" -b "${d}/fe_hb.s.bed" 2>/dev/null | awk 'NR==2 {print $3}')"
    fi
  fi
  echo -e "high_score_ppois_jaccard\t${hj_p}" | tee -a "${sum_p}"
  echo -e "high_score_fe_jaccard\t${hj_f}" | tee -a "${sum_f}"

  local ok_p ok_f
  ok_p="$(awk -F'\t' '$1=="bdg_rows_identical_sorted"{print $2}' "${sum_p}")"
  ok_f="$(awk -F'\t' '$1=="bdg_rows_identical_sorted"{print $2}' "${sum_f}")"

  local pear_p pear_f maxp maxf rows_pa rows_pb rows_fa rows_fb
  pear_p="$(awk -F'\t' '$1=="pearson_bp_weighted_either_nonempty"{print $2}' "${sum_p}")"
  pear_f="$(awk -F'\t' '$1=="pearson_bp_weighted_either_nonempty"{print $2}' "${sum_f}")"
  maxp="$(awk -F'\t' '$1=="max_abs_diff_segment_either_nonempty"{print $2}' "${sum_p}")"
  maxf="$(awk -F'\t' '$1=="max_abs_diff_segment_either_nonempty"{print $2}' "${sum_f}")"
  rows_pa="$(awk -F'\t' '$1=="bdg_row_count_a"{print $2}' "${sum_p}")"
  rows_pb="$(awk -F'\t' '$1=="bdg_row_count_b"{print $2}' "${sum_p}")"
  rows_fa="$(awk -F'\t' '$1=="bdg_row_count_a"{print $2}' "${sum_f}")"
  rows_fb="$(awk -F'\t' '$1=="bdg_row_count_b"{print $2}' "${sum_f}")"

  local pass_p pass_f
  pass_p=0
  pass_f=0
  if [[ "${ok_p}" == "True" ]]; then
    pass_p=1
  else
    local rd
    rd=$((rows_pb > rows_pa ? rows_pb - rows_pa : rows_pa - rows_pb))
    if awk -v p="${pear_p}" -v m="${maxp}" -v r="${rd}" 'BEGIN {
      exit !((p+0)>=0.99999 && (m+0)<=0.0005 && (r+0)<=80)
    }'; then
      pass_p=1
    fi
  fi
  if [[ "${ok_f}" == "True" ]]; then
    pass_f=1
  else
    local rdf
    rdf=$((rows_fb > rows_fa ? rows_fb - rows_fa : rows_fa - rows_fb))
    if awk -v p="${pear_f}" -v m="${maxf}" -v r="${rdf}" 'BEGIN {
      exit !((p+0)>=0.99999 && (m+0)<=0.02 && (r+0)<=30)
    }'; then
      pass_f=1
    fi
  fi

  local summary="${d}/phase4_score_summary.tsv"
  {
    echo "metric	value"
    echo "phase	4_bdgcmp_score_parity"
    echo "ppois_compare_tsv	${sum_p}"
    echo "fe_compare_tsv	${sum_f}"
    echo "bdg_ppois_rows_match	${ok_p}"
    echo "bdg_fe_rows_match	${ok_f}"
    echo "ppois_pass_engineering	${pass_p}"
    echo "fe_pass_engineering	${pass_f}"
    echo "conclusion	see harness stdout"
  } | tee "${summary}"

  {
    echo "PEAK_100K_RUN_ROOT=${PEAK_100K_RUN_ROOT:-}"
    echo "FRAGMENTS_TSV_GZ=${fr}"
    echo "CMD_MACS3_CALLPEAK_B=${cmd_macs}"
    echo "CMD_MACS3_BDGCM=${cmd_bdgcmp}"
    echo "CMD_CHROMAP_SCORES=${cmd_cpp}"
  } >"${OUTDIR}/phase4_commands.txt"

  log_msg "Done. Summary TSV: ${summary}"
  log_msg "Output directory: ${OUTDIR}"

  if [[ "${pass_p}" != "1" ]]; then
    echo "FAIL: C++ ppois vs MACS3 bdgcmp outside tolerance (see ${sum_p})" >&2
    exit 1
  fi
  if [[ "${pass_f}" != "1" ]]; then
    echo "FAIL: C++ FE vs MACS3 bdgcmp outside tolerance (see ${sum_f})" >&2
    exit 1
  fi
  if [[ "${ok_p}" == "True" && "${ok_f}" == "True" ]]; then
    log_msg "Score bedGraphs: exact sorted row match to MACS3 bdgcmp (ppois and FE)."
  else
    log_msg "Score bedGraphs: engineering parity (high Pearson, low MAE/max-abs, small row delta). See ${summary}."
  fi
  exit 0
}

main "$@"
