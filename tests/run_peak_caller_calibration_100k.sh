#!/usr/bin/env bash
# Grid search / calibration for chromap_callpeaks on the 100K benchmark.
# Reuses the same input discovery as run_peak_caller_100k.sh (source peak_caller_100k_common.sh).
# Runs MACS3 defaults + shiftTAG once, then many internal-caller invocations.
# Machine-readable TSV: OUTDIR/calibration_summary.tsv
#
# CALIB_MODE=full|reduced|minimal  — preset grids (default: full for this script; make uses reduced).
# Or override with space-separated lists, e.g.:
#   CALIB_BIN_WIDTHS="50 100" CALIB_EXT_SIZES="100 150" ...
# Optional: CALIB_MIN_PEAK_BP, CALIB_MIN_SUMMIT_CUTS passed to every run.
# RUN_MACS3=0 skips MACS3 (Jaccard/frac columns become NA where applicable).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

BENCH_ROOT="${CHROMAP_100K_BENCH:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k}"
CHROMAP_PEAK_RUN_ROOT="${CHROMAP_PEAK_RUN_ROOT:-}"
FRAGMENTS_TSV_GZ="${FRAGMENTS_TSV_GZ:-}"
ATAC_BAM="${ATAC_BAM:-}"
ARC_PEAKS_BED="${ARC_PEAKS_BED:-}"
OUTDIR="${OUTDIR:-}"
MACS3_RUNNER="${MACS3_RUNNER:-/mnt/pikachu/multiomic-atac-scrna/scripts/run_macs3_peak_call_from_bam.sh}"
RUN_MACS3="${RUN_MACS3:-1}"
CALIB_MODE="${CALIB_MODE:-full}"
CALIB_MIN_PEAK_BP="${CALIB_MIN_PEAK_BP:-0}"
CALIB_MIN_SUMMIT_CUTS="${CALIB_MIN_SUMMIT_CUTS:-0}"

# shellcheck source=peak_caller_100k_common.sh
source "${SCRIPT_DIR}/peak_caller_100k_common.sh"

# shellcheck disable=SC2034
tsv_one_line() {
  # Replace tabs and newlines so one TSV field stays on one line
  printf '%s' "$1" | tr '\t\n' '  '
}

load_calib_grid() {
  case "${CALIB_MODE}" in
  full)
    CALIB_BIN_WIDTHS="${CALIB_BIN_WIDTHS:-50 100}"
    CALIB_EXT_SIZES="${CALIB_EXT_SIZES:-100 150}"
    CALIB_LOCAL_WINDOWS="${CALIB_LOCAL_WINDOWS:-200 500 1000}"
    CALIB_P_VALUES="${CALIB_P_VALUES:-0.01 0.001 0.0001}"
    CALIB_FDRS="${CALIB_FDRS:-0.05 0.01}"
    CALIB_MERGE_GAPS="${CALIB_MERGE_GAPS:-0 100 200}"
    ;;
  reduced)
    # ~32 runs: good for make test-peak-calibration-100k
    CALIB_BIN_WIDTHS="${CALIB_BIN_WIDTHS:-50 100}"
    CALIB_EXT_SIZES="${CALIB_EXT_SIZES:-100 150}"
    CALIB_LOCAL_WINDOWS="${CALIB_LOCAL_WINDOWS:-200 500}"
    CALIB_P_VALUES="${CALIB_P_VALUES:-0.01 0.001}"
    CALIB_FDRS="${CALIB_FDRS:-0.05}"
    CALIB_MERGE_GAPS="${CALIB_MERGE_GAPS:-0 100}"
    ;;
  minimal)
    CALIB_BIN_WIDTHS="${CALIB_BIN_WIDTHS:-50}"
    CALIB_EXT_SIZES="${CALIB_EXT_SIZES:-150}"
    CALIB_LOCAL_WINDOWS="${CALIB_LOCAL_WINDOWS:-200 1000}"
    CALIB_P_VALUES="${CALIB_P_VALUES:-0.01 0.001}"
    CALIB_FDRS="${CALIB_FDRS:-0.05}"
    CALIB_MERGE_GAPS="${CALIB_MERGE_GAPS:-0 200}"
    ;;
  *)
    echo "Unknown CALIB_MODE=${CALIB_MODE} (use full, reduced, or minimal)" >&2
    return 1
    ;;
  esac
}

p_label() {
  # Short label for p-value in a name, e.g. 0.0001 -> p1e-4
  awk -v p="$1" 'BEGIN {
    if (p >= 0.1) { printf "p%.0e", p; exit }
    for (e=0; e>=-20; e--) {
      t = 10^e
      if (t <= p) { printf "p1e%d", e; exit }
    }
    printf "p%g", p
  }'
}

run_macs3_reference() {
  local bam=$1
  local wdir=$2
  mkdir -p "${wdir}/macs3"
  if [[ "${RUN_MACS3}" == "0" ]]; then
    return 0
  fi
  if ! command -v macs3 &>/dev/null || ! command -v samtools &>/dev/null || ! command -v bedtools &>/dev/null; then
    return 0
  fi
  if [[ -z "$bam" || ! -f "$bam" || ! -f "${MACS3_RUNNER}" || ! -x "${MACS3_RUNNER}" ]]; then
    return 0
  fi
  "${MACS3_RUNNER}" --input-bam "${bam}" --outdir "${wdir}/macs3" --prefix atac_100k --mode both
  if [[ -f "${wdir}/macs3/atac_100k.macs3_defaults_peaks.narrowPeak" ]]; then
    awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' "${wdir}/macs3/atac_100k.macs3_defaults_peaks.narrowPeak" > "${wdir}/ref_macs3_defaults.3.bed"
  fi
  if [[ -f "${wdir}/macs3/atac_100k.macs3_shiftTAG_peaks.narrowPeak" ]]; then
    awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' "${wdir}/macs3/atac_100k.macs3_shiftTAG_peaks.narrowPeak" > "${wdir}/ref_macs3_shiftTAG.3.bed"
  fi
}

prepare_arc_bed() {
  local wdir=$1
  local arcbed
  if ! arcbed="$(discover_arc_peaks 2>/dev/null)"; then
    return 0
  fi
  if ! command -v bedtools &>/dev/null; then
    return 0
  fi
  awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=3 { s=$2+0; e=$3+0; if(s<e) print $1,s,e }' "${arcbed}" > "${wdir}/arc.3raw.bed" || true
  if [[ -s "${wdir}/arc.3raw.bed" ]]; then
    bedtools sort -i "${wdir}/arc.3raw.bed" > "${wdir}/arc.3.bed"
  fi
}

metrics_one_run() {
  local internal_np=$1
  local wdir=$2
  local run_id="${3:-run}"
  local internal3="${wdir}/tmp_internal_${run_id}.3.bed"
  if [[ ! -s "${internal_np}" ]]; then
    printf '%s' "0" $'\t' "0" $'\t' "0" $'\t' "NA" $'\t' "NA" $'\t' "NA" $'\t' "NA" $'\n'
    return 0
  fi
  awk 'BEGIN{FS=OFS="\t"} !/^#/{print $1,$2,$3}' "${internal_np}" > "${internal3}"
  local n_peaks total med jdef jst arc_frac int_arc
  n_peaks=$(wc -l < "${internal3}" | tr -d ' ')
  total=$(awk 'BEGIN{FS=OFS="\t"; t=0} {t+=$3-$2} END{print t}' "${internal3}")
  med=$(median_width "${internal3}")
  jdef="NA"
  jst="NA"
  if [[ -f "${wdir}/ref_macs3_defaults.3.bed" && -s "${internal3}" ]]; then
    jdef="$(jaccard_bed3 "${internal3}" "${wdir}/ref_macs3_defaults.3.bed" || echo NA)"
  fi
  if [[ -f "${wdir}/ref_macs3_shiftTAG.3.bed" && -s "${internal3}" ]]; then
    jst="$(jaccard_bed3 "${internal3}" "${wdir}/ref_macs3_shiftTAG.3.bed" || echo NA)"
  fi
  arc_frac="NA"
  int_arc="NA"
  if [[ -f "${wdir}/arc.3.bed" && -s "${wdir}/arc.3.bed" ]]; then
    arc_n=$(wc -l < "${wdir}/arc.3.bed" | tr -d ' ')
    int_n="${n_peaks}"
    if [[ "${arc_n}" -gt 0 && "${int_n}" -gt 0 ]]; then
      hit_a=$(bedtools intersect -a "${wdir}/arc.3.bed" -b "${internal3}" -u 2>/dev/null | wc -l | tr -d ' ')
      hit_i=$(bedtools intersect -a "${internal3}" -b "${wdir}/arc.3.bed" -u 2>/dev/null | wc -l | tr -d ' ')
      arc_frac=$(awk -v h="${hit_a}" -v n="${arc_n}" 'BEGIN{ if(n<=0) print "NA"; else printf "%.6f", h/n }')
      int_arc=$(awk -v h="${hit_i}" -v n="${int_n}" 'BEGIN{ if(n<=0) print "NA"; else printf "%.6f", h/n }')
    fi
  fi
  printf '%s' "${n_peaks}" $'\t' "${total}" $'\t' "${med}" $'\t' "${jdef}" $'\t' "${jst}" $'\t' "${arc_frac}" $'\t' "${int_arc}" $'\n'
}

main() {
  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/chromap_peaks_calib_100k.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}"
  local log="${OUTDIR}/calibration.log"
  log_msg() { echo "$@" | tee -a "${log}"; }

  if ! load_calib_grid; then
    exit 1
  fi

  if ! peak_100k_resolve_inputs; then
    exit 1
  fi
  local fr bam
  fr="${PEAK_100K_FRAG}"
  bam="${PEAK_100K_BAM:-}"

  log_msg "[1] make chromap_callpeaks"
  (cd "${REPO_ROOT}" && make chromap_callpeaks) >> "${log}" 2>&1
  local caller="${REPO_ROOT}/chromap_callpeaks"

  log_msg "[2] reference MACS3 (once) and ARC BEDs under ${OUTDIR}"
  run_macs3_reference "${bam}" "${OUTDIR}"
  prepare_arc_bed "${OUTDIR}"

  local summary="${OUTDIR}/calibration_summary.tsv"
  {
    printf '%s' "param_set"
    printf '\t%s' bin_width ext_size local_window p_value fdr merge_gap
    printf '\t%s' min_peak_bp min_summit_cuts
    printf '\t%s' internal_peak_count internal_total_bp internal_median_width_bp
    printf '\t%s' jaccard_vs_macs3_defaults jaccard_vs_macs3_shiftTAG
    printf '\t%s' frac_arc_peaks_with_any_internal_overlap frac_internal_peaks_with_any_arc_overlap
    printf '\t%s' command
    echo
  } > "${summary}.partial"
  log_msg "Header -> ${summary}.partial"

  local ncomb
  ncomb=0
  for bw in ${CALIB_BIN_WIDTHS}; do
    for ex in ${CALIB_EXT_SIZES}; do
      for lw in ${CALIB_LOCAL_WINDOWS}; do
        for pv in ${CALIB_P_VALUES}; do
          for fdr in ${CALIB_FDRS}; do
            for mg in ${CALIB_MERGE_GAPS}; do
              ncomb=$((ncomb + 1))
            done
          done
        done
      done
    done
  done
  log_msg "CALIB_MODE=${CALIB_MODE}  grid_size=${ncomb}  OUTDIR=${OUTDIR}"

  local r=0
  for bw in ${CALIB_BIN_WIDTHS}; do
    for ex in ${CALIB_EXT_SIZES}; do
      for lw in ${CALIB_LOCAL_WINDOWS}; do
        for pv in ${CALIB_P_VALUES}; do
          for fdr in ${CALIB_FDRS}; do
            for mg in ${CALIB_MERGE_GAPS}; do
              r=$((r + 1))
              pl="$(p_label "$pv")"
              pname="bw${bw}_e${ex}_lw${lw}_${pl}_f${fdr}_mg${mg}_mp${CALIB_MIN_PEAK_BP}_ms${CALIB_MIN_SUMMIT_CUTS}"
              sub="${OUTDIR}/grid/${pname}"
              mkdir -p "${sub}"
              prefix="${sub}/internal"
              cmd=("${caller}" -i "${fr}" --out-prefix "${prefix}" --bin-width "${bw}" --ext-size "${ex}" \
                --local-window "${lw}" --fdr "${fdr}" --p-value "${pv}" --merge-gap "${mg}" \
                --min-peak-bp "${CALIB_MIN_PEAK_BP}" --min-summit-cuts "${CALIB_MIN_SUMMIT_CUTS}")
              local cmdstr
              cmdstr=$(printf '%q ' "${cmd[@]}")
              log_msg "[${r}/${ncomb}] ${pname}"
              if ! "${cmd[@]}" >> "${log}" 2>&1; then
                {
                  printf '%s' "${pname}"
                  printf '\t%s' "${bw}" "${ex}" "${lw}" "${pv}" "${fdr}" "${mg}"
                  printf '\t%s' "${CALIB_MIN_PEAK_BP}" "${CALIB_MIN_SUMMIT_CUTS}"
                  printf '\tNA\tNA\tNA\tNA\tNA\tNA\tNA'
                  printf '\t%s' "$(tsv_one_line "FAILED ${cmdstr}")"
                  echo
                } >> "${summary}.partial"
                continue
              fi
              IFS=$'\t' read -r n_peaks total med jdef jst arc_f int_f < <(metrics_one_run "${prefix}.narrowPeak" "${OUTDIR}" "${r}")
              {
                printf '%s' "${pname}"
                printf '\t%s' "${bw}" "${ex}" "${lw}" "${pv}" "${fdr}" "${mg}"
                printf '\t%s' "${CALIB_MIN_PEAK_BP}" "${CALIB_MIN_SUMMIT_CUTS}"
                printf '\t%s' "${n_peaks}" "${total}" "${med}" "${jdef}" "${jst}" "${arc_f}" "${int_f}"
                printf '\t%s' "$(tsv_one_line "${cmdstr}")"
                echo
              } >> "${summary}.partial"
            done
          done
        done
      done
    done
  done

  mv -f "${summary}.partial" "${summary}"
  log_msg "Wrote ${summary} (${ncomb} parameter sets)"
  # Columns: 1=param_set ... 10=internal_peak_count 11=internal_total_bp 12=internal_median_width_bp
  # 13=jaccard_defaults 14=jaccard_shiftTAG ...
  log_msg "Top by jaccard_vs_macs3_defaults (non-NA, head):"
  if command -v sort &>/dev/null; then
    tail -n +2 "${summary}" 2>/dev/null | awk -F'\t' '$13!="NA" && $13!="" {print $13"\t"$0}' | sort -t$'\t' -gr | head -5 | cut -f2- | tee -a "${log}" || true
  fi
  log_msg "Top by jaccard_vs_macs3_shiftTAG (non-NA, head):"
  if command -v sort &>/dev/null; then
    tail -n +2 "${summary}" 2>/dev/null | awk -F'\t' '$14!="NA" && $14!="" {print $14"\t"$0}' | sort -t$'\t' -gr | head -5 | cut -f2- | tee -a "${log}" || true
  fi
}

main "$@"
