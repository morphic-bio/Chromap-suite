#!/usr/bin/env bash
# Validation harness: internal fragment peak caller on the 100K multiome ATAC
# fixture. Optional: MACS3 (coordination script) and ARC peak comparisons.
# Next phase: direct Chromap/STAR integration is out of scope for this script.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

BENCH_ROOT="${CHROMAP_100K_BENCH:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k}"
FRAGMENTS_TSV_GZ="${FRAGMENTS_TSV_GZ:-}"
ATAC_BAM="${ATAC_BAM:-}"
ARC_PEAKS_BED="${ARC_PEAKS_BED:-}"
OUTDIR="${OUTDIR:-}"
MACS3_RUNNER="${MACS3_RUNNER:-/mnt/pikachu/multiomic-atac-scrna/scripts/run_macs3_peak_call_from_bam.sh}"
# Set RUN_MACS3=0 to skip external MACS3 (faster; Jaccard columns stay NA).
RUN_MACS3="${RUN_MACS3:-1}"

discover_fragments() {
  if [[ -n "${FRAGMENTS_TSV_GZ}" && -f "${FRAGMENTS_TSV_GZ}" ]]; then
    printf '%s' "${FRAGMENTS_TSV_GZ}"
    return 0
  fi
  if [[ -d "${BENCH_ROOT}" ]]; then
    local f
    f="$(find "${BENCH_ROOT}" \( -path '*/dual/fragments.tsv.gz' -o -path '*/fragments.tsv.gz' \) -type f 2>/dev/null | sort | head -1 || true)"
    if [[ -n "${f}" && -f "${f}" ]]; then
      printf '%s' "${f}"
      return 0
    fi
  fi
  return 1
}

discover_bam() {
  if [[ -n "${ATAC_BAM}" && -f "${ATAC_BAM}" ]]; then
    printf '%s' "${ATAC_BAM}"
    return 0
  fi
  local b
  b="$(find "${BENCH_ROOT}" -name 'atac_possorted_bam.bam' -type f 2>/dev/null | sort | head -1 || true)"
  if [[ -n "${b}" && -f "${b}" ]]; then
    printf '%s' "${b}"
    return 0
  fi
  return 1
}

discover_arc_peaks() {
  if [[ -n "${ARC_PEAKS_BED}" && -f "${ARC_PEAKS_BED}" ]]; then
    printf '%s' "${ARC_PEAKS_BED}"
    return 0
  fi
  local p
  p="${BENCH_ROOT}/pbmc_unsorted_3k_100k_arc/outs/atac_peaks.bed"
  if [[ -f "${p}" ]]; then
    printf '%s' "${p}"
    return 0
  fi
  return 1
}

validate_bed3_sorted() {
  local f=$1
  local label=$2
  [[ -f "${f}" ]] || { echo "Missing ${label}: ${f}" >&2; return 1; }
  if [[ ! -s "${f}" ]]; then
    echo "Empty BED: ${f}" >&2
    return 1
  fi
  awk 'BEGIN{FS=OFS="\t"} !/^#/ && !/^track/ { if(NF<3) exit 2; s=$2+0; e=$3+0; if(s>=e) exit 3; print $1, s, e }' "${f}" > "${f}.3clean" || {
    echo "Invalid intervals in ${label} (${f})" >&2
    return 1
  }
  LC_ALL=C sort -c -k1,1 -k2,2n -k3,3n "${f}.3clean" 2>/dev/null || {
    echo "Not sorted: ${f}" >&2
    return 1
  }
  return 0
}

jaccard_bed3() {
  local a=$1
  local b=$2
  if ! command -v bedtools &>/dev/null; then
    echo "NA"
    return 0
  fi
  bedtools sort -i "${a}" > "${a}.s" 2>/dev/null
  bedtools sort -i "${b}" > "${b}.s" 2>/dev/null
  bedtools jaccard -a "${a}.s" -b "${b}.s" 2>/dev/null | awk 'NR==2 {print $3}'
}

median_width() {
  local f=$1
  awk '{print $3-$2}' "${f}" | sort -n | awk '{
    a[NR]=$0
  } END {
    if(NR==0) { print 0; exit }
    if(NR%2) print a[(NR+1)/2]
    else print (a[NR/2]+a[NR/2+1])/2
  }'
}

main() {
  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/chromap_peaks_100k.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}"
  local log="${OUTDIR}/harness.log"
  log_msg() { echo "$@" | tee -a "${log}"; }

  log_msg "[1] make chromap_callpeaks"
  (cd "${REPO_ROOT}" && make chromap_callpeaks) >> "${log}" 2>&1

  if ! fr="$(discover_fragments)"; then
    echo "Set FRAGMENTS_TSV_GZ or point CHROMAP_100K_BENCH to a directory containing fragments.tsv(.gz)" >&2
    exit 1
  fi
  local caller="${REPO_ROOT}/chromap_callpeaks"
  local prefix="${OUTDIR}/internal"
  {
    echo "REPO_ROOT=${REPO_ROOT}"
    echo "FRAGMENTS_TSV_GZ=${fr}"
    echo "BENCH_ROOT=${BENCH_ROOT}"
  } >> "${log}"

  log_msg "CMD_INTERNAL=${caller} -i <frag> --out-prefix <prefix> ..."
  "${caller}" -i "${fr}" --out-prefix "${prefix}" --bin-width 50 --ext-size 150 \
    --local-window 200 --fdr 0.05 --p-value 0.01 --merge-gap 0 >> "${log}" 2>&1

  if [[ ! -s "${prefix}.narrowPeak" ]]; then
    echo "Empty internal narrowPeak" >&2
    exit 1
  fi
  awk 'BEGIN{FS=OFS="\t"} NF>=3 && $0!~/^#/ { if($2+0 >= $3+0) { print "Invalid interval line "NR": "$0 > "/dev/stderr"; exit 1 } }' "${prefix}.narrowPeak"
  if ! LC_ALL=C sort -c -k1,1 -k2,2n -k3,3n "${prefix}.narrowPeak" 2>/dev/null; then
    echo "internal narrowPeak not sorted (C locale, chrom then start)" >&2
    exit 1
  fi

  awk 'BEGIN{FS=OFS="\t"} {if(NF<10) exit 1} END{exit 0}' "${prefix}.narrowPeak" || {
    echo "expected >=10 column narrowPeak" >&2
    exit 1
  }
  if ! validate_bed3_sorted "${prefix}.narrowPeak" "internal_narrowPeak"; then
    exit 1
  fi

  local n_peaks
  n_peaks=$(wc -l < "${prefix}.narrowPeak" | tr -d ' ')
  local internal3="${OUTDIR}/internal.3.bed"
  awk 'BEGIN{FS=OFS="\t"} !/^#/{print $1,$2,$3}' "${prefix}.narrowPeak" > "${internal3}"
  local total_bp med
  total_bp=$(awk 'BEGIN{FS=OFS="\t"; t=0} {t+=$3-$2} END{print t}' "${internal3}")
  med=$(median_width "${internal3}")

  JMAC_DEF="NA"
  JMAC_ST="NA"
  local bam
  bam="$(discover_bam || true)"
  if [[ "${RUN_MACS3}" == "0" ]]; then
    log_msg "Skipping MACS3: RUN_MACS3=0"
  elif [[ -z "${bam}" || ! -f "${MACS3_RUNNER}" || ! -x "${MACS3_RUNNER}" ]]; then
    log_msg "Skipping MACS3: no atac BAM under BENCH_ROOT, or ${MACS3_RUNNER} missing or not executable"
  else
    log_msg "[2] optional MACS3 (defaults + shiftTAG) via ${MACS3_RUNNER}"
    echo "ATAC_BAM=${bam}" >> "${log}"
    if command -v macs3 &>/dev/null && command -v samtools &>/dev/null && command -v bedtools &>/dev/null; then
      mkdir -p "${OUTDIR}/macs3"
      "${MACS3_RUNNER}" --input-bam "${bam}" --outdir "${OUTDIR}/macs3" --prefix atac_100k --mode both >> "${log}" 2>&1
      for suf in "macs3_defaults" "macs3_shiftTAG"; do
        if [[ -f "${OUTDIR}/macs3/atac_100k.${suf}_peaks.narrowPeak" ]]; then
          validate_bed3_sorted "${OUTDIR}/macs3/atac_100k.${suf}_peaks.narrowPeak" "macs3_${suf}" || true
        fi
      done
      if [[ -f "${OUTDIR}/macs3/atac_100k.macs3_defaults_peaks.narrowPeak" ]]; then
        awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' "${OUTDIR}/macs3/atac_100k.macs3_defaults_peaks.narrowPeak" > "${OUTDIR}/mdef3.bed"
        JMAC_DEF="$(jaccard_bed3 "${internal3}" "${OUTDIR}/mdef3.bed" || echo NA)"
      fi
      if [[ -f "${OUTDIR}/macs3/atac_100k.macs3_shiftTAG_peaks.narrowPeak" ]]; then
        awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3}' "${OUTDIR}/macs3/atac_100k.macs3_shiftTAG_peaks.narrowPeak" > "${OUTDIR}/mst3.bed"
        JMAC_ST="$(jaccard_bed3 "${internal3}" "${OUTDIR}/mst3.bed" || echo NA)"
      fi
    else
      log_msg "Skipping MACS3: need macs3, samtools, bedtools on PATH"
    fi
  fi

  ARC_FRA_INT="NA"
  INT_FRA_ARC="NA"
  local arcbed
  if arcbed="$(discover_arc_peaks 2>/dev/null)"; then
    if command -v bedtools &>/dev/null; then
      awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=3 { s=$2+0; e=$3+0; if(s<e) print $1,s,e }' "${arcbed}" > "${OUTDIR}/arc.3raw.bed"
      if [[ -s "${OUTDIR}/arc.3raw.bed" ]]; then
        bedtools sort -i "${OUTDIR}/arc.3raw.bed" > "${OUTDIR}/arc.3.bed"
        arc_n=$(wc -l < "${OUTDIR}/arc.3.bed" | tr -d ' ')
        int_n=$(wc -l < "${internal3}" | tr -d ' ')
        if [[ "${arc_n}" -gt 0 && "${int_n}" -gt 0 ]]; then
          hit_a=$(bedtools intersect -a "${OUTDIR}/arc.3.bed" -b "${internal3}" -u 2>/dev/null | wc -l | tr -d ' ')
          hit_i=$(bedtools intersect -a "${internal3}" -b "${OUTDIR}/arc.3.bed" -u 2>/dev/null | wc -l | tr -d ' ')
          ARC_FRA_INT=$(awk -v h="${hit_a}" -v n="${arc_n}" 'BEGIN{ if(n<=0) print "NA"; else printf "%.6f", h/n }')
          INT_FRA_ARC=$(awk -v h="${hit_i}" -v n="${int_n}" 'BEGIN{ if(n<=0) print "NA"; else printf "%.6f", h/n }')
        fi
      fi
    fi
  fi

  local summary="${OUTDIR}/benchmark_summary.tsv"
  {
    echo -e "metric\tvalue"
    echo -e "internal_peak_count\t${n_peaks}"
    echo -e "internal_total_merged_bp\t${total_bp}"
    echo -e "internal_median_width_bp\t${med}"
    echo -e "jaccard_1bp_internal_vs_macs3_defaults\t${JMAC_DEF}"
    echo -e "jaccard_1bp_internal_vs_macs3_shiftTAG\t${JMAC_ST}"
    echo -e "frac_arc_peaks_with_any_internal_overlap\t${ARC_FRA_INT}"
    echo -e "frac_internal_peaks_with_any_arc_overlap\t${INT_FRA_ARC}"
  } | tee "${summary}"
  {
    echo "FRAGMENTS_TSV_GZ=${fr}"
    echo "internal_prefix=${prefix}"
    echo "make_target=chromap_callpeaks"
  } > "${OUTDIR}/commands.txt"
  log_msg "Done. Output: ${OUTDIR}"
}

main "$@"
