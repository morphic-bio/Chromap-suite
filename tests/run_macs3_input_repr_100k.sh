#!/usr/bin/env bash
# Phase 1: MACS3 input-representation parity (diagnostic). Compares peak sets from:
#   - BAM-derived tagAlign + callpeak -f BED (MACS3 defaults on tagAlign)
#   - fragments.tsv.gz + callpeak -f FRAG
#   - fragments.tsv.gz + callpeak -f FRAG --max-count 1
# Writes phase1_input_repr_summary.tsv (two columns: metric, value) and phase1_commands.txt.
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
GENOME_SIZE="${GENOME_SIZE:-hs}"
TAGALIGN_THREADS="${TAGALIGN_THREADS:-4}"

# shellcheck source=peak_caller_100k_common.sh
source "${SCRIPT_DIR}/peak_caller_100k_common.sh"

RUN_MACS3=1

bed3_stats() {
  local f=$1
  [[ -s "${f}" ]] || {
    echo -e "0\t0\t0"
    return 0
  }
  awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=3 {
    w=$3-$2; if(w>0){n++; t+=w; a[n]=w}
  } END {
    if(n<=0){print 0,0,0; exit}
    asort(a)
    med=(n%2)?a[(n+1)/2]:(a[n/2]+a[n/2+1])/2
    print n,t,med
  }' "${f}"
}

sort_bed3() {
  local in=$1
  local out=$2
  awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=3 { s=$2+0; e=$3+0; if(s<e) print $1,s,e }' "${in}" |
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n >"${out}"
}

emit_pair_metrics() {
  local label_a=$1
  local bed_a=$2
  local label_b=$3
  local bed_b=$4
  local sum=$5
  local tmp="${OUTDIR}/.bedtmp_${RANDOM}"
  if ! command -v bedtools &>/dev/null; then
    {
      echo -e "pair_jaccard_${label_a}__${label_b}\tNA"
      echo -e "pair_frac_${label_a}_with_any_overlap_in_${label_b}\tNA"
      echo -e "pair_frac_${label_b}_with_any_overlap_in_${label_a}\tNA"
    } >>"${sum}"
    return 0
  fi
  bedtools sort -i "${bed_a}" >"${tmp}.a.s"
  bedtools sort -i "${bed_b}" >"${tmp}.b.s"
  local jac
  jac="$(bedtools jaccard -a "${tmp}.a.s" -b "${tmp}.b.s" 2>/dev/null | awk 'NR==2 {print $3}')"
  [[ -n "${jac}" ]] || jac="NA"
  local na nb ha hb
  na=$(wc -l <"${bed_a}" | tr -d ' ')
  nb=$(wc -l <"${bed_b}" | tr -d ' ')
  ha=$(bedtools intersect -a "${bed_a}" -b "${bed_b}" -u 2>/dev/null | wc -l | tr -d ' ')
  hb=$(bedtools intersect -a "${bed_b}" -b "${bed_a}" -u 2>/dev/null | wc -l | tr -d ' ')
  local fa fb
  fa=$(awk -v h="${ha}" -v n="${na}" 'BEGIN{ if(n<=0) print "NA"; else printf "%.8f", h/n }')
  fb=$(awk -v h="${hb}" -v n="${nb}" 'BEGIN{ if(n<=0) print "NA"; else printf "%.8f", h/n }')
  {
    echo -e "pair_jaccard_${label_a}__${label_b}\t${jac}"
    echo -e "pair_frac_${label_a}_with_any_overlap_in_${label_b}\t${fa}"
    echo -e "pair_frac_${label_b}_with_any_overlap_in_${label_a}\t${fb}"
  } >>"${sum}"
  rm -f "${tmp}.a.s" "${tmp}.b.s"
}

main() {
  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/chromap_macs3_input_repr.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}"
  local log="${OUTDIR}/phase1_harness.log"
  log_msg() { echo "$@" | tee -a "${log}"; }

  for x in "${MACS3_BIN}" samtools bedtools; do
    if ! command -v "${x}" &>/dev/null; then
      echo "Missing ${x} on PATH" >&2
      exit 2
    fi
  done

  if ! peak_100k_resolve_inputs; then
    exit 1
  fi
  local fr bam
  fr="${PEAK_100K_FRAG}"
  bam="${PEAK_100K_BAM}"

  local d="${OUTDIR}/phase1_input_repr"
  mkdir -p "${d}/callpeak_tagalign" "${d}/callpeak_frag" "${d}/callpeak_frag_max1"

  local tag="${d}/tagAlign_from_bam.gz"
  log_msg "[1] tagAlign from BAM (paired with fragments)"
  peak_100k_tagalign_from_bam "${bam}" "${tag}" "${TAGALIGN_THREADS}" >>"${log}" 2>&1

  log_msg "[2] macs3 callpeak tagAlign -f BED (defaults)"
  local cmd_tag
  cmd_tag="${MACS3_BIN} callpeak -t ${tag} -f BED -g ${GENOME_SIZE} -n macs3_tagalign_bed_defaults --outdir ${d}/callpeak_tagalign"
  echo "${cmd_tag}" >>"${log}"
  ${MACS3_BIN} callpeak -t "${tag}" -f BED -g "${GENOME_SIZE}" -n macs3_tagalign_bed_defaults \
    --outdir "${d}/callpeak_tagalign" >>"${log}" 2>&1

  log_msg "[3] macs3 callpeak fragments -f FRAG"
  local cmd_frag
  cmd_frag="${MACS3_BIN} callpeak -t ${fr} -f FRAG -g ${GENOME_SIZE} -n macs3_frag --outdir ${d}/callpeak_frag"
  echo "${cmd_frag}" >>"${log}"
  ${MACS3_BIN} callpeak -t "${fr}" -f FRAG -g "${GENOME_SIZE}" -n macs3_frag \
    --outdir "${d}/callpeak_frag" >>"${log}" 2>&1

  log_msg "[4] macs3 callpeak fragments -f FRAG --max-count 1"
  local cmd_m1
  cmd_m1="${MACS3_BIN} callpeak -t ${fr} -f FRAG --max-count 1 -g ${GENOME_SIZE} -n macs3_frag_max1 --outdir ${d}/callpeak_frag_max1"
  echo "${cmd_m1}" >>"${log}"
  ${MACS3_BIN} callpeak -t "${fr}" -f FRAG --max-count 1 -g "${GENOME_SIZE}" -n macs3_frag_max1 \
    --outdir "${d}/callpeak_frag_max1" >>"${log}" 2>&1

  local np_t np_f np_m
  np_t="${d}/callpeak_tagalign/macs3_tagalign_bed_defaults_peaks.narrowPeak"
  np_f="${d}/callpeak_frag/macs3_frag_peaks.narrowPeak"
  np_m="${d}/callpeak_frag_max1/macs3_frag_max1_peaks.narrowPeak"

  sort_bed3 "${np_t}" "${d}/tagalign.bed3"
  sort_bed3 "${np_f}" "${d}/frag.bed3"
  sort_bed3 "${np_m}" "${d}/frag_max1.bed3"

  local nt tt mtw nfw tfw mfw nmw tmw mmw
  IFS=$'\t' read -r nt tt mtw < <(bed3_stats "${d}/tagalign.bed3")
  IFS=$'\t' read -r nfw tfw mfw < <(bed3_stats "${d}/frag.bed3")
  IFS=$'\t' read -r nmw tmw mmw < <(bed3_stats "${d}/frag_max1.bed3")

  local sum="${OUTDIR}/phase1_input_repr_summary.tsv"
  {
    echo -e "metric\tvalue"
    echo -e "phase\tphase1_input_repr"
    echo -e "PEAK_100K_RUN_ROOT\t${PEAK_100K_RUN_ROOT:-}"
    echo -e "FRAGMENTS_TSV_GZ\t${fr}"
    echo -e "ATAC_BAM\t${bam}"
    echo -e "tagAlign_gz\t${tag}"
    echo -e "macs3_tagalign_peak_count\t${nt}"
    echo -e "macs3_tagalign_total_bp\t${tt}"
    echo -e "macs3_tagalign_median_width_bp\t${mtw}"
    echo -e "macs3_frag_peak_count\t${nfw}"
    echo -e "macs3_frag_total_bp\t${tfw}"
    echo -e "macs3_frag_median_width_bp\t${mfw}"
    echo -e "macs3_frag_max1_peak_count\t${nmw}"
    echo -e "macs3_frag_max1_total_bp\t${tmw}"
    echo -e "macs3_frag_max1_median_width_bp\t${mmw}"
  } | tee "${sum}"

  emit_pair_metrics "tagalign_bed_defaults" "${d}/tagalign.bed3" "frag_FRAG" "${d}/frag.bed3" "${sum}"
  emit_pair_metrics "tagalign_bed_defaults" "${d}/tagalign.bed3" "frag_FRAG_max1" "${d}/frag_max1.bed3" "${sum}"
  emit_pair_metrics "frag_FRAG" "${d}/frag.bed3" "frag_FRAG_max1" "${d}/frag_max1.bed3" "${sum}"

  {
    echo "PEAK_100K_RUN_ROOT=${PEAK_100K_RUN_ROOT:-}"
    echo "FRAGMENTS_TSV_GZ=${fr}"
    echo "ATAC_BAM=${bam}"
    echo "CMD_TAGALIGN=peak_100k_tagalign_from_bam (see peak_caller_100k_common.sh)"
    echo "CMD_CALLPEAK_TAGALIGN=${cmd_tag}"
    echo "CMD_CALLPEAK_FRAG=${cmd_frag}"
    echo "CMD_CALLPEAK_FRAG_MAX1=${cmd_m1}"
  } >"${OUTDIR}/phase1_commands.txt"

  log_msg "Done. Summary: ${sum}  Commands: ${OUTDIR}/phase1_commands.txt"
}

main "$@"
