#!/usr/bin/env bash
# Integration: chromap --call-macs3-frag-peaks vs chromap_callpeaks + MACS3 BED3 (100K fixture).
# - Verifies sorted fragment lines match dual run without peak calling.
# - Verifies integrated narrowPeak/summits match standalone chromap_callpeaks on the same fragments.
# - Optional (RUN_MACS3=1): BED3 geometry vs MACS3 callpeak -f FRAG -p 1e-5.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

CHROMAP="${CHROMAP:-${REPO_ROOT}/chromap}"
CALLPEAKS="${CALLPEAKS:-${REPO_ROOT}/chromap_callpeaks}"

FIXTURE_ATAC="${FIXTURE_ATAC:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac}"
INDEX="${INDEX:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index}"
REF="${REF:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa}"
WHITELIST="${WHITELIST:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"
OUTDIR="${OUTDIR:-}"
THREADS="${THREADS:-8}"
MACS3_BIN="${MACS3_BIN:-macs3}"
export RUN_MACS3="${RUN_MACS3:-0}"
GENOME_SIZE="${GENOME_SIZE:-hs}"
MINLEN="${MINLEN:-200}"
MAXGAP="${MAXGAP:-30}"

R1="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R1_001.fastq.gz"
R2="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R3_001.fastq.gz"
BC="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R2_001.fastq.gz"

chromap_common_args=(
  -t "${THREADS}" -x "${INDEX}" -r "${REF}" -1 "${R1}" -2 "${R2}" -b "${BC}"
  --barcode-whitelist "${WHITELIST}" -l 2000 --trim-adapters --remove-pcr-duplicates
  --remove-pcr-duplicates-at-cell-level --Tn5-shift --BAM --sort-bam
)

bed3_from_narrow() {
  local in_np=$1
  local out_bed=$2
  awk 'BEGIN{FS=OFS="\t"} !/^track/ && !/^#/ && NF>=3 { print $1, $2+0, $3+0 }' "${in_np}" |
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n >"${out_bed}"
}

main() {
  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/chromap_peak_integration_100k.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}/baseline" "${OUTDIR}/integrated" "${OUTDIR}/standalone_ref" "${OUTDIR}/logs"
  local log="${OUTDIR}/integration_harness.log"
  log_msg() { echo "$@" | tee -a "${log}"; }

  if [[ ! -f "${CHROMAP}" || ! -f "${CALLPEAKS}" ]]; then
    log_msg "Build chromap and chromap_callpeaks at repo root first."
    exit 2
  fi
  if [[ ! -d "${FIXTURE_ATAC}" || ! -f "${INDEX}" || ! -f "${REF}" || ! -f "${WHITELIST}" ]]; then
    log_msg "Missing fixture paths (FIXTURE_ATAC, INDEX, REF, WHITELIST). See script header."
    exit 2
  fi

  (cd "${REPO_ROOT}" && make chromap chromap_callpeaks) >>"${log}" 2>&1

  log_msg "[1] Dual ATAC baseline (no peak calling)"
  "${CHROMAP}" "${chromap_common_args[@]}" \
    --atac-fragments "${OUTDIR}/baseline/fragments.tsv.gz" \
    --summary "${OUTDIR}/baseline/summary.tsv" \
    -o "${OUTDIR}/baseline/possorted_bam.bam" \
    >>"${log}" 2>&1

  log_msg "[2] Dual ATAC + opt-in MACS3 FRAG peaks"
  "${CHROMAP}" "${chromap_common_args[@]}" \
    --atac-fragments "${OUTDIR}/integrated/fragments.tsv.gz" \
    --summary "${OUTDIR}/integrated/summary.tsv" \
    -o "${OUTDIR}/integrated/possorted_bam.bam" \
    --call-macs3-frag-peaks \
    --macs3-frag-peaks-output "${OUTDIR}/integrated/chromap_macs3_frag.narrowPeak" \
    --macs3-frag-summits-output "${OUTDIR}/integrated/chromap_macs3_frag_summits.bed" \
    --macs3-frag-pvalue 1e-5 \
    --macs3-frag-min-length "${MINLEN}" \
    --macs3-frag-max-gap "${MAXGAP}" \
    >>"${log}" 2>&1

  log_msg "[3] Compare sorted fragment multiset (baseline vs integrated)"
  if ! cmp -s <(zcat "${OUTDIR}/baseline/fragments.tsv.gz" | LC_ALL=C sort) \
            <(zcat "${OUTDIR}/integrated/fragments.tsv.gz" | LC_ALL=C sort); then
    echo "FAIL: fragment lines differ with vs without --call-macs3-frag-peaks" >&2
    exit 1
  fi

  log_msg "[4] Standalone chromap_callpeaks on baseline fragments (reference C++ outputs)"
  "${CALLPEAKS}" -i "${OUTDIR}/baseline/fragments.tsv.gz" \
    --frag-pileup-macs3-uint8-counts \
    --macs3-frag-narrowpeak "${OUTDIR}/standalone_ref/cpp_macs3_frag.narrowPeak" \
    --macs3-frag-summits "${OUTDIR}/standalone_ref/cpp_macs3_frag_summits.bed" \
    --bdgpeakcall-cutoff 5 --bdgpeakcall-min-len "${MINLEN}" --bdgpeakcall-max-gap "${MAXGAP}" \
    --frag-score-pseudocount 0 \
    >>"${log}" 2>&1

  log_msg "[5] Integrated narrowPeak/summits vs standalone (expect byte identity)"
  if ! cmp -s "${OUTDIR}/integrated/chromap_macs3_frag.narrowPeak" \
              "${OUTDIR}/standalone_ref/cpp_macs3_frag.narrowPeak"; then
    echo "FAIL: integrated narrowPeak differs from chromap_callpeaks" >&2
    exit 1
  fi
  if ! cmp -s "${OUTDIR}/integrated/chromap_macs3_frag_summits.bed" \
              "${OUTDIR}/standalone_ref/cpp_macs3_frag_summits.bed"; then
    echo "FAIL: integrated summits differ from chromap_callpeaks" >&2
    exit 1
  fi

  local summary_tsv="${OUTDIR}/integration_summary.tsv"
  local macs3_note="skipped_RUN_MACS3_0"
  local bed3_id="na"
  local jac="na"
  local narrowpeak_cmp="${OUTDIR}/narrowpeak_compare.tsv"
  local bed3_cmp="${OUTDIR}/bed3_compare.tsv"

  if [[ "${RUN_MACS3}" != "0" ]]; then
    if ! command -v "${MACS3_BIN}" &>/dev/null; then
      log_msg "RUN_MACS3!=0 but ${MACS3_BIN} not on PATH; skipping MACS3 BED3 check"
      macs3_note="skipped_no_macs3"
    else
      log_msg "[6] MACS3 callpeak -f FRAG -p 1e-5 (BED3 reference)"
      mkdir -p "${OUTDIR}/macs3_callpeak"
      ${MACS3_BIN} callpeak -t "${OUTDIR}/baseline/fragments.tsv.gz" -f FRAG -g "${GENOME_SIZE}" \
        -n macs3_frag -p 1e-5 --min-length "${MINLEN}" --max-gap "${MAXGAP}" \
        --outdir "${OUTDIR}/macs3_callpeak" >>"${log}" 2>&1
      local macs_np="${OUTDIR}/macs3_callpeak/macs3_frag_peaks.narrowPeak"
      local int_bed3="${OUTDIR}/integrated/chromap_macs3_frag.bed3"
      local mac_bed3="${OUTDIR}/macs3_callpeak/macs3_frag.bed3"
      bed3_from_narrow "${OUTDIR}/integrated/chromap_macs3_frag.narrowPeak" "${int_bed3}"
      bed3_from_narrow "${macs_np}" "${mac_bed3}"
      python3 "${SCRIPT_DIR}/compare_bdgpeakcall_regions.py" "${mac_bed3}" "${int_bed3}" \
        --label-a macs3_narrow_bed3 --label-b chromap_integrated_bed3 \
        | tee "${bed3_cmp}" | tee -a "${log}"
      jac="$(awk -F'\t' '$1=="jaccard_interval_bp"{print $2; exit}' "${bed3_cmp}" 2>/dev/null || echo na)"
      if cmp -s "${mac_bed3}" "${int_bed3}"; then
        bed3_id="true"
      else
        bed3_id="false"
      fi
      python3 "${SCRIPT_DIR}/compare_narrowpeak_parity.py" \
        --macs3-np "${macs_np}" \
        --cpp-np "${OUTDIR}/integrated/chromap_macs3_frag.narrowPeak" \
        --out-tsv "${narrowpeak_cmp}" >>"${log}" 2>&1 || true
      macs3_note="ran"
    fi
  fi

  {
    echo "harness	chromap_peak_integration_100k"
    echo "outdir	${OUTDIR}"
    echo "fragments_unchanged_vs_baseline	true"
    echo "integrated_matches_standalone_cpp_narrowPeak	true"
    echo "integrated_matches_standalone_cpp_summits	true"
    echo "macs3_bed3_check	${macs3_note}"
    echo "bed3_sorted_identity	${bed3_id}"
    echo "jaccard_bed3	${jac}"
    if [[ "${macs3_note}" == "ran" ]]; then
      echo "bed3_compare_tsv	${bed3_cmp}"
      echo "narrowpeak_parity_tsv	${narrowpeak_cmp}"
    fi
    echo "sidecar_tsv	${OUTDIR}/integrated/summary.tsv.macs3_frag_peaks.tsv"
    echo "log	${log}"
    echo "note	with_RUN_MACS3=1_see_narrowpeak_parity_tsv_for_p_q_signal_summaries"
  } | tee "${summary_tsv}"

  log_msg "OK. Summary: ${summary_tsv}"
}

main "$@"
