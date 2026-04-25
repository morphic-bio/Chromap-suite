#!/usr/bin/env bash
# Compare --macs3-frag-peaks-source file vs memory on the 100K fixture.
# Expects identical fragment lines, narrowPeak, summits; records storage mode from logs.
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

main() {
  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/chromap_peak_memory_source_100k.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}/file" "${OUTDIR}/memory" "${OUTDIR}/logs"
  local log="${OUTDIR}/memory_source_100k.log"
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

  log_msg "[1] file-source integrated peaks"
  /usr/bin/time -f 'elapsed_sec %e max_rss_kb %M' -o "${OUTDIR}/logs/time_file.txt" \
    "${CHROMAP}" "${chromap_common_args[@]}" \
    --atac-fragments "${OUTDIR}/file/fragments.tsv.gz" \
    --summary "${OUTDIR}/file/summary.tsv" \
    -o "${OUTDIR}/file/possorted_bam.bam" \
    --call-macs3-frag-peaks \
    --macs3-frag-peaks-source file \
    --macs3-frag-peaks-output "${OUTDIR}/file/chromap_macs3_frag.narrowPeak" \
    --macs3-frag-summits-output "${OUTDIR}/file/chromap_macs3_frag_summits.bed" \
    --macs3-frag-pvalue 1e-5 \
    --macs3-frag-min-length "${MINLEN}" \
    --macs3-frag-max-gap "${MAXGAP}" \
    >>"${log}" 2>&1

  log_msg "[2] memory-source integrated peaks"
  /usr/bin/time -f 'elapsed_sec %e max_rss_kb %M' -o "${OUTDIR}/logs/time_memory.txt" \
    "${CHROMAP}" "${chromap_common_args[@]}" \
    --atac-fragments "${OUTDIR}/memory/fragments.tsv.gz" \
    --summary "${OUTDIR}/memory/summary.tsv" \
    -o "${OUTDIR}/memory/possorted_bam.bam" \
    --call-macs3-frag-peaks \
    --macs3-frag-peaks-source memory \
    --macs3-frag-peaks-output "${OUTDIR}/memory/chromap_macs3_frag.narrowPeak" \
    --macs3-frag-summits-output "${OUTDIR}/memory/chromap_macs3_frag_summits.bed" \
    --macs3-frag-pvalue 1e-5 \
    --macs3-frag-min-length "${MINLEN}" \
    --macs3-frag-max-gap "${MAXGAP}" \
    >>"${log}" 2>&1

  log_msg "[3] fragment multiset identity (sorted lines)"
  if ! cmp -s <(zcat "${OUTDIR}/file/fragments.tsv.gz" | LC_ALL=C sort) \
            <(zcat "${OUTDIR}/memory/fragments.tsv.gz" | LC_ALL=C sort); then
    echo "FAIL: file vs memory fragment lines differ" >&2
    exit 1
  fi

  log_msg "[4] narrowPeak + summits identity"
  if ! cmp -s "${OUTDIR}/file/chromap_macs3_frag.narrowPeak" \
              "${OUTDIR}/memory/chromap_macs3_frag.narrowPeak"; then
    echo "FAIL: narrowPeak differs (file vs memory source)" >&2
    exit 1
  fi
  if ! cmp -s "${OUTDIR}/file/chromap_macs3_frag_summits.bed" \
              "${OUTDIR}/memory/chromap_macs3_frag_summits.bed"; then
    echo "FAIL: summits differ (file vs memory source)" >&2
    exit 1
  fi

  log_msg "[5] reference chromap_callpeaks on file fragments (parity check)"
  "${CALLPEAKS}" -i "${OUTDIR}/file/fragments.tsv.gz" \
    --frag-pileup-macs3-uint8-counts \
    --macs3-frag-narrowpeak "${OUTDIR}/ref_cpp_macs3_frag.narrowPeak" \
    --macs3-frag-summits "${OUTDIR}/ref_cpp_macs3_frag_summits.bed" \
    --bdgpeakcall-cutoff 5 --bdgpeakcall-min-len "${MINLEN}" --bdgpeakcall-max-gap "${MAXGAP}" \
    --frag-score-pseudocount 0 \
    >>"${log}" 2>&1
  if ! cmp -s "${OUTDIR}/file/chromap_macs3_frag.narrowPeak" "${OUTDIR}/ref_cpp_macs3_frag.narrowPeak"; then
    echo "FAIL: file-source narrowPeak differs from chromap_callpeaks" >&2
    exit 1
  fi

  local storage_mode
  storage_mode="$(awk '/MACS3 FRAG memory storage mode:/{print $NF; exit}' "${log}" || true)"
  if [[ -z "${storage_mode}" ]]; then
    storage_mode="n/a"
  fi

  local tfile tmem
  tfile="$(cat "${OUTDIR}/logs/time_file.txt" 2>/dev/null || true)"
  tmem="$(cat "${OUTDIR}/logs/time_memory.txt" 2>/dev/null || true)"

  local summary_tsv="${OUTDIR}/memory_source_100k_summary.tsv"
  {
    echo "harness	chromap_peak_memory_source_100k"
    echo "outdir	${OUTDIR}"
    echo "memory_storage_mode	${storage_mode}"
    echo "file_run_cpu_time	${tfile:-na}"
    echo "memory_run_cpu_time	${tmem:-na}"
    echo "fragments_unchanged_file_vs_memory	true"
    echo "narrowpeak_identical_file_vs_memory	true"
    echo "summits_identical_file_vs_memory	true"
    echo "file_matches_standalone_chromap_callpeaks_narrowpeak	true"
    echo "log	${log}"
  } | tee "${summary_tsv}"

  log_msg "OK. Summary: ${summary_tsv}"
}

main "$@"
