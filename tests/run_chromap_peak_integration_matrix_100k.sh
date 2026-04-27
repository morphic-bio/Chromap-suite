#!/usr/bin/env bash
# Integration matrix for chromap --call-macs3-frag-peaks (100K fixture).
# Validates that all three reachable peak-pipeline configurations produce
# byte-identical narrowPeak, and that all match the standalone
# chromap_callpeaks reference. The matrix axes are:
#
#   macs3-frag-peaks-source : memory | file
#   --macs3-frag-low-mem    : on (sweep workspace) | off (events workspace)
#
# Effective cells (chromap --low-mem is mutually exclusive with
# --atac-fragments, so it's not a usable axis):
#
#   #1: kMemory + events  (default)
#   #2: kMemory + sweep   (--macs3-frag-low-mem)
#   #3: kFile   + events  (--macs3-frag-peaks-source file)
#
# BED-only peak-calling (no BAM, no atac-fragments) is currently rejected
# by chromap_driver.cc validation; that path is a deferred follow-up.
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
OUTROOT="${OUTROOT:-$(mktemp -d /tmp/chromap_peak_matrix_100k.XXXXXX)}"
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

mkdir -p "${OUTROOT}"

run_cell() {
  local label=$1
  shift
  local cell_dir="${OUTROOT}/${label}"
  mkdir -p "${cell_dir}"
  echo "[matrix] ${label}: chromap $*" >&2
  "${CHROMAP}" "${chromap_common_args[@]}" \
    --atac-fragments "${cell_dir}/fragments.tsv.gz" \
    --summary "${cell_dir}/summary.tsv" \
    -o "${cell_dir}/possorted_bam.bam" \
    --call-macs3-frag-peaks \
    --macs3-frag-peaks-output "${cell_dir}/chromap_macs3_frag.narrowPeak" \
    --macs3-frag-summits-output "${cell_dir}/chromap_macs3_frag_summits.bed" \
    --macs3-frag-pvalue 1e-5 \
    --macs3-frag-min-length "${MINLEN}" \
    --macs3-frag-max-gap "${MAXGAP}" \
    "$@" \
    >>"${cell_dir}/run.log" 2>&1
}

# Standalone reference (chromap_callpeaks on the same fragments).
echo "[matrix] standalone reference (chromap_callpeaks on baseline fragments)" >&2
mkdir -p "${OUTROOT}/standalone_ref"
run_cell baseline_no_peaks
"${CALLPEAKS}" -i "${OUTROOT}/baseline_no_peaks/fragments.tsv.gz" \
  --frag-pileup-macs3-uint8-counts \
  --macs3-frag-narrowpeak "${OUTROOT}/standalone_ref/cpp_macs3_frag.narrowPeak" \
  --macs3-frag-summits "${OUTROOT}/standalone_ref/cpp_macs3_frag_summits.bed" \
  --bdgpeakcall-cutoff 5 --bdgpeakcall-min-len "${MINLEN}" --bdgpeakcall-max-gap "${MAXGAP}" \
  --frag-score-pseudocount 0 \
  >>"${OUTROOT}/standalone_ref/run.log" 2>&1

REF_NP="${OUTROOT}/standalone_ref/cpp_macs3_frag.narrowPeak"
REF_SM="${OUTROOT}/standalone_ref/cpp_macs3_frag_summits.bed"

run_cell cell1_kMemory_events --macs3-frag-peaks-source memory
run_cell cell2_kMemory_sweep  --macs3-frag-peaks-source memory --macs3-frag-low-mem
run_cell cell3_kFile_events   --macs3-frag-peaks-source file

# Validate each cell produces byte-identical narrowPeak vs the standalone reference,
# and that all cells are byte-identical to each other (transitive via #1).
echo "[matrix] validating byte-identity..." >&2
fail=0
for cell in cell1_kMemory_events cell2_kMemory_sweep cell3_kFile_events ; do
  np="${OUTROOT}/${cell}/chromap_macs3_frag.narrowPeak"
  sm="${OUTROOT}/${cell}/chromap_macs3_frag_summits.bed"
  if ! cmp -s "${np}" "${REF_NP}" ; then
    echo "FAIL: ${cell} narrowPeak diverges from standalone ref" >&2
    fail=1
  fi
  if ! cmp -s "${sm}" "${REF_SM}" ; then
    echo "FAIL: ${cell} summits diverge from standalone ref" >&2
    fail=1
  fi
done

if [[ "${fail}" -ne 0 ]] ; then
  echo "MATRIX FAIL — outputs in ${OUTROOT}" >&2
  exit 1
fi

NP_MD5=$(md5sum "${REF_NP}" | cut -d' ' -f1)
SM_MD5=$(md5sum "${REF_SM}" | cut -d' ' -f1)
N_PEAKS=$(wc -l < "${REF_NP}")
echo "MATRIX PASS: 3 cells × byte-identical narrowPeak (md5=${NP_MD5}, peaks=${N_PEAKS}, summits md5=${SM_MD5})"
echo "outputs: ${OUTROOT}"
