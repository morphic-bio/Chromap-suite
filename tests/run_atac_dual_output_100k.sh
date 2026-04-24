#!/usr/bin/env bash
# Regression: ATAC dual output (fragments + BAM) on 100K PBMC fixture.
#
# Contract (intended):
#   - Fragment lines from dual mode match fragment-only BED (same mapping options).
#   - Dual BAM record count == 2 * (fragment line count). Not comparable to BAM-only.
#
# Prerequisites: samtools, python3, zcat; Chromap built at repo root (./chromap).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CHROMAP="${CHROMAP:-${REPO_ROOT}/chromap}"

FIXTURE_ATAC="${FIXTURE_ATAC:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac}"
INDEX="${INDEX:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index}"
REF="${REF:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa}"
WHITELIST="${WHITELIST:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"
BENCH_ROOT="${BENCH_ROOT:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k}"

TS="$(date +%Y%m%d_%H%M%S)"
OUT="${OUT:-${BENCH_ROOT}/chromap_atac_dual_regress_${TS}}"
THREADS="${THREADS:-8}"

R1="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R1_001.fastq.gz"
R2="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R3_001.fastq.gz"
BC="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R2_001.fastq.gz"

mkdir -p "${OUT}/logs" "${OUT}/fragment_only" "${OUT}/dual"

if [[ ! -f "${CHROMAP}" ]]; then
  echo "ERROR: Chromap binary not found: ${CHROMAP}" >&2
  exit 2
fi

summarize() {
  echo "==== $1 ====" >> "${OUT}/REGRESSION_SUMMARY.txt"
  cat >> "${OUT}/REGRESSION_SUMMARY.txt"
  echo >> "${OUT}/REGRESSION_SUMMARY.txt"
}

{
  echo "artifact_dir ${OUT}"
  echo "chromap ${CHROMAP}"
  echo "threads ${THREADS}"
  echo "fixture_atac ${FIXTURE_ATAC}"
} | summarize "run_config"

echo "Running fragment-only baseline..."
"${CHROMAP}" -t "${THREADS}" -x "${INDEX}" -r "${REF}" -1 "${R1}" -2 "${R2}" -b "${BC}" \
  --barcode-whitelist "${WHITELIST}" -l 2000 --trim-adapters --remove-pcr-duplicates \
  --remove-pcr-duplicates-at-cell-level --Tn5-shift --BED \
  --summary "${OUT}/logs/summary_fragment.tsv" \
  -o "${OUT}/fragment_only/fragments.tsv" \
  2>&1 | tee "${OUT}/logs/fragment_only.log"
gzip -f "${OUT}/fragment_only/fragments.tsv"

echo "Running dual fragments + BAM..."
"${CHROMAP}" -t "${THREADS}" -x "${INDEX}" -r "${REF}" -1 "${R1}" -2 "${R2}" -b "${BC}" \
  --barcode-whitelist "${WHITELIST}" -l 2000 --trim-adapters --remove-pcr-duplicates \
  --remove-pcr-duplicates-at-cell-level --Tn5-shift --BAM --sort-bam \
  --atac-fragments "${OUT}/dual/fragments.tsv.gz" \
  --summary "${OUT}/logs/summary_dual.tsv" \
  -o "${OUT}/dual/possorted_bam.bam" \
  2>&1 | tee "${OUT}/logs/dual.log"

echo "Optional: BAM-only (informational counts only)..."
mkdir -p "${OUT}/bam_only"
set +e
"${CHROMAP}" -t "${THREADS}" -x "${INDEX}" -r "${REF}" -1 "${R1}" -2 "${R2}" -b "${BC}" \
  --barcode-whitelist "${WHITELIST}" -l 2000 --trim-adapters --remove-pcr-duplicates \
  --remove-pcr-duplicates-at-cell-level --Tn5-shift --BAM --sort-bam \
  --summary "${OUT}/logs/summary_bam_only.tsv" \
  -o "${OUT}/bam_only/possorted_bam.bam" \
  2>&1 | tee "${OUT}/logs/bam_only.log"
set -e
if [[ -f "${OUT}/bam_only/possorted_bam.bam" ]]; then
  BAM_ONLY_REC="$(samtools view -c "${OUT}/bam_only/possorted_bam.bam")"
  echo "BAM-only records (informational, not compared to dual): ${BAM_ONLY_REC}" | tee -a "${OUT}/REGRESSION_SUMMARY.txt"
fi

if ! cmp -s <(zcat "${OUT}/fragment_only/fragments.tsv.gz" | sort) \
            <(zcat "${OUT}/dual/fragments.tsv.gz" | sort); then
  echo "ERROR: sorted fragment lines differ between fragment-only and dual" >&2
  exit 1
fi
echo "fragment tuple multiset: EXACT match (sorted lines)" | summarize "fragments"

samtools quickcheck "${OUT}/dual/possorted_bam.bam"
echo "samtools quickcheck dual BAM: OK" | summarize "dual_bam"

FRAG_LINES="$(zcat "${OUT}/dual/fragments.tsv.gz" | wc -l)"
BAM_REC="$(samtools view -c "${OUT}/dual/possorted_bam.bam")"
EXP=$((FRAG_LINES * 2))
echo "dual fragment lines: ${FRAG_LINES}" | summarize "counts"
echo "dual BAM records: ${BAM_REC} (expected 2 * fragments = ${EXP})" | summarize "counts"

if [[ "${BAM_REC}" -ne "${EXP}" ]]; then
  echo "ERROR: dual BAM records ${BAM_REC} != 2 * fragment lines (${EXP})" >&2
  exit 1
fi

python3 "${SCRIPT_DIR}/check_dual_bam_mate_pairs.py" "${OUT}/dual/possorted_bam.bam" | summarize "mate_pair_check"

echo "PASS" | tee -a "${OUT}/REGRESSION_SUMMARY.txt"
echo "All checks passed. Artifacts: ${OUT}"
