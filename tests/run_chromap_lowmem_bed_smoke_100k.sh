#!/usr/bin/env bash
# Smoke: chromap --low-mem produces BED output identical to non-low-mem mode
# (sorted byte-identical). Regression test for the empty-overflow-stub bug
# that silently dropped BED output under --low-mem on master and earlier
# relink-libmacs3 commits (fixed in 9bec184 by adding WriteToFile /
# LoadFromFile / SerializedSize to the four BED mapping types and
# replacing the empty-body specializations with explicit template
# instantiations).
#
# Three variants:
#   - paired-end + barcode (scATAC)
#   - paired-end + no barcode (bulk)
#   - single-end + barcode (single-end scATAC)
#
# 100K fixture; ~30s wall total.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

CHROMAP="${CHROMAP:-${REPO_ROOT}/chromap}"
FIXTURE_ATAC="${FIXTURE_ATAC:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac}"
INDEX="${INDEX:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index}"
REF="${REF:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa}"
WHITELIST="${WHITELIST:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"
THREADS="${THREADS:-4}"
OUT="${OUT:-$(mktemp -d /tmp/chromap_lowmem_bed_smoke.XXXXXX)}"

R1="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R1_001.fastq.gz"
R2="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R3_001.fastq.gz"
BC="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R2_001.fastq.gz"

if [[ ! -x "${CHROMAP}" ]] ; then
  echo "ERROR: ${CHROMAP} not built" >&2
  exit 2
fi

mkdir -p "${OUT}"
cd "${OUT}"

run_pair() {
  local label=$1
  shift
  echo "[lowmem-bed] ${label}: pair (no --low-mem vs --low-mem)" >&2
  "${CHROMAP}" "$@" -o "${label}.A.bed" 2>"${label}.A.log"
  "${CHROMAP}" "$@" --low-mem -o "${label}.B.bed" 2>"${label}.B.log"
  local a=$(wc -l <"${label}.A.bed")
  local b=$(wc -l <"${label}.B.bed")
  if [[ "${a}" -eq 0 ]] ; then
    echo "FAIL: ${label} baseline produced 0 lines" >&2
    return 1
  fi
  if [[ "${a}" -ne "${b}" ]] ; then
    echo "FAIL: ${label} line count mismatch: A=${a} B=${b}" >&2
    return 1
  fi
  if ! diff -q <(LC_ALL=C sort "${label}.A.bed") <(LC_ALL=C sort "${label}.B.bed") >/dev/null ; then
    echo "FAIL: ${label} sorted BED differs" >&2
    return 1
  fi
  echo "  ${label}: ${a} lines, sorted byte-identical" >&2
}

# Paired-end + barcode (scATAC fragments)
run_pair pe_bc \
  -t "${THREADS}" -x "${INDEX}" -r "${REF}" -1 "${R1}" -2 "${R2}" -b "${BC}" \
  --barcode-whitelist "${WHITELIST}" -l 2000 --trim-adapters --remove-pcr-duplicates \
  --remove-pcr-duplicates-at-cell-level --Tn5-shift

# Paired-end no barcode (bulk ATAC fragments)
run_pair pe_nobc \
  -t "${THREADS}" -x "${INDEX}" -r "${REF}" -1 "${R1}" -2 "${R2}" \
  -l 2000 --trim-adapters --remove-pcr-duplicates --Tn5-shift

# Single-end + barcode
run_pair se_bc \
  -t "${THREADS}" -x "${INDEX}" -r "${REF}" -1 "${R1}" -b "${BC}" \
  --barcode-whitelist "${WHITELIST}" -l 2000 --trim-adapters --remove-pcr-duplicates \
  --remove-pcr-duplicates-at-cell-level

echo "PASS: chromap --low-mem BED output (3 variants × byte-identical sorted vs baseline) — outputs in ${OUT}"
