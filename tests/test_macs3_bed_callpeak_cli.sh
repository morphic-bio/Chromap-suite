#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CALLPEAKS="${CALLPEAKS:-${REPO_ROOT}/chromap_callpeaks}"
MACS3_BIN="${MACS3_BIN:-macs3}"
ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
OUTDIR="${ARTIFACT_ROOT}/macs3_bed_callpeak_cli"
rm -rf "${OUTDIR}"
mkdir -p "${OUTDIR}/cpp" "${OUTDIR}/macs3"

if "${CALLPEAKS}" \
  --macs-bed-input "${OUTDIR}/missing.bed" \
  --macs-input-format BED \
  --nomodel --extsize 200 --shift -100 \
  --effective-genome-size 2700000000 \
  --qvalue 0.05 --pvalue 1e-5 \
  --macs-bed-peaks-out "${OUTDIR}/bad.narrowPeak" \
  --macs-bed-summits-out "${OUTDIR}/bad.summits.bed" \
  >"${OUTDIR}/mutual_exclusion.stdout" \
  2>"${OUTDIR}/mutual_exclusion.stderr"; then
  echo "FAIL: BED q/p mutual exclusion command unexpectedly succeeded" >&2
  exit 1
fi
if ! grep -q -- "--qvalue and --pvalue are mutually exclusive" \
  "${OUTDIR}/mutual_exclusion.stderr"; then
  echo "FAIL: expected BED q/p mutual exclusion error was not reported" >&2
  cat "${OUTDIR}/mutual_exclusion.stderr" >&2
  exit 1
fi

fragments="${OUTDIR}/fragments.tsv"
cells="${OUTDIR}/cells.tsv"
bed="${OUTDIR}/input.bed"
for i in $(seq 0 199); do
  printf 'chr1\t%d\t%d\tBC1\t1\n' "$((1000 + i))" "$((1050 + i))"
done >"${fragments}"
for i in $(seq 0 19); do
  printf 'chr1\t%d\t%d\tBC2\t1\n' "$((5000 + i * 20))" "$((5050 + i * 20))"
done >>"${fragments}"
printf 'BC1\nBC2\n' >"${cells}"
awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,".",".","+"}' "${fragments}" >"${bed}"

"${CALLPEAKS}" \
  --input-fragments "${fragments}" \
  --cells "${cells}" \
  --macs-profile signac-atac \
  --macs-bed-out "${OUTDIR}/projected.bed" \
  --name fixture \
  --outdir "${OUTDIR}/cpp" \
  --macs-bed-peaks-out "${OUTDIR}/cpp/fixture_profile_peaks.narrowPeak" \
  --macs-bed-summits-out "${OUTDIR}/cpp/fixture_profile_summits.bed" \
  --macs-bed-summary "${OUTDIR}/cpp/fixture_profile_summary.tsv"

if [[ "$(wc -l < "${OUTDIR}/projected.bed" | tr -d ' ')" != "220" ]]; then
  echo "FAIL: projected BED row count mismatch" >&2
  exit 1
fi
for expected in \
  $'macs_profile\tsignac-atac' \
  $'macs_input_format\tBED' \
  $'nomodel\ttrue' \
  $'extsize\t200' \
  $'shift\t-100' \
  $'effective_genome_size\t2700000000' \
  $'threshold_mode\tqvalue' \
  $'qvalue\t0.05' \
  $'llocal\t10000' \
  $'projection_input_fragment_rows\t220' \
  $'projection_written_bed_rows\t220'; do
  if ! grep -Fqx "${expected}" "${OUTDIR}/cpp/fixture_profile_summary.tsv"; then
    echo "FAIL: missing summary row: ${expected}" >&2
    cat "${OUTDIR}/cpp/fixture_profile_summary.tsv" >&2
    exit 1
  fi
done

"${CALLPEAKS}" \
  --macs-bed-input "${bed}" \
  --macs-input-format BED \
  --nomodel --extsize 200 --shift -100 \
  --effective-genome-size 2700000000 \
  --qvalue 0.05 \
  --keep-dup all \
  --llocal 10000 \
  --name fixture \
  --macs-bed-peaks-out "${OUTDIR}/cpp/fixture_peaks.narrowPeak" \
  --macs-bed-summits-out "${OUTDIR}/cpp/fixture_summits.bed" \
  --macs-bed-summary "${OUTDIR}/cpp/fixture_summary.tsv"

if command -v "${MACS3_BIN}" >/dev/null 2>&1; then
  "${MACS3_BIN}" callpeak \
    -t "${bed}" \
    -g 2700000000 \
    -f BED \
    --nomodel \
    --extsize 200 \
    --shift -100 \
    -q 0.05 \
    --keep-dup all \
    --llocal 10000 \
    -n fixture \
    --outdir "${OUTDIR}/macs3" \
    >"${OUTDIR}/macs3.stdout" \
    2>"${OUTDIR}/macs3.stderr"

  cmp -s "${OUTDIR}/macs3/fixture_peaks.narrowPeak" \
         "${OUTDIR}/cpp/fixture_peaks.narrowPeak" || {
    echo "FAIL: BED callpeak narrowPeak differs from MACS3" >&2
    diff -u "${OUTDIR}/macs3/fixture_peaks.narrowPeak" \
            "${OUTDIR}/cpp/fixture_peaks.narrowPeak" >&2 || true
    exit 1
  }
  cmp -s "${OUTDIR}/macs3/fixture_summits.bed" \
         "${OUTDIR}/cpp/fixture_summits.bed" || {
    echo "FAIL: BED callpeak summits differ from MACS3" >&2
    diff -u "${OUTDIR}/macs3/fixture_summits.bed" \
            "${OUTDIR}/cpp/fixture_summits.bed" >&2 || true
    exit 1
  }
fi

echo "PASS: MACS3 BED callpeak CLI/profile/projection"
