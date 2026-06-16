#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CHROMAP_BIN="${CHROMAP_BIN:-${REPO_ROOT}/chromap}"
ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
OUTDIR="${ARTIFACT_ROOT}/macs3_frag_qvalue_cli"
mkdir -p "${OUTDIR}"

if "${CHROMAP_BIN}" \
  --preset atac \
  --index "${OUTDIR}/missing.index" \
  --ref "${OUTDIR}/missing.fa" \
  --read1 "${OUTDIR}/missing_R1.fastq.gz" \
  --read2 "${OUTDIR}/missing_R2.fastq.gz" \
  --output "${OUTDIR}/out.bed" \
  --macs3-frag-peaks-source memory \
  --call-macs3-frag-peaks \
  --macs3-frag-pvalue 1e-5 \
  --macs3-frag-qvalue 0.05 \
  >"${OUTDIR}/mutual_exclusion.stdout" \
  2>"${OUTDIR}/mutual_exclusion.stderr"; then
  echo "FAIL: p/q mutual exclusion command unexpectedly succeeded" >&2
  exit 1
fi

if ! grep -q -- "--macs3-frag-pvalue and --macs3-frag-qvalue are mutually exclusive" \
  "${OUTDIR}/mutual_exclusion.stderr"; then
  echo "FAIL: expected p/q mutual exclusion error was not reported" >&2
  cat "${OUTDIR}/mutual_exclusion.stderr" >&2
  exit 1
fi

echo "PASS: --macs3-frag-pvalue/--macs3-frag-qvalue mutual exclusion"
