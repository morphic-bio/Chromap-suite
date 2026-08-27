#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
artifact_root=${CHROMAP_ARTIFACT_ROOT:-${repo_root}/plans/artifacts/atac_mergeable_spill_materializer}
run_root=${CHROMAP_ATAC_SPILL_RUN_ROOT:-${artifact_root}/$(date -u +%Y%m%dT%H%M%SZ)-$$}

[[ ! -e ${run_root} ]] || {
  echo "FATAL: CHROMAP_ATAC_SPILL_RUN_ROOT must be fresh: ${run_root}" >&2
  exit 2
}
mkdir -p "${run_root}"
"${repo_root}/tests/test_atac_mergeable_spill_materializer" "${run_root}"
samtools quickcheck "${run_root}/materialized.bam" \
  "${run_root}/materialized.noY.bam" "${run_root}/materialized.Y.bam" \
  "${run_root}/materialized.cram"
test "$(samtools view -c "${run_root}/materialized.bam")" -eq 6
test "$(samtools view -c "${run_root}/materialized.noY.bam")" -eq 4
test "$(samtools view -c "${run_root}/materialized.Y.bam")" -eq 2
test "$(samtools view -c -T "${run_root}/reference.fa" \
  "${run_root}/materialized.cram")" -eq 6
test "$(dd if="${run_root}/materialized.aev1" bs=1 count=4 status=none)" = "AEV1"
