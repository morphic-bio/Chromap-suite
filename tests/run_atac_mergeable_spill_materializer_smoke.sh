#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
artifact_root=${CHROMAP_ARTIFACT_ROOT:-${repo_root}/plans/artifacts/atac_mergeable_spill_materializer}
run_root=${artifact_root}/unit_v3

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
