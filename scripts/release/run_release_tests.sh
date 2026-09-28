#!/usr/bin/env bash
# Hermetic release gate: no GRCh38 or other out-of-tree biological fixtures.
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
artifact_root="$(realpath -m "${1:-$repo_root/plans/artifacts/release-tests}")"
mkdir -p "$artifact_root"
cd "$repo_root"
for tool in python3 samtools bgzip; do command -v "$tool" >/dev/null; done
test -x /usr/bin/time
export CHROMAP_ARTIFACT_ROOT="$artifact_root" BUILD=0 THREADS=4
printf 'target\tstatus\n' > "$artifact_root/tests.tsv"
targets=(
  test-unit test-barcode-sampling test-materialized-reference
  test-frag-compact-store test-macs3-fragment-buckets test-macs3-frag-qvalue-cli
  test-atac-spill-record-roundtrip test-atac-mergeable-spill-materializer
  test-fastq-intake-smoke test-input-format-smoke test-atac-sidecar-only-smoke
  test-libchromap-core-smoke test-cbq-atac-smoke test-cbq-modality-matrix
)
for target in "${targets[@]}"; do
  echo "[release-tests] $target"
  if make -j1 "$target" > "$artifact_root/$target.log" 2>&1 &&
     ! grep -Eq '\] SKIP:' "$artifact_root/$target.log"; then
    printf '%s\tPASS\n' "$target" >> "$artifact_root/tests.tsv"
  else
    printf '%s\tFAIL\n' "$target" >> "$artifact_root/tests.tsv"
    tail -n 80 "$artifact_root/$target.log" >&2
    exit 1
  fi
done
echo "[release-tests] PASS: ${#targets[@]} targets"
