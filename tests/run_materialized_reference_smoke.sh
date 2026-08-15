#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
CHROMAP="${CHROMAP:-${REPO_ROOT}/chromap}"
ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
OUT_ROOT="${OUT_ROOT:-${ARTIFACT_ROOT}/materialized_reference_smoke/$(date -u +%Y%m%dT%H%M%SZ)}"
THREADS="${THREADS:-4}"

mkdir -p "${OUT_ROOT}/logs"

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  exit 1
}

expect_failure() {
  local label="$1"
  local expected="$2"
  shift 2
  if "$@" >"${OUT_ROOT}/logs/${label}.stdout" \
      2>"${OUT_ROOT}/logs/${label}.stderr"; then
    fail "${label} unexpectedly succeeded"
  fi
  grep -Fq -- "${expected}" "${OUT_ROOT}/logs/${label}.stderr" || {
    sed -n '1,160p' "${OUT_ROOT}/logs/${label}.stderr" >&2
    fail "${label} did not report: ${expected}"
  }
}

[[ -x "${CHROMAP}" ]] || fail "Chromap binary is not executable: ${CHROMAP}"

python3 - "${OUT_ROOT}" <<'PY'
import random
import sys
from pathlib import Path

root = Path(sys.argv[1])
rng = random.Random(1701)

def random_sequence(length):
    return "".join(rng.choice("ACGT") for _ in range(length))

chr1 = random_sequence(50000)
chr2 = random_sequence(35000)
chr3 = "NNNN" + random_sequence(11992) + "acgtnn"

with (root / "reference.fa").open("wt", encoding="ascii") as out:
    for name, sequence in (("chr1", chr1), ("chr2", chr2), ("chrTiny", chr3)):
        out.write(f">{name}\n")
        for start in range(0, len(sequence), 61):
            out.write(sequence[start:start + 61] + "\n")

def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]

records = []
for i in range(120):
    sequence = chr1 if i % 2 == 0 else chr2
    start = 200 + i * 127
    read = sequence[start:start + 80]
    if i % 3 == 0:
        read = reverse_complement(read)
    records.append((f"read{i:04d}", read))

with (root / "reads.fastq").open("wt", encoding="ascii") as out:
    for name, sequence in records:
        out.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")

alternate = list(chr1)
alternate[1234] = "A" if alternate[1234] != "A" else "C"
with (root / "alternate.fa").open("wt", encoding="ascii") as out:
    out.write(">chr1\n")
    sequence = "".join(alternate)
    for start in range(0, len(sequence), 61):
        out.write(sequence[start:start + 61] + "\n")
    out.write(">chr2\n")
    for start in range(0, len(chr2), 61):
        out.write(chr2[start:start + 61] + "\n")
    out.write(">chrTiny\n")
    for start in range(0, len(chr3), 61):
        out.write(chr3[start:start + 61] + "\n")
PY

"${CHROMAP}" --build-index \
  --ref "${OUT_ROOT}/reference.fa" \
  --output "${OUT_ROOT}/bound.index" \
  --reference-sidecar "${OUT_ROOT}/reference.chromapref" \
  --kmer 11 --window 5 --num-threads "${THREADS}" \
  >"${OUT_ROOT}/logs/build_bound.stdout" \
  2>"${OUT_ROOT}/logs/build_bound.stderr"

[[ -s "${OUT_ROOT}/bound.index" ]] || fail "bound index was not created"
[[ -s "${OUT_ROOT}/reference.chromapref" ]] || \
  fail "materialized reference was not created"
grep -Fq "Saved materialized reference successfully" \
  "${OUT_ROOT}/logs/build_bound.stderr" || \
  fail "index generation did not report sidecar creation"

COMMON_ARGS=(
  --index "${OUT_ROOT}/bound.index"
  --read1 "${OUT_ROOT}/reads.fastq"
  --BED
  --min-num-seeds 1
  --error-threshold 2
  --MAPQ-threshold 0
  --num-threads "${THREADS}"
)

/usr/bin/time -f 'wall_seconds=%e\nuser_seconds=%U\nsys_seconds=%S\nmax_rss_kb=%M' \
  -o "${OUT_ROOT}/fasta.time.txt" \
  "${CHROMAP}" "${COMMON_ARGS[@]}" \
  --ref "${OUT_ROOT}/reference.fa" \
  --output "${OUT_ROOT}/fasta.bed" \
  >"${OUT_ROOT}/logs/fasta.stdout" 2>"${OUT_ROOT}/logs/fasta.stderr"

/usr/bin/time -f 'wall_seconds=%e\nuser_seconds=%U\nsys_seconds=%S\nmax_rss_kb=%M' \
  -o "${OUT_ROOT}/sidecar.time.txt" \
  "${CHROMAP}" "${COMMON_ARGS[@]}" \
  --reference-sidecar "${OUT_ROOT}/reference.chromapref" \
  --output "${OUT_ROOT}/sidecar.bed" \
  >"${OUT_ROOT}/logs/sidecar.stdout" 2>"${OUT_ROOT}/logs/sidecar.stderr"

cmp "${OUT_ROOT}/fasta.bed" "${OUT_ROOT}/sidecar.bed" || \
  fail "FASTA and materialized-reference mappings differ"
grep -Fq "Loaded all sequences successfully" \
  "${OUT_ROOT}/logs/fasta.stderr" || fail "FASTA path was not exercised"
grep -Fq "Loaded materialized reference successfully" \
  "${OUT_ROOT}/logs/sidecar.stderr" || fail "sidecar path was not exercised"
grep -Fq \
  "Materialized reference input path: parallel direct I/O positioned reader" \
  "${OUT_ROOT}/logs/sidecar.stderr" || \
  fail "parallel direct-I/O sidecar loader was not exercised"
if grep -Fq "Loaded all sequences successfully" \
    "${OUT_ROOT}/logs/sidecar.stderr"; then
  fail "sidecar mapping unexpectedly parsed the FASTA"
fi

"${CHROMAP}" --build-index \
  --ref "${OUT_ROOT}/reference.fa" \
  --output "${OUT_ROOT}/legacy.index" \
  --kmer 11 --window 5 \
  >"${OUT_ROOT}/logs/build_legacy.stdout" \
  2>"${OUT_ROOT}/logs/build_legacy.stderr"

"${CHROMAP}" --index "${OUT_ROOT}/legacy.index" \
  --ref "${OUT_ROOT}/reference.fa" \
  --read1 "${OUT_ROOT}/reads.fastq" --BED \
  --output "${OUT_ROOT}/legacy-fasta.bed" \
  --min-num-seeds 1 --error-threshold 2 --MAPQ-threshold 0 \
  --num-threads "${THREADS}" \
  >"${OUT_ROOT}/logs/legacy-fasta.stdout" \
  2>"${OUT_ROOT}/logs/legacy-fasta.stderr"
cmp "${OUT_ROOT}/fasta.bed" "${OUT_ROOT}/legacy-fasta.bed" || \
  fail "legacy FASTA/index fallback differs from the bound index"

expect_failure "legacy_index_rejects_sidecar" \
  "does not declare a materialized reference" \
  "${CHROMAP}" --index "${OUT_ROOT}/legacy.index" \
  --reference-sidecar "${OUT_ROOT}/reference.chromapref" \
  --read1 "${OUT_ROOT}/reads.fastq" --BED \
  --output "${OUT_ROOT}/legacy-sidecar.bed"

"${CHROMAP}" --build-index \
  --ref "${OUT_ROOT}/alternate.fa" \
  --output "${OUT_ROOT}/alternate.index" \
  --reference-sidecar "${OUT_ROOT}/alternate.chromapref" \
  --kmer 11 --window 5 --num-threads "${THREADS}" \
  >"${OUT_ROOT}/logs/build_alternate.stdout" \
  2>"${OUT_ROOT}/logs/build_alternate.stderr"

expect_failure "mismatched_sidecar" \
  "does not match the Chromap index" \
  "${CHROMAP}" --index "${OUT_ROOT}/bound.index" \
  --reference-sidecar "${OUT_ROOT}/alternate.chromapref" \
  --read1 "${OUT_ROOT}/reads.fastq" --BED \
  --output "${OUT_ROOT}/mismatched.bed" --num-threads "${THREADS}"

python3 - "${OUT_ROOT}/reference.chromapref" \
  "${OUT_ROOT}/truncated.chromapref" <<'PY'
import sys
from pathlib import Path

source = Path(sys.argv[1]).read_bytes()
Path(sys.argv[2]).write_bytes(source[:-1])
PY

expect_failure "truncated_sidecar" \
  "invalid materialized reference dimensions" \
  "${CHROMAP}" --index "${OUT_ROOT}/bound.index" \
  --reference-sidecar "${OUT_ROOT}/truncated.chromapref" \
  --read1 "${OUT_ROOT}/reads.fastq" --BED \
  --output "${OUT_ROOT}/truncated.bed" --num-threads "${THREADS}"

{
  printf 'status=passed\n'
  printf 'threads=%s\n' "${THREADS}"
  printf 'reference_bases=97002\n'
  printf 'mapping_records=%s\n' "$(wc -l < "${OUT_ROOT}/sidecar.bed")"
  printf 'fasta_output_sha256=%s\n' \
    "$(sha256sum "${OUT_ROOT}/fasta.bed" | awk '{print $1}')"
  printf 'sidecar_output_sha256=%s\n' \
    "$(sha256sum "${OUT_ROOT}/sidecar.bed" | awk '{print $1}')"
} >"${OUT_ROOT}/RESULTS.txt"

printf 'PASS: materialized reference smoke\n'
printf 'Artifacts: %s\n' "${OUT_ROOT}"
