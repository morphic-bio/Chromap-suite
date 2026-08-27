#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
image="${1:-local/chromap-suite:v1.0.1}"
artifact_root="${CHROMAP_ARTIFACT_ROOT:-${repo_root}/plans/artifacts}"
run_id="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
out="${OUT:-${artifact_root}/container_smoke/${run_id}}"

mkdir -p "${out}"

fail() {
  echo "[container-smoke] FAIL: $*" >&2
  exit 1
}

run_container() {
  docker run --rm \
    --user "$(id -u):$(id -g)" \
    --mount "type=bind,src=${out},dst=/work" \
    --workdir /work \
    "${image}" "$@"
}

docker image inspect "${image}" >/dev/null

version="$(docker run --rm "${image}" chromap --version 2>&1)"
[[ "${version}" == "1.0.1" ]] || fail "chromap --version=${version}"

revision="$(docker image inspect \
  --format '{{ index .Config.Labels "org.opencontainers.image.revision" }}' \
  "${image}")"
[[ "${revision}" == "98a4da086f81b7cb159d8fe44efff2fb168e0785" ]] || \
  fail "unexpected source revision label: ${revision}"

docker run --rm "${image}" sh -ceu '
  for binary in chromap rapidmacs chromap_callpeaks chromap_lib_runner chromap_atac_spill_materializer; do
    command -v "${binary}" >/dev/null
    if ldd "$(command -v "${binary}")" | grep -q "not found"; then
      echo "missing shared library for ${binary}" >&2
      exit 1
    fi
  done
'

python3 - "${out}" <<'PY'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
rng = random.Random(1701)
genome = "".join(rng.choice("ACGT") for _ in range(12000))

def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

(out / "ref.fa").write_text(">chrSynthetic\n" + genome + "\n", encoding="ascii")
with (out / "read1.fastq").open("wt", encoding="ascii") as r1, \
     (out / "read2.fastq").open("wt", encoding="ascii") as r2:
    for i in range(16):
        start = 500 + i * 311
        end = start + 220
        seq1 = genome[start:start + 90]
        seq2 = rc(genome[end - 90:end])
        r1.write(f"@pair{i}/1\n{seq1}\n+\n{'I' * len(seq1)}\n")
        r2.write(f"@pair{i}/2\n{seq2}\n+\n{'I' * len(seq2)}\n")
PY

run_container chromap --build-index \
  -r /work/ref.fa -o /work/ref.index -k 11 -w 5 \
  >"${out}/index.stdout" 2>"${out}/index.stderr"

common=(
  -x /work/ref.index
  -r /work/ref.fa
  -1 /work/read1.fastq
  -2 /work/read2.fastq
  -t 1
  --BED
)

run_container chromap "${common[@]}" -o /work/cli.bed \
  >"${out}/cli.stdout" 2>"${out}/cli.stderr"
run_container chromap_lib_runner "${common[@]}" -o /work/lib.bed \
  >"${out}/lib.stdout" 2>"${out}/lib.stderr"

[[ -s "${out}/cli.bed" ]] || fail "Chromap produced no alignments"
[[ -s "${out}/lib.bed" ]] || fail "libchromap runner produced no alignments"
LC_ALL=C sort "${out}/cli.bed" >"${out}/cli.sorted.bed"
LC_ALL=C sort "${out}/lib.bed" >"${out}/lib.sorted.bed"
cmp "${out}/cli.sorted.bed" "${out}/lib.sorted.bed" || \
  fail "container CLI and libchromap outputs differ"

{
  printf 'image\t%s\n' "${image}"
  printf 'source_revision\t%s\n' "${revision}"
  printf 'chromap_suite_version\t%s\n' "${version}"
  printf 'alignment_rows\t%s\n' "$(wc -l < "${out}/cli.bed")"
  printf 'result\tPASS\n'
} >"${out}/summary.tsv"

echo "[container-smoke] PASS: ${out}"
cat "${out}/summary.tsv"
