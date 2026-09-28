#!/usr/bin/env bash
# Sidecar-only ATAC output (--atac-sidecar-only) on a synthetic scATAC fixture.
#
# For each case the dual mode (--BAM --atac-fragments
# --atac-fragment-binary-output) and the sidecar-only mode run on the same
# input and options, through both the chromap CLI and chromap_lib_runner.
# Required, byte for byte:
#   - the AEV1 sidecar and <sidecar>.chroms.tsv,
#   - the --summary table.
# Also required: the sidecar-only run leaves no BAM, fragments text, primary
# output or temporary sidecar file; its records match the dual fragment rows
# (chrom, start, end, count); the invalid combinations are rejected.
#
# The fixture is generated here (about 60 kb of reference, a few thousand
# read pairs); no production genome or index is loaded.
#
# Environment: CHROMAP_BIN, LIBRUNNER_BIN, THREADS (default 2), BUILD (1 builds
# chromap and chromap_lib_runner with MAKE_JOBS, default 2), OUTROOT.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)-$$}"
OUTROOT="${OUTROOT:-${ARTIFACT_ROOT}/atac_sidecar_only_smoke/${RUN_ID}}"
CHROMAP_BIN="${CHROMAP_BIN:-${REPO_ROOT}/chromap}"
LIBRUNNER_BIN="${LIBRUNNER_BIN:-${REPO_ROOT}/chromap_lib_runner}"
THREADS="${THREADS:-2}"
BUILD="${BUILD:-1}"

mkdir -p "${OUTROOT}"/{fixture,index,cases,logs}
SUMMARY="${OUTROOT}/summary.tsv"
printf 'case\tstatus\tdetail\n' >"${SUMMARY}"

log() { printf '[sidecar-only-smoke] %s\n' "$*" >&2; }
fail() { printf '[sidecar-only-smoke] FAIL: %s\n' "$*" >&2; exit 1; }
record() { printf '%s\t%s\t%s\n' "$1" "$2" "$3" >>"${SUMMARY}"; }

run_cmd() {
  local id="$1"; shift
  log "${id}: $*"
  "$@" >"${OUTROOT}/logs/${id}.stdout" 2>"${OUTROOT}/logs/${id}.stderr" ||
    fail "${id} exited $? (see ${OUTROOT}/logs/${id}.stderr)"
}

expect_reject() {
  local id="$1"; shift
  local want="$1"; shift
  log "${id} (expect rejection): $*"
  if "$@" >"${OUTROOT}/logs/${id}.stdout" 2>"${OUTROOT}/logs/${id}.stderr"; then
    fail "${id}: invalid combination was accepted"
  fi
  grep -q -- "${want}" "${OUTROOT}/logs/${id}.stderr" ||
    fail "${id}: rejection message lacks '${want}'"
  record "${id}" "PASS" "rejected: ${want}"
}

generate_fixture() {
  python3 - "$1" <<'PY'
import os
import random
import sys

out = sys.argv[1]
rng = random.Random(20260928)
read_len = 50
qual = "I" * read_len


def randseq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def rc(s):
    return s.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


chroms = {"chr1": randseq(30000), "chr2": randseq(20000), "chrY": randseq(10000)}
# A repeated block makes some fragments multi-map (low MAPQ, filtered).
block = chroms["chr1"][5000:5600]
chroms["chr2"] = chroms["chr2"][:8000] + block + chroms["chr2"][8600:]

with open(os.path.join(out, "ref.fa"), "w") as f:
    for name, seq in chroms.items():
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i + 80] + "\n")

whitelist = sorted({randseq(16) for _ in range(40)})[:32]
with open(os.path.join(out, "whitelist.txt"), "w") as f:
    f.write("\n".join(whitelist) + "\n")

sites = []
for name, seq in chroms.items():
    for _ in range(len(seq) // 150):
        length = rng.randint(120, 700)
        start = rng.randint(0, len(seq) - length - 1)
        sites.append((name, start, length))
weights = [1.0 / (1 + i) ** 0.7 for i in range(len(sites))]


def mutate(s):
    i = rng.randrange(len(s))
    return s[:i] + rng.choice([b for b in "ACGT" if b != s[i]]) + s[i + 1:]


with open(os.path.join(out, "read1.fq"), "w") as r1, \
     open(os.path.join(out, "read2.fq"), "w") as r2, \
     open(os.path.join(out, "barcode.fq"), "w") as bc:
    for i in range(6000):
        name, start, length = rng.choices(sites, weights)[0]
        seq = chroms[name]
        left = seq[start:start + read_len]
        right = rc(seq[start + length - read_len:start + length])
        a, b = (left, right) if rng.random() < 0.5 else (right, left)
        if rng.random() < 0.1:
            a = mutate(a)
        barcode = rng.choice(whitelist[:8] if rng.random() < 0.6 else whitelist)
        roll = rng.random()
        if roll < 0.05:
            barcode = mutate(barcode)       # correctable
        elif roll < 0.07:
            barcode = randseq(16)           # usually not in the whitelist
        r1.write(f"@pair{i}/1\n{a}\n+\n{qual}\n")
        r2.write(f"@pair{i}/2\n{b}\n+\n{qual}\n")
        bc.write(f"@pair{i}/3\n{barcode}\n+\n{'I' * 16}\n")
PY
}

sidecar_rows() {  # AEV1 sidecar -> "chrom\tstart\tend\tcount" rows
  python3 - "$1" <<'PY'
import struct
import sys

path = sys.argv[1]
names = {}
with open(path + ".chroms.tsv") as f:
    for line in f:
        idx, name = line.rstrip("\n").split("\t")
        names[int(idx)] = name
with open(path, "rb") as f:
    data = f.read()
magic, version, rsize, bclen, nchrom, flags, nrec = struct.unpack_from("<4sIIIIIQ", data, 0)
if magic != b"AEV1" or version != 1 or rsize != 24:
    sys.exit(f"bad AEV1 header in {path}")
if len(data) != 32 + nrec * 24:
    sys.exit(f"{path}: {len(data)} bytes for {nrec} records")
if nrec == 0:
    sys.exit(f"{path}: no records")
for i in range(nrec):
    chrom, start, end, count, _key = struct.unpack_from("<iiiIQ", data, 32 + 24 * i)
    print(f"{names[chrom]}\t{start}\t{end}\t{count}")
PY
}

run_case() {  # run_case ID TOOL EXTRA...
  local id="$1"; shift
  local tool="$1"; shift
  local dir="${OUTROOT}/cases/${id}"
  mkdir -p "${dir}/dual" "${dir}/sidecar" "${dir}/tmp_dual" "${dir}/tmp_sidecar"
  run_cmd "${id}.dual" "${tool}" "${COMMON[@]}" "$@" \
    --temp-dir "${dir}/tmp_dual" --summary "${dir}/dual/summary.tsv" \
    --BAM -o "${dir}/dual/atac.bam" --atac-fragments "${dir}/dual/fragments.tsv" \
    --atac-fragment-binary-output "${dir}/dual/atac_fragments.bin"
  run_cmd "${id}.sidecar" "${tool}" "${COMMON[@]}" "$@" \
    --temp-dir "${dir}/tmp_sidecar" --summary "${dir}/sidecar/summary.tsv" \
    --atac-sidecar-only \
    --atac-fragment-binary-output "${dir}/sidecar/atac_fragments.bin"

  for f in atac_fragments.bin atac_fragments.bin.chroms.tsv summary.tsv; do
    [[ -s "${dir}/dual/${f}" ]] || fail "${id}: dual ${f} missing"
    cmp "${dir}/dual/${f}" "${dir}/sidecar/${f}" ||
      fail "${id}: ${f} differs between dual and sidecar-only"
  done
  local extra
  extra="$(cd "${dir}/sidecar" && find . -mindepth 1 \
    ! -name atac_fragments.bin ! -name atac_fragments.bin.chroms.tsv \
    ! -name summary.tsv -print)"
  [[ -z "${extra}" ]] || fail "${id}: sidecar-only left extra files: ${extra}"

  sidecar_rows "${dir}/sidecar/atac_fragments.bin" >"${dir}/sidecar.rows"
  cut -f1,2,3,5 "${dir}/dual/fragments.tsv" >"${dir}/dual.rows"
  cmp "${dir}/sidecar.rows" "${dir}/dual.rows" ||
    fail "${id}: sidecar records differ from dual fragment rows"
  local n
  n="$(wc -l <"${dir}/sidecar.rows")"
  record "${id}" "PASS" "sidecar,chroms,summary byte-identical; ${n} records; no BAM/text"
}

main() {
  if [[ "${BUILD}" == "1" ]]; then
    log "building chromap and chromap_lib_runner (-j${MAKE_JOBS:-2})"
    (cd "${REPO_ROOT}" && make -j"${MAKE_JOBS:-2}" chromap chromap_lib_runner) \
      >"${OUTROOT}/logs/make.stdout" 2>"${OUTROOT}/logs/make.stderr"
  fi
  [[ -x "${CHROMAP_BIN}" ]] || fail "missing ${CHROMAP_BIN}"
  [[ -x "${LIBRUNNER_BIN}" ]] || fail "missing ${LIBRUNNER_BIN}"

  local fx="${OUTROOT}/fixture"
  generate_fixture "${fx}"
  local index="${OUTROOT}/index/ref.index"
  run_cmd index "${CHROMAP_BIN}" -i -r "${fx}/ref.fa" -o "${index}"

  COMMON=(-t "${THREADS}" -x "${index}" -r "${fx}/ref.fa"
    -1 "${fx}/read1.fq" -2 "${fx}/read2.fq" -b "${fx}/barcode.fq"
    --barcode-whitelist "${fx}/whitelist.txt" -l 2000 --trim-adapters
    --remove-pcr-duplicates --Tn5-shift)

  for tool_id in cli lib; do
    local tool="${CHROMAP_BIN}"
    [[ "${tool_id}" == "lib" ]] && tool="${LIBRUNNER_BIN}"
    run_case "S01_${tool_id}_cell_inmem" "${tool}" --remove-pcr-duplicates-at-cell-level
    run_case "S02_${tool_id}_cell_lowmem_spill" "${tool}" \
      --remove-pcr-duplicates-at-cell-level --low-mem --low-mem-ram 16K
    run_case "S03_${tool_id}_cell_lowmem_drain" "${tool}" \
      --remove-pcr-duplicates-at-cell-level --low-mem
    run_case "S04_${tool_id}_bulk_lowmem_spill" "${tool}" --low-mem --low-mem-ram 16K
    run_case "S05_${tool_id}_bulk_inmem" "${tool}"
  done
  # A spill-forcing case must actually have spilled.
  grep -q "overflow files for k-way merge" \
    "${OUTROOT}/logs/S02_cli_cell_lowmem_spill.sidecar.stderr" ||
    fail "S02: sidecar-only low-mem run did not spill"

  local s="${OUTROOT}/cases/S01_cli_cell_inmem/sidecar/x.bin"
  expect_reject R01_with_output "do not pass -o" "${CHROMAP_BIN}" "${COMMON[@]}" \
    --atac-sidecar-only --atac-fragment-binary-output "${s}" -o "${OUTROOT}/cases/x.bed"
  expect_reject R02_with_bam "uses the BED fragment path" "${CHROMAP_BIN}" "${COMMON[@]}" \
    --atac-sidecar-only --BAM --atac-fragment-binary-output "${s}"
  expect_reject R03_no_sidecar_path "requires --atac-fragment-binary-output" \
    "${CHROMAP_BIN}" "${COMMON[@]}" --atac-sidecar-only
  expect_reject R04_with_fragments "writes no fragment text rows" "${CHROMAP_BIN}" \
    "${COMMON[@]}" --atac-sidecar-only --atac-fragments "${OUTROOT}/cases/x.tsv" \
    --atac-fragment-binary-output "${s}"
  expect_reject R05_single_end "requires paired-end reads" "${CHROMAP_BIN}" -t 1 \
    -x "${index}" -r "${fx}/ref.fa" -1 "${fx}/read1.fq" --atac-sidecar-only \
    --atac-fragment-binary-output "${s}"
  expect_reject R06_lib_with_bam "uses the BED fragment path" "${LIBRUNNER_BIN}" \
    "${COMMON[@]}" --atac-sidecar-only --BAM --atac-fragment-binary-output "${s}"
  expect_reject R07_file_peaks "need the memory source" "${CHROMAP_BIN}" "${COMMON[@]}" \
    --atac-sidecar-only --atac-fragment-binary-output "${s}" --call-macs3-frag-peaks \
    --macs3-frag-peaks-output "${OUTROOT}/cases/x.np" \
    --macs3-frag-summits-output "${OUTROOT}/cases/x.summits"
  [[ ! -e "${s}" && ! -e "${s}.tmp" ]] || fail "a rejected run wrote ${s}"

  log "PASS: ${SUMMARY}"
  cat "${SUMMARY}" >&2
}

main "$@"
