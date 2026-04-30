#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
OUTROOT="${OUTROOT:-${ARTIFACT_ROOT}/chromap_core_smoke/${RUN_ID}}"
CHROMAP_BIN="${CHROMAP_BIN:-${REPO_ROOT}/chromap}"
LIBRUNNER_BIN="${LIBRUNNER_BIN:-${REPO_ROOT}/chromap_lib_runner}"
THREADS="${THREADS:-1}"
BUILD="${BUILD:-1}"

mkdir -p "${OUTROOT}"/{fixture,index,cases,logs}
SUMMARY="${OUTROOT}/summary.tsv"

log() {
  printf '[libchromap-smoke] %s\n' "$*" >&2
}

fail() {
  printf '[libchromap-smoke] FAIL: %s\n' "$*" >&2
  exit 1
}

record() {
  local id="$1"; shift
  local status="$1"; shift
  local detail="$*"
  printf '%s\t%s\t%s\n' "${id}" "${status}" "${detail}" >>"${SUMMARY}"
}

require_file_nonempty() {
  local path="$1"
  [[ -s "${path}" ]] || fail "missing or empty file: ${path}"
}

normalize_text_rows() {
  local in="$1"
  LC_ALL=C grep -v '^#' "${in}" | sed '/^[[:space:]]*$/d' | LC_ALL=C sort || true
}

assert_same_rows() {
  local id="$1"
  local a="$2"
  local b="$3"
  local work="${OUTROOT}/cases/${id}"
  mkdir -p "${work}"
  normalize_text_rows "${a}" >"${work}/a.rows"
  normalize_text_rows "${b}" >"${work}/b.rows"
  cmp -s "${work}/a.rows" "${work}/b.rows" || {
    diff -u "${work}/a.rows" "${work}/b.rows" | sed -n '1,80p' >&2 || true
    fail "${id}: normalized rows differ"
  }
}

assert_bam_semantic_equal() {
  local id="$1"
  local a="$2"
  local b="$3"
  local work="${OUTROOT}/cases/${id}"
  mkdir -p "${work}"
  if ! command -v samtools >/dev/null 2>&1; then
    log "${id}: samtools missing; checking BAM existence only"
    require_file_nonempty "${a}"
    require_file_nonempty "${b}"
    return
  fi
  samtools view "${a}" | LC_ALL=C sort >"${work}/a.samrows"
  samtools view "${b}" | LC_ALL=C sort >"${work}/b.samrows"
  cmp -s "${work}/a.samrows" "${work}/b.samrows" || {
    diff -u "${work}/a.samrows" "${work}/b.samrows" | sed -n '1,80p' >&2 || true
    fail "${id}: BAM alignment rows differ"
  }
}

assert_index_exists() {
  local bam="$1"
  [[ -s "${bam}.bai" || -s "${bam}.csi" ]] || fail "missing BAM index for ${bam}"
}

run_cmd() {
  local id="$1"; shift
  log "${id}: $*"
  "$@" >"${OUTROOT}/logs/${id}.stdout" 2>"${OUTROOT}/logs/${id}.stderr"
}

run_cli_and_lib_text_case() {
  local id="$1"; shift
  local extra=("$@")
  local case_dir="${OUTROOT}/cases/${id}"
  mkdir -p "${case_dir}"
  local cli_out="${case_dir}/cli.out"
  local lib_out="${case_dir}/lib.out"

  run_cmd "${id}.cli" "${CHROMAP_BIN}" "${COMMON[@]}" "${extra[@]}" -o "${cli_out}"
  run_cmd "${id}.lib" "${LIBRUNNER_BIN}" "${COMMON[@]}" "${extra[@]}" -o "${lib_out}"
  require_file_nonempty "${cli_out}"
  require_file_nonempty "${lib_out}"
  assert_same_rows "${id}" "${cli_out}" "${lib_out}"
  record "${id}" "PASS" "cli_vs_lib_text_rows_equal"
}

generate_fixture() {
  local fixture_dir="$1"
  python3 - "${fixture_dir}" <<'PY'
import os
import random
import sys

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
rng = random.Random(20260430)
read_len = 76
qual = "I" * read_len

def randseq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))

def rc(s):
    table = str.maketrans("ACGTN", "TGCAN")
    return s.translate(table)[::-1]

chroms = {
    "chr1": randseq(3000),
    "chr2": randseq(3000),
    "chrY": randseq(3000),
}

with open(os.path.join(out, "ref.fa"), "w") as f:
    for name, seq in chroms.items():
        f.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

records = []
def add_pair(name, chrom1, pos1, chrom2=None, pos2=None):
    chrom2 = chrom2 or chrom1
    pos2 = pos2 if pos2 is not None else pos1 + 180
    records.append((name, chrom1, pos1, chrom2, pos2))

for i in range(12):
    add_pair(f"chr1_{i}", "chr1", 100 + i * 37)
for i in range(8):
    add_pair(f"chr2_{i}", "chr2", 120 + i * 41)
for i in range(6):
    add_pair(f"chry_{i}", "chrY", 140 + i * 43)
for i in range(4):
    add_pair(f"trans_{i}", "chr1", 900 + i * 29, "chr2", 1100 + i * 31)

barcodes = ["ACGTACGTACGTACGT", "TGCATGCATGCATGCA"]
with open(os.path.join(out, "read1.fq"), "w") as r1, \
     open(os.path.join(out, "read2.fq"), "w") as r2, \
     open(os.path.join(out, "barcode.fq"), "w") as bc:
    for i, (name, c1, p1, c2, p2) in enumerate(records):
        s1 = chroms[c1][p1:p1+read_len]
        s2 = rc(chroms[c2][p2:p2+read_len])
        barcode = barcodes[i % len(barcodes)]
        r1.write(f"@{name}/1\n{s1}\n+\n{qual}\n")
        r2.write(f"@{name}/2\n{s2}\n+\n{qual}\n")
        bc.write(f"@{name}/bc\n{barcode}\n+\n{'I' * len(barcode)}\n")

with open(os.path.join(out, "whitelist.txt"), "w") as f:
    for barcode in barcodes:
        f.write(barcode + "\n")

with open(os.path.join(out, "chrom_order.txt"), "w") as f:
    f.write("chr1\nchr2\nchrY\n")
PY
}

main() {
  printf 'case\tstatus\tdetail\n' >"${SUMMARY}"

  if [[ "${BUILD}" == "1" ]]; then
    log "building chromap and chromap_lib_runner"
    (cd "${REPO_ROOT}" && make -j"${MAKE_JOBS:-4}" chromap chromap_lib_runner) \
      >"${OUTROOT}/logs/make.stdout" 2>"${OUTROOT}/logs/make.stderr"
  fi

  [[ -x "${CHROMAP_BIN}" ]] || fail "missing executable ${CHROMAP_BIN}"
  [[ -x "${LIBRUNNER_BIN}" ]] || fail "missing executable ${LIBRUNNER_BIN}"

  generate_fixture "${OUTROOT}/fixture"
  REF="${OUTROOT}/fixture/ref.fa"
  R1="${OUTROOT}/fixture/read1.fq"
  R2="${OUTROOT}/fixture/read2.fq"
  BC="${OUTROOT}/fixture/barcode.fq"
  WHITELIST="${OUTROOT}/fixture/whitelist.txt"
  INDEX="${OUTROOT}/index/ref.index"

  run_cmd "C01_index_build" "${CHROMAP_BIN}" -i -r "${REF}" -o "${INDEX}"
  require_file_nonempty "${INDEX}"
  record "C01_index_build" "PASS" "index=${INDEX}"

  COMMON=(-x "${INDEX}" -r "${REF}" -1 "${R1}" -2 "${R2}" -t "${THREADS}")

  run_cli_and_lib_text_case "C02_basic_bed" --BED
  run_cli_and_lib_text_case "C03_chip_tagalign" --preset chip --TagAlign
  run_cli_and_lib_text_case "C04_atac_bed" --preset atac --Tn5-shift-mode symmetric
  run_cli_and_lib_text_case "C07_lowmem_bed" --BED --low-mem --low-mem-ram 1M
  run_cli_and_lib_text_case "C10_hic_pairs" --split-alignment --pairs --MAPQ-threshold 1

  local scatac_dir="${OUTROOT}/cases/C05_scatac_barcode_summary"
  mkdir -p "${scatac_dir}"
  run_cmd "C05.cli" "${CHROMAP_BIN}" "${COMMON[@]}" --preset atac \
    -b "${BC}" --barcode-whitelist "${WHITELIST}" --skip-barcode-check \
    --summary "${scatac_dir}/cli.summary.tsv" -o "${scatac_dir}/cli.bed"
  run_cmd "C05.lib" "${LIBRUNNER_BIN}" "${COMMON[@]}" --preset atac \
    -b "${BC}" --barcode-whitelist "${WHITELIST}" --skip-barcode-check \
    --summary "${scatac_dir}/lib.summary.tsv" -o "${scatac_dir}/lib.bed"
  require_file_nonempty "${scatac_dir}/cli.bed"
  require_file_nonempty "${scatac_dir}/lib.bed"
  require_file_nonempty "${scatac_dir}/cli.summary.tsv"
  require_file_nonempty "${scatac_dir}/lib.summary.tsv"
  assert_same_rows "C05_scatac_barcode_summary.bed" \
    "${scatac_dir}/cli.bed" "${scatac_dir}/lib.bed"
  record "C05_scatac_barcode_summary" "PASS" "bed_rows_equal;summary_present"

  local bam_dir="${OUTROOT}/cases/C06_sorted_bam"
  mkdir -p "${bam_dir}"
  run_cmd "C06.cli" "${CHROMAP_BIN}" "${COMMON[@]}" --BAM --sort-bam --write-index \
    -o "${bam_dir}/cli.bam"
  run_cmd "C06.lib" "${LIBRUNNER_BIN}" "${COMMON[@]}" --BAM --sort-bam --write-index \
    -o "${bam_dir}/lib.bam"
  require_file_nonempty "${bam_dir}/cli.bam"
  require_file_nonempty "${bam_dir}/lib.bam"
  assert_index_exists "${bam_dir}/cli.bam"
  assert_index_exists "${bam_dir}/lib.bam"
  assert_bam_semantic_equal "C06_sorted_bam" "${bam_dir}/cli.bam" "${bam_dir}/lib.bam"
  record "C06_sorted_bam" "PASS" "bam_rows_equal;index_present"

  local dual_dir="${OUTROOT}/cases/C08_atac_bam_fragments"
  mkdir -p "${dual_dir}"
  run_cmd "C08.cli" "${CHROMAP_BIN}" "${COMMON[@]}" --preset atac \
    -b "${BC}" --barcode-whitelist "${WHITELIST}" --skip-barcode-check \
    --BAM --sort-bam --atac-fragments "${dual_dir}/cli.fragments.tsv.gz" \
    --summary "${dual_dir}/cli.summary.tsv" -o "${dual_dir}/cli.bam"
  run_cmd "C08.lib" "${LIBRUNNER_BIN}" "${COMMON[@]}" --preset atac \
    -b "${BC}" --barcode-whitelist "${WHITELIST}" --skip-barcode-check \
    --BAM --sort-bam --atac-fragments "${dual_dir}/lib.fragments.tsv.gz" \
    --summary "${dual_dir}/lib.summary.tsv" -o "${dual_dir}/lib.bam"
  require_file_nonempty "${dual_dir}/cli.bam"
  require_file_nonempty "${dual_dir}/lib.bam"
  require_file_nonempty "${dual_dir}/cli.fragments.tsv.gz"
  require_file_nonempty "${dual_dir}/lib.fragments.tsv.gz"
  assert_bam_semantic_equal "C08_atac_bam_fragments.bam" \
    "${dual_dir}/cli.bam" "${dual_dir}/lib.bam"
  zcat "${dual_dir}/cli.fragments.tsv.gz" | LC_ALL=C sort >"${dual_dir}/cli.fragments.rows"
  zcat "${dual_dir}/lib.fragments.tsv.gz" | LC_ALL=C sort >"${dual_dir}/lib.fragments.rows"
  cmp -s "${dual_dir}/cli.fragments.rows" "${dual_dir}/lib.fragments.rows" || \
    fail "C08_atac_bam_fragments: fragment rows differ"
  record "C08_atac_bam_fragments" "PASS" "bam_rows_equal;fragment_rows_equal"

  local y_dir="${OUTROOT}/cases/C11_y_noy_split"
  mkdir -p "${y_dir}"
  run_cmd "C11.cli" "${CHROMAP_BIN}" "${COMMON[@]}" --SAM \
    --emit-Y-bam --emit-noY-bam --emit-Y-read-names \
    --Y-read-names-output "${y_dir}/cli.Y.names.txt" \
    -o "${y_dir}/cli.sam"
  run_cmd "C11.lib" "${LIBRUNNER_BIN}" "${COMMON[@]}" --SAM \
    --emit-Y-bam --emit-noY-bam --emit-Y-read-names \
    --Y-read-names-output "${y_dir}/lib.Y.names.txt" \
    -o "${y_dir}/lib.sam"
  require_file_nonempty "${y_dir}/cli.sam"
  require_file_nonempty "${y_dir}/lib.sam"
  require_file_nonempty "${y_dir}/cli.noY.sam"
  require_file_nonempty "${y_dir}/cli.Y.sam"
  require_file_nonempty "${y_dir}/lib.noY.sam"
  require_file_nonempty "${y_dir}/lib.Y.sam"
  require_file_nonempty "${y_dir}/cli.Y.names.txt"
  require_file_nonempty "${y_dir}/lib.Y.names.txt"
  assert_same_rows "C11_y_noy_split.all" "${y_dir}/cli.sam" "${y_dir}/lib.sam"
  assert_same_rows "C11_y_noy_split.noY" "${y_dir}/cli.noY.sam" "${y_dir}/lib.noY.sam"
  assert_same_rows "C11_y_noy_split.Y" "${y_dir}/cli.Y.sam" "${y_dir}/lib.Y.sam"
  cmp -s <(LC_ALL=C sort "${y_dir}/cli.Y.names.txt") \
         <(LC_ALL=C sort "${y_dir}/lib.Y.names.txt") || \
    fail "C11_y_noy_split: Y read names differ"
  record "C11_y_noy_split" "PASS" "all/Y/noY rows equal;names_equal"

  log "PASS summary: ${SUMMARY}"
  cat "${SUMMARY}"
}

main "$@"
