#!/usr/bin/env bash
# FASTQ intake smoke: the zlib (kseq) reader and the BGZF reader mirrored from
# STAR Suite deliver the same records, and Chromap's outputs do not depend on
# which reader, how many inflate workers, or how the input is fed.
#
# Record level (tests/fastq_intake_harness):
#   - dumps of every (name, comment, sequence, quality) record are identical
#     for kseq-serial and kseq-threads on gzip and BGZF, and for the BGZF
#     provider with 0, 1, 3 and 17 inflate workers; the fixture has more than
#     two 500,000-pair batches, header comments, and a record with an empty
#     sequence in every file (kseq skips it; so must the BGZF reader);
#   - a CRLF fixture reads the same through both readers;
#   - auto selects BGZF only when every file is BGZF;
#   - the BGZF provider rejects mispaired names and unequal record counts.
# Chromap level (BED fragments and --summary, byte for byte):
#   - --input-bgzf-mode off on gzip = auto on BGZF = on on BGZF with explicit
#     reader threads = auto on gzip, for the CLI and chromap_lib_runner;
#   - a lane mix (one BGZF lane, one gzip lane) = both lanes through zlib;
#   - --atac-sidecar-only sidecars match across the readers;
#   - FIFO inputs fed by a producer that interleaves records across the
#     FIFOs (bulk R1/R2, and R1/barcode/R2 for a mergeable spill worker) give
#     the same output as regular files, at 2 and 16 threads, within a timeout;
#   - invalid settings are rejected.
#
# Environment: CHROMAP_BIN, LIBRUNNER_BIN, HARNESS_BIN, BGZIP (default bgzip),
# THREADS (default 4), PAIRS (default 1100000), BUILD (1 builds with
# MAKE_JOBS, default 2), OUTROOT.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)-$$}"
OUTROOT="${OUTROOT:-${ARTIFACT_ROOT}/fastq_intake_smoke/${RUN_ID}}"
CHROMAP_BIN="${CHROMAP_BIN:-${REPO_ROOT}/chromap}"
LIBRUNNER_BIN="${LIBRUNNER_BIN:-${REPO_ROOT}/chromap_lib_runner}"
HARNESS_BIN="${HARNESS_BIN:-${REPO_ROOT}/tests/fastq_intake_harness}"
BGZIP="${BGZIP:-bgzip}"
THREADS="${THREADS:-4}"
PAIRS="${PAIRS:-1100000}"
BUILD="${BUILD:-1}"
FIFO_TIMEOUT="${FIFO_TIMEOUT:-600}"

mkdir -p "${OUTROOT}"/{fixture,index,cases,logs,fifo}
SUMMARY="${OUTROOT}/summary.tsv"
printf 'case\tstatus\tdetail\n' >"${SUMMARY}"

log() { printf '[fastq-intake-smoke] %s\n' "$*" >&2; }
fail() { printf '[fastq-intake-smoke] FAIL: %s\n' "$*" >&2; exit 1; }
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
    fail "${id}: invalid input or setting was accepted"
  fi
  grep -q -- "${want}" "${OUTROOT}/logs/${id}.stderr" ||
    fail "${id}: rejection message lacks '${want}'"
  record "${id}" "PASS" "rejected: ${want}"
}

generate_fixture() {
  python3 - "$1" "$2" <<'PY'
import os
import random
import sys

out, n_pairs = sys.argv[1], int(sys.argv[2])
rng = random.Random(20260928)
read_len = 50


def randseq(n):
    return "".join(rng.choice("ACGT") for _ in range(n))


comp = str.maketrans("ACGTN", "TGCAN")
chroms = {"chr1": randseq(30000), "chr2": randseq(20000), "chrY": randseq(10000)}
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
picks = rng.choices(range(len(sites)), weights, k=n_pairs)
quals = ["".join(rng.choice("F:,#") for _ in range(read_len)) for _ in range(64)]
empty = {1000, n_pairs * 2 // 3}  # empty sequence in every file

with open(os.path.join(out, "read1.fq"), "w") as r1, \
     open(os.path.join(out, "read2.fq"), "w") as r2, \
     open(os.path.join(out, "barcode.fq"), "w") as bc:
    for i, pick in enumerate(picks):
        name, start, length = sites[pick]
        seq = chroms[name]
        left = seq[start:start + read_len]
        right = seq[start + length - read_len:start + length].translate(comp)[::-1]
        a, b = (left, right) if i % 2 else (right, left)
        barcode = whitelist[i % 8] if i % 5 else whitelist[(i * 7) % 32]
        head = f"@SRRX.{i + 1} F00:1:{i % 97}:{i}"
        q = quals[i % 64]
        if i in empty:
            a = b = barcode = ""
            q = ""
        r1.write(f"{head} length={len(a)}\n{a}\n+\n{q}\n")
        r2.write(f"{head} length={len(b)}\n{b}\n+\n{q}\n")
        bc.write(f"{head} length={len(barcode)}\n{barcode}\n+\n{q[:len(barcode)]}\n")
PY
}

# The harness dump of one lane with one reader; prints its sha256.
dump_sha() {
  local id="$1"; shift
  "${HARNESS_BIN}" "$@" --dump - 2>"${OUTROOT}/logs/${id}.stderr" | sha256sum |
    cut -d' ' -f1
}

check_same() {  # check_same ID EXPECTED ACTUAL DETAIL
  [[ "$2" == "$3" ]] || fail "$1: $4 differs ($2 vs $3)"
  record "$1" "PASS" "$4"
}

cmp_outputs() {  # cmp_outputs ID DIR_A DIR_B FILE...
  local id="$1" a="$2" b="$3"; shift 3
  local f
  for f in "$@"; do
    [[ -s "${a}/${f}" ]] || fail "${id}: ${a}/${f} missing or empty"
    cmp "${a}/${f}" "${b}/${f}" || fail "${id}: ${f} differs (${a} vs ${b})"
  done
}

fifo_producer() {  # fifo_producer FIFO... -- FILE...: interleave records
  python3 - "$@" <<'PY'
import sys

args = sys.argv[1:]
cut = args.index("--")
fifos, files = args[:cut], args[cut + 1:]
readers = [open(path, "rb") for path in files]
writers = []
for path in fifos:  # the same order Chromap opens them
    writers.append(open(path, "wb", buffering=0))
while True:
    records = [b"".join(r.readline() for _ in range(4)) for r in readers]
    if not any(records):
        break
    for w, rec in zip(writers, records):
        w.write(rec)
for w in writers:
    w.close()
PY
}

main() {
  if [[ "${BUILD}" == "1" ]]; then
    log "building chromap, chromap_lib_runner and the harness (-j${MAKE_JOBS:-2})"
    (cd "${REPO_ROOT}" && make -j"${MAKE_JOBS:-2}" chromap chromap_lib_runner \
      tests/fastq_intake_harness) \
      >"${OUTROOT}/logs/make.stdout" 2>"${OUTROOT}/logs/make.stderr"
  fi
  local bin
  for bin in "${CHROMAP_BIN}" "${LIBRUNNER_BIN}" "${HARNESS_BIN}"; do
    [[ -x "${bin}" ]] || fail "missing ${bin}"
  done
  command -v "${BGZIP}" >/dev/null || fail "bgzip not found (set BGZIP)"

  local fx="${OUTROOT}/fixture"
  log "generating ${PAIRS} read pairs"
  generate_fixture "${fx}" "${PAIRS}"
  local r
  for r in read1 read2 barcode; do
    gzip -c "${fx}/${r}.fq" >"${fx}/${r}.fq.gz"
    "${BGZIP}" -c "${fx}/${r}.fq" >"${fx}/${r}.fq.bgz"
    # Two lanes: the first 600,000 records and the rest.
    head -n 2400000 "${fx}/${r}.fq" | "${BGZIP}" -c >"${fx}/${r}.laneA.fq.bgz"
    tail -n +2400001 "${fx}/${r}.fq" | gzip -c >"${fx}/${r}.laneB.fq.gz"
    head -n 2400000 "${fx}/${r}.fq" | gzip -c >"${fx}/${r}.laneA.fq.gz"
  done
  # Mispaired read 2 (one renamed record) and a read 2 one record short.
  awk 'NR==4*300000+1 {sub(/^@SRRX\.[0-9]+/, "@SRRX.other")} {print}' \
    "${fx}/read2.fq" | "${BGZIP}" -c >"${fx}/read2.misnamed.fq.bgz"
  head -n -4 "${fx}/read2.fq" | "${BGZIP}" -c >"${fx}/read2.short.fq.bgz"
  printf '@c1 x\r\nACGT\r\n+\r\nIIII\r\n@c2\r\nAC\r\n+c2\r\nII\r\n@c3 y z\r\nGGT\r\n+\r\nIII\r\n' \
    >"${fx}/crlf1.fq"
  cp "${fx}/crlf1.fq" "${fx}/crlf2.fq"
  for r in crlf1 crlf2; do
    gzip -c "${fx}/${r}.fq" >"${fx}/${r}.fq.gz"
    "${BGZIP}" -c "${fx}/${r}.fq" >"${fx}/${r}.fq.bgz"
  done

  # ---- record level ----
  local gz=("${fx}/read1.fq.gz" "${fx}/read2.fq.gz" "${fx}/barcode.fq.gz")
  local bgz=("${fx}/read1.fq.bgz" "${fx}/read2.fq.bgz" "${fx}/barcode.fq.bgz")
  local ref_sha
  ref_sha="$(dump_sha H00 --reader kseq-serial "${gz[@]}")"
  grep -q "pairs=$((PAIRS - 2)) " "${OUTROOT}/logs/H00.stderr" ||
    fail "H00: expected $((PAIRS - 2)) pairs (two empty records skipped)"
  check_same H01_kseq_threads_gz "${ref_sha}" \
    "$(dump_sha H01 --reader kseq-threads "${gz[@]}")" "record dump, kseq threads on gzip"
  check_same H02_kseq_serial_bgzf "${ref_sha}" \
    "$(dump_sha H02 --reader kseq-serial "${bgz[@]}")" "record dump, kseq on BGZF"
  local workers
  for workers in 0 1 3 17; do
    check_same "H03_bgzf_workers_${workers}" "${ref_sha}" \
      "$(dump_sha "H03_${workers}" --reader bgzf --threads 1 \
        --bgzf-threads "${workers}" "${bgz[@]}")" \
      "record dump, BGZF provider, ${workers} inflate workers"
  done
  check_same H04_bgzf_bulk "$(dump_sha H04a --reader kseq-serial "${gz[@]:0:2}")" \
    "$(dump_sha H04b --reader bgzf --threads 8 "${bgz[@]:0:2}")" \
    "record dump, two files (bulk)"
  check_same H05_crlf "$(dump_sha H05a --reader kseq-serial \
      "${fx}/crlf1.fq.gz" "${fx}/crlf2.fq.gz")" \
    "$(dump_sha H05b --reader bgzf --threads 4 "${fx}/crlf1.fq.bgz" "${fx}/crlf2.fq.bgz")" \
    "record dump, CRLF line endings"
  dump_sha H06 --reader auto --threads 8 "${bgz[@]}" >/dev/null
  grep -q "auto selected bgzf" "${OUTROOT}/logs/H06.stderr" || fail "H06: auto did not select BGZF"
  dump_sha H07 --reader auto --threads 8 "${bgz[0]}" "${gz[1]}" "${bgz[2]}" >/dev/null
  grep -q "auto selected kseq-threads (.*read2.fq.gz: not BGZF)" "${OUTROOT}/logs/H07.stderr" ||
    fail "H07: auto did not fall back for a gzip file"
  record H06_H07_auto PASS "auto: BGZF when all files are BGZF, zlib otherwise"
  expect_reject H08_misnamed "do not pair at pair 299999" "${HARNESS_BIN}" --reader bgzf \
    --threads 8 "${bgz[0]}" "${fx}/read2.misnamed.fq.bgz" "${bgz[2]}"
  expect_reject H09_short "record counts differ" "${HARNESS_BIN}" --reader bgzf \
    --threads 8 "${bgz[0]}" "${fx}/read2.short.fq.bgz" "${bgz[2]}"

  # ---- Chromap level ----
  local index="${OUTROOT}/index/ref.index"
  run_cmd index "${CHROMAP_BIN}" -i -r "${fx}/ref.fa" -o "${index}"
  local common=(-x "${index}" -r "${fx}/ref.fa" --barcode-whitelist "${fx}/whitelist.txt"
    -l 2000 --trim-adapters --remove-pcr-duplicates --remove-pcr-duplicates-at-cell-level
    --Tn5-shift)
  chromap_case() {  # chromap_case ID TOOL THREADS R1 R2 BC EXTRA...
    local id="$1" tool="$2" t="$3" r1="$4" r2="$5" b="$6"; shift 6
    local dir="${OUTROOT}/cases/${id}"
    mkdir -p "${dir}/tmp"
    run_cmd "${id}" "${tool}" -t "${t}" "${common[@]}" -1 "${r1}" -2 "${r2}" -b "${b}" \
      --temp-dir "${dir}/tmp" --summary "${dir}/summary.tsv" -o "${dir}/fragments.bed" "$@"
  }
  chromap_case C00_off_gz "${CHROMAP_BIN}" "${THREADS}" "${gz[@]}" --input-bgzf-mode off
  grep -q "zlib reader (--input-bgzf-mode off)" "${OUTROOT}/logs/C00_off_gz.stderr" ||
    fail "C00: zlib reader not reported"
  chromap_case C01_auto_bgzf "${CHROMAP_BIN}" "${THREADS}" "${bgz[@]}"
  grep -q "FASTQ lane 1: BGZF reader" "${OUTROOT}/logs/C01_auto_bgzf.stderr" ||
    fail "C01: BGZF reader not used"
  grep -q "Barcode abundance input 1: BGZF reader" "${OUTROOT}/logs/C01_auto_bgzf.stderr" ||
    fail "C01: barcode abundance pass did not use the BGZF reader"
  chromap_case C02_on_bgzf_w5 "${CHROMAP_BIN}" "${THREADS}" "${bgz[@]}" \
    --input-bgzf-mode on --input-bgzf-reader-threads 5
  chromap_case C03_auto_gz "${CHROMAP_BIN}" "${THREADS}" "${gz[@]}"
  chromap_case C04_lib_auto_bgzf "${LIBRUNNER_BIN}" "${THREADS}" "${bgz[@]}"
  local c
  for c in C01_auto_bgzf C02_on_bgzf_w5 C03_auto_gz C04_lib_auto_bgzf; do
    cmp_outputs "${c}" "${OUTROOT}/cases/C00_off_gz" "${OUTROOT}/cases/${c}" \
      fragments.bed summary.tsv
    record "${c}" PASS "fragments.bed and summary byte-identical to C00 (off, gzip)"
  done
  chromap_case C05_lanes_off "${CHROMAP_BIN}" "${THREADS}" \
    "${fx}/read1.laneA.fq.gz,${fx}/read1.laneB.fq.gz" \
    "${fx}/read2.laneA.fq.gz,${fx}/read2.laneB.fq.gz" \
    "${fx}/barcode.laneA.fq.gz,${fx}/barcode.laneB.fq.gz" --input-bgzf-mode off
  chromap_case C06_lanes_mixed "${CHROMAP_BIN}" "${THREADS}" \
    "${fx}/read1.laneA.fq.bgz,${fx}/read1.laneB.fq.gz" \
    "${fx}/read2.laneA.fq.bgz,${fx}/read2.laneB.fq.gz" \
    "${fx}/barcode.laneA.fq.bgz,${fx}/barcode.laneB.fq.gz"
  grep -q "FASTQ lane 1: BGZF reader" "${OUTROOT}/logs/C06_lanes_mixed.stderr" &&
    grep -q "FASTQ lane 2: zlib reader" "${OUTROOT}/logs/C06_lanes_mixed.stderr" ||
    fail "C06: expected BGZF for lane 1 and zlib for lane 2"
  cmp_outputs C06 "${OUTROOT}/cases/C05_lanes_off" "${OUTROOT}/cases/C06_lanes_mixed" \
    fragments.bed summary.tsv
  record C06_lanes_mixed PASS "BGZF lane + gzip lane byte-identical to both lanes through zlib"
  local sc
  for sc in off:gz auto:bgz; do
    local mode="${sc%%:*}" ext="${sc##*:}" dir="${OUTROOT}/cases/C07_sidecar_${sc%%:*}"
    mkdir -p "${dir}/tmp"
    run_cmd "C07_sidecar_${mode}" "${CHROMAP_BIN}" -t "${THREADS}" "${common[@]}" \
      -1 "${fx}/read1.fq.${ext}" -2 "${fx}/read2.fq.${ext}" -b "${fx}/barcode.fq.${ext}" \
      --temp-dir "${dir}/tmp" --summary "${dir}/summary.tsv" --input-bgzf-mode "${mode}" \
      --atac-sidecar-only --atac-fragment-binary-output "${dir}/atac_fragments.bin"
  done
  cmp_outputs C07 "${OUTROOT}/cases/C07_sidecar_off" "${OUTROOT}/cases/C07_sidecar_auto" \
    atac_fragments.bin atac_fragments.bin.chroms.tsv summary.tsv
  record C07_sidecar_only PASS "AEV1 sidecar, chroms and summary byte-identical"
  # Y/noY FASTQ output carries read comments. Same file names in both runs
  # (the output names derive from the input names), different compression.
  local ymode
  for ymode in off auto; do
    local ydir="${OUTROOT}/cases/C08_ynoy_${ymode}"
    mkdir -p "${ydir}/in" "${ydir}/tmp"
    for r in read1 read2; do
      if [[ "${ymode}" == off ]]; then
        head -n 800000 "${fx}/${r}.fq" | gzip -c >"${ydir}/in/${r}.fq.gz"
      else
        head -n 800000 "${fx}/${r}.fq" | "${BGZIP}" -c >"${ydir}/in/${r}.fq.gz"
      fi
    done
    run_cmd "C08_ynoy_${ymode}" "${CHROMAP_BIN}" -t "${THREADS}" --SAM -x "${index}" \
      -r "${fx}/ref.fa" -1 "${ydir}/in/read1.fq.gz" -2 "${ydir}/in/read2.fq.gz" \
      --emit-Y-read-names --emit-Y-noY-fastq --emit-Y-noY-fastq-compression none \
      --temp-dir "${ydir}/tmp" --input-bgzf-mode "${ymode}" -o "${ydir}/out.sam"
    grep -v '^@PG' "${ydir}/out.sam" >"${ydir}/out.body.sam"
  done
  grep -q "FASTQ lane 1: BGZF reader" "${OUTROOT}/logs/C08_ynoy_auto.stderr" ||
    fail "C08: BGZF reader not used"
  local yfiles
  yfiles="$(cd "${OUTROOT}/cases/C08_ynoy_off/y_separated" && ls)"
  [[ -n "${yfiles}" ]] || fail "C08: no Y/noY FASTQ files"
  # shellcheck disable=SC2086
  cmp_outputs C08 "${OUTROOT}/cases/C08_ynoy_off/y_separated" \
    "${OUTROOT}/cases/C08_ynoy_auto/y_separated" ${yfiles}
  cmp_outputs C08 "${OUTROOT}/cases/C08_ynoy_off" "${OUTROOT}/cases/C08_ynoy_auto" out.body.sam
  record C08_ynoy_fastq PASS "Y/noY FASTQ (with comments) and SAM records byte-identical"

  # ---- FIFO inputs with an interleaving producer ----
  local plain=("${fx}/read1.fq" "${fx}/read2.fq" "${fx}/barcode.fq")
  local t
  for t in 2 16; do
    local ref_dir="${OUTROOT}/cases/F01_bulk_files_t${t}" fifo_dir="${OUTROOT}/cases/F01_bulk_fifo_t${t}"
    mkdir -p "${ref_dir}/tmp" "${fifo_dir}/tmp"
    run_cmd "F01_bulk_files_t${t}" "${CHROMAP_BIN}" -t "${t}" -x "${index}" -r "${fx}/ref.fa" \
      -l 2000 --remove-pcr-duplicates -1 "${plain[0]}" -2 "${plain[1]}" \
      --temp-dir "${ref_dir}/tmp" --summary "${ref_dir}/summary.tsv" -o "${ref_dir}/fragments.bed"
    local f1="${OUTROOT}/fifo/b1_t${t}" f2="${OUTROOT}/fifo/b2_t${t}"
    rm -f "${f1}" "${f2}"; mkfifo "${f1}" "${f2}"
    fifo_producer "${f1}" "${f2}" -- "${plain[0]}" "${plain[1]}" &
    local producer=$!
    run_cmd "F01_bulk_fifo_t${t}" timeout "${FIFO_TIMEOUT}" "${CHROMAP_BIN}" -t "${t}" \
      -x "${index}" -r "${fx}/ref.fa" -l 2000 --remove-pcr-duplicates -1 "${f1}" -2 "${f2}" \
      --temp-dir "${fifo_dir}/tmp" --summary "${fifo_dir}/summary.tsv" -o "${fifo_dir}/fragments.bed"
    wait "${producer}" || fail "F01 t${t}: FIFO producer failed"
    grep -q "zlib reader (.*not a regular file)" "${OUTROOT}/logs/F01_bulk_fifo_t${t}.stderr" ||
      fail "F01 t${t}: FIFO lane did not use the zlib reader"
    cmp_outputs "F01_t${t}" "${ref_dir}" "${fifo_dir}" fragments.bed summary.tsv
    record "F01_bulk_fifo_t${t}" PASS "interleaved R1/R2 FIFOs = regular files, ${t} threads"

    local spill=(--create-mergeable-spill-record SPILL --mergeable-spill-sample-id s
      --mergeable-spill-input-id i --mergeable-spill-shard-ordinal 0
      --mergeable-spill-shard-count 1)
    local sref="${OUTROOT}/cases/F02_spill_files_t${t}" sfifo="${OUTROOT}/cases/F02_spill_fifo_t${t}"
    mkdir -p "${sref}/tmp" "${sfifo}/tmp"
    run_cmd "F02_spill_files_t${t}" "${CHROMAP_BIN}" -t "${t}" "${common[@]}" \
      -1 "${plain[0]}" -2 "${plain[1]}" -b "${plain[2]}" --temp-dir "${sref}/tmp" \
      "${spill[@]/SPILL/${sref}/shard.atacms}"
    local g1="${OUTROOT}/fifo/s1_t${t}" g2="${OUTROOT}/fifo/s2_t${t}" g3="${OUTROOT}/fifo/s3_t${t}"
    rm -f "${g1}" "${g2}" "${g3}"; mkfifo "${g1}" "${g2}" "${g3}"
    # Chromap opens read 1, read 2, then the barcode file; the producer writes
    # each record to read 1, the barcode and then read 2.
    python3 - "${g1}" "${g2}" "${g3}" "${plain[@]}" <<'PY' &
import sys

f1, f2, f3, p1, p2, p3 = sys.argv[1:]
readers = [open(p, "rb") for p in (p1, p3, p2)]
w1 = open(f1, "wb", buffering=0)
w2 = open(f2, "wb", buffering=0)
w3 = open(f3, "wb", buffering=0)
writers = [w1, w3, w2]  # read 1, barcode, read 2
while True:
    records = [b"".join(r.readline() for _ in range(4)) for r in readers]
    if not any(records):
        break
    for w, rec in zip(writers, records):
        w.write(rec)
for w in writers:
    w.close()
PY
    producer=$!
    run_cmd "F02_spill_fifo_t${t}" timeout "${FIFO_TIMEOUT}" "${CHROMAP_BIN}" -t "${t}" \
      "${common[@]}" -1 "${g1}" -2 "${g2}" -b "${g3}" --temp-dir "${sfifo}/tmp" \
      "${spill[@]/SPILL/${sfifo}/shard.atacms}"
    wait "${producer}" || fail "F02 t${t}: FIFO producer failed"
    cmp_outputs "F02_t${t}" "${sref}" "${sfifo}" shard.atacms
    record "F02_spill_fifo_t${t}" PASS "interleaved R1/barcode/R2 FIFOs = regular files, ${t} threads"
  done

  # ---- rejections ----
  expect_reject R01_on_gzip "requires regular BGZF FASTQ input" "${CHROMAP_BIN}" \
    -t 2 "${common[@]}" -1 "${gz[0]}" -2 "${gz[1]}" -b "${gz[2]}" \
    --input-bgzf-mode on -o "${OUTROOT}/cases/x.bed"
  expect_reject R02_bad_mode "must be auto, on or off" "${CHROMAP_BIN}" -t 2 \
    "${common[@]}" -1 "${bgz[0]}" -2 "${bgz[1]}" -b "${bgz[2]}" --input-bgzf-mode yes \
    -o "${OUTROOT}/cases/x.bed"
  expect_reject R03_bad_threads "must be >= 0" "${CHROMAP_BIN}" -t 2 "${common[@]}" \
    -1 "${bgz[0]}" -2 "${bgz[1]}" -b "${bgz[2]}" --input-bgzf-reader-threads -1 \
    -o "${OUTROOT}/cases/x.bed"
  expect_reject R04_lib_on_gzip "requires regular BGZF FASTQ input" "${LIBRUNNER_BIN}" \
    -t 2 "${common[@]}" -1 "${gz[0]}" -2 "${gz[1]}" -b "${gz[2]}" \
    --input-bgzf-mode on -o "${OUTROOT}/cases/x.bed"
  expect_reject R05_misnamed "do not pair" "${CHROMAP_BIN}" -t 2 "${common[@]}" \
    -1 "${bgz[0]}" -2 "${fx}/read2.misnamed.fq.bgz" -b "${bgz[2]}" -o "${OUTROOT}/cases/x.bed"

  log "PASS: ${SUMMARY}"
  cat "${SUMMARY}" >&2
}

main "$@"
