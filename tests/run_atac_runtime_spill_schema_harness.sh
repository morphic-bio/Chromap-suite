#!/usr/bin/env bash
# Harness: ATAC runtime spill schema + low-memory parity (Chromap-suite).
#
# Environment (see docs/atac_runtime_spill_schema_runbook.md):
#   CHROMAP_ARTIFACT_ROOT  (default: REPO_ROOT/plans/artifacts)
#   CHROMAP, THREADS, LOW_MEM_RAM, FIXTURE_ATAC, REF, INDEX, WHITELIST
#
# Requires: samtools, python3, zcat, /usr/bin/time (optional for RSS lines).
# Skips non-destructively if the 100K PBMC fixture paths are missing.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

CHROMAP="${CHROMAP:-${REPO_ROOT}/chromap}"
CHROMAP_ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
THREADS="${THREADS:-8}"
FIXTURE_ATAC="${FIXTURE_ATAC:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac}"
INDEX="${INDEX:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index}"
REF="${REF:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa}"
WHITELIST="${WHITELIST:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"

# Force mid-run spills (tune if fixture grows).
LOW_MEM_FORCE="${LOW_MEM_FORCE:-65536}"
# Final-drain-only: large enough buffer to avoid mid-run flush on 100K.
LOW_MEM_DRAIN="${LOW_MEM_DRAIN:-$((512 * 1024 * 1024))}"

R1="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R1_001.fastq.gz"
R2="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R3_001.fastq.gz"
BC="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R2_001.fastq.gz"

BASE_FLAGS=( -t "${THREADS}" -x "${INDEX}" -r "${REF}" -1 "${R1}" -2 "${R2}" -b "${BC}"
  --barcode-whitelist "${WHITELIST}" -l 2000 --trim-adapters --remove-pcr-duplicates
  --remove-pcr-duplicates-at-cell-level --Tn5-shift )

TS="$(date +%Y%m%d_%H%M%S)"
RUN_ROOT="${CHROMAP_ARTIFACT_ROOT}/atac_runtime_spill_schema/${TS}"
SUMMARY="${RUN_ROOT}/HARNESS_SUMMARY.txt"

log() { echo "$*" | tee -a "${SUMMARY}"; }

check_aev1_sidecar() {
  local bin_path="$1"
  python3 - <<'PY' "${bin_path}"
import struct, sys
path = sys.argv[1]
with open(path, "rb") as f:
    b = f.read(32)
if len(b) != 32:
    sys.exit("header short")
magic, fv, rs, bl, nc, fl, nr = struct.unpack("<4sIIIIIQ", b)
if magic != b"AEV1":
    sys.exit("bad magic %r" % (magic,))
if rs != 24:
    sys.exit("bad record_size")
if nr <= 0:
    sys.exit("empty num_records")
sz = __import__("os").path.getsize(path)
if (sz - 32) % 24 != 0:
    sys.exit("size/record parity")
if (sz - 32) // 24 != nr:
    sys.exit("record_count mismatch header vs file")
PY
}

# Decode AEV1 sidecar to full 5-column fragment tuples and compare to BED:
# chrom,start,end,barcode,count.
check_aev1_tuple_parity_vs_bed() {
  local bin_path="$1"
  local bed_sorted="$2"
  python3 - <<'PY' "${bin_path}" "${bed_sorted}"
import struct, sys
bin_path, bed_path = sys.argv[1], sys.argv[2]
BASES = "ACGT"
def seed_to_sequence(seed, length):
    return "".join(BASES[(seed >> ((length - 1 - i) * 2)) & 3] for i in range(length))
chrom_lines = open(bin_path + ".chroms.tsv", "r", encoding="utf-8", errors="replace").read().splitlines()
id_to_name = {}
for ln in chrom_lines:
    parts = ln.split("\t", 1)
    if len(parts) < 2:
        continue
    id_to_name[int(parts[0])] = parts[1]
with open(bin_path, "rb") as f:
    hdr = f.read(32)
    magic, fv, rs, bl, nc, fl, nr = struct.unpack("<4sIIIIIQ", hdr)
    body = f.read()
if magic != b"AEV1" or fv != 1 or rs != 24:
    sys.exit("bad sidecar header")
if bl <= 0 or bl > 32:
    sys.exit("unsupported barcode_length %d" % bl)
if len(body) != nr * 24:
    sys.exit("sidecar body size mismatch")
rows = []
for i in range(nr):
    chunk = body[i * 24 : (i + 1) * 24]
    cid, start, end, count, bc = struct.unpack("<iiiIQ", chunk)
    chrom = id_to_name.get(cid)
    if chrom is None:
        sys.exit("unknown chrom_id %d" % cid)
    rows.append((chrom, int(start), int(end), seed_to_sequence(bc, bl), int(count)))
rows.sort()
bed_rows = []
with open(bed_path, "r", encoding="utf-8", errors="replace") as bf:
    for ln in bf:
        c = ln.rstrip("\n").split("\t")
        if len(c) < 5:
            sys.exit("bad bed line")
        bed_rows.append((c[0], int(c[1]), int(c[2]), c[3], int(c[4])))
bed_rows.sort()
if rows != bed_rows:
    sys.exit("sidecar decoded fragment multiset != BED fragment multiset")
PY
}

mid_batch_flush_count() {
  local logf="$1"
  # Emitted by ProcessAndOutputMappingsInLowMemoryFromOverflow
  grep -oE 'Low-memory overflow: mid-batch flush count: [0-9]+' "${logf}" | head -1 | grep -oE '[0-9]+$' || echo "0"
}

check_bam_coordinate_sorted_records() {
  local bam="$1"
  python3 - <<'PY' "${bam}"
import subprocess, sys
bam = sys.argv[1]
hdr = subprocess.check_output(["samtools", "view", "-H", bam], text=True)
order = []
for ln in hdr.splitlines():
    if not ln.startswith("@SQ"):
        continue
    sn = None
    for f in ln.split("\t")[1:]:
        if f.startswith("SN:"):
            sn = f[3:]
            break
    if sn is None:
        sys.exit("bad @SQ line")
    order.append(sn)
name_i = {n: i for i, n in enumerate(order)}
p = subprocess.Popen(["samtools", "view", bam], stdout=subprocess.PIPE, text=True)
assert p.stdout is not None
prev = None
for line in p.stdout:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 4:
        continue
    flag = int(parts[1])
    if flag & 0x4:
        continue
    chrom, pos = parts[2], int(parts[3])
    if chrom == "*":
        continue
    if chrom not in name_i:
        sys.exit("unknown chrom %r" % chrom)
    key = (name_i[chrom], pos)
    if prev is not None and key < prev:
        sys.exit("BAM records not in non-decreasing (ref_index, pos) order")
    prev = key
p.stdout.close()
if p.wait() != 0:
    sys.exit("samtools view failed")
PY
}

check_chroms() {
  local bin_path="$1"
  local ref_fa="$2"
  local chroms="${bin_path}.chroms.tsv"
  [[ -f "${chroms}" ]] || { echo "missing ${chroms}" >&2; return 1; }
  local nref
  nref=$(grep -c "^>" "${ref_fa}" || true)
  local nlines
  nlines=$(wc -l <"${chroms}")
  [[ "${nlines}" -eq "${nref}" ]] || {
    echo "chroms lines ${nlines} != ref seqs ${nref}" >&2
    return 1
  }
}

check_sorted_so() {
  samtools view -H "$1" | grep -E "^@HD" | grep -q "SO:coordinate"
}

run_with_time() {
  local logf="$1"
  shift
  if [[ -x /usr/bin/time ]]; then
    /usr/bin/time -v "$@" 2>&1 | tee "${logf}"
  else
    "$@" 2>&1 | tee "${logf}"
  fi
}

if [[ ! -x "${CHROMAP}" ]]; then
  echo "ERROR: build chromap first (${CHROMAP})" >&2
  exit 2
fi

if [[ ! -f "${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R1_001.fastq.gz" ]]; then
  echo "SKIP: fixture not found at FIXTURE_ATAC=${FIXTURE_ATAC}" >&2
  exit 0
fi

mkdir -p "${RUN_ROOT}/logs"
: >"${SUMMARY}"
log "run_root ${RUN_ROOT}"
log "chromap ${CHROMAP}"
log "threads ${THREADS}"
log "low_mem_force_bytes ${LOW_MEM_FORCE}"
log "low_mem_drain_bytes ${LOW_MEM_DRAIN}"

# --- Case 1: fragment-only baseline (no --low-mem) ---
log "=== case1_fragment_baseline ==="
run_with_time "${RUN_ROOT}/logs/case1.log" "${CHROMAP}" "${BASE_FLAGS[@]}" --BED \
  -o "${RUN_ROOT}/case1_frag.bed" \
  --summary "${RUN_ROOT}/logs/summary_case1.tsv"
SORTED_BASE="${RUN_ROOT}/case1_frag.sorted.bed"
LC_ALL=C sort "${RUN_ROOT}/case1_frag.bed" -o "${SORTED_BASE}"

# --- Case 2: fragment-only + low-mem forced spill ---
log "=== case2_fragment_lowmem_forced ==="
run_with_time "${RUN_ROOT}/logs/case2.log" "${CHROMAP}" "${BASE_FLAGS[@]}" --BED --low-mem \
  --low-mem-ram "${LOW_MEM_FORCE}" \
  -o "${RUN_ROOT}/case2_frag.bed" \
  --summary "${RUN_ROOT}/logs/summary_case2.tsv"
if grep -q "Processing .* overflow files" "${RUN_ROOT}/logs/case2.log"; then
  MB2="$(mid_batch_flush_count "${RUN_ROOT}/logs/case2.log")"
  [[ "${MB2}" -ge 1 ]] || {
    echo "ERROR: case2 mid-batch flush count ${MB2} (expected >=1 when overflow merge runs)" >&2
    exit 1
  }
  log "case2 mid-batch flush count ${MB2} (>=1): OK"
fi
if ! cmp -s "${SORTED_BASE}" <(LC_ALL=C sort "${RUN_ROOT}/case2_frag.bed"); then
  echo "ERROR: case1 vs case2 sorted fragments differ" >&2
  exit 1
fi
log "case1 vs case2 fragment multiset: OK"

# --- Case 3: dual BAM + gz fragments + sidecar, low-mem forced ---
log "=== case3_dual_lowmem_forced ==="
BIN3="${RUN_ROOT}/case3_fragments.bin"
run_with_time "${RUN_ROOT}/logs/case3.log" "${CHROMAP}" "${BASE_FLAGS[@]}" --BAM \
  --atac-fragments "${RUN_ROOT}/case3_fragments.tsv.gz" \
  --atac-fragment-binary-output "${BIN3}" \
  --low-mem --low-mem-ram "${LOW_MEM_FORCE}" \
  -o "${RUN_ROOT}/case3.bam" \
  --summary "${RUN_ROOT}/logs/summary_case3.tsv"
if grep -q "Processing .* overflow files" "${RUN_ROOT}/logs/case3.log"; then
  MB3="$(mid_batch_flush_count "${RUN_ROOT}/logs/case3.log")"
  [[ "${MB3}" -ge 1 ]] || {
    echo "ERROR: case3 mid-batch flush count ${MB3} (expected >=1 when overflow merge runs)" >&2
    exit 1
  }
  log "case3 mid-batch flush count ${MB3} (>=1): OK"
fi
samtools quickcheck "${RUN_ROOT}/case3.bam"
FRAG_LINES="$(zcat "${RUN_ROOT}/case3_fragments.tsv.gz" | wc -l)"
BAM_REC="$(samtools view -c "${RUN_ROOT}/case3.bam")"
EXP=$((FRAG_LINES * 2))
[[ "${BAM_REC}" -eq "${EXP}" ]] || {
  echo "ERROR: case3 BAM rec ${BAM_REC} != 2*frag ${EXP}" >&2
  exit 1
}
python3 "${SCRIPT_DIR}/check_dual_bam_mate_pairs.py" "${RUN_ROOT}/case3.bam" | tee -a "${SUMMARY}"
check_aev1_sidecar "${BIN3}"
check_chroms "${BIN3}" "${REF}"
check_aev1_tuple_parity_vs_bed "${BIN3}" "${SORTED_BASE}"
if ! cmp -s "${SORTED_BASE}" <(LC_ALL=C sort <(zcat "${RUN_ROOT}/case3_fragments.tsv.gz")); then
  echo "ERROR: case1 vs case3 dual fragments differ" >&2
  exit 1
fi
log "case1 vs case3 dual fragments: OK"

# --- Case 4: dual low-mem final-drain-only (high RAM threshold) ---
log "=== case4_dual_lowmem_drain_only ==="
BIN4="${RUN_ROOT}/case4_fragments.bin"
run_with_time "${RUN_ROOT}/logs/case4.log" "${CHROMAP}" "${BASE_FLAGS[@]}" --BAM \
  --atac-fragments "${RUN_ROOT}/case4_fragments.tsv.gz" \
  --atac-fragment-binary-output "${BIN4}" \
  --low-mem --low-mem-ram "${LOW_MEM_DRAIN}" \
  -o "${RUN_ROOT}/case4.bam" \
  --summary "${RUN_ROOT}/logs/summary_case4.tsv"
MB4="$(mid_batch_flush_count "${RUN_ROOT}/logs/case4.log")"
[[ "${MB4}" -eq 0 ]] || {
  echo "ERROR: case4 expected mid-batch flush count 0 (drain-only), got ${MB4}" >&2
  exit 1
}
log "case4 mid-batch flush count 0 (drain-only): OK"
if grep -q "Processing .* overflow files" "${RUN_ROOT}/logs/case4.log"; then
  log "case4: k-way merge of final-drain overflow files (expected)"
else
  log "case4: no overflow merge message (unexpected for --low-mem)"
fi
samtools quickcheck "${RUN_ROOT}/case4.bam"
check_aev1_sidecar "${BIN4}"
check_chroms "${BIN4}" "${REF}"
check_aev1_tuple_parity_vs_bed "${BIN4}" "${SORTED_BASE}"
if ! cmp -s "${SORTED_BASE}" <(LC_ALL=C sort <(zcat "${RUN_ROOT}/case4_fragments.tsv.gz")); then
  echo "ERROR: case1 vs case4 dual fragments differ" >&2
  exit 1
fi
log "case4 fragments vs baseline: OK"

# --- Case 5: dual + coordinate sort ---
log "=== case5_dual_sort_bam ==="
BIN5="${RUN_ROOT}/case5_fragments.bin"
run_with_time "${RUN_ROOT}/logs/case5.log" "${CHROMAP}" "${BASE_FLAGS[@]}" --BAM --sort-bam \
  --atac-fragments "${RUN_ROOT}/case5_fragments.tsv.gz" \
  --atac-fragment-binary-output "${BIN5}" \
  --low-mem --low-mem-ram "${LOW_MEM_FORCE}" \
  -o "${RUN_ROOT}/case5_sorted.bam" \
  --summary "${RUN_ROOT}/logs/summary_case5.tsv"
samtools quickcheck "${RUN_ROOT}/case5_sorted.bam"
check_sorted_so "${RUN_ROOT}/case5_sorted.bam" || {
  echo "ERROR: case5 BAM @HD not coordinate sorted" >&2
  exit 1
}
check_bam_coordinate_sorted_records "${RUN_ROOT}/case5_sorted.bam" || {
  echo "ERROR: case5 BAM records not coordinate-sorted" >&2
  exit 1
}
check_aev1_sidecar "${BIN5}"
check_chroms "${BIN5}" "${REF}"
check_aev1_tuple_parity_vs_bed "${BIN5}" "${SORTED_BASE}"
FRAG5="$(zcat "${RUN_ROOT}/case5_fragments.tsv.gz" | wc -l)"
BAM5="$(samtools view -c "${RUN_ROOT}/case5_sorted.bam")"
[[ "$((FRAG5 * 2))" -eq "${BAM5}" ]] || {
  echo "ERROR: case5 BAM count" >&2
  exit 1
}
log "case5 sort-bam: OK"

# --- Case 6 (optional): Y/noY dual — requires emit flags; skip by default ---
if [[ "${RUN_ATAC_Y_HARNESS:-0}" == "1" ]]; then
  log "=== case6_dual_y (RUN_ATAC_Y_HARNESS=1) ==="
  echo "TODO: extend with emit-noY-bam / emit-Y-bam when a small ATAC Y fixture is wired" | tee -a "${SUMMARY}"
else
  log "=== case6_dual_y SKIPPED (set RUN_ATAC_Y_HARNESS=1 to attempt) ==="
fi

log "PASS — artifacts under ${RUN_ROOT}"
echo "Full log: ${SUMMARY}"
