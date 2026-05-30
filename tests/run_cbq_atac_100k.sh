#!/usr/bin/env bash
# 100K PBMC ATAC CBQ parity gate.
#
# Compares barcoded paired-end ATAC fragment output produced from native CBQ
# input against the FASTQ baseline on the 100K PBMC fixture:
#
#   - FASTQ baseline through `chromap`,
#   - native CBQ through `chromap`,
#   - native CBQ through `chromap_lib_runner`,
#   - canonical comparison: LC_ALL=C sorted fragment (BED) rows.
#
# Encoder choice (important):
#   The barcoded ATAC CBQ path consumes a paired-read CBQ lane and a separate
#   barcode CBQ lane, and requires record i of the read lane to correspond to
#   record i of the barcode lane (the same original read). chromap's CBQ reader
#   enforces this with a read/barcode name-match guard.
#
#   Stock `bqtools encode` does NOT preserve input record order at multi-block
#   scale: the paired-read lane and the barcode lane are independently
#   reordered, so they no longer align and chromap rejects them. (The 4-record
#   synthetic smoke happens to fit one block and survives; 100K does not.)
#   This gate therefore uses the order-preserving `cbq_ordered_encoder`, which
#   emits records in input FASTQ order so the two lanes stay aligned. The zstd
#   decompression path is covered separately by the synthetic smoke, which
#   feeds chromap a compressed bqtools CBQ.
#
# No intermediate FASTQ is materialized in the production CBQ mapping path; the
# only FASTQ-to-CBQ conversion is the test-only encode step below. The encoder
# is an optional, test-only external dependency; released Chromap-suite does
# not depend on STAR-suite or any STAR adapter.
#
# Benchmarks run serially. The manifest records command lines, git state,
# fixture/index/reference paths, thread count, output directory, and
# wall/user/sys/max-RSS for each chromap invocation (GNU /usr/bin/time when
# available).
#
# Prerequisites: chromap + chromap_lib_runner built at repo root, the
# cbq_ordered_encoder (set CBQ_ORDERED_ENCODER=/path or place on PATH), and the
# 100K fixture, 12 GB index, GRCh38 reference, and 10x ATAC whitelist. A missing
# encoder or fixture input is reported as a SKIP, not a failure.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

CHROMAP="${CHROMAP:-${REPO_ROOT}/chromap}"
CHROMAP_LIB_RUNNER="${CHROMAP_LIB_RUNNER:-${REPO_ROOT}/chromap_lib_runner}"

FIXTURE_ATAC="${FIXTURE_ATAC:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac}"
INDEX="${INDEX:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index}"
REF="${REF:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa}"
WHITELIST="${WHITELIST:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt}"
# Lane numbers to include (space separated). Default exercises the multi-lane
# CBQ reader loop. Override e.g. LANES="1" for a faster single-lane gate.
LANES="${LANES:-1 2 3 4}"
THREADS="${THREADS:-8}"
# Zstd compression level for the ordered encoder (0 = store).
CBQ_COMPRESSION_LEVEL="${CBQ_COMPRESSION_LEVEL:-0}"

ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
OUT_ROOT="${OUT_ROOT:-${ARTIFACT_ROOT}/cbq_atac_100k/$(date -u +%Y%m%dT%H%M%SZ)}"
CBQ_DIR="${OUT_ROOT}/cbq"
RUN_DIR="${OUT_ROOT}/runs"
mkdir -p "${CBQ_DIR}" "${RUN_DIR}"

MANIFEST="${OUT_ROOT}/MANIFEST.txt"

log() {
  printf '[cbq-atac-100k] %s\n' "$*"
}

skip() {
  log "SKIP: $*"
  printf 'status=skipped\nreason=%s\n' "$*" > "${OUT_ROOT}/SKIPPED.txt"
  exit 0
}

resolve_encoder() {
  if [[ -n "${CBQ_ORDERED_ENCODER:-}" ]]; then
    if [[ -x "${CBQ_ORDERED_ENCODER}" ]]; then
      printf '%s\n' "${CBQ_ORDERED_ENCODER}"
      return 0
    fi
    if command -v "${CBQ_ORDERED_ENCODER}" >/dev/null 2>&1; then
      command -v "${CBQ_ORDERED_ENCODER}"
      return 0
    fi
    return 1
  fi
  if [[ -x /mnt/pikachu/STAR-suite/core/legacy/source/cbq_ordered_encoder ]]; then
    printf '%s\n' /mnt/pikachu/STAR-suite/core/legacy/source/cbq_ordered_encoder
    return 0
  fi
  command -v cbq_ordered_encoder 2>/dev/null
}

encode_pair_cbq() {
  local enc="$1" r1="$2" r2="$3" out="$4"
  rm -f "${out}"
  "${enc}" --readFilesIn "${r1}" "${r2}" --outFile "${out}" \
      -l "${CBQ_COMPRESSION_LEVEL}" \
      > "${out}.encode.stdout" 2> "${out}.encode.stderr"
  [[ -s "${out}" ]]
}

encode_single_cbq() {
  local enc="$1" r1="$2" out="$3"
  rm -f "${out}"
  "${enc}" --readFilesIn "${r1}" --outFile "${out}" \
      -l "${CBQ_COMPRESSION_LEVEL}" \
      > "${out}.encode.stdout" 2> "${out}.encode.stderr"
  [[ -s "${out}" ]]
}

# GNU time wrapper that appends wall/user/sys/max-RSS for a labeled command to
# the manifest. Falls back to plain execution when /usr/bin/time is absent.
TIME_BIN=""
if [[ -x /usr/bin/time ]]; then
  TIME_BIN=/usr/bin/time
fi

timed_run() {
  local label="$1"
  shift
  local time_log="${RUN_DIR}/${label}.time.txt"
  log "${label}: $*"
  {
    echo "==== ${label} ===="
    printf '  cmd:'
    printf ' %q' "$@"
    printf '\n'
  } >> "${MANIFEST}"
  if [[ -n "${TIME_BIN}" ]]; then
    "${TIME_BIN}" -v -o "${time_log}" "$@"
    {
      echo "  wall(elapsed): $(grep -m1 'wall clock' "${time_log}" | sed 's/.*: //')"
      echo "  user(s): $(grep -m1 'User time' "${time_log}" | sed 's/.*: //')"
      echo "  sys(s): $(grep -m1 'System time' "${time_log}" | sed 's/.*: //')"
      echo "  max_rss(kb): $(grep -m1 'Maximum resident' "${time_log}" | sed 's/.*: //')"
    } >> "${MANIFEST}"
  else
    "$@"
  fi
}

normalize_fragments() {
  LC_ALL=C sort "$1" > "$2"
}

assert_fragment_parity() {
  local expected="$1" observed="$2" label="$3"
  local expected_norm="${expected}.norm" observed_norm="${observed}.norm"
  normalize_fragments "${expected}" "${expected_norm}"
  normalize_fragments "${observed}" "${observed_norm}"
  if ! cmp -s "${expected_norm}" "${observed_norm}"; then
    echo "ERROR: fragment parity failed for ${label}" >&2
    diff -u "${expected_norm}" "${observed_norm}" | head -80 >&2 || true
    exit 1
  fi
  log "PASS: ${label} ($(wc -l < "${expected_norm}") rows)"
  echo "  parity ${label}: PASS ($(wc -l < "${expected_norm}") sorted rows)" >> "${MANIFEST}"
}

[[ -x "${CHROMAP}" ]] || { echo "ERROR: chromap not built (${CHROMAP})" >&2; exit 2; }
[[ -x "${CHROMAP_LIB_RUNNER}" ]] || { echo "ERROR: chromap_lib_runner not built (${CHROMAP_LIB_RUNNER})" >&2; exit 2; }

encoder_bin="$(resolve_encoder)" || skip "cbq_ordered_encoder not found; set CBQ_ORDERED_ENCODER=/path/to/cbq_ordered_encoder"
[[ -d "${FIXTURE_ATAC}" ]] || skip "fixture not found: ${FIXTURE_ATAC}"
[[ -f "${INDEX}" ]] || skip "index not found: ${INDEX}"
[[ -f "${REF}" ]] || skip "reference not found: ${REF}"
[[ -f "${WHITELIST}" ]] || skip "whitelist not found: ${WHITELIST}"

# Assemble per-lane input lists (FASTQ: R1=read1, R3=read2, R2=barcode).
declare -a R1_LIST R3_LIST BC_LIST PAIR_CBQ_LIST BC_CBQ_LIST
for lane in ${LANES}; do
  printf -v lane3 'L%03d' "${lane}"
  r1="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_${lane3}_R1_001.fastq.gz"
  r3="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_${lane3}_R3_001.fastq.gz"
  bc="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_${lane3}_R2_001.fastq.gz"
  for f in "${r1}" "${r3}" "${bc}"; do
    [[ -f "${f}" ]] || skip "missing fixture lane file: ${f}"
  done
  R1_LIST+=("${r1}")
  R3_LIST+=("${r3}")
  BC_LIST+=("${bc}")
  PAIR_CBQ_LIST+=("${CBQ_DIR}/${lane3}.reads_pair.cbq")
  BC_CBQ_LIST+=("${CBQ_DIR}/${lane3}.barcodes.cbq")
done

join_csv() {
  local IFS=','
  printf '%s' "$*"
}

R1_CSV="$(join_csv "${R1_LIST[@]}")"
R3_CSV="$(join_csv "${R3_LIST[@]}")"
BC_CSV="$(join_csv "${BC_LIST[@]}")"
PAIR_CBQ_CSV="$(join_csv "${PAIR_CBQ_LIST[@]}")"
BC_CBQ_CSV="$(join_csv "${BC_CBQ_LIST[@]}")"

GIT_STATE="$(git -C "${REPO_ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
# --porcelain reports staged, unstaged, AND untracked changes, so an untracked
# input or script that affected the run still marks the SHA dirty.
if [[ -n "$(git -C "${REPO_ROOT}" status --porcelain 2>/dev/null)" ]]; then
  GIT_STATE="${GIT_STATE}-dirty"
fi

{
  echo "cbq_atac_100k manifest"
  echo "out_root ${OUT_ROOT}"
  echo "git_state ${GIT_STATE}"
  echo "chromap ${CHROMAP}"
  echo "chromap_lib_runner ${CHROMAP_LIB_RUNNER}"
  echo "cbq_ordered_encoder ${encoder_bin}"
  echo "cbq_compression_level ${CBQ_COMPRESSION_LEVEL}"
  echo "fixture_atac ${FIXTURE_ATAC}"
  echo "index ${INDEX}"
  echo "reference ${REF}"
  echo "whitelist ${WHITELIST}"
  echo "lanes ${LANES}"
  echo "threads ${THREADS}"
  echo "canonicalization LC_ALL=C sort over fragment BED rows"
  echo
} > "${MANIFEST}"

log "outputs: ${OUT_ROOT}"
log "encoder: ${encoder_bin}"
log "lanes: ${LANES}  threads: ${THREADS}  git: ${GIT_STATE}"

# Encode CBQ inputs (test-only, order-preserving FASTQ -> CBQ conversion).
for i in "${!PAIR_CBQ_LIST[@]}"; do
  log "encode lane $((i + 1))/${#PAIR_CBQ_LIST[@]} -> CBQ"
  encode_pair_cbq "${encoder_bin}" "${R1_LIST[$i]}" "${R3_LIST[$i]}" "${PAIR_CBQ_LIST[$i]}"
  encode_single_cbq "${encoder_bin}" "${BC_LIST[$i]}" "${BC_CBQ_LIST[$i]}"
done

# Common mapping options. Identical across FASTQ and CBQ runs so the only
# difference is the input loader.
common_args=(
  -t "${THREADS}"
  -x "${INDEX}"
  -r "${REF}"
  --barcode-whitelist "${WHITELIST}"
  -l 2000
  --trim-adapters
  --remove-pcr-duplicates
  --remove-pcr-duplicates-at-cell-level
  --Tn5-shift
  --BED
)

fastq_bed="${RUN_DIR}/fastq.fragments.bed"
cbq_bed="${RUN_DIR}/cbq.fragments.bed"
cbq_lib_bed="${RUN_DIR}/cbq.lib.fragments.bed"

timed_run fastq_baseline \
  "${CHROMAP}" "${common_args[@]}" \
  -1 "${R1_CSV}" -2 "${R3_CSV}" -b "${BC_CSV}" \
  --summary "${RUN_DIR}/fastq.summary.tsv" \
  -o "${fastq_bed}"

timed_run cbq_cli \
  "${CHROMAP}" "${common_args[@]}" \
  --input-format cbq \
  --read-pair-cbq "${PAIR_CBQ_CSV}" \
  --barcode-cbq "${BC_CBQ_CSV}" \
  --summary "${RUN_DIR}/cbq.summary.tsv" \
  -o "${cbq_bed}"

timed_run cbq_lib_runner \
  "${CHROMAP_LIB_RUNNER}" "${common_args[@]}" \
  --input-format cbq \
  --read-pair-cbq "${PAIR_CBQ_CSV}" \
  --barcode-cbq "${BC_CBQ_CSV}" \
  --summary "${RUN_DIR}/cbq.lib.summary.tsv" \
  -o "${cbq_lib_bed}"

for f in "${fastq_bed}" "${cbq_bed}" "${cbq_lib_bed}"; do
  [[ -s "${f}" ]] || { echo "ERROR: empty fragment output: ${f}" >&2; exit 1; }
done

assert_fragment_parity "${fastq_bed}" "${cbq_bed}" "CLI FASTQ vs CLI CBQ"
assert_fragment_parity "${fastq_bed}" "${cbq_lib_bed}" "CLI FASTQ vs libchromap CBQ"

{
  echo
  echo "status=pass"
  echo "fastq_rows=$(wc -l < "${fastq_bed}")"
} >> "${MANIFEST}"

log "PASS: CBQ ATAC 100K parity completed at ${OUT_ROOT}"
