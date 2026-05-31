#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
CACHE_ROOT="${ENCODE_CROSS_ASSAY_CACHE:-${ARTIFACT_ROOT}/encode_cross_assay_cache}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
OUTROOT="${OUTROOT:-${ARTIFACT_ROOT}/encode_cbq_cross_assay_smoke/${RUN_ID}}"
GENERATED_MANIFEST="${ENCODE_CROSS_ASSAY_GENERATED_MANIFEST:-${CACHE_ROOT}/manifest.generated.tsv}"
CHROMAP_BIN="${CHROMAP_BIN:-${REPO_ROOT}/chromap}"
LIBRUNNER_BIN="${LIBRUNNER_BIN:-${REPO_ROOT}/chromap_lib_runner}"
THREADS="${THREADS:-1}"
BUILD="${BUILD:-1}"
ASSAYS="${ENCODE_CBQ_ASSAYS:-${ENCODE_ASSAYS:-all}}"
CBQ_COMPRESSION_LEVEL="${ENCODE_CBQ_COMPRESSION_LEVEL:-3}"

log() {
  printf '[encode-cbq-cross-smoke] %s\n' "$*" >&2
}

fail() {
  printf '[encode-cbq-cross-smoke] FAIL: %s\n' "$*" >&2
  exit 1
}

require_file() {
  [[ -s "$1" ]] || fail "missing or empty file: $1"
}

want_assay() {
  local assay="$1"
  [[ "${ASSAYS}" == "all" ]] && return 0
  [[ ",${ASSAYS}," == *",${assay},"* ]]
}

parse_tsv_line() {
  python3 - "$1" <<'PY'
import csv
import sys

row = next(csv.reader([sys.argv[1]], delimiter="\t"))
for field in row:
    print(field)
PY
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

git_state() {
  local state
  state="$(git -C "${REPO_ROOT}" rev-parse --short HEAD 2>/dev/null || echo unknown)"
  if [[ -n "$(git -C "${REPO_ROOT}" status --porcelain 2>/dev/null)" ]]; then
    state="${state}-dirty"
  fi
  printf '%s\n' "${state}"
}

normalize_text_rows() {
  local in="$1"
  awk 'length($0) && substr($0, 1, 1) != "#" {print}' "${in}" | LC_ALL=C sort
}

normalize_summary_rows() {
  local in="$1"
  tail -n +2 "${in}" | LC_ALL=C sort
}

assert_same_rows() {
  local id="$1"
  local a="$2"
  local b="$3"
  local label="$4"
  local work="${OUTROOT}/compare/${id}"
  mkdir -p "${work}"
  normalize_text_rows "${a}" >"${work}/expected.rows"
  normalize_text_rows "${b}" >"${work}/observed.rows"
  cmp -s "${work}/expected.rows" "${work}/observed.rows" || {
    diff -u "${work}/expected.rows" "${work}/observed.rows" | sed -n '1,80p' >&2 || true
    fail "${label}: rows differ"
  }
  [[ -s "${work}/expected.rows" ]] || fail "${label}: no non-header rows emitted"
}

assert_same_summary() {
  local id="$1"
  local a="$2"
  local b="$3"
  local label="$4"
  local work="${OUTROOT}/compare/${id}"
  mkdir -p "${work}"
  normalize_summary_rows "${a}" >"${work}/expected.summary.rows"
  normalize_summary_rows "${b}" >"${work}/observed.summary.rows"
  cmp -s "${work}/expected.summary.rows" "${work}/observed.summary.rows" || {
    diff -u "${work}/expected.summary.rows" "${work}/observed.summary.rows" | sed -n '1,80p' >&2 || true
    fail "${label}: summary rows differ"
  }
  [[ -s "${work}/expected.summary.rows" ]] || fail "${label}: no summary rows emitted"
}

run_cmd() {
  local id="$1"
  shift
  log "${id}: $*"
  {
    echo "==== ${id} ===="
    printf 'cmd:'
    printf ' %q' "$@"
    printf '\n'
  } >>"${OUTROOT}/RUN_MANIFEST.txt"
  /usr/bin/time -v -o "${OUTROOT}/logs/${id}.time.txt" \
    "$@" >"${OUTROOT}/logs/${id}.stdout" 2>"${OUTROOT}/logs/${id}.stderr"
}

encode_pair_cbq() {
  local id="$1"
  local r1="$2"
  local r2="$3"
  local out="$4"
  if [[ -s "${out}" ]]; then
    return
  fi
  run_cmd "${id}" "${ENCODER_BIN}" --readFilesIn "${r1}" "${r2}" \
    --outFile "${out}" -l "${CBQ_COMPRESSION_LEVEL}"
  require_file "${out}"
}

encode_single_cbq() {
  local id="$1"
  local r1="$2"
  local out="$3"
  if [[ -s "${out}" ]]; then
    return
  fi
  run_cmd "${id}" "${ENCODER_BIN}" --readFilesIn "${r1}" \
    --outFile "${out}" -l "${CBQ_COMPRESSION_LEVEL}"
  require_file "${out}"
}

run_case() {
  local assay="$1"
  local case_id="$2"
  local experiment="$3"
  local layout="$4"
  local r1="$5"
  local r2="$6"
  local barcode="$7"
  local whitelist="$8"
  local chromap_args="$9"
  local case_dir="${OUTROOT}/cases/${case_id}_${assay}"
  local cbq_dir="${OUTROOT}/cbq/${case_id}_${assay}"
  mkdir -p "${case_dir}" "${cbq_dir}"

  local -a args
  read -r -a args <<<"${chromap_args}"

  local read_pair_cbq="${cbq_dir}/reads.cbq"
  local barcode_cbq="${cbq_dir}/barcode.cbq"
  encode_pair_cbq "${case_id}.encode.reads" "${r1}" "${r2}" "${read_pair_cbq}"

  local -a fastq_common cbq_common lib_cbq_common
  fastq_common=(-x "${CHROMAP_GRCH38_INDEX}" -r "${CHROMAP_GRCH38_REF}" -t "${THREADS}")
  cbq_common=(-x "${CHROMAP_GRCH38_INDEX}" -r "${CHROMAP_GRCH38_REF}" -t "${THREADS}" \
    --input-format cbq --read-pair-cbq "${read_pair_cbq}")
  lib_cbq_common=("${cbq_common[@]}")

  if [[ "${layout}" == "paired" ]]; then
    fastq_common+=(-1 "${r1}" -2 "${r2}")
  elif [[ "${layout}" == "scatac_10x_atac" ]]; then
    require_file "${barcode}"
    encode_single_cbq "${case_id}.encode.barcode" "${barcode}" "${barcode_cbq}"
    fastq_common+=(-1 "${r1}" -2 "${r2}" -b "${barcode}" \
      --summary "${case_dir}/fastq.summary.tsv")
    cbq_common+=(--barcode-cbq "${barcode_cbq}" --summary "${case_dir}/cbq.summary.tsv")
    lib_cbq_common+=(--barcode-cbq "${barcode_cbq}" --summary "${case_dir}/cbq.lib.summary.tsv")
    if [[ -n "${whitelist}" ]]; then
      require_file "${whitelist}"
      fastq_common+=(--barcode-whitelist "${whitelist}")
      cbq_common+=(--barcode-whitelist "${whitelist}")
      lib_cbq_common+=(--barcode-whitelist "${whitelist}")
    fi
  else
    fail "${case_id}: unsupported layout ${layout}"
  fi

  run_cmd "${case_id}.fastq" "${CHROMAP_BIN}" "${fastq_common[@]}" "${args[@]}" \
    -o "${case_dir}/fastq.out"
  run_cmd "${case_id}.cbq.cli" env CHROMAP_REQUIRE_CBQ_INDEX=1 \
    "${CHROMAP_BIN}" "${cbq_common[@]}" "${args[@]}" -o "${case_dir}/cbq.out"
  run_cmd "${case_id}.cbq.lib" env CHROMAP_REQUIRE_CBQ_INDEX=1 \
    "${LIBRUNNER_BIN}" "${lib_cbq_common[@]}" "${args[@]}" \
    -o "${case_dir}/cbq.lib.out"

  require_file "${case_dir}/fastq.out"
  require_file "${case_dir}/cbq.out"
  require_file "${case_dir}/cbq.lib.out"
  for log_file in "${OUTROOT}/logs/${case_id}.cbq.cli.stderr" \
                  "${OUTROOT}/logs/${case_id}.cbq.lib.stderr"; do
    if ! grep -q "Using indexed CBQ range producer" "${log_file}"; then
      fail "${case_id}: indexed CBQ range producer was not reported in ${log_file}"
    fi
  done
  assert_same_rows "${case_id}_${assay}_cli" "${case_dir}/fastq.out" \
    "${case_dir}/cbq.out" "${case_id}: FASTQ vs CBQ CLI"
  assert_same_rows "${case_id}_${assay}_lib" "${case_dir}/fastq.out" \
    "${case_dir}/cbq.lib.out" "${case_id}: FASTQ vs CBQ lib"

  if [[ "${layout}" == "scatac_10x_atac" ]]; then
    assert_same_summary "${case_id}_${assay}_summary_cli" \
      "${case_dir}/fastq.summary.tsv" "${case_dir}/cbq.summary.tsv" \
      "${case_id}: FASTQ vs CBQ CLI summary"
    assert_same_summary "${case_id}_${assay}_summary_lib" \
      "${case_dir}/fastq.summary.tsv" "${case_dir}/cbq.lib.summary.tsv" \
      "${case_id}: FASTQ vs CBQ lib summary"
  fi

  local row_count
  row_count="$(wc -l <"${OUTROOT}/compare/${case_id}_${assay}_cli/expected.rows")"
  printf '%s\tPASS\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${case_id}" "${assay}" "${layout}" "${experiment}" "${row_count}" \
    "${chromap_args}" "${OUTROOT}/compare/${case_id}_${assay}_cli/expected.rows" \
    >>"${OUTROOT}/summary.tsv"
}

main() {
  : "${CHROMAP_GRCH38_REF:?set CHROMAP_GRCH38_REF to a GRCh38 FASTA}"
  : "${CHROMAP_GRCH38_INDEX:?set CHROMAP_GRCH38_INDEX to a matching Chromap index}"
  require_file "${CHROMAP_GRCH38_REF}"
  require_file "${CHROMAP_GRCH38_INDEX}"

  ENCODER_BIN="$(resolve_encoder)" || fail "cbq_ordered_encoder not found; set CBQ_ORDERED_ENCODER"
  export ENCODER_BIN

  mkdir -p "${OUTROOT}"/{cases,cbq,compare,logs}
  {
    echo "encode_cbq_cross_assay_smoke manifest"
    echo "outroot ${OUTROOT}"
    echo "git_state $(git_state)"
    echo "generated_manifest ${GENERATED_MANIFEST}"
    echo "reference ${CHROMAP_GRCH38_REF}"
    echo "index ${CHROMAP_GRCH38_INDEX}"
    echo "encoder ${ENCODER_BIN}"
    echo "cbq_compression_level ${CBQ_COMPRESSION_LEVEL}"
    echo "threads ${THREADS}"
    echo "assays ${ASSAYS}"
    echo "canonicalization non-comment non-empty rows + LC_ALL=C sort"
    echo
  } >"${OUTROOT}/RUN_MANIFEST.txt"
  printf 'case\tstatus\tassay\tlayout\texperiment\tnon_header_rows\tchromap_args\tcanonical_rows\n' \
    >"${OUTROOT}/summary.tsv"

  if [[ "${BUILD}" == "1" ]]; then
    log "building chromap and chromap_lib_runner"
    (cd "${REPO_ROOT}" && make -j"${MAKE_JOBS:-4}" chromap chromap_lib_runner) \
      >"${OUTROOT}/logs/make.stdout" 2>"${OUTROOT}/logs/make.stderr"
  fi
  require_file "${CHROMAP_BIN}"
  require_file "${LIBRUNNER_BIN}"

  if [[ "${ENCODE_SKIP_PREPARE:-0}" != "1" ]]; then
    ENCODE_CROSS_ASSAY_GENERATED_MANIFEST="${GENERATED_MANIFEST}" \
      "${REPO_ROOT}/tests/prepare_encode_cross_assay_fixtures.sh" \
      >"${OUTROOT}/logs/prepare.stdout" \
      2>"${OUTROOT}/logs/prepare.stderr"
  fi
  require_file "${GENERATED_MANIFEST}"

  local rows=0
  local line
  while IFS= read -r line || [[ -n "${line}" ]]; do
    [[ -z "${line}" ]] && continue
    mapfile -t fields < <(parse_tsv_line "${line}")
    local assay="${fields[0]:-}"
    [[ "${assay}" == "assay" ]] && continue
    [[ -n "${assay}" ]] || continue
    want_assay "${assay}" || continue
    [[ "${#fields[@]}" -ge 19 ]] || fail "malformed generated manifest row for ${assay}: expected 19 columns"

    local case_id="${fields[1]}"
    local experiment="${fields[2]}"
    local layout="${fields[3]}"
    local r1_fastq="${fields[5]}"
    local r2_fastq="${fields[6]}"
    local barcode_fastq="${fields[7]}"
    local whitelist="${fields[8]}"
    local chromap_args="${fields[13]}"

    require_file "${r1_fastq}"
    require_file "${r2_fastq}"
    rows=$((rows + 1))
    run_case "${assay}" "${case_id}" "${experiment}" "${layout}" \
      "${r1_fastq}" "${r2_fastq}" "${barcode_fastq}" "${whitelist}" \
      "${chromap_args}"
  done <"${GENERATED_MANIFEST}"

  [[ "${rows}" -gt 0 ]] || fail "no assays selected by ENCODE_CBQ_ASSAYS=${ASSAYS}"
  log "PASS summary: ${OUTROOT}/summary.tsv"
  cat "${OUTROOT}/summary.tsv"
}

main "$@"
