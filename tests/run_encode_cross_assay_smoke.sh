#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
CACHE_ROOT="${ENCODE_CROSS_ASSAY_CACHE:-${ARTIFACT_ROOT}/encode_cross_assay_cache}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
OUTROOT="${OUTROOT:-${ARTIFACT_ROOT}/encode_cross_assay_smoke/${RUN_ID}}"
GENERATED_MANIFEST="${ENCODE_CROSS_ASSAY_GENERATED_MANIFEST:-${CACHE_ROOT}/manifest.generated.tsv}"
CHROMAP_BIN="${CHROMAP_BIN:-${REPO_ROOT}/chromap}"
LIBRUNNER_BIN="${LIBRUNNER_BIN:-${REPO_ROOT}/chromap_lib_runner}"
THREADS="${THREADS:-1}"
BUILD="${BUILD:-1}"

log() {
  printf '[encode-cross-smoke] %s\n' "$*" >&2
}

fail() {
  printf '[encode-cross-smoke] FAIL: %s\n' "$*" >&2
  exit 1
}

require_file() {
  [[ -s "$1" ]] || fail "missing or empty file: $1"
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

assert_same_rows() {
  local id="$1"
  local a="$2"
  local b="$3"
  local work="${OUTROOT}/compare/${id}"
  mkdir -p "${work}"
  normalize_text_rows "${a}" >"${work}/cli.rows"
  normalize_text_rows "${b}" >"${work}/lib.rows"
  cmp -s "${work}/cli.rows" "${work}/lib.rows" || {
    diff -u "${work}/cli.rows" "${work}/lib.rows" | sed -n '1,80p' >&2 || true
    fail "${id}: CLI and libchromap rows differ"
  }
  [[ -s "${work}/cli.rows" ]] || fail "${id}: no non-header rows emitted"
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
  "$@" >"${OUTROOT}/logs/${id}.stdout" 2>"${OUTROOT}/logs/${id}.stderr"
}

assert_summary_has_totals() {
  local path="$1"
  require_file "${path}"
  awk -F, 'NR > 1 && $2 + 0 > 0 {found=1} END {exit found ? 0 : 1}' "${path}" || \
    fail "summary has no barcode rows with total > 0: ${path}"
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
  mkdir -p "${case_dir}"

  local -a args
  read -r -a args <<<"${chromap_args}"

  local -a cli_common lib_common
  cli_common=(-x "${CHROMAP_GRCH38_INDEX}" -r "${CHROMAP_GRCH38_REF}" -t "${THREADS}")
  lib_common=(-x "${CHROMAP_GRCH38_INDEX}" -r "${CHROMAP_GRCH38_REF}" -t "${THREADS}")

  if [[ "${layout}" == "paired" ]]; then
    cli_common+=(-1 "${r1}" -2 "${r2}")
    lib_common+=(-1 "${r1}" -2 "${r2}")
  elif [[ "${layout}" == "scatac_10x_atac" ]]; then
    require_file "${barcode}"
    cli_common+=(-1 "${r1}" -2 "${r2}" -b "${barcode}" \
      --summary "${case_dir}/cli.summary.tsv")
    lib_common+=(-1 "${r1}" -2 "${r2}" -b "${barcode}" \
      --summary "${case_dir}/lib.summary.tsv")
    if [[ -n "${whitelist}" ]]; then
      require_file "${whitelist}"
      cli_common+=(--barcode-whitelist "${whitelist}")
      lib_common+=(--barcode-whitelist "${whitelist}")
    fi
  else
    fail "${case_id}: unsupported layout ${layout}"
  fi

  run_cmd "${case_id}.cli" "${CHROMAP_BIN}" "${cli_common[@]}" "${args[@]}" \
    -o "${case_dir}/cli.out"
  run_cmd "${case_id}.lib" "${LIBRUNNER_BIN}" "${lib_common[@]}" "${args[@]}" \
    -o "${case_dir}/lib.out"

  require_file "${case_dir}/cli.out"
  require_file "${case_dir}/lib.out"
  assert_same_rows "${case_id}_${assay}" "${case_dir}/cli.out" "${case_dir}/lib.out"
  if [[ "${layout}" == "scatac_10x_atac" ]]; then
    assert_summary_has_totals "${case_dir}/cli.summary.tsv"
    assert_summary_has_totals "${case_dir}/lib.summary.tsv"
  fi

  local row_count
  row_count="$(wc -l <"${OUTROOT}/compare/${case_id}_${assay}/cli.rows")"
  printf '%s\tPASS\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "${case_id}" "${assay}" "${layout}" "${experiment}" "${row_count}" \
    "${chromap_args}" "${OUTROOT}/compare/${case_id}_${assay}/cli.rows" \
    >>"${OUTROOT}/summary.tsv"
}

main() {
  : "${CHROMAP_GRCH38_REF:?set CHROMAP_GRCH38_REF to a GRCh38 FASTA}"
  : "${CHROMAP_GRCH38_INDEX:?set CHROMAP_GRCH38_INDEX to a matching Chromap index}"
  require_file "${CHROMAP_GRCH38_REF}"
  require_file "${CHROMAP_GRCH38_INDEX}"

  mkdir -p "${OUTROOT}"/{cases,compare,logs}
  {
    echo "encode_cross_assay_smoke manifest"
    echo "outroot ${OUTROOT}"
    echo "git_state $(git_state)"
    echo "source_manifest ${ENCODE_CROSS_ASSAY_SOURCE_MANIFEST:-${REPO_ROOT}/tests/encode_cross_assay_manifest.tsv}"
    echo "generated_manifest ${GENERATED_MANIFEST}"
    echo "reference ${CHROMAP_GRCH38_REF}"
    echo "index ${CHROMAP_GRCH38_INDEX}"
    echo "threads ${THREADS}"
    echo "assays ${ENCODE_ASSAYS:-all}"
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

  local line
  while IFS= read -r line || [[ -n "${line}" ]]; do
    [[ -z "${line}" ]] && continue
    mapfile -t fields < <(parse_tsv_line "${line}")
    local assay="${fields[0]:-}"
    [[ "${assay}" == "assay" ]] && continue
    [[ -n "${assay}" ]] || continue
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
    run_case "${assay}" "${case_id}" "${experiment}" "${layout}" \
      "${r1_fastq}" "${r2_fastq}" "${barcode_fastq}" "${whitelist}" \
      "${chromap_args}"
  done <"${GENERATED_MANIFEST}"

  log "PASS summary: ${OUTROOT}/summary.tsv"
  cat "${OUTROOT}/summary.tsv"
}

main "$@"
