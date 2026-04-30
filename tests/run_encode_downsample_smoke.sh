#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
CACHE_ROOT="${ENCODE_FIXTURE_CACHE:-${ARTIFACT_ROOT}/encode_fixture_cache}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
OUTROOT="${OUTROOT:-${ARTIFACT_ROOT}/encode_downsample_smoke/${RUN_ID}}"
GENERATED_MANIFEST="${ENCODE_GENERATED_MANIFEST:-${CACHE_ROOT}/manifest.generated.tsv}"
CHROMAP_BIN="${CHROMAP_BIN:-${REPO_ROOT}/chromap}"
LIBRUNNER_BIN="${LIBRUNNER_BIN:-${REPO_ROOT}/chromap_lib_runner}"
THREADS="${THREADS:-1}"
BUILD="${BUILD:-1}"

log() {
  printf '[encode-smoke] %s\n' "$*" >&2
}

fail() {
  printf '[encode-smoke] FAIL: %s\n' "$*" >&2
  exit 1
}

require_file() {
  [[ -s "$1" ]] || fail "missing or empty file: $1"
}

normalize_text_rows() {
  local in="$1"
  LC_ALL=C grep -v '^#' "${in}" | sed '/^[[:space:]]*$/d' | LC_ALL=C sort || true
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
  local id="$1"; shift
  log "${id}: $*"
  "$@" >"${OUTROOT}/logs/${id}.stdout" 2>"${OUTROOT}/logs/${id}.stderr"
}

run_case() {
  local assay="$1"
  local case_id="$2"
  local r1="$3"
  local r2="$4"
  local chromap_args="$5"
  local case_dir="${OUTROOT}/cases/${case_id}_${assay}"
  mkdir -p "${case_dir}"

  local -a args
  read -r -a args <<<"${chromap_args}"
  local common=(-x "${CHROMAP_GRCH38_INDEX}" -r "${CHROMAP_GRCH38_REF}" \
    -1 "${r1}" -2 "${r2}" -t "${THREADS}")

  run_cmd "${case_id}.cli" "${CHROMAP_BIN}" "${common[@]}" "${args[@]}" \
    -o "${case_dir}/cli.out"
  run_cmd "${case_id}.lib" "${LIBRUNNER_BIN}" "${common[@]}" "${args[@]}" \
    -o "${case_dir}/lib.out"
  require_file "${case_dir}/cli.out"
  require_file "${case_dir}/lib.out"
  assert_same_rows "${case_id}_${assay}" "${case_dir}/cli.out" "${case_dir}/lib.out"
  printf '%s\tPASS\t%s\t%s\t%s\n' "${case_id}" "${assay}" \
    "$(wc -l <"${OUTROOT}/compare/${case_id}_${assay}/cli.rows")" \
    "${chromap_args}" >>"${OUTROOT}/summary.tsv"
}

main() {
  : "${CHROMAP_GRCH38_REF:?set CHROMAP_GRCH38_REF to a GRCh38 FASTA}"
  : "${CHROMAP_GRCH38_INDEX:?set CHROMAP_GRCH38_INDEX to a matching Chromap index}"
  require_file "${CHROMAP_GRCH38_REF}"
  require_file "${CHROMAP_GRCH38_INDEX}"

  mkdir -p "${OUTROOT}"/{cases,compare,logs}
  printf 'case\tstatus\tassay\tnon_header_rows\tchromap_args\n' >"${OUTROOT}/summary.tsv"

  if [[ "${BUILD}" == "1" ]]; then
    log "building chromap and chromap_lib_runner"
    (cd "${REPO_ROOT}" && make -j"${MAKE_JOBS:-4}" chromap chromap_lib_runner) \
      >"${OUTROOT}/logs/make.stdout" 2>"${OUTROOT}/logs/make.stderr"
  fi

  if [[ "${ENCODE_SKIP_PREPARE:-0}" != "1" ]]; then
    ENCODE_GENERATED_MANIFEST="${GENERATED_MANIFEST}" \
      "${REPO_ROOT}/tests/prepare_encode_downsample_fixtures.sh" \
      >"${OUTROOT}/logs/prepare.stdout" \
      2>"${OUTROOT}/logs/prepare.stderr"
  fi
  require_file "${GENERATED_MANIFEST}"

  local line
  while IFS=$'\t' read -r assay case_id experiment r1_acc r2_acc read_count r1_fastq r2_fastq r1_md5 r2_md5 chromap_args description; do
    [[ "${assay}" == "assay" ]] && continue
    [[ -n "${assay}" ]] || continue
    require_file "${r1_fastq}"
    require_file "${r2_fastq}"
    run_case "${assay}" "${case_id}" "${r1_fastq}" "${r2_fastq}" "${chromap_args}"
  done <"${GENERATED_MANIFEST}"

  log "PASS summary: ${OUTROOT}/summary.tsv"
  cat "${OUTROOT}/summary.tsv"
}

main "$@"
