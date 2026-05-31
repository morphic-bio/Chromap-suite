#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

HARNESS="${CBQ_RANGE_READER_HARNESS:-${REPO_ROOT}/tests/cbq_range_reader_harness}"
ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
OUT_ROOT="${OUT_ROOT:-${ARTIFACT_ROOT}/cbq_range_reader_smoke/$(date -u +%Y%m%dT%H%M%SZ)}"
THREADS="${THREADS:-4}"
BATCH_SIZE="${CBQ_RANGE_BATCH_SIZE:-131072}"
EXISTING_PAIR_CBQ="${EXISTING_PAIR_CBQ:-${REPO_ROOT}/plans/artifacts/cbq_producer_downsample_2m_20260530T160612Z/cbq/L001.reads_pair.cbq}"

mkdir -p "${OUT_ROOT}/logs" "${OUT_ROOT}/synthetic"

log() {
  printf '[cbq-range-reader-smoke] %s\n' "$*"
}

skip() {
  log "SKIP: $*"
  printf 'status=skipped\nreason=%s\n' "$*" > "${OUT_ROOT}/SKIPPED.txt"
  exit 0
}

resolve_bqtools() {
  if [[ -n "${BQTOOLS:-}" ]]; then
    if [[ -x "${BQTOOLS}" ]]; then
      printf '%s\n' "${BQTOOLS}"
      return 0
    fi
    if command -v "${BQTOOLS}" >/dev/null 2>&1; then
      command -v "${BQTOOLS}"
      return 0
    fi
    return 1
  fi
  if [[ -x /tmp/star_suite_bqtools/bin/bqtools ]]; then
    printf '%s\n' /tmp/star_suite_bqtools/bin/bqtools
    return 0
  fi
  command -v bqtools 2>/dev/null
}

encode_pair_cbq() {
  local bqtools="$1"
  local r1="$2"
  local r2="$3"
  local out="$4"
  local level="${5:-}"
  local -a level_args=()
  if [[ -n "${level}" ]]; then
    level_args=(-l "${level}")
  fi

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" "${r2}" --mode cbq "${level_args[@]}" \
      -o "${out}" -T 2 > "${out}.encode.stdout" 2> "${out}.encode.stderr"; then
    [[ -s "${out}" ]] && return 0
  fi

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" "${r2}" -o "${out}" --mode cbq \
      "${level_args[@]}" -T 2 >> "${out}.encode.stdout" \
      2>> "${out}.encode.stderr"; then
    [[ -s "${out}" ]] && return 0
  fi

  rm -f "${out}"
  "${bqtools}" encode "${r1}" "${r2}" -o "${out}" "${level_args[@]}" -T 2 \
      >> "${out}.encode.stdout" 2>> "${out}.encode.stderr"
  [[ -s "${out}" ]]
}

encode_single_cbq() {
  local bqtools="$1"
  local r1="$2"
  local out="$3"
  local level="${4:-}"
  local -a level_args=()
  if [[ -n "${level}" ]]; then
    level_args=(-l "${level}")
  fi

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" --mode cbq "${level_args[@]}" \
      -o "${out}" -T 2 > "${out}.encode.stdout" 2> "${out}.encode.stderr"; then
    [[ -s "${out}" ]] && return 0
  fi

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" -o "${out}" --mode cbq "${level_args[@]}" \
      -T 2 >> "${out}.encode.stdout" 2>> "${out}.encode.stderr"; then
    [[ -s "${out}" ]] && return 0
  fi

  rm -f "${out}"
  "${bqtools}" encode "${r1}" -o "${out}" "${level_args[@]}" -T 2 \
      >> "${out}.encode.stdout" 2>> "${out}.encode.stderr"
  [[ -s "${out}" ]]
}

run_harness() {
  local label="$1"
  local cbq="$2"
  local mate_count="$3"

  log "run ${label}: ${cbq}"
  "${HARNESS}" --cbq "${cbq}" --mate-count "${mate_count}" \
    --threads "${THREADS}" --batch-size "${BATCH_SIZE}" \
    --materialize none > "${OUT_ROOT}/logs/${label}.none.stdout" \
    2> "${OUT_ROOT}/logs/${label}.none.stderr"
  "${HARNESS}" --cbq "${cbq}" --mate-count "${mate_count}" \
    --threads "${THREADS}" --batch-size "${BATCH_SIZE}" \
    --materialize sequence > "${OUT_ROOT}/logs/${label}.sequence.stdout" \
    2> "${OUT_ROOT}/logs/${label}.sequence.stderr"
}

[[ -x "${HARNESS}" ]] || {
  echo "ERROR: range reader harness not found: ${HARNESS}" >&2
  exit 1
}

log "outputs: ${OUT_ROOT}"
log "threads: ${THREADS}; batch_size: ${BATCH_SIZE}"

ran_any=0
if [[ -s "${EXISTING_PAIR_CBQ}" ]]; then
  run_harness "existing_pair" "${EXISTING_PAIR_CBQ}" 2
  ran_any=1
fi

if bqtools_bin="$(resolve_bqtools)"; then
  log "bqtools: ${bqtools_bin}"
  python3 - "${OUT_ROOT}/synthetic" <<'PY'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
rng = random.Random(53131)
alphabet = "ACGT"
records = 20000
read_len = 72
barcode_len = 16

with (out / "R1.fastq").open("wt", encoding="ascii") as r1, \
     (out / "R2.fastq").open("wt", encoding="ascii") as r2, \
     (out / "BC.fastq").open("wt", encoding="ascii") as bc:
    for i in range(records):
        name = f"range{i:07d}"
        s1 = "".join(rng.choice(alphabet) for _ in range(read_len))
        s2 = "".join(rng.choice(alphabet) for _ in range(read_len))
        sb = "".join(rng.choice(alphabet) for _ in range(barcode_len))
        q1 = "I" * read_len
        q2 = "H" * read_len
        qb = "J" * barcode_len
        r1.write(f"@{name}/1 1:N:0:ACGT\n{s1}\n+\n{q1}\n")
        r2.write(f"@{name}/2 2:N:0:ACGT\n{s2}\n+\n{q2}\n")
        bc.write(f"@{name}/3 3:N:0:ACGT\n{sb}\n+\n{qb}\n")
PY
  encode_pair_cbq "${bqtools_bin}" "${OUT_ROOT}/synthetic/R1.fastq" \
    "${OUT_ROOT}/synthetic/R2.fastq" "${OUT_ROOT}/synthetic/pair.cbq"
  encode_pair_cbq "${bqtools_bin}" "${OUT_ROOT}/synthetic/R1.fastq" \
    "${OUT_ROOT}/synthetic/R2.fastq" "${OUT_ROOT}/synthetic/pair_level0.cbq" 0
  encode_single_cbq "${bqtools_bin}" "${OUT_ROOT}/synthetic/BC.fastq" \
    "${OUT_ROOT}/synthetic/barcode.cbq"
  encode_single_cbq "${bqtools_bin}" "${OUT_ROOT}/synthetic/BC.fastq" \
    "${OUT_ROOT}/synthetic/barcode_level0.cbq" 0

  run_harness "synthetic_pair" "${OUT_ROOT}/synthetic/pair.cbq" 2
  run_harness "synthetic_pair_level0" "${OUT_ROOT}/synthetic/pair_level0.cbq" 2
  run_harness "synthetic_barcode" "${OUT_ROOT}/synthetic/barcode.cbq" 1
  run_harness "synthetic_barcode_level0" \
    "${OUT_ROOT}/synthetic/barcode_level0.cbq" 1
  ran_any=1
else
  log "bqtools not found; synthetic compressed/level0 cases skipped"
fi

if [[ "${ran_any}" -eq 0 ]]; then
  skip "no existing CBQ fixture and bqtools not available"
fi

{
  echo "cbq_range_reader_smoke manifest"
  echo "out_root ${OUT_ROOT}"
  echo "threads ${THREADS}"
  echo "batch_size ${BATCH_SIZE}"
  echo "existing_pair_cbq ${EXISTING_PAIR_CBQ}"
} > "${OUT_ROOT}/MANIFEST.txt"

log "PASS: CBQ range reader smoke completed at ${OUT_ROOT}"
