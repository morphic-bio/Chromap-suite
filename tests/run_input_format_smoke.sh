#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

CHROMAP="${CHROMAP:-${REPO_ROOT}/chromap}"
if [[ ! -x "${CHROMAP}" ]]; then
  if command -v chromap >/dev/null 2>&1; then
    CHROMAP="$(command -v chromap)"
  else
    echo "ERROR: chromap binary not found; run make chromap or set CHROMAP=/path/to/chromap" >&2
    exit 1
  fi
fi

ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
OUT_ROOT="${OUT_ROOT:-${ARTIFACT_ROOT}/input_format_smoke/$(date -u +%Y%m%dT%H%M%SZ)}"
DATA_DIR="${OUT_ROOT}/data"
RUN_DIR="${OUT_ROOT}/runs"
BINSEQ_DIR="${OUT_ROOT}/binseq"
mkdir -p "${DATA_DIR}" "${RUN_DIR}" "${BINSEQ_DIR}"

log() {
  printf '[input-format-smoke] %s\n' "$*"
}

skip_binseq() {
  log "SKIP BINSEQ: $*"
  printf 'binseq_status=skipped\nreason=%s\n' "$*" > "${OUT_ROOT}/BINSEQ_SKIPPED.txt"
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
  command -v bqtools 2>/dev/null
}

normalize_sam_alignments() {
  python3 - "$1" "$2" <<'PY'
from collections import Counter
import sys

sam_path, out_path = sys.argv[1:3]
rows = Counter()
with open(sam_path, "rt", encoding="utf-8") as handle:
    for line in handle:
        if line.startswith("@"):
            continue
        rows[line.rstrip("\n")] += 1
with open(out_path, "wt", encoding="utf-8") as out:
    for row, count in sorted(rows.items()):
        out.write(f"{count}\t{row}\n")
PY
}

assert_sam_parity() {
  local expected="$1"
  local observed="$2"
  local label="$3"
  local expected_norm="${expected}.norm.tsv"
  local observed_norm="${observed}.norm.tsv"
  normalize_sam_alignments "${expected}" "${expected_norm}"
  normalize_sam_alignments "${observed}" "${observed_norm}"
  if ! cmp -s "${expected_norm}" "${observed_norm}"; then
    echo "ERROR: SAM alignment parity failed for ${label}" >&2
    diff -u "${expected_norm}" "${observed_norm}" | head -80 >&2 || true
    exit 1
  fi
  log "PASS: ${label}"
}

run_pe_sam() {
  local label="$1"
  local r1="$2"
  local r2="$3"
  local out="${RUN_DIR}/${label}.sam"
  "${CHROMAP}" --SAM \
    -r "${DATA_DIR}/test_ref.fa" -x "${DATA_DIR}/test_ref.idx" \
    -1 "${r1}" -2 "${r2}" -o "${out}" \
    --min-num-seeds 1 --error-threshold 10 --max-insert-size 1000 -t 1 \
    > "${RUN_DIR}/${label}.stdout" 2> "${RUN_DIR}/${label}.stderr"
  printf '%s\n' "${out}"
}

run_se_sam() {
  local label="$1"
  local r1="$2"
  local out="${RUN_DIR}/${label}.sam"
  "${CHROMAP}" --SAM \
    -r "${DATA_DIR}/test_ref.fa" -x "${DATA_DIR}/test_ref.idx" \
    -1 "${r1}" -o "${out}" \
    --min-num-seeds 1 --error-threshold 10 -t 1 \
    > "${RUN_DIR}/${label}.stdout" 2> "${RUN_DIR}/${label}.stderr"
  printf '%s\n' "${out}"
}

count_fastq_reads() {
  local path="$1"
  case "${path}" in
    *.gz) gzip -cd "${path}" | awk 'END {print NR / 4}' ;;
    *) awk 'END {print NR / 4}' "${path}" ;;
  esac
}

encode_pair_cbq() {
  local bqtools="$1"
  local r1="$2"
  local r2="$3"
  local out="$4"
  local level="${5:-}"
  local level_args=()
  [[ -n "${level}" ]] && level_args=(-l "${level}")

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" "${r2}" --mode cbq "${level_args[@]}" -o "${out}" -T 2 \
      > "${out}.encode.stdout" 2> "${out}.encode.stderr"; then
    [[ -s "${out}" ]] && return 0
  fi

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" "${r2}" -o "${out}" --mode cbq "${level_args[@]}" -T 2 \
      >> "${out}.encode.stdout" 2>> "${out}.encode.stderr"; then
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
  local level_args=()
  [[ -n "${level}" ]] && level_args=(-l "${level}")

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" --mode cbq "${level_args[@]}" -o "${out}" -T 2 \
      > "${out}.encode.stdout" 2> "${out}.encode.stderr"; then
    [[ -s "${out}" ]] && return 0
  fi

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" -o "${out}" --mode cbq "${level_args[@]}" -T 2 \
      >> "${out}.encode.stdout" 2>> "${out}.encode.stderr"; then
    [[ -s "${out}" ]] && return 0
  fi

  rm -f "${out}"
  "${bqtools}" encode "${r1}" -o "${out}" "${level_args[@]}" -T 2 \
      >> "${out}.encode.stdout" 2>> "${out}.encode.stderr"
  [[ -s "${out}" ]]
}

decode_pair_cbq() {
  local bqtools="$1"
  local cbq="$2"
  local prefix="$3"
  "${bqtools}" decode "${cbq}" -f q --prefix "${prefix}" \
    > "${prefix}.decode.stdout" 2> "${prefix}.decode.stderr"
  local r1=""
  local r2=""
  for candidate in "${prefix}_R1.fastq" "${prefix}.R1.fastq" "${prefix}_1.fastq" "${prefix}.1.fastq" \
                   "${prefix}_R1.fq" "${prefix}.R1.fq" "${prefix}_1.fq" "${prefix}.1.fq"; do
    [[ -f "${candidate}" ]] && { r1="${candidate}"; break; }
  done
  for candidate in "${prefix}_R2.fastq" "${prefix}.R2.fastq" "${prefix}_2.fastq" "${prefix}.2.fastq" \
                   "${prefix}_R2.fq" "${prefix}.R2.fq" "${prefix}_2.fq" "${prefix}.2.fq"; do
    [[ -f "${candidate}" ]] && { r2="${candidate}"; break; }
  done
  [[ -n "${r1}" && -n "${r2}" ]] || {
    echo "ERROR: bqtools decode did not produce paired FASTQ for ${cbq}" >&2
    exit 1
  }
  printf '%s\t%s\n' "${r1}" "${r2}"
}

decode_single_cbq() {
  local bqtools="$1"
  local cbq="$2"
  local out="$3"
  "${bqtools}" decode "${cbq}" -f q -o "${out}" \
    > "${out}.decode.stdout" 2> "${out}.decode.stderr"
  [[ -s "${out}" ]] || {
    echo "ERROR: bqtools decode did not produce single-end FASTQ for ${cbq}" >&2
    exit 1
  }
  printf '%s\n' "${out}"
}

log "outputs: ${OUT_ROOT}"
log "chromap: ${CHROMAP}"

python3 "${SCRIPT_DIR}/data/generate_test_data.py" -o "${DATA_DIR}" > "${OUT_ROOT}/generate_test_data.log"
"${CHROMAP}" --build-index -r "${DATA_DIR}/test_ref.fa" -o "${DATA_DIR}/test_ref.idx" \
  > "${OUT_ROOT}/build_index.stdout" 2> "${OUT_ROOT}/build_index.stderr"

gzip -c "${DATA_DIR}/test_pe_R1.fq" > "${DATA_DIR}/test_pe_R1.fq.gz"
gzip -c "${DATA_DIR}/test_pe_R2.fq" > "${DATA_DIR}/test_pe_R2.fq.gz"
gzip -c "${DATA_DIR}/test_se.fq" > "${DATA_DIR}/test_se.fq.gz"

pe_plain="$(run_pe_sam pe_fastq_plain "${DATA_DIR}/test_pe_R1.fq" "${DATA_DIR}/test_pe_R2.fq")"
pe_gzip="$(run_pe_sam pe_fastq_gzip "${DATA_DIR}/test_pe_R1.fq.gz" "${DATA_DIR}/test_pe_R2.fq.gz")"
assert_sam_parity "${pe_plain}" "${pe_gzip}" "paired FASTQ plain vs gzip input"

se_plain="$(run_se_sam se_fastq_plain "${DATA_DIR}/test_se.fq")"
se_gzip="$(run_se_sam se_fastq_gzip "${DATA_DIR}/test_se.fq.gz")"
assert_sam_parity "${se_plain}" "${se_gzip}" "single-end FASTQ plain vs gzip input"

"${CHROMAP}" --SAM --emit-Y-read-names --emit-Y-noY-fastq \
  --emit-Y-noY-fastq-compression none \
  -r "${DATA_DIR}/test_ref.fa" -x "${DATA_DIR}/test_ref.idx" \
  -1 "${DATA_DIR}/test_pe_R1.fq.gz" -2 "${DATA_DIR}/test_pe_R2.fq.gz" \
  -o "${RUN_DIR}/pe_split_none.sam" \
  --min-num-seeds 1 --error-threshold 10 --max-insert-size 1000 -t 1 \
  > "${RUN_DIR}/pe_split_none.stdout" 2> "${RUN_DIR}/pe_split_none.stderr"

"${CHROMAP}" --SAM --emit-Y-read-names --emit-Y-noY-fastq \
  --emit-Y-noY-fastq-compression gz \
  -r "${DATA_DIR}/test_ref.fa" -x "${DATA_DIR}/test_ref.idx" \
  -1 "${DATA_DIR}/test_pe_R1.fq.gz" -2 "${DATA_DIR}/test_pe_R2.fq.gz" \
  -o "${RUN_DIR}/pe_split_gz.sam" \
  --min-num-seeds 1 --error-threshold 10 --max-insert-size 1000 -t 1 \
  > "${RUN_DIR}/pe_split_gz.stdout" 2> "${RUN_DIR}/pe_split_gz.stderr"

for stream in Y noY; do
  for mate in R1 R2; do
    none_path="${RUN_DIR}/test_pe_${stream}_${mate}.fq"
    gz_path="${RUN_DIR}/test_pe_${stream}_${mate}.fq.gz"
    [[ -f "${none_path}" && -f "${gz_path}" ]] || {
      echo "ERROR: missing Y/noY FASTQ split output ${none_path} or ${gz_path}" >&2
      exit 1
    }
    none_count="$(count_fastq_reads "${none_path}")"
    gz_count="$(count_fastq_reads "${gz_path}")"
    [[ "${none_count}" == "${gz_count}" ]] || {
      echo "ERROR: split FASTQ count mismatch for ${stream}_${mate}: ${none_count} vs ${gz_count}" >&2
      exit 1
    }
  done
done
log "PASS: Y/noY FASTQ split compression none vs gz"

if bqtools_bin="$(resolve_bqtools)"; then
  log "bqtools: ${bqtools_bin}"
  pair_cbq="${BINSEQ_DIR}/test_pe.cbq"
  pair_cbq_uncompressed="${BINSEQ_DIR}/test_pe_uncompressed.cbq"
  single_cbq="${BINSEQ_DIR}/test_se.cbq"
  single_cbq_uncompressed="${BINSEQ_DIR}/test_se_uncompressed.cbq"

  encode_pair_cbq "${bqtools_bin}" "${DATA_DIR}/test_pe_R1.fq" "${DATA_DIR}/test_pe_R2.fq" "${pair_cbq}"
  encode_pair_cbq "${bqtools_bin}" "${DATA_DIR}/test_pe_R1.fq" "${DATA_DIR}/test_pe_R2.fq" "${pair_cbq_uncompressed}" 0
  encode_single_cbq "${bqtools_bin}" "${DATA_DIR}/test_se.fq" "${single_cbq}"
  encode_single_cbq "${bqtools_bin}" "${DATA_DIR}/test_se.fq" "${single_cbq_uncompressed}" 0

  "${bqtools_bin}" info --json "${pair_cbq}" > "${pair_cbq}.info.json" 2> "${pair_cbq}.info.stderr" || \
    "${bqtools_bin}" info "${pair_cbq}" > "${pair_cbq}.info.txt" 2> "${pair_cbq}.info.stderr"
  "${bqtools_bin}" info --json "${pair_cbq_uncompressed}" > "${pair_cbq_uncompressed}.info.json" 2> "${pair_cbq_uncompressed}.info.stderr" || \
    "${bqtools_bin}" info "${pair_cbq_uncompressed}" > "${pair_cbq_uncompressed}.info.txt" 2> "${pair_cbq_uncompressed}.info.stderr"
  "${bqtools_bin}" info --json "${single_cbq}" > "${single_cbq}.info.json" 2> "${single_cbq}.info.stderr" || \
    "${bqtools_bin}" info "${single_cbq}" > "${single_cbq}.info.txt" 2> "${single_cbq}.info.stderr"
  "${bqtools_bin}" info --json "${single_cbq_uncompressed}" > "${single_cbq_uncompressed}.info.json" 2> "${single_cbq_uncompressed}.info.stderr" || \
    "${bqtools_bin}" info "${single_cbq_uncompressed}" > "${single_cbq_uncompressed}.info.txt" 2> "${single_cbq_uncompressed}.info.stderr"

  IFS=$'\t' read -r pair_dec_r1 pair_dec_r2 < <(decode_pair_cbq "${bqtools_bin}" "${pair_cbq}" "${BINSEQ_DIR}/decoded_pair")
  IFS=$'\t' read -r pair_uc_dec_r1 pair_uc_dec_r2 < <(decode_pair_cbq "${bqtools_bin}" "${pair_cbq_uncompressed}" "${BINSEQ_DIR}/decoded_pair_uncompressed")
  single_dec="$(decode_single_cbq "${bqtools_bin}" "${single_cbq}" "${BINSEQ_DIR}/decoded_single.fastq")"
  single_uc_dec="$(decode_single_cbq "${bqtools_bin}" "${single_cbq_uncompressed}" "${BINSEQ_DIR}/decoded_single_uncompressed.fastq")"

  pe_binseq="$(run_pe_sam pe_binseq_decoded "${pair_dec_r1}" "${pair_dec_r2}")"
  pe_binseq_uc="$(run_pe_sam pe_binseq_uncompressed_decoded "${pair_uc_dec_r1}" "${pair_uc_dec_r2}")"
  assert_sam_parity "${pe_plain}" "${pe_binseq}" "paired FASTQ vs default-compressed CBQ decoded input"
  assert_sam_parity "${pe_plain}" "${pe_binseq_uc}" "paired FASTQ vs uncompressed CBQ decoded input"

  se_binseq="$(run_se_sam se_binseq_decoded "${single_dec}")"
  se_binseq_uc="$(run_se_sam se_binseq_uncompressed_decoded "${single_uc_dec}")"
  assert_sam_parity "${se_plain}" "${se_binseq}" "single-end FASTQ vs default-compressed CBQ decoded input"
  assert_sam_parity "${se_plain}" "${se_binseq_uc}" "single-end FASTQ vs uncompressed CBQ decoded input"
else
  skip_binseq "bqtools not found; set BQTOOLS=/path/to/bqtools to exercise CBQ cases"
fi

cat > "${OUT_ROOT}/SUMMARY.txt" <<EOF
chromap=${CHROMAP}
out_root=${OUT_ROOT}
fastq_plain_gzip=pass
y_noy_fastq_compression=pass
binseq_cbq=$([[ -f "${OUT_ROOT}/BINSEQ_SKIPPED.txt" ]] && echo skipped || echo pass)
EOF

log "PASS: Chromap input format smoke completed at ${OUT_ROOT}"
