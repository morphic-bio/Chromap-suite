#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"

CHROMAP="${CHROMAP:-${REPO_ROOT}/chromap}"
CHROMAP_LIB_RUNNER="${CHROMAP_LIB_RUNNER:-${REPO_ROOT}/chromap_lib_runner}"
ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
OUT_ROOT="${OUT_ROOT:-${ARTIFACT_ROOT}/cbq_atac_smoke/$(date -u +%Y%m%dT%H%M%SZ)}"
DATA_DIR="${OUT_ROOT}/data"
RUN_DIR="${OUT_ROOT}/runs"
CBQ_DIR="${OUT_ROOT}/cbq"

mkdir -p "${DATA_DIR}" "${RUN_DIR}" "${CBQ_DIR}"

log() {
  printf '[cbq-atac-smoke] %s\n' "$*"
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

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" "${r2}" --mode cbq -o "${out}" -T 2 \
      > "${out}.encode.stdout" 2> "${out}.encode.stderr"; then
    [[ -s "${out}" ]] && return 0
  fi

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" "${r2}" -o "${out}" --mode cbq -T 2 \
      >> "${out}.encode.stdout" 2>> "${out}.encode.stderr"; then
    [[ -s "${out}" ]] && return 0
  fi

  rm -f "${out}"
  "${bqtools}" encode "${r1}" "${r2}" -o "${out}" -T 2 \
      >> "${out}.encode.stdout" 2>> "${out}.encode.stderr"
  [[ -s "${out}" ]]
}

encode_single_cbq() {
  local bqtools="$1"
  local r1="$2"
  local out="$3"

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" --mode cbq -o "${out}" -T 2 \
      > "${out}.encode.stdout" 2> "${out}.encode.stderr"; then
    [[ -s "${out}" ]] && return 0
  fi

  rm -f "${out}"
  if "${bqtools}" encode "${r1}" -o "${out}" --mode cbq -T 2 \
      >> "${out}.encode.stdout" 2>> "${out}.encode.stderr"; then
    [[ -s "${out}" ]] && return 0
  fi

  rm -f "${out}"
  "${bqtools}" encode "${r1}" -o "${out}" -T 2 \
      >> "${out}.encode.stdout" 2>> "${out}.encode.stderr"
  [[ -s "${out}" ]]
}

normalize_fragments() {
  local input="$1"
  local output="$2"
  LC_ALL=C sort "${input}" > "${output}"
}

assert_fragment_parity() {
  local expected="$1"
  local observed="$2"
  local label="$3"
  local expected_norm="${expected}.norm"
  local observed_norm="${observed}.norm"
  normalize_fragments "${expected}" "${expected_norm}"
  normalize_fragments "${observed}" "${observed_norm}"
  if ! cmp -s "${expected_norm}" "${observed_norm}"; then
    echo "ERROR: fragment parity failed for ${label}" >&2
    diff -u "${expected_norm}" "${observed_norm}" | head -80 >&2 || true
    exit 1
  fi
  log "PASS: ${label}"
}

count_fastq_reads() {
  local path="$1"
  if [[ ! -f "${path}" ]]; then
    printf '0\n'
    return
  fi
  grep -c '^@' "${path}" || true
}

[[ -x "${CHROMAP}" ]] || {
  echo "ERROR: chromap binary not found; run make chromap or set CHROMAP" >&2
  exit 1
}
[[ -x "${CHROMAP_LIB_RUNNER}" ]] || {
  echo "ERROR: chromap_lib_runner not found; run make chromap_lib_runner or set CHROMAP_LIB_RUNNER" >&2
  exit 1
}

bqtools_bin="$(resolve_bqtools)" || skip "bqtools not found; set BQTOOLS=/path/to/bqtools"

log "outputs: ${OUT_ROOT}"
log "bqtools: ${bqtools_bin}"

python3 - "${DATA_DIR}" <<'PY'
import random
import sys
from pathlib import Path

out = Path(sys.argv[1])
rng = random.Random(74219)
alphabet = "ACGT"
genome = "".join(rng.choice(alphabet) for _ in range(5000))

def rc(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

barcodes = [
    "ACGTACGTACGTACGT",
    "TGCATGCATGCATGCA",
    "GGGGAAAACCCCTTTT",
    "AAAACCCCGGGGTTTT",
]
fragments = [
    ("read001", 200, 340, barcodes[0], "I"),
    ("read002", 700, 850, barcodes[1], "H"),
    ("read003", 1500, 1660, barcodes[2], "J"),
    ("read004", 2600, 2765, barcodes[3], "I"),
]

(out / "genome.fa").write_text(">chrSynthetic\n" + genome + "\n", encoding="ascii")
(out / "whitelist.txt").write_text("\n".join(barcodes) + "\n", encoding="ascii")

with (out / "reads_R1.fastq").open("wt", encoding="ascii") as r1, \
     (out / "reads_R2.fastq").open("wt", encoding="ascii") as r2, \
     (out / "barcodes.fastq").open("wt", encoding="ascii") as bc:
    for name, start, end, barcode, qchar in fragments:
        seq1 = genome[start:start + 70]
        seq2 = rc(genome[end - 70:end])
        r1.write(f"@{name}/1 1:N:0:ACGT\n{seq1}\n+\n{qchar * len(seq1)}\n")
        r2.write(f"@{name}/2 2:N:0:ACGT\n{seq2}\n+\n{qchar * len(seq2)}\n")
        bc.write(f"@{name}/3 3:N:0:ACGT\n{barcode}\n+\n{qchar * len(barcode)}\n")
PY

"${CHROMAP}" --build-index -r "${DATA_DIR}/genome.fa" -o "${DATA_DIR}/genome.idx" \
  -k 11 -w 5 > "${RUN_DIR}/index.stdout" 2> "${RUN_DIR}/index.stderr"

read_pair_cbq="${CBQ_DIR}/reads_pair.cbq"
barcode_cbq="${CBQ_DIR}/barcodes.cbq"
encode_pair_cbq "${bqtools_bin}" "${DATA_DIR}/reads_R1.fastq" \
  "${DATA_DIR}/reads_R2.fastq" "${read_pair_cbq}"
encode_single_cbq "${bqtools_bin}" "${DATA_DIR}/barcodes.fastq" "${barcode_cbq}"

common_args=(
  --preset atac
  -r "${DATA_DIR}/genome.fa"
  -x "${DATA_DIR}/genome.idx"
  --barcode-whitelist "${DATA_DIR}/whitelist.txt"
  --read-format "bc:0:-1"
  --BED
  -t 1
  --min-read-length 20
  --min-num-seeds 1
  -q 0
)

"${CHROMAP}" "${common_args[@]}" \
  -1 "${DATA_DIR}/reads_R1.fastq" \
  -2 "${DATA_DIR}/reads_R2.fastq" \
  -b "${DATA_DIR}/barcodes.fastq" \
  -o "${RUN_DIR}/fastq.fragments.bed" \
  > "${RUN_DIR}/fastq.stdout" 2> "${RUN_DIR}/fastq.stderr"

"${CHROMAP}" "${common_args[@]}" \
  --input-format cbq \
  --read-pair-cbq "${read_pair_cbq}" \
  --barcode-cbq "${barcode_cbq}" \
  -o "${RUN_DIR}/cbq.fragments.bed" \
  > "${RUN_DIR}/cbq.stdout" 2> "${RUN_DIR}/cbq.stderr"

"${CHROMAP_LIB_RUNNER}" "${common_args[@]}" \
  --input-format cbq \
  --read-pair-cbq "${read_pair_cbq}" \
  --barcode-cbq "${barcode_cbq}" \
  -o "${RUN_DIR}/cbq.lib.fragments.bed" \
  > "${RUN_DIR}/cbq.lib.stdout" 2> "${RUN_DIR}/cbq.lib.stderr"

"${CHROMAP}" "${common_args[@]}" \
  --input-format cbq \
  --read-pair-cbq "${read_pair_cbq}" \
  --barcode-cbq "${barcode_cbq}" \
  --emit-Y-noY-fastq \
  --emit-Y-noY-fastq-compression none \
  -o "${RUN_DIR}/cbq.split.fragments.bed" \
  > "${RUN_DIR}/cbq.split.stdout" 2> "${RUN_DIR}/cbq.split.stderr"

[[ -s "${RUN_DIR}/fastq.fragments.bed" ]] || {
  echo "ERROR: FASTQ baseline produced no fragments" >&2
  exit 1
}
[[ -s "${RUN_DIR}/cbq.fragments.bed" ]] || {
  echo "ERROR: CBQ CLI produced no fragments" >&2
  exit 1
}
[[ -s "${RUN_DIR}/cbq.lib.fragments.bed" ]] || {
  echo "ERROR: CBQ lib runner produced no fragments" >&2
  exit 1
}
[[ -s "${RUN_DIR}/cbq.split.fragments.bed" ]] || {
  echo "ERROR: CBQ CLI Y/noY FASTQ run produced no fragments" >&2
  exit 1
}

assert_fragment_parity "${RUN_DIR}/fastq.fragments.bed" \
  "${RUN_DIR}/cbq.fragments.bed" "CLI FASTQ vs CLI CBQ"
assert_fragment_parity "${RUN_DIR}/fastq.fragments.bed" \
  "${RUN_DIR}/cbq.lib.fragments.bed" "CLI FASTQ vs libchromap CBQ"
assert_fragment_parity "${RUN_DIR}/fastq.fragments.bed" \
  "${RUN_DIR}/cbq.split.fragments.bed" "CLI FASTQ vs CLI CBQ Y/noY sidecar"

for mate in 1 2; do
  y_fastq="${RUN_DIR}/y_separated/Y_reads.mate${mate}.fastq"
  noy_fastq="${RUN_DIR}/y_separated/noY_reads.mate${mate}.fastq"
  [[ -f "${y_fastq}" && -f "${noy_fastq}" ]] || {
    echo "ERROR: missing CBQ Y/noY FASTQ sidecars for mate ${mate}" >&2
    exit 1
  }
  [[ "$(count_fastq_reads "${y_fastq}")" == "0" ]] || {
    echo "ERROR: expected empty Y FASTQ sidecar for no-Y synthetic reference: ${y_fastq}" >&2
    exit 1
  }
  [[ "$(count_fastq_reads "${noy_fastq}")" == "4" ]] || {
    echo "ERROR: expected four noY FASTQ records for mate ${mate}: ${noy_fastq}" >&2
    exit 1
  }
done
log "PASS: CBQ Y/noY FASTQ sidecars"

cat > "${OUT_ROOT}/SUMMARY.txt" <<EOF
chromap=${CHROMAP}
chromap_lib_runner=${CHROMAP_LIB_RUNNER}
bqtools=${bqtools_bin}
out_root=${OUT_ROOT}
cbq_cli=pass
cbq_lib_runner=pass
cbq_y_noy_fastq=pass
EOF

log "PASS: CBQ ATAC smoke completed at ${OUT_ROOT}"
