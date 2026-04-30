#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
CACHE_ROOT="${ENCODE_FIXTURE_CACHE:-${ARTIFACT_ROOT}/encode_fixture_cache}"
SOURCE_MANIFEST="${ENCODE_FIXTURE_SOURCE_MANIFEST:-${REPO_ROOT}/tests/encode_downsample_manifest.tsv}"
READ_COUNT="${ENCODE_DOWNSAMPLE_READS:-}"
ALLOW_DOWNLOAD="${ENCODE_ALLOW_DOWNLOAD:-0}"
ASSAYS="${ENCODE_ASSAYS:-all}"
GENERATED_MANIFEST="${ENCODE_GENERATED_MANIFEST:-${CACHE_ROOT}/manifest.generated.tsv}"

log() {
  printf '[encode-fixtures] %s\n' "$*" >&2
}

fail() {
  printf '[encode-fixtures] FAIL: %s\n' "$*" >&2
  exit 1
}

want_assay() {
  local assay="$1"
  [[ "${ASSAYS}" == "all" ]] && return 0
  [[ ",${ASSAYS}," == *",${assay},"* ]]
}

download_if_needed() {
  local url="$1"
  local out="$2"
  if [[ -s "${out}" ]]; then
    return
  fi
  if [[ "${ALLOW_DOWNLOAD}" != "1" ]]; then
    fail "missing ${out}; set ENCODE_ALLOW_DOWNLOAD=1 to download ${url}"
  fi
  mkdir -p "$(dirname "${out}")"
  log "downloading ${url}"
  curl -L --fail --retry 3 --connect-timeout 30 -o "${out}.tmp" "${url}"
  mv "${out}.tmp" "${out}"
}

downsample_fastq() {
  local in="$1"
  local out="$2"
  local reads="$3"
  if [[ -s "${out}" ]]; then
    return
  fi
  log "downsampling $(basename "${in}") to ${reads} read(s)"
  python3 - "${in}" "${reads}" <<'PY' | gzip -n >"${out}.tmp"
import gzip
import sys

path = sys.argv[1]
limit = int(sys.argv[2])
emitted = 0
with gzip.open(path, "rt", encoding="ascii", errors="replace") as handle:
    while emitted < limit:
        record = [handle.readline() for _ in range(4)]
        if not record[0]:
            break
        if any(line == "" for line in record):
            raise SystemExit(f"truncated FASTQ record in {path}")
        sys.stdout.write("".join(record))
        emitted += 1
if emitted == 0:
    raise SystemExit(f"no FASTQ records found in {path}")
PY
  mv "${out}.tmp" "${out}"
}

fastq_count() {
  local path="$1"
  python3 - "${path}" <<'PY'
import gzip
import sys

count = 0
with gzip.open(sys.argv[1], "rt", encoding="ascii", errors="replace") as handle:
    for count, _ in enumerate(handle, 1):
        pass
if count % 4:
    raise SystemExit(f"FASTQ line count is not divisible by four: {sys.argv[1]}")
print(count // 4)
PY
}

verify_pairing() {
  local r1="$1"
  local r2="$2"
  python3 - "${r1}" "${r2}" <<'PY'
import gzip
import sys

def norm(name):
    name = name.strip()
    if name.startswith("@"):
        name = name[1:]
    name = name.split()[0]
    for suffix in ("/1", "/2"):
        if name.endswith(suffix):
            name = name[:-2]
    return name

with gzip.open(sys.argv[1], "rt", encoding="ascii", errors="replace") as a, \
     gzip.open(sys.argv[2], "rt", encoding="ascii", errors="replace") as b:
    i = 0
    while True:
        r1 = [a.readline() for _ in range(4)]
        r2 = [b.readline() for _ in range(4)]
        if not r1[0] and not r2[0]:
            break
        if not r1[0] or not r2[0]:
            raise SystemExit("paired FASTQs have different record counts")
        i += 1
        if norm(r1[0]) != norm(r2[0]):
            raise SystemExit(
                f"read name mismatch at record {i}: {r1[0].strip()} vs {r2[0].strip()}"
            )
print(i)
PY
}

md5_file() {
  md5sum "$1" | awk '{print $1}'
}

main() {
  [[ -s "${SOURCE_MANIFEST}" ]] || fail "missing source manifest ${SOURCE_MANIFEST}"
  mkdir -p "${CACHE_ROOT}"/{downloads,downsampled}
  printf 'assay\tcase_id\texperiment\tr1_accession\tr2_accession\tread_count\tr1_fastq\tr2_fastq\tr1_md5\tr2_md5\tchromap_args\tdescription\n' \
    >"${GENERATED_MANIFEST}"

  local rows=0
  local line
  while IFS=$'\t' read -r assay case_id experiment r1_acc r2_acc r1_url r2_url default_reads chromap_args description; do
    [[ "${assay}" == "assay" ]] && continue
    [[ -n "${assay}" ]] || continue
    want_assay "${assay}" || continue
    rows=$((rows + 1))

    local reads="${READ_COUNT:-${default_reads}}"
    local pair_dir="${CACHE_ROOT}/downsampled/${assay}_${experiment}_${reads}"
    local full_r1="${CACHE_ROOT}/downloads/${r1_acc}.fastq.gz"
    local full_r2="${CACHE_ROOT}/downloads/${r2_acc}.fastq.gz"
    local ds_r1="${pair_dir}/${r1_acc}.${reads}.fastq.gz"
    local ds_r2="${pair_dir}/${r2_acc}.${reads}.fastq.gz"

    download_if_needed "${r1_url}" "${full_r1}"
    download_if_needed "${r2_url}" "${full_r2}"
    mkdir -p "${pair_dir}"
    downsample_fastq "${full_r1}" "${ds_r1}" "${reads}"
    downsample_fastq "${full_r2}" "${ds_r2}" "${reads}"

    local observed_reads
    observed_reads="$(verify_pairing "${ds_r1}" "${ds_r2}")"
    [[ "${observed_reads}" == "${reads}" ]] || \
      fail "${case_id}: expected ${reads} read pairs, observed ${observed_reads}"

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${assay}" "${case_id}" "${experiment}" "${r1_acc}" "${r2_acc}" \
      "${observed_reads}" "${ds_r1}" "${ds_r2}" \
      "$(md5_file "${ds_r1}")" "$(md5_file "${ds_r2}")" \
      "${chromap_args}" "${description}" >>"${GENERATED_MANIFEST}"
  done <"${SOURCE_MANIFEST}"

  [[ "${rows}" -gt 0 ]] || fail "no assays selected by ENCODE_ASSAYS=${ASSAYS}"
  log "wrote ${GENERATED_MANIFEST}"
  cat "${GENERATED_MANIFEST}"
}

main "$@"
