#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
CACHE_ROOT="${ENCODE_CROSS_ASSAY_CACHE:-${ARTIFACT_ROOT}/encode_cross_assay_cache}"
SOURCE_MANIFEST="${ENCODE_CROSS_ASSAY_SOURCE_MANIFEST:-${REPO_ROOT}/tests/encode_cross_assay_manifest.tsv}"
GENERATED_MANIFEST="${ENCODE_CROSS_ASSAY_GENERATED_MANIFEST:-${CACHE_ROOT}/manifest.generated.tsv}"
READ_COUNT="${ENCODE_DOWNSAMPLE_READS:-}"
ALLOW_DOWNLOAD="${ENCODE_ALLOW_DOWNLOAD:-0}"
ASSAYS="${ENCODE_ASSAYS:-all}"
FETCH_METADATA="${ENCODE_FETCH_METADATA:-1}"

log() {
  printf '[encode-cross-fixtures] %s\n' "$*" >&2
}

fail() {
  printf '[encode-cross-fixtures] FAIL: %s\n' "$*" >&2
  exit 1
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

safe_id() {
  printf '%s' "$1" | tr -c 'A-Za-z0-9._-' '_'
}

download_if_needed() {
  local url="$1"
  local out="$2"
  if [[ -z "${url}" ]]; then
    fail "download URL is empty for ${out}"
  fi
  if [[ -s "${out}" ]]; then
    return
  fi
  if [[ "${ALLOW_DOWNLOAD}" != "1" ]]; then
    fail "missing ${out}; set ENCODE_ALLOW_DOWNLOAD=1 to download ${url}"
  fi
  mkdir -p "$(dirname "${out}")"
  log "downloading ${url}"
  rm -f "${out}.tmp"
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
  mkdir -p "$(dirname "${out}")"
  log "downsampling $(basename "${in}") to ${reads} read(s)"
  rm -f "${out}.tmp"
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

verify_fastq_alignment() {
  python3 - "$@" <<'PY'
import gzip
import sys

paths = sys.argv[1:]
if len(paths) < 2:
    raise SystemExit("need at least two FASTQs to verify alignment")

def norm(header):
    name = header.strip()
    if name.startswith("@"):
        name = name[1:]
    name = name.split()[0]
    for suffix in ("/1", "/2", "/3"):
        if name.endswith(suffix):
            name = name[:-2]
    return name

handles = [gzip.open(path, "rt", encoding="ascii", errors="replace") for path in paths]
try:
    count = 0
    while True:
        records = [[handle.readline() for _ in range(4)] for handle in handles]
        present = [bool(record[0]) for record in records]
        if not any(present):
            break
        if not all(present):
            raise SystemExit("FASTQs have different record counts")
        count += 1
        for path, record in zip(paths, records):
            if any(line == "" for line in record):
                raise SystemExit(f"truncated FASTQ record in {path}")
        names = [norm(record[0]) for record in records]
        if any(name != names[0] for name in names[1:]):
            raise SystemExit(
                f"read name mismatch at record {count}: " + " vs ".join(names)
            )
finally:
    for handle in handles:
        handle.close()
print(count)
PY
}

md5_file() {
  if [[ -n "$1" && -s "$1" ]]; then
    md5sum "$1" | awk '{print $1}'
  fi
}

fetch_metadata() {
  if [[ "${FETCH_METADATA}" != "1" ]]; then
    printf 'not_fetched\tnot_fetched\tnot_fetched\tnot_fetched\n'
    return
  fi
  python3 - "$@" <<'PY'
import json
import sys
import urllib.request

statuses = []
formats = []
outputs = []
read_lengths = []
for accession in sys.argv[1:]:
    if not accession:
        statuses.append("")
        formats.append("")
        outputs.append("")
        read_lengths.append("")
        continue
    try:
        req = urllib.request.Request(
            f"https://www.encodeproject.org/files/{accession}/?format=json",
            headers={"Accept": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=20) as response:
            obj = json.load(response)
        statuses.append(str(obj.get("status", "unknown")))
        formats.append(str(obj.get("file_format", "unknown")))
        outputs.append(str(obj.get("output_type", "unknown")))
        read_lengths.append(str(obj.get("read_length", "")))
    except Exception:
        statuses.append("unknown")
        formats.append("unknown")
        outputs.append("unknown")
        read_lengths.append("")

print(
    "\t".join(
        [
            "|".join(statuses),
            "|".join(formats),
            "|".join(outputs),
            "|".join(read_lengths),
        ]
    )
)
PY
}

prepare_whitelist() {
  local whitelist_id="$1"
  local whitelist_url="$2"
  if [[ -n "${SCATAC_WHITELIST:-}" ]]; then
    [[ -s "${SCATAC_WHITELIST}" ]] || fail "SCATAC_WHITELIST does not exist: ${SCATAC_WHITELIST}"
    printf '%s\n' "${SCATAC_WHITELIST}"
    return
  fi
  if [[ -z "${whitelist_id}" || -z "${whitelist_url}" ]]; then
    printf '\n'
    return
  fi

  local safe
  safe="$(safe_id "${whitelist_id}")"
  local raw="${CACHE_ROOT}/downloads/${safe}.txt.gz"
  local out="${CACHE_ROOT}/references/${safe}.txt"
  download_if_needed "${whitelist_url}" "${raw}"
  if [[ ! -s "${out}" ]]; then
    mkdir -p "$(dirname "${out}")"
    log "decompressing whitelist ${whitelist_id}"
    rm -f "${out}.tmp"
    gzip -cd "${raw}" >"${out}.tmp"
    mv "${out}.tmp" "${out}"
  fi
  printf '%s\n' "${out}"
}

main() {
  [[ -s "${SOURCE_MANIFEST}" ]] || fail "missing source manifest ${SOURCE_MANIFEST}"
  mkdir -p "${CACHE_ROOT}"/{downloads,downsampled,references}
  printf 'assay\tcase_id\texperiment\tlayout\tread_count\tr1_fastq\tr2_fastq\tbarcode_fastq\twhitelist\tr1_md5\tr2_md5\tbarcode_md5\twhitelist_md5\tchromap_args\tsource_status\tsource_file_formats\tsource_output_types\tsource_read_lengths\tdescription\tr1_accession\tr2_accession\tbarcode_accession\twhitelist_id\tr1_url\tr2_url\tbarcode_url\twhitelist_url\n' \
    >"${GENERATED_MANIFEST}"

  local rows=0
  local line
  while IFS= read -r line || [[ -n "${line}" ]]; do
    [[ -z "${line}" ]] && continue
    mapfile -t fields < <(parse_tsv_line "${line}")
    local assay="${fields[0]:-}"
    [[ "${assay}" == "assay" ]] && continue
    [[ -n "${assay}" ]] || continue
    want_assay "${assay}" || continue
    [[ "${#fields[@]}" -ge 15 ]] || fail "malformed manifest row for ${assay}: expected 15 columns"

    local case_id="${fields[1]}"
    local experiment="${fields[2]}"
    local layout="${fields[3]}"
    local r1_acc="${fields[4]}"
    local r2_acc="${fields[5]}"
    local barcode_acc="${fields[6]}"
    local whitelist_id="${fields[7]}"
    local r1_url="${fields[8]}"
    local r2_url="${fields[9]}"
    local barcode_url="${fields[10]}"
    local whitelist_url="${fields[11]}"
    local default_reads="${fields[12]}"
    local chromap_args="${fields[13]}"
    local description="${fields[14]}"

    rows=$((rows + 1))
    local reads="${READ_COUNT:-${default_reads}}"
    local case_dir="${CACHE_ROOT}/downsampled/${case_id}_${reads}"
    local full_r1="${CACHE_ROOT}/downloads/${r1_acc}.fastq.gz"
    local full_r2="${CACHE_ROOT}/downloads/${r2_acc}.fastq.gz"
    local full_barcode=""
    local ds_r1="${case_dir}/${r1_acc}.${reads}.fastq.gz"
    local ds_r2="${case_dir}/${r2_acc}.${reads}.fastq.gz"
    local ds_barcode=""
    local whitelist=""

    download_if_needed "${r1_url}" "${full_r1}"
    download_if_needed "${r2_url}" "${full_r2}"
    downsample_fastq "${full_r1}" "${ds_r1}" "${reads}"
    downsample_fastq "${full_r2}" "${ds_r2}" "${reads}"

    local observed_reads
    if [[ "${layout}" == "paired" ]]; then
      observed_reads="$(verify_fastq_alignment "${ds_r1}" "${ds_r2}")"
    elif [[ "${layout}" == "scatac_10x_atac" ]]; then
      [[ -n "${barcode_acc}" ]] || fail "${case_id}: scatac layout requires barcode_accession"
      [[ -n "${barcode_url}" ]] || fail "${case_id}: scatac layout requires barcode_url"
      full_barcode="${CACHE_ROOT}/downloads/${barcode_acc}.fastq.gz"
      ds_barcode="${case_dir}/${barcode_acc}.${reads}.fastq.gz"
      download_if_needed "${barcode_url}" "${full_barcode}"
      downsample_fastq "${full_barcode}" "${ds_barcode}" "${reads}"
      whitelist="$(prepare_whitelist "${whitelist_id}" "${whitelist_url}")"
      observed_reads="$(verify_fastq_alignment "${ds_r1}" "${ds_r2}" "${ds_barcode}")"
    else
      fail "${case_id}: unsupported layout ${layout}"
    fi
    [[ "${observed_reads}" == "${reads}" ]] || \
      fail "${case_id}: expected ${reads} records, observed ${observed_reads}"

    local source_status source_formats source_outputs source_read_lengths
    IFS=$'\t' read -r source_status source_formats source_outputs source_read_lengths < <(
      fetch_metadata "${r1_acc}" "${r2_acc}" "${barcode_acc}"
    )

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "${assay}" "${case_id}" "${experiment}" "${layout}" "${observed_reads}" \
      "${ds_r1}" "${ds_r2}" "${ds_barcode}" "${whitelist}" \
      "$(md5_file "${ds_r1}")" "$(md5_file "${ds_r2}")" \
      "$(md5_file "${ds_barcode}")" "$(md5_file "${whitelist}")" \
      "${chromap_args}" "${source_status}" "${source_formats}" \
      "${source_outputs}" "${source_read_lengths}" "${description}" \
      "${r1_acc}" "${r2_acc}" "${barcode_acc}" "${whitelist_id}" \
      "${r1_url}" "${r2_url}" "${barcode_url}" "${whitelist_url}" \
      >>"${GENERATED_MANIFEST}"
  done <"${SOURCE_MANIFEST}"

  [[ "${rows}" -gt 0 ]] || fail "no assays selected by ENCODE_ASSAYS=${ASSAYS}"
  log "wrote ${GENERATED_MANIFEST}"
  cat "${GENERATED_MANIFEST}"
}

main "$@"
