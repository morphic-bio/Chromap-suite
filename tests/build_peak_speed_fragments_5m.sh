#!/usr/bin/env bash
# Build a fixed fragment subset for MACS3/C++ speed benchmarks.
# Does not commit outputs. Requires the source fragments TSV (plain, not .gz).
# shellcheck shell=bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

N_LINES="${N_LINES:-5000000}"
SRC="${SRC:-/mnt/pikachu/atac-seq/10xMultiome/pbmc_unsorted_3k/star_chromap_concurrent_full_20260424_174426/chromap_fragments.tsv}"
OUT_DIR="${OUT_DIR:-${REPO_ROOT}/out/peak_speed_parity_20260425}"
if (( N_LINES % 1000000 == 0 )); then
  LABEL="$((N_LINES / 1000000))m"
elif (( N_LINES % 1000 == 0 )); then
  LABEL="$((N_LINES / 1000))k"
else
  LABEL="${N_LINES}"
fi
OUT_GZ="${OUT_GZ:-${OUT_DIR}/fragments_${LABEL}.tsv.gz}"
META="${META:-${OUT_DIR}/fragments_${LABEL}.metadata.tsv}"

main() {
  [[ -f "${SRC}" ]] || {
    echo "Missing source fragments TSV: ${SRC}" >&2
    echo "Override with SRC=/path/to/chromap_fragments.tsv" >&2
    return 1
  }
  mkdir -p "${OUT_DIR}"
  echo "[subset] first ${N_LINES} lines of ${SRC} -> ${OUT_GZ}" >&2
  head -n "${N_LINES}" "${SRC}" | gzip -nc >"${OUT_GZ}.tmp"
  mv -f "${OUT_GZ}.tmp" "${OUT_GZ}"

  local md5 line_count
  md5="$(md5sum "${OUT_GZ}" | awk '{print $1}')"
  line_count="$(gzip -dc "${OUT_GZ}" | wc -l | tr -d ' ')"

  {
    echo -e "metric\tvalue"
    echo -e "source_tsv\t${SRC}"
    echo -e "output_gz\t${OUT_GZ}"
    echo -e "row_count\t${line_count}"
    echo -e "expected_rows\t${N_LINES}"
    echo -e "md5_gz\t${md5}"
    echo -e "built_at\t$(date -Iseconds 2>/dev/null || date)"
  } | tee "${META}"

  echo "Wrote ${META}" >&2
  return 0
}

main "$@"
