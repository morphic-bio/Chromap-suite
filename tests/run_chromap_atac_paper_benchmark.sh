#!/usr/bin/env bash
# One-dataset Chromap Suite paper benchmark.
#
# Required environment:
#   DATASET_ID, ASSAY_TYPE=(bulk|scatac), R1_CSV, R2_CSV, OUTDIR
#   BARCODE_CSV and WHITELIST are additionally required for scatac.
#
# The harness deliberately maps twice: once without peak calling and once with
# integrated RapidMACS.  The two decompressed fragment streams must compare
# byte-for-byte before standalone MACS3 or any derived metric is accepted.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

DATASET_ID="${DATASET_ID:?set DATASET_ID}"
ASSAY_TYPE="${ASSAY_TYPE:?set ASSAY_TYPE to bulk or scatac}"
R1_CSV="${R1_CSV:?set R1_CSV}"
R2_CSV="${R2_CSV:?set R2_CSV}"
OUTDIR="${OUTDIR:?set OUTDIR to a fresh path}"
BARCODE_CSV="${BARCODE_CSV:-}"
WHITELIST="${WHITELIST:-}"

CHROMAP="${CHROMAP:-${REPO_ROOT}/chromap}"
MACS3="${MACS3:-macs3}"
INDEX="${INDEX:?set INDEX to the Chromap reference index}"
REF="${REF:?set REF to the matching reference FASTA}"
THREADS="${THREADS:-16}"
MACS3_PVALUE="${MACS3_PVALUE:-0.01}"
MACS3_MIN_LENGTH="${MACS3_MIN_LENGTH:-200}"
MACS3_MAX_GAP="${MACS3_MAX_GAP:-30}"
LOW_MEM_RAM="${LOW_MEM_RAM:-1G}"
TEMP_SAMPLE_SECONDS="${TEMP_SAMPLE_SECONDS:-5}"
HASH_INPUTS="${HASH_INPUTS:-1}"

fail() {
  printf 'FAIL: %s\n' "$*" >&2
  exit 1
}

on_exit() {
  local rc=$?
  if ((rc != 0)) && [[ -n "${OUTDIR:-}" && -d "${OUTDIR}" && ! -e "${OUTDIR}/REJECTED.tsv" ]]; then
    {
      printf 'phase\tharness\n'
      printf 'reason\tunhandled_failure\n'
      printf 'exit_status\t%s\n' "${rc}"
      printf 'utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    } >"${OUTDIR}/REJECTED.tsv"
  fi
}
trap on_exit EXIT

require_file() {
  [[ -s "$1" ]] || fail "missing or empty file: $1"
}

check_csv_files() {
  local csv="$1"
  local old_ifs="${IFS}"
  IFS=','
  read -r -a paths <<<"${csv}"
  IFS="${old_ifs}"
  ((${#paths[@]} > 0)) || fail "empty input list"
  local path
  for path in "${paths[@]}"; do
    require_file "${path}"
  done
  printf '%s\n' "${#paths[@]}"
}

csv_sha256() {
  local role="$1"
  local csv="$2"
  local old_ifs="${IFS}"
  IFS=','
  read -r -a paths <<<"${csv}"
  IFS="${old_ifs}"
  local path digest
  for path in "${paths[@]}"; do
    digest="$(sha256sum "${path}" | awk '{print $1}')"
    printf '%s\t%s\t%s\t%s\n' "${role}" "${digest}" "$(stat -c '%s' "${path}")" "${path}"
  done
}

quote_command() {
  printf '%q ' "$@"
  printf '\n'
}

run_timed() {
  local label="$1"
  local temp_dir="$2"
  shift 2
  local time_file="${OUTDIR}/timing/${label}.tsv"
  local stdout_file="${OUTDIR}/logs/${label}.stdout"
  local stderr_file="${OUTDIR}/logs/${label}.stderr"
  local temp_file="${OUTDIR}/timing/${label}.peak_temp_bytes"
  local rc_file="${OUTDIR}/timing/${label}.exit_status"

  mkdir -p "${temp_dir}"
  {
    printf '[%s]\n' "${label}"
    quote_command "$@"
  } >>"${OUTDIR}/COMMANDS.txt"

  /usr/bin/time -f $'wall_seconds\t%e\nuser_seconds\t%U\nsystem_seconds\t%S\nmax_rss_kbytes\t%M\nexit_status\t%x' \
    -o "${time_file}" "$@" >"${stdout_file}" 2>"${stderr_file}" &
  local pid=$!
  local peak=0 current=0
  while kill -0 "${pid}" 2>/dev/null; do
    current="$(du -sb "${temp_dir}" 2>/dev/null | awk '{print $1+0}')"
    if ((current > peak)); then
      peak="${current}"
    fi
    sleep "${TEMP_SAMPLE_SECONDS}"
  done
  set +e
  wait "${pid}"
  local rc=$?
  set -e
  current="$(du -sb "${temp_dir}" 2>/dev/null | awk '{print $1+0}')"
  if ((current > peak)); then
    peak="${current}"
  fi
  printf '%s\n' "${peak}" >"${temp_file}"
  printf '%s\n' "${rc}" >"${rc_file}"
  return "${rc}"
}

metric_from_time() {
  local label="$1"
  local metric="$2"
  awk -F'\t' -v key="${metric}" '$1 == key {print $2; exit}' "${OUTDIR}/timing/${label}.tsv"
}

stream_fragments() {
  local fragments="$1"
  local magic
  magic="$(od -An -tx1 -N2 "${fragments}" | tr -d '[:space:]')"
  if [[ "${magic}" == "1f8b" ]]; then
    gzip -cd "${fragments}"
  else
    command cat "${fragments}"
  fi
}

fragment_stats() {
  local fragments="$1"
  stream_fragments "${fragments}" | awk '
    BEGIN {n=0; weighted=0}
    NF >= 3 && substr($1,1,1) != "#" {
      count=1
      if (NF == 5 || NF == 7) count=$NF+0
      n++
      weighted+=count
    }
    END {printf "%d\t%.0f\n", n, weighted}'
}

frip_stats() {
  local label="$1"
  local fragments="$2"
  local peaks="$3"
  local total_records total_weight overlap_records overlap_weight
  read -r total_records total_weight < <(fragment_stats "${fragments}")
  read -r overlap_records overlap_weight < <(
    bedtools intersect -u -a <(stream_fragments "${fragments}") -b "${peaks}" 2>/dev/null | awk '
      BEGIN {n=0; weighted=0}
      NF >= 3 && substr($1,1,1) != "#" {
        count=1
        if (NF == 5 || NF == 7) count=$NF+0
        n++
        weighted+=count
      }
      END {printf "%d\t%.0f\n", n, weighted}'
  )
  awk -v label="${label}" \
      -v tr="${total_records}" -v tw="${total_weight}" \
      -v ovrec="${overlap_records}" -v ovweight="${overlap_weight}" \
      'BEGIN {
         printf "%s_total_fragment_records\t%d\n", label, tr
         printf "%s_total_fragment_weight\t%.0f\n", label, tw
         printf "%s_peak_overlap_fragment_records\t%d\n", label, ovrec
         printf "%s_peak_overlap_fragment_weight\t%.0f\n", label, ovweight
         printf "%s_frip_unique\t%.12g\n", label, tr ? ovrec/tr : 0
         printf "%s_frip_count_weighted\t%.12g\n", label, tw ? ovweight/tw : 0
       }'
}

main() {
  [[ "${ASSAY_TYPE}" == "bulk" || "${ASSAY_TYPE}" == "scatac" ]] || \
    fail "ASSAY_TYPE must be bulk or scatac"
  [[ ! -e "${OUTDIR}" ]] || fail "OUTDIR already exists; benchmarks require a fresh directory: ${OUTDIR}"
  require_file "${CHROMAP}"
  require_file "${INDEX}"
  require_file "${REF}"
  command -v "${MACS3}" >/dev/null || fail "MACS3 not found: ${MACS3}"
  command -v bedtools >/dev/null || fail "bedtools not found"
  command -v /usr/bin/time >/dev/null || fail "GNU time not found"

  local r1_count r2_count barcode_count=0
  r1_count="$(check_csv_files "${R1_CSV}")"
  r2_count="$(check_csv_files "${R2_CSV}")"
  [[ "${r1_count}" == "${r2_count}" ]] || fail "R1/R2 lane counts differ"
  if [[ "${ASSAY_TYPE}" == "scatac" ]]; then
    [[ -n "${BARCODE_CSV}" ]] || fail "BARCODE_CSV is required for scatac"
    require_file "${WHITELIST}"
    barcode_count="$(check_csv_files "${BARCODE_CSV}")"
    [[ "${barcode_count}" == "${r1_count}" ]] || fail "barcode/genomic lane counts differ"
  fi

  mkdir -p "${OUTDIR}"/{mapping,integrated,standalone,logs,timing,compare,tmp/mapping,tmp/integrated}
  : >"${OUTDIR}/COMMANDS.txt"

  local head rapidmacs_head source_diff_sha source_tree_sha
  local rapidmacs_source_tree_sha chromap_sha
  head="$(git -C "${REPO_ROOT}" rev-parse HEAD)"
  rapidmacs_head="$(git -C "${REPO_ROOT}/third_party/rapidmacs" rev-parse HEAD)"
  source_diff_sha="$(git -C "${REPO_ROOT}" diff -- src | sha256sum | awk '{print $1}')"
  git -C "${REPO_ROOT}" ls-files --cached --others --exclude-standard -z -- src | \
    LC_ALL=C sort -z | \
    (cd "${REPO_ROOT}" && xargs -0 sha256sum) >"${OUTDIR}/SOURCE_MANIFEST.sha256"
  source_tree_sha="$(sha256sum "${OUTDIR}/SOURCE_MANIFEST.sha256" | awk '{print $1}')"
  git -C "${REPO_ROOT}/third_party/rapidmacs" \
    ls-files --cached --others --exclude-standard -z -- include src cli | \
    LC_ALL=C sort -z | \
    (cd "${REPO_ROOT}/third_party/rapidmacs" && xargs -0 sha256sum) \
    >"${OUTDIR}/RAPIDMACS_SOURCE_MANIFEST.sha256"
  rapidmacs_source_tree_sha="$(sha256sum "${OUTDIR}/RAPIDMACS_SOURCE_MANIFEST.sha256" | awk '{print $1}')"
  chromap_sha="$(sha256sum "${CHROMAP}" | awk '{print $1}')"
  {
    printf 'schema\tchromap-suite-atac-paper-benchmark-v1\n'
    printf 'dataset_id\t%s\n' "${DATASET_ID}"
    printf 'assay_type\t%s\n' "${ASSAY_TYPE}"
    printf 'utc_started\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'host\t%s\n' "$(hostname)"
    printf 'kernel\t%s\n' "$(uname -srmo)"
    printf 'cpu\t%s\n' "$(lscpu | awk -F: '$1 ~ /^Model name/ {sub(/^[[:space:]]+/,"",$2); print $2; exit}')"
    printf 'logical_cpus\t%s\n' "$(nproc)"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'repo_head\t%s\n' "${head}"
    printf 'repo_source_diff_sha256\t%s\n' "${source_diff_sha}"
    printf 'repo_source_tree_sha256\t%s\n' "${source_tree_sha}"
    printf 'rapidmacs_head\t%s\n' "${rapidmacs_head}"
    printf 'rapidmacs_source_tree_sha256\t%s\n' "${rapidmacs_source_tree_sha}"
    printf 'chromap_sha256\t%s\n' "${chromap_sha}"
    printf 'chromap_version\t%s\n' "$("${CHROMAP}" --version 2>&1 | tr '\n' ' ')"
    printf 'macs3_version\t%s\n' "$("${MACS3}" --version 2>&1 | tr '\n' ' ')"
    printf 'bedtools_version\t%s\n' "$(bedtools --version 2>&1 | tr '\n' ' ')"
    printf 'reference\t%s\n' "${REF}"
    printf 'index\t%s\n' "${INDEX}"
    printf 'r1_csv\t%s\n' "${R1_CSV}"
    printf 'r2_csv\t%s\n' "${R2_CSV}"
    printf 'barcode_csv\t%s\n' "${BARCODE_CSV}"
    printf 'whitelist\t%s\n' "${WHITELIST}"
    printf 'macs3_pvalue\t%s\n' "${MACS3_PVALUE}"
    printf 'macs3_min_length\t%s\n' "${MACS3_MIN_LENGTH}"
    printf 'macs3_max_gap\t%s\n' "${MACS3_MAX_GAP}"
    printf 'tn5_shift_mode\tsymmetric\n'
    printf 'chromap_low_mem_ram\t%s\n' "${LOW_MEM_RAM}"
    printf 'integrated_peak_workspace\tlow_mem_sweep\n'
  } >"${OUTDIR}/PROVENANCE.tsv"
  git -C "${REPO_ROOT}" status --porcelain=v1 >"${OUTDIR}/REPO_STATUS.txt"
  git -C "${REPO_ROOT}/third_party/rapidmacs" status --porcelain=v1 \
    >"${OUTDIR}/RAPIDMACS_REPO_STATUS.txt"
  ldd "${CHROMAP}" >"${OUTDIR}/CHROMAP_LDD.txt"

  if [[ "${HASH_INPUTS}" == "1" ]]; then
    {
      printf 'role\tsha256\tbytes\tpath\n'
      csv_sha256 r1 "${R1_CSV}"
      csv_sha256 r2 "${R2_CSV}"
      if [[ "${ASSAY_TYPE}" == "scatac" ]]; then
        csv_sha256 barcode "${BARCODE_CSV}"
        printf 'whitelist\t%s\t%s\t%s\n' \
          "$(sha256sum "${WHITELIST}" | awk '{print $1}')" \
          "$(stat -c '%s' "${WHITELIST}")" "${WHITELIST}"
      fi
    } >"${OUTDIR}/INPUTS.sha256.tsv"
  else
    printf 'hashing\tdisabled\n' >"${OUTDIR}/INPUTS.sha256.tsv"
  fi

  local -a common
  common=(
    -t "${THREADS}" -x "${INDEX}" -r "${REF}"
    -1 "${R1_CSV}" -2 "${R2_CSV}"
    --preset atac -l 2000 --trim-adapters --remove-pcr-duplicates
    --Tn5-shift-mode symmetric --low-mem --low-mem-ram "${LOW_MEM_RAM}" --BED
  )
  if [[ "${ASSAY_TYPE}" == "scatac" ]]; then
    common+=(
      -b "${BARCODE_CSV}" --barcode-whitelist "${WHITELIST}"
      --remove-pcr-duplicates-at-cell-level
    )
  else
    common+=(--remove-pcr-duplicates-at-bulk-level)
  fi

  if ! run_timed mapping "${OUTDIR}/tmp/mapping" \
      "${CHROMAP}" "${common[@]}" \
      --temp-dir "${OUTDIR}/tmp/mapping" \
      --summary "${OUTDIR}/mapping/summary.tsv" \
      -o "${OUTDIR}/mapping/fragments.tsv"; then
    printf 'phase\tmapping\nreason\tcommand_failed\n' >"${OUTDIR}/REJECTED.tsv"
    fail "Chromap mapping-only run failed"
  fi
  require_file "${OUTDIR}/mapping/fragments.tsv"

  if ! run_timed integrated "${OUTDIR}/tmp/integrated" \
      "${CHROMAP}" "${common[@]}" \
      --temp-dir "${OUTDIR}/tmp/integrated" \
      --summary "${OUTDIR}/integrated/summary.tsv" \
      -o "${OUTDIR}/integrated/fragments.tsv" \
      --call-macs3-frag-peaks \
      --macs3-frag-peaks-source memory \
      --macs3-frag-low-mem \
      --macs3-frag-peaks-output "${OUTDIR}/integrated/rapidmacs_peaks.narrowPeak" \
      --macs3-frag-summits-output "${OUTDIR}/integrated/rapidmacs_summits.bed" \
      --macs3-frag-pvalue "${MACS3_PVALUE}" \
      --macs3-frag-min-length "${MACS3_MIN_LENGTH}" \
      --macs3-frag-max-gap "${MACS3_MAX_GAP}"; then
    printf 'phase\tintegrated\nreason\tcommand_failed\n' >"${OUTDIR}/REJECTED.tsv"
    fail "integrated Chromap + RapidMACS run failed"
  fi
  require_file "${OUTDIR}/integrated/fragments.tsv"
  require_file "${OUTDIR}/integrated/rapidmacs_peaks.narrowPeak"
  require_file "${OUTDIR}/integrated/rapidmacs_summits.bed"

  if ! cmp -s <(stream_fragments "${OUTDIR}/mapping/fragments.tsv") \
              <(stream_fragments "${OUTDIR}/integrated/fragments.tsv"); then
    printf 'phase\tfragment_parity\nreason\tdecompressed_streams_differ\n' >"${OUTDIR}/REJECTED.tsv"
    fail "mapping-only and integrated fragments differ"
  fi
  {
    printf 'comparison\tPASS\n'
    printf 'mapping_file_sha256\t%s\n' \
      "$(sha256sum "${OUTDIR}/mapping/fragments.tsv" | awk '{print $1}')"
    printf 'integrated_file_sha256\t%s\n' \
      "$(sha256sum "${OUTDIR}/integrated/fragments.tsv" | awk '{print $1}')"
    printf 'fragment_stream_sha256\t%s\n' \
      "$(stream_fragments "${OUTDIR}/mapping/fragments.tsv" | sha256sum | awk '{print $1}')"
  } >"${OUTDIR}/compare/FRAGMENT_PARITY.tsv"

  local macs_input="${OUTDIR}/mapping/fragments.tsv"
  if [[ "${ASSAY_TYPE}" == "bulk" ]]; then
    macs_input="${OUTDIR}/standalone/mapping_fragments.frag.tsv.gz"
    stream_fragments "${OUTDIR}/mapping/fragments.tsv" | \
      awk 'BEGIN {OFS="\t"} NF >= 7 {print $1,$2,$3,"bulk",$7}' | \
      gzip -n >"${macs_input}"
    require_file "${macs_input}"
  fi

  if ! run_timed standalone_macs3 "${OUTDIR}/standalone" \
      "${MACS3}" callpeak -t "${macs_input}" -f FRAG -g hs \
      -n NA -p "${MACS3_PVALUE}" \
      --min-length "${MACS3_MIN_LENGTH}" --max-gap "${MACS3_MAX_GAP}" \
      --outdir "${OUTDIR}/standalone"; then
    printf 'phase\tstandalone_macs3\nreason\tcommand_failed\n' >"${OUTDIR}/REJECTED.tsv"
    fail "standalone MACS3 failed"
  fi
  local macs_peaks="${OUTDIR}/standalone/NA_peaks.narrowPeak"
  local macs_summits="${OUTDIR}/standalone/NA_summits.bed"
  require_file "${macs_peaks}"
  require_file "${macs_summits}"

  local narrowpeak_exact=0 summits_exact=0
  cmp -s "${macs_peaks}" \
    "${OUTDIR}/integrated/rapidmacs_peaks.narrowPeak" && narrowpeak_exact=1
  cmp -s "${macs_summits}" \
    "${OUTDIR}/integrated/rapidmacs_summits.bed" && summits_exact=1
  {
    printf 'artifact\texact\tmacs3_sha256\trapidmacs_sha256\n'
    printf 'narrowPeak\t%s\t%s\t%s\n' \
      "${narrowpeak_exact}" \
      "$(sha256sum "${macs_peaks}" | awk '{print $1}')" \
      "$(sha256sum "${OUTDIR}/integrated/rapidmacs_peaks.narrowPeak" | awk '{print $1}')"
    printf 'summits\t%s\t%s\t%s\n' \
      "${summits_exact}" \
      "$(sha256sum "${macs_summits}" | awk '{print $1}')" \
      "$(sha256sum "${OUTDIR}/integrated/rapidmacs_summits.bed" | awk '{print $1}')"
  } >"${OUTDIR}/compare/PEAK_FILE_PARITY.tsv"

  python3 "${SCRIPT_DIR}/compare_narrowpeak_parity.py" \
    --macs3-np "${macs_peaks}" \
    --cpp-np "${OUTDIR}/integrated/rapidmacs_peaks.narrowPeak" \
    --macs3-summits "${macs_summits}" \
    --cpp-summits "${OUTDIR}/integrated/rapidmacs_summits.bed" \
    --out-tsv "${OUTDIR}/compare/PEAK_SUMMIT_PARITY.tsv" \
    >"${OUTDIR}/logs/compare_narrowpeak.stdout"

  awk 'BEGIN {OFS="\t"} NF >= 3 && substr($1,1,1) != "#" {print $1,$2,$3}' \
    "${macs_peaks}" | LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
    >"${OUTDIR}/compare/macs3.peaks.bed3"
  awk 'BEGIN {OFS="\t"} NF >= 3 && substr($1,1,1) != "#" {print $1,$2,$3}' \
    "${OUTDIR}/integrated/rapidmacs_peaks.narrowPeak" | \
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
    >"${OUTDIR}/compare/rapidmacs.peaks.bed3"
  python3 "${SCRIPT_DIR}/compare_bdgpeakcall_regions.py" \
    "${OUTDIR}/compare/macs3.peaks.bed3" \
    "${OUTDIR}/compare/rapidmacs.peaks.bed3" \
    --label-a macs3 --label-b rapidmacs \
    >"${OUTDIR}/compare/PEAK_OVERLAP.tsv"

  local n_macs n_rapid r50_macs r50_rapid
  n_macs="$(wc -l <"${OUTDIR}/compare/macs3.peaks.bed3")"
  n_rapid="$(wc -l <"${OUTDIR}/compare/rapidmacs.peaks.bed3")"
  r50_macs="$(bedtools intersect -u -f 0.5 -r \
    -a "${OUTDIR}/compare/macs3.peaks.bed3" \
    -b "${OUTDIR}/compare/rapidmacs.peaks.bed3" 2>/dev/null | wc -l)"
  r50_rapid="$(bedtools intersect -u -f 0.5 -r \
    -a "${OUTDIR}/compare/rapidmacs.peaks.bed3" \
    -b "${OUTDIR}/compare/macs3.peaks.bed3" 2>/dev/null | wc -l)"
  {
    printf 'macs3_reciprocal_50pct_overlap_count\t%s\n' "${r50_macs}"
    printf 'macs3_reciprocal_50pct_overlap_fraction\t%s\n' \
      "$(awk -v n="${n_macs}" -v k="${r50_macs}" 'BEGIN {printf "%.12g", n ? k/n : 0}')"
    printf 'rapidmacs_reciprocal_50pct_overlap_count\t%s\n' "${r50_rapid}"
    printf 'rapidmacs_reciprocal_50pct_overlap_fraction\t%s\n' \
      "$(awk -v n="${n_rapid}" -v k="${r50_rapid}" 'BEGIN {printf "%.12g", n ? k/n : 0}')"
  } >"${OUTDIR}/compare/RECIPROCAL_50_OVERLAP.tsv"

  {
    frip_stats macs3 "${OUTDIR}/mapping/fragments.tsv" "${macs_peaks}"
    frip_stats rapidmacs "${OUTDIR}/mapping/fragments.tsv" \
      "${OUTDIR}/integrated/rapidmacs_peaks.narrowPeak"
  } >"${OUTDIR}/compare/FRIP.tsv"

  local map_wall int_wall macs_wall
  map_wall="$(metric_from_time mapping wall_seconds)"
  int_wall="$(metric_from_time integrated wall_seconds)"
  macs_wall="$(metric_from_time standalone_macs3 wall_seconds)"
  {
    printf 'metric\tvalue\n'
    printf 'dataset_id\t%s\n' "${DATASET_ID}"
    printf 'assay_type\t%s\n' "${ASSAY_TYPE}"
    printf 'status\tPASS\n'
    printf 'fragment_exact_decompressed\t1\n'
    printf 'narrowpeak_file_exact\t%s\n' "${narrowpeak_exact}"
    printf 'summits_file_exact\t%s\n' "${summits_exact}"
    printf 'mapping_wall_seconds\t%s\n' "${map_wall}"
    printf 'mapping_user_seconds\t%s\n' "$(metric_from_time mapping user_seconds)"
    printf 'mapping_system_seconds\t%s\n' "$(metric_from_time mapping system_seconds)"
    printf 'mapping_max_rss_kbytes\t%s\n' "$(metric_from_time mapping max_rss_kbytes)"
    printf 'integrated_wall_seconds\t%s\n' "${int_wall}"
    printf 'integrated_user_seconds\t%s\n' "$(metric_from_time integrated user_seconds)"
    printf 'integrated_system_seconds\t%s\n' "$(metric_from_time integrated system_seconds)"
    printf 'integrated_max_rss_kbytes\t%s\n' "$(metric_from_time integrated max_rss_kbytes)"
    printf 'standalone_macs3_wall_seconds\t%s\n' "${macs_wall}"
    printf 'standalone_macs3_user_seconds\t%s\n' "$(metric_from_time standalone_macs3 user_seconds)"
    printf 'standalone_macs3_system_seconds\t%s\n' "$(metric_from_time standalone_macs3 system_seconds)"
    printf 'standalone_macs3_max_rss_kbytes\t%s\n' "$(metric_from_time standalone_macs3 max_rss_kbytes)"
    printf 'standalone_workflow_wall_seconds\t%s\n' \
      "$(awk -v a="${map_wall}" -v b="${macs_wall}" 'BEGIN {printf "%.6f", a+b}')"
    printf 'integrated_incremental_peak_wall_seconds\t%s\n' \
      "$(awk -v a="${int_wall}" -v b="${map_wall}" 'BEGIN {printf "%.6f", a-b}')"
    printf 'mapping_peak_temp_bytes_sampled\t%s\n' "$(<"${OUTDIR}/timing/mapping.peak_temp_bytes")"
    printf 'integrated_peak_temp_bytes_sampled\t%s\n' "$(<"${OUTDIR}/timing/integrated.peak_temp_bytes")"
    printf 'mapping_output_bytes\t%s\n' "$(du -sb "${OUTDIR}/mapping" | awk '{print $1}')"
    printf 'integrated_output_bytes\t%s\n' "$(du -sb "${OUTDIR}/integrated" | awk '{print $1}')"
    printf 'standalone_output_bytes\t%s\n' "$(du -sb "${OUTDIR}/standalone" | awk '{print $1}')"
    printf 'mapping_fragment_bytes\t%s\n' "$(stat -c '%s' "${OUTDIR}/mapping/fragments.tsv")"
    printf 'integrated_fragment_bytes\t%s\n' "$(stat -c '%s' "${OUTDIR}/integrated/fragments.tsv")"
    awk 'NR > 1 {
      printf "peak_parity_%s", $1
      for (i=2; i<=NF; i++) printf "\t%s", $i
      printf "\n"
    }' "${OUTDIR}/compare/PEAK_SUMMIT_PARITY.tsv"
    cat "${OUTDIR}/compare/RECIPROCAL_50_OVERLAP.tsv"
    cat "${OUTDIR}/compare/FRIP.tsv"
    printf 'utc_completed\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  } >"${OUTDIR}/RESULTS.tsv"
  sha256sum "${OUTDIR}/RESULTS.tsv" "${OUTDIR}/PROVENANCE.tsv" \
    "${OUTDIR}/SOURCE_MANIFEST.sha256" \
    "${OUTDIR}/RAPIDMACS_SOURCE_MANIFEST.sha256" \
    "${OUTDIR}/REPO_STATUS.txt" \
    "${OUTDIR}/RAPIDMACS_REPO_STATUS.txt" \
    "${OUTDIR}/INPUTS.sha256.tsv" \
    "${OUTDIR}/compare/FRAGMENT_PARITY.tsv" \
    "${OUTDIR}/compare/PEAK_FILE_PARITY.tsv" \
    "${OUTDIR}/compare/PEAK_SUMMIT_PARITY.tsv" \
    "${OUTDIR}/compare/RECIPROCAL_50_OVERLAP.tsv" \
    "${OUTDIR}/compare/FRIP.tsv" \
    >"${OUTDIR}/ACCEPTED.sha256"
  printf 'PASS: %s\n' "${OUTDIR}"
  cat "${OUTDIR}/RESULTS.tsv"
}

main "$@"
