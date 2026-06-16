#!/usr/bin/env bash
# Full PBMC 3K ATAC benchmark for the Chromap-suite libMACS3 integration.
#
# Baseline:
#   unfixed upstream Chromap 0.3.3 emits SAM to stdout, piped through samtools
#   view/sort/index for a coordinate-sorted BAM. The same upstream Chromap
#   binary emits BED/FRAG rows for standalone MACS3 callpeak.
#
# Integrated:
#   Chromap-suite emits sorted BAM, indexed BAM, fragments, and libMACS3 FRAG
#   peaks in one command.
#
# The runner executes modes serially. RUN_MODES=normal runs without --low-mem,
# RUN_MODES=lowmem adds --low-mem to both baseline and integrated Chromap
# commands, and RUN_MODES=both runs normal first, then lowmem.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# shellcheck source=scripts/benchmarks/libmacs3_chromap_common.sh
source "${SCRIPT_DIR}/libmacs3_chromap_common.sh"

export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

BENCH_ROOT="${BENCH_ROOT:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k}"
FIXTURE_ATAC="${FIXTURE_ATAC:-${BENCH_ROOT}/extracted/pbmc_unsorted_3k/atac}"
INDEX="${INDEX:-${BENCH_ROOT}/chromap_index/genome.index}"
REF="${REF:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa}"
WHITELIST="${WHITELIST:-${BENCH_ROOT}/chromap_index/737K-arc-v1_atac.txt}"

DEFAULT_UNFIXED_CHROMAP="${REPO_ROOT}/plans/artifacts/upstream_chromap_unfixed/c9d8ae058bfdf04d45bc5f99a164a460f759a6a7/chromap"
CHROMAP_ORIGINAL="${CHROMAP_ORIGINAL:-${DEFAULT_UNFIXED_CHROMAP}}"
CHROMAP_SUITE="${CHROMAP_SUITE:-${REPO_ROOT}/chromap}"
MACS3_BIN="${MACS3_BIN:-macs3}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

THREADS="${THREADS:-16}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-16}"
RUN_MODES="${RUN_MODES:-both}"
SKIP_MAKE="${SKIP_MAKE:-0}"
BENCHMARK_NOTE="${BENCHMARK_NOTE:-}"
RESUME_EXISTING="${RESUME_EXISTING:-0}"
RESUME_COMPARISONS="${RESUME_COMPARISONS:-0}"

MACS3_PVALUE="${MACS3_PVALUE:-1e-5}"
MACS3_QVALUE="${MACS3_QVALUE:-}"
GENOME_SIZE="${GENOME_SIZE:-hs}"
MACS3_MIN_LENGTH="${MACS3_MIN_LENGTH:-200}"
MACS3_MAX_GAP="${MACS3_MAX_GAP:-30}"

MACS3_THRESHOLD_MODE="pvalue"
macs3_threshold_args=(-p "${MACS3_PVALUE}")
chromap_suite_threshold_args=(--macs3-frag-pvalue "${MACS3_PVALUE}")
if [[ -n "${MACS3_QVALUE}" ]]; then
  MACS3_THRESHOLD_MODE="qvalue"
  macs3_threshold_args=(-q "${MACS3_QVALUE}")
  chromap_suite_threshold_args=(--macs3-frag-qvalue "${MACS3_QVALUE}")
fi

ARTIFACT_ROOT="${ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts/libmacs3_chromap_atac_panel}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
OUTDIR="${OUTDIR:-${ARTIFACT_ROOT}/${RUN_ID}/full_pbmc3k}"

R1="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R1_001.fastq.gz"
R2="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R3_001.fastq.gz"
BC="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R2_001.fastq.gz"

mkdir -p "${OUTDIR}/logs" "${OUTDIR}/times"
COMMANDS="${OUTDIR}/commands.sh"
RUN_LOG="${OUTDIR}/logs/run.log"
VERSIONS="${OUTDIR}/versions.txt"
MANIFEST="${OUTDIR}/manifest.tsv"
SUMMARY="${OUTDIR}/full_summary.tsv"

{
  printf '#!/usr/bin/env bash\n'
  printf 'set -euo pipefail\n\n'
  printf '# Re-run from repo root: %s\n' "${REPO_ROOT}"
  printf '# Output directory from this run: %s\n' "${OUTDIR}"
} >"${COMMANDS}"
chmod +x "${COMMANDS}"

log_msg() {
  printf '%s\n' "$*" | tee -a "${RUN_LOG}"
}

quote_cmd_inline() {
  printf '%q ' "$@"
}

shell_quote_one() {
  printf '%q' "$1"
}

sort_fragments_5col_cmd() {
  local input=$1
  local output=$2
  printf 'source %s; bench_sort_fragments_5col %s %s' \
    "$(shell_quote_one "${SCRIPT_DIR}/libmacs3_chromap_common.sh")" \
    "$(shell_quote_one "${input}")" \
    "$(shell_quote_one "${output}")"
}

run_timed_try() {
  local label=$1
  local time_file=$2
  local log_file=$3
  local commands_file=$4
  shift 4
  {
    printf '\n# [%s]\n' "${label}"
    bench_quote_cmd "$@"
  } >>"${commands_file}"
  set +e
  /usr/bin/time -v -o "${time_file}" "$@" >>"${log_file}" 2>&1
  local rc=$?
  set -e
  printf 'Codex benchmark exit status: %s\n' "${rc}" >>"${time_file}"
  return "${rc}"
}

run_timed_shell_try() {
  local label=$1
  local time_file=$2
  local log_file=$3
  local commands_file=$4
  local command_text=$5
  {
    printf '\n# [%s]\n' "${label}"
    printf '%s\n' "${command_text}"
  } >>"${commands_file}"
  set +e
  /usr/bin/time -v -o "${time_file}" bash -o pipefail -c "${command_text}" >>"${log_file}" 2>&1
  local rc=$?
  set -e
  printf 'Codex benchmark exit status: %s\n' "${rc}" >>"${time_file}"
  return "${rc}"
}

time_field() {
  local input=$1
  local key=$2
  if [[ ! -s "${input}" ]]; then
    printf 'na'
    return
  fi
  awk -v key="${key}" '
    {
      line = $0
      sub(/^[ \t]+/, "", line)
    }
    index(line, key ":") == 1 {
      value = substr(line, length(key) + 2)
      sub(/^[ \t]+/, "", value)
      print value
      found = 1
      exit
    }
    END {
      if (!found) {
        print "na"
      }
    }
  ' "${input}"
}

write_time_summary() {
  local mode_dir=$1
  local output=$2
  local time_file step
  {
    printf 'step\telapsed\tuser_seconds\tsys_seconds\tmax_rss_kbytes\tfs_inputs\tfs_outputs\texit_status\n'
    for time_file in "${mode_dir}/times/"*.time.txt; do
      [[ -f "${time_file}" ]] || continue
      step="$(basename "${time_file}" .time.txt)"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${step}" \
        "$(time_field "${time_file}" "Elapsed (wall clock) time (h:mm:ss or m:ss)")" \
        "$(time_field "${time_file}" "User time (seconds)")" \
        "$(time_field "${time_file}" "System time (seconds)")" \
        "$(time_field "${time_file}" "Maximum resident set size (kbytes)")" \
        "$(time_field "${time_file}" "File system inputs")" \
        "$(time_field "${time_file}" "File system outputs")" \
        "$(time_field "${time_file}" "Exit status")"
    done
  } >"${output}"
}

count_lines_or_na() {
  local input=$1
  if [[ -s "${input}" ]]; then
    wc -l <"${input}" | awk '{print $1}'
  else
    printf 'na'
  fi
}

count_noncomment_or_na() {
  local input=$1
  if [[ -s "${input}" ]]; then
    bench_count_noncomment_lines "${input}"
  else
    printf 'na'
  fi
}

md5_or_na() {
  local input=$1
  if [[ -s "${input}" ]]; then
    bench_md5_file "${input}"
  else
    printf 'na'
  fi
}

cmp_bool_or_na() {
  local a=$1
  local b=$2
  if [[ ! -s "${a}" || ! -s "${b}" ]]; then
    printf 'na'
  elif cmp -s "${a}" "${b}"; then
    printf 'true'
  else
    printf 'false'
  fi
}

line_set_counts_sorted_or_na() {
  local a=$1
  local b=$2
  if [[ ! -s "${a}" || ! -s "${b}" ]]; then
    printf 'na\tna\tna\n'
    return
  fi
  local a_only b_only common
  a_only="$(comm -23 "${a}" "${b}" | wc -l | awk '{print $1}')"
  b_only="$(comm -13 "${a}" "${b}" | wc -l | awk '{print $1}')"
  common="$(comm -12 "${a}" "${b}" | wc -l | awk '{print $1}')"
  printf '%s\t%s\t%s\n' "${a_only}" "${b_only}" "${common}"
}

bam_records_or_na() {
  local bam=$1
  if [[ -s "${bam}" ]]; then
    bench_bam_total_records "${SAMTOOLS_BIN}" "${bam}"
  else
    printf 'na'
  fi
}

summit_distance_stats_or_na() {
  local a=$1
  local b=$2
  if [[ ! -s "${a}" || ! -s "${b}" ]] || ! command -v bedtools >/dev/null 2>&1; then
    printf 'na\tna\tna\tna\n'
    return
  fi
  bedtools closest -d -a "${a}" -b "${b}" 2>/dev/null |
    awk '
      {
        d = $NF + 0
        n++
        sum += d
        if (d == 0) {
          zero++
        }
        if (n == 1 || d > max) {
          max = d
        }
      }
      END {
        if (n == 0) {
          print "na\tna\tna\tna"
        } else {
          printf "%.6f\t%d\t%d\t%d\n", sum / n, max, zero + 0, n
        }
      }'
}

write_static_metadata() {
  {
    printf 'run_id\t%s\n' "${RUN_ID}"
    printf 'outdir\t%s\n' "${OUTDIR}"
    printf 'bench_root\t%s\n' "${BENCH_ROOT}"
    printf 'fixture_atac\t%s\n' "${FIXTURE_ATAC}"
    printf 'index\t%s\n' "${INDEX}"
    printf 'ref\t%s\n' "${REF}"
    printf 'whitelist\t%s\n' "${WHITELIST}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'samtools_threads\t%s\n' "${SAMTOOLS_THREADS}"
    printf 'run_modes\t%s\n' "${RUN_MODES}"
    printf 'benchmark_note\t%s\n' "${BENCHMARK_NOTE}"
    printf 'resume_existing\t%s\n' "${RESUME_EXISTING}"
    printf 'resume_comparisons\t%s\n' "${RESUME_COMPARISONS}"
    printf 'macs3_threshold_mode\t%s\n' "${MACS3_THRESHOLD_MODE}"
    printf 'macs3_pvalue\t%s\n' "${MACS3_PVALUE}"
    printf 'macs3_qvalue\t%s\n' "${MACS3_QVALUE:-na}"
    printf 'macs3_genome_size\t%s\n' "${GENOME_SIZE}"
    printf 'macs3_min_length\t%s\n' "${MACS3_MIN_LENGTH}"
    printf 'macs3_max_gap\t%s\n' "${MACS3_MAX_GAP}"
    printf 'mapping_flags\t%s\n' "explicit_atac_flags_no_preset"
    printf 'chromap_original_unfixed\t%s\n' "${CHROMAP_ORIGINAL}"
    printf 'chromap_suite\t%s\n' "${CHROMAP_SUITE}"
    printf 'macs3\t%s\n' "$(command -v "${MACS3_BIN}")"
    printf 'samtools\t%s\n' "$(command -v "${SAMTOOLS_BIN}")"
  } >"${MANIFEST}"

  {
    printf 'Original unfixed Chromap version: '
    "${CHROMAP_ORIGINAL}" --version 2>&1 || true
    printf 'Chromap-suite version: '
    "${CHROMAP_SUITE}" --version 2>&1 || true
    "${MACS3_BIN}" --version 2>&1 || true
    "${SAMTOOLS_BIN}" --version 2>&1 | sed -n '1,3p' || true
    printf 'Original unfixed Chromap sha256: %s\n' "$(sha256sum "${CHROMAP_ORIGINAL}" | awk '{print $1}')"
    printf 'Chromap-suite sha256: %s\n' "$(sha256sum "${CHROMAP_SUITE}" | awk '{print $1}')"
    printf '\n'
  } >"${VERSIONS}"
  bench_write_git_state "${VERSIONS}" "original unfixed Chromap repo" "$(cd "$(dirname "${CHROMAP_ORIGINAL}")" && pwd)"
  bench_write_git_state "${VERSIONS}" "Chromap-suite" "${REPO_ROOT}"
}

prepare_inputs_and_tools() {
  for path in "${INDEX}" "${REF}" "${WHITELIST}"; do
    bench_require_file "${path}" "benchmark input"
  done
  for lane in L001 L002 L003 L004; do
    bench_require_file "${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_${lane}_R1_001.fastq.gz" "R1 ${lane}"
    bench_require_file "${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_${lane}_R2_001.fastq.gz" "barcode ${lane}"
    bench_require_file "${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_${lane}_R3_001.fastq.gz" "R2 ${lane}"
  done
  if [[ ! -x "${CHROMAP_ORIGINAL}" && "${CHROMAP_ORIGINAL}" == "${DEFAULT_UNFIXED_CHROMAP}" ]]; then
    log_msg "[full-benchmark] preparing unfixed upstream Chromap baseline"
    CHROMAP_ORIGINAL="$("${SCRIPT_DIR}/prepare_unfixed_chromap.sh")"
  fi
  bench_require_exe "${CHROMAP_ORIGINAL}" "unfixed upstream Chromap"
  bench_require_exe "${MACS3_BIN}" "MACS3"
  bench_require_exe "${SAMTOOLS_BIN}" "samtools"

  if [[ "${SKIP_MAKE}" != "1" ]]; then
    log_msg "[full-benchmark] building Chromap-suite binaries"
    if ! run_timed_try "build Chromap-suite chromap + chromap_callpeaks" \
      "${OUTDIR}/times/build.time.txt" "${RUN_LOG}" "${COMMANDS}" \
      make -C "${REPO_ROOT}" chromap chromap_callpeaks; then
      echo "ERROR: build failed; see ${RUN_LOG}" >&2
      exit 1
    fi
  fi
  bench_require_exe "${CHROMAP_SUITE}" "Chromap-suite"
}

run_mode() {
  local mode=$1
  local low_mem_flag=$2
  local mode_dir="${OUTDIR}/${mode}"
  local original_dir="${mode_dir}/original_chromap"
  local suite_dir="${mode_dir}/chromap_suite"
  local log="${mode_dir}/logs/run.log"
  local status="PASS"
  local failures=()

  mkdir -p \
    "${original_dir}/macs3" \
    "${original_dir}/compare" \
    "${suite_dir}/compare" \
    "${mode_dir}/logs" \
    "${mode_dir}/times" \
    "${mode_dir}/tmp" \
    "${mode_dir}/tmp/chromap_suite"

  export TMPDIR="${mode_dir}/tmp"

  local low_mem_args=()
  if [[ "${low_mem_flag}" == "1" ]]; then
    low_mem_args=(--low-mem)
  fi

  local chromap_common_args=(
    -x "${INDEX}"
    -r "${REF}"
    -1 "${R1}"
    -2 "${R2}"
    -b "${BC}"
    --barcode-whitelist "${WHITELIST}"
    -l 2000
    --trim-adapters
    --remove-pcr-duplicates
    --remove-pcr-duplicates-at-cell-level
    --Tn5-shift
  )

  log_msg "[full-benchmark:${mode}] original unfixed Chromap baseline: SAM -> samtools sorted BAM"
  local original_sam_cmd
  original_sam_cmd="$(quote_cmd_inline "${CHROMAP_ORIGINAL}" -t "${THREADS}" "${chromap_common_args[@]}" "${low_mem_args[@]}" --SAM --summary "${original_dir}/summary.sam.tsv" -o /dev/stdout)"
  original_sam_cmd+=" | $(quote_cmd_inline "${SAMTOOLS_BIN}" view -@ "${SAMTOOLS_THREADS}" -b -)"
  original_sam_cmd+=" | $(quote_cmd_inline "${SAMTOOLS_BIN}" sort -@ "${SAMTOOLS_THREADS}" -o "${original_dir}/possorted.bam" -)"
  if [[ "${RESUME_EXISTING}" == "1" && -s "${original_dir}/possorted.bam" ]]; then
    log_msg "[full-benchmark:${mode}] reusing existing original baseline BAM"
  else
    rm -f "${original_dir}/possorted.bam" "${original_dir}/possorted.bam.bai"
    if ! run_timed_shell_try "original unfixed Chromap ${mode} SAM piped to samtools sorted BAM" \
      "${mode_dir}/times/original_sam_to_bam.time.txt" "${log}" "${COMMANDS}" "${original_sam_cmd}"; then
      status="FAIL"
      failures+=("original_sam_to_bam")
    fi
  fi
  if [[ -s "${original_dir}/possorted.bam" && ! -s "${original_dir}/possorted.bam.bai" ]]; then
    if ! run_timed_try "index original ${mode} baseline BAM" \
      "${mode_dir}/times/original_bam_index.time.txt" "${log}" "${COMMANDS}" \
      "${SAMTOOLS_BIN}" index "${original_dir}/possorted.bam"; then
      status="FAIL"
      failures+=("original_bam_index")
    fi
  fi

  log_msg "[full-benchmark:${mode}] original unfixed Chromap baseline: BED/FRAG rows"
  if [[ "${RESUME_EXISTING}" == "1" && -s "${original_dir}/fragments.tsv" ]]; then
    log_msg "[full-benchmark:${mode}] reusing existing original baseline fragments"
  else
    rm -f "${original_dir}/fragments.tsv"
    if ! run_timed_try "original unfixed Chromap ${mode} BED/FRAG rows" \
      "${mode_dir}/times/original_fragments.time.txt" "${log}" "${COMMANDS}" \
      "${CHROMAP_ORIGINAL}" -t "${THREADS}" "${chromap_common_args[@]}" "${low_mem_args[@]}" \
      --BED --summary "${original_dir}/summary.fragments.tsv" \
      -o "${original_dir}/fragments.tsv"; then
      status="FAIL"
      failures+=("original_fragments")
    fi
  fi

  if [[ -s "${original_dir}/fragments.tsv" ]]; then
    log_msg "[full-benchmark:${mode}] original unfixed Chromap baseline: MACS3 FRAG peaks"
    if [[ "${RESUME_EXISTING}" == "1" &&
          -s "${original_dir}/macs3/original_chromap_macs3_peaks.narrowPeak" &&
          -s "${original_dir}/macs3/original_chromap_macs3_summits.bed" ]]; then
      log_msg "[full-benchmark:${mode}] reusing existing original baseline MACS3 peaks"
    else
      rm -f "${original_dir}/macs3"/original_chromap_macs3_*
      if ! run_timed_try "MACS3 callpeak on original ${mode} Chromap fragments" \
        "${mode_dir}/times/original_macs3.time.txt" "${log}" "${COMMANDS}" \
        "${MACS3_BIN}" callpeak \
        -t "${original_dir}/fragments.tsv" \
        -f FRAG -g "${GENOME_SIZE}" \
        -n original_chromap_macs3 \
        "${macs3_threshold_args[@]}" \
        --min-length "${MACS3_MIN_LENGTH}" \
        --max-gap "${MACS3_MAX_GAP}" \
        --outdir "${original_dir}/macs3"; then
        status="FAIL"
        failures+=("original_macs3")
      fi
    fi
  fi

  log_msg "[full-benchmark:${mode}] Chromap-suite integrated command"
  local suite_cmd=(
    "${CHROMAP_SUITE}" -t "${THREADS}" "${chromap_common_args[@]}" "${low_mem_args[@]}"
    --BAM --sort-bam --write-index
    --atac-fragments "${suite_dir}/fragments.tsv.gz"
    --summary "${suite_dir}/summary.tsv"
    --call-macs3-frag-peaks
    --macs3-frag-peaks-output "${suite_dir}/chromap_suite_libmacs3_peaks.narrowPeak"
    --macs3-frag-summits-output "${suite_dir}/chromap_suite_libmacs3_summits.bed"
    "${chromap_suite_threshold_args[@]}"
    -o "${suite_dir}/possorted.bam"
  )
  if "${CHROMAP_SUITE}" --help 2>&1 | grep -q -- '--temp-dir'; then
    suite_cmd+=(--temp-dir "${mode_dir}/tmp/chromap_suite")
  fi
  if [[ "${RESUME_EXISTING}" == "1" &&
        -s "${suite_dir}/possorted.bam" &&
        -s "${suite_dir}/possorted.bam.bai" &&
        -s "${suite_dir}/fragments.tsv.gz" &&
        -s "${suite_dir}/chromap_suite_libmacs3_peaks.narrowPeak" &&
        -s "${suite_dir}/chromap_suite_libmacs3_summits.bed" ]]; then
    log_msg "[full-benchmark:${mode}] reusing existing Chromap-suite integrated outputs"
  else
    rm -f \
      "${suite_dir}/possorted.bam" \
      "${suite_dir}/possorted.bam.bai" \
      "${suite_dir}/fragments.tsv.gz" \
      "${suite_dir}/chromap_suite_libmacs3_peaks.narrowPeak" \
      "${suite_dir}/chromap_suite_libmacs3_summits.bed"
    rm -f "${mode_dir}/tmp/chromap_suite"/chromap_sort_spill_*.dat 2>/dev/null || true
    if ! run_timed_try "Chromap-suite ${mode} integrated BAM + fragments + libMACS3 peaks" \
      "${mode_dir}/times/chromap_suite_integrated.time.txt" "${log}" "${COMMANDS}" \
      "${suite_cmd[@]}"; then
      status="FAIL"
      failures+=("chromap_suite_integrated")
    fi
  fi

  if [[ -s "${original_dir}/fragments.tsv" ]]; then
    if [[ "${RESUME_COMPARISONS}" == "1" && -s "${original_dir}/compare/fragments.5col.sorted.tsv" ]]; then
      log_msg "[full-benchmark:${mode}] reusing sorted original fragments for identity check"
    else
      log_msg "[full-benchmark:${mode}] sorting original fragments for identity check"
      if ! run_timed_shell_try "sort original ${mode} fragments 5-col" \
        "${mode_dir}/times/sort_original_fragments.time.txt" "${log}" "${COMMANDS}" \
        "$(sort_fragments_5col_cmd "${original_dir}/fragments.tsv" "${original_dir}/compare/fragments.5col.sorted.tsv")"; then
        status="FAIL"
        failures+=("sort_original_fragments")
      fi
    fi
  fi
  if [[ -s "${suite_dir}/fragments.tsv.gz" ]]; then
    if [[ "${RESUME_COMPARISONS}" == "1" && -s "${suite_dir}/compare/fragments.5col.sorted.tsv" ]]; then
      log_msg "[full-benchmark:${mode}] reusing sorted Chromap-suite fragments for identity check"
    else
      log_msg "[full-benchmark:${mode}] sorting Chromap-suite fragments for identity check"
      if ! run_timed_shell_try "sort Chromap-suite ${mode} fragments 5-col" \
        "${mode_dir}/times/sort_suite_fragments.time.txt" "${log}" "${COMMANDS}" \
        "$(sort_fragments_5col_cmd "${suite_dir}/fragments.tsv.gz" "${suite_dir}/compare/fragments.5col.sorted.tsv")"; then
        status="FAIL"
        failures+=("sort_suite_fragments")
      fi
    fi
  fi

  if [[ -s "${original_dir}/macs3/original_chromap_macs3_peaks.narrowPeak" ]]; then
    bench_bed3_from_narrowpeak "${original_dir}/macs3/original_chromap_macs3_peaks.narrowPeak" \
      "${original_dir}/macs3/original_chromap_macs3_peaks.bed3"
  fi
  if [[ -s "${original_dir}/macs3/original_chromap_macs3_summits.bed" ]]; then
    bench_bed3_from_narrowpeak "${original_dir}/macs3/original_chromap_macs3_summits.bed" \
      "${original_dir}/macs3/original_chromap_macs3_summits.bed3"
  fi
  if [[ -s "${suite_dir}/chromap_suite_libmacs3_peaks.narrowPeak" ]]; then
    bench_bed3_from_narrowpeak "${suite_dir}/chromap_suite_libmacs3_peaks.narrowPeak" \
      "${suite_dir}/chromap_suite_libmacs3_peaks.bed3"
  fi
  if [[ -s "${suite_dir}/chromap_suite_libmacs3_summits.bed" ]]; then
    bench_bed3_from_narrowpeak "${suite_dir}/chromap_suite_libmacs3_summits.bed" \
      "${suite_dir}/chromap_suite_libmacs3_summits.bed3"
  fi

  local fragments_identical peaks_bed3_identical summits_bed3_identical
  fragments_identical="$(cmp_bool_or_na "${original_dir}/compare/fragments.5col.sorted.tsv" "${suite_dir}/compare/fragments.5col.sorted.tsv")"
  peaks_bed3_identical="$(cmp_bool_or_na "${original_dir}/macs3/original_chromap_macs3_peaks.bed3" "${suite_dir}/chromap_suite_libmacs3_peaks.bed3")"
  summits_bed3_identical="$(cmp_bool_or_na "${original_dir}/macs3/original_chromap_macs3_summits.bed3" "${suite_dir}/chromap_suite_libmacs3_summits.bed3")"
  local original_only suite_only common_fragments
  read -r original_only suite_only common_fragments < <(line_set_counts_sorted_or_na \
    "${original_dir}/compare/fragments.5col.sorted.tsv" \
    "${suite_dir}/compare/fragments.5col.sorted.tsv")
  local peak_jaccard summit_jaccard
  peak_jaccard="$(bench_jaccard_or_na \
    "${original_dir}/macs3/original_chromap_macs3_peaks.bed3" \
    "${suite_dir}/chromap_suite_libmacs3_peaks.bed3")"
  summit_jaccard="$(bench_jaccard_or_na \
    "${original_dir}/macs3/original_chromap_macs3_summits.bed3" \
    "${suite_dir}/chromap_suite_libmacs3_summits.bed3")"
  local summit_distance_mean summit_distance_max summit_distance_zero summit_distance_n
  read -r summit_distance_mean summit_distance_max summit_distance_zero summit_distance_n < <(summit_distance_stats_or_na \
    "${suite_dir}/chromap_suite_libmacs3_summits.bed3" \
    "${original_dir}/macs3/original_chromap_macs3_summits.bed3")

  if [[ "${fragments_identical}" == "false" || "${peaks_bed3_identical}" == "false" ]]; then
    status="FAIL"
  fi

  local mode_summary="${mode_dir}/summary.tsv"
  {
    printf 'metric\tvalue\n'
    printf 'mode\t%s\n' "${mode}"
    printf 'status\t%s\n' "${status}"
    printf 'failed_steps\t%s\n' "$(IFS=,; printf '%s' "${failures[*]:-}")"
    printf 'low_mem_enabled\t%s\n' "${low_mem_flag}"
    printf 'threads\t%s\n' "${THREADS}"
    printf 'samtools_threads\t%s\n' "${SAMTOOLS_THREADS}"
    printf 'original_fragments_count\t%s\n' "$(count_lines_or_na "${original_dir}/compare/fragments.5col.sorted.tsv")"
    printf 'suite_fragments_count\t%s\n' "$(count_lines_or_na "${suite_dir}/compare/fragments.5col.sorted.tsv")"
    printf 'original_fragments_sorted_md5\t%s\n' "$(md5_or_na "${original_dir}/compare/fragments.5col.sorted.tsv")"
    printf 'suite_fragments_sorted_md5\t%s\n' "$(md5_or_na "${suite_dir}/compare/fragments.5col.sorted.tsv")"
    printf 'suite_vs_original_fragments_identical\t%s\n' "${fragments_identical}"
    printf 'suite_vs_original_fragments_original_only\t%s\n' "${original_only}"
    printf 'suite_vs_original_fragments_suite_only\t%s\n' "${suite_only}"
    printf 'suite_vs_original_fragments_common\t%s\n' "${common_fragments}"
    printf 'original_bam_records\t%s\n' "$(bam_records_or_na "${original_dir}/possorted.bam")"
    printf 'suite_bam_records\t%s\n' "$(bam_records_or_na "${suite_dir}/possorted.bam")"
    printf 'macs3_peak_count\t%s\n' "$(count_noncomment_or_na "${original_dir}/macs3/original_chromap_macs3_peaks.narrowPeak")"
    printf 'chromap_suite_peak_count\t%s\n' "$(count_noncomment_or_na "${suite_dir}/chromap_suite_libmacs3_peaks.narrowPeak")"
    printf 'suite_vs_macs3_bed3_identical\t%s\n' "${peaks_bed3_identical}"
    printf 'suite_vs_macs3_bed3_jaccard\t%s\n' "${peak_jaccard}"
    printf 'suite_vs_macs3_summits_bed3_identical\t%s\n' "${summits_bed3_identical}"
    printf 'suite_vs_macs3_summits_jaccard\t%s\n' "${summit_jaccard}"
    printf 'suite_to_macs3_summit_distance_mean_bp\t%s\n' "${summit_distance_mean}"
    printf 'suite_to_macs3_summit_distance_max_bp\t%s\n' "${summit_distance_max}"
    printf 'suite_to_macs3_summit_distance_exact_zero_count\t%s\n' "${summit_distance_zero}"
    printf 'suite_to_macs3_summit_distance_n\t%s\n' "${summit_distance_n}"
    printf 'commands\t%s\n' "${COMMANDS}"
    printf 'versions\t%s\n' "${VERSIONS}"
    printf 'manifest\t%s\n' "${MANIFEST}"
    printf 'time_summary\t%s\n' "${mode_dir}/time_summary.tsv"
    printf 'log\t%s\n' "${log}"
    printf 'note_bam_compare\t%s\n' "standalone SAM/samtools BAM and integrated dual-fragment BAM can vary in non-biological read-level fields, especially read names/tie-selected duplicate representatives; fragment coordinates, barcode, duplicate count, and peak geometry are the mapping/peak parity gates"
  } >"${mode_summary}"

  write_time_summary "${mode_dir}" "${mode_dir}/time_summary.tsv"
  log_msg "[full-benchmark:${mode}] wrote ${mode_summary}"
}

main() {
  prepare_inputs_and_tools
  write_static_metadata

  case "${RUN_MODES}" in
    normal)
      run_mode normal 0
      ;;
    lowmem)
      run_mode lowmem 1
      ;;
    both)
      run_mode normal 0
      run_mode lowmem 1
      ;;
    *)
      echo "ERROR: RUN_MODES must be normal, lowmem, or both; got ${RUN_MODES}" >&2
      exit 2
      ;;
  esac

  {
    printf 'mode\tstatus\tfragments_identical\toriginal_fragments\tsuite_fragments\tmacs3_peaks\tsuite_peaks\tpeak_jaccard\ttime_summary\n'
    for mode in normal lowmem; do
      local mode_summary="${OUTDIR}/${mode}/summary.tsv"
      [[ -f "${mode_summary}" ]] || continue
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${mode}" \
        "$(awk -F'\t' '$1=="status"{print $2}' "${mode_summary}")" \
        "$(awk -F'\t' '$1=="suite_vs_original_fragments_identical"{print $2}' "${mode_summary}")" \
        "$(awk -F'\t' '$1=="original_fragments_count"{print $2}' "${mode_summary}")" \
        "$(awk -F'\t' '$1=="suite_fragments_count"{print $2}' "${mode_summary}")" \
        "$(awk -F'\t' '$1=="macs3_peak_count"{print $2}' "${mode_summary}")" \
        "$(awk -F'\t' '$1=="chromap_suite_peak_count"{print $2}' "${mode_summary}")" \
        "$(awk -F'\t' '$1=="suite_vs_macs3_bed3_jaccard"{print $2}' "${mode_summary}")" \
        "$(awk -F'\t' '$1=="time_summary"{print $2}' "${mode_summary}")"
    done
  } >"${SUMMARY}"

  log_msg "[full-benchmark] wrote ${SUMMARY}"
  cat "${SUMMARY}"
}

main "$@"
