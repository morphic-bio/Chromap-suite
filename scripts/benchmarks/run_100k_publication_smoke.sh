#!/usr/bin/env bash
# Publication smoke for the Chromap-suite libMACS3 integration.
#
# Baseline:
#   unfixed upstream Chromap 0.3.3 in normal-memory mode emits SAM to stdout,
#   which is piped through samtools view/sort/index for sorted BAM. The same
#   original upstream Chromap binary emits BED/FRAG rows for standalone MACS3.
#
# Integrated:
#   Chromap-suite emits sorted BAM, fragments, and libMACS3 FRAG peaks in one
#   command.
#
# Important: do not use `--preset atac` here. Upstream Chromap's ATAC preset
# enables low-memory mode. This smoke uses explicit ATAC-equivalent flags in
# normal-memory mode for both the baseline and integrated runs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"
# shellcheck source=scripts/benchmarks/libmacs3_chromap_common.sh
source "${SCRIPT_DIR}/libmacs3_chromap_common.sh"

export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

BENCH_ROOT="${BENCH_ROOT:-/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k}"
FIXTURE_ATAC="${FIXTURE_ATAC:-${BENCH_ROOT}/fixture/atac}"
INDEX="${INDEX:-${BENCH_ROOT}/chromap_index/genome.index}"
REF="${REF:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa}"
WHITELIST="${WHITELIST:-${BENCH_ROOT}/chromap_index/737K-arc-v1_atac.txt}"

DEFAULT_UNFIXED_CHROMAP="${REPO_ROOT}/plans/artifacts/upstream_chromap_unfixed/c9d8ae058bfdf04d45bc5f99a164a460f759a6a7/chromap"
CHROMAP_ORIGINAL="${CHROMAP_ORIGINAL:-${DEFAULT_UNFIXED_CHROMAP}}"
CHROMAP_SUITE="${CHROMAP_SUITE:-${REPO_ROOT}/chromap}"
MACS3_BIN="${MACS3_BIN:-macs3}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

THREADS="${THREADS:-1}"
SAMTOOLS_THREADS="${SAMTOOLS_THREADS:-1}"
MACS3_PVALUE="${MACS3_PVALUE:-1e-5}"
GENOME_SIZE="${GENOME_SIZE:-hs}"
MACS3_MIN_LENGTH="${MACS3_MIN_LENGTH:-200}"
MACS3_MAX_GAP="${MACS3_MAX_GAP:-30}"

ARTIFACT_ROOT="${ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts/libmacs3_chromap_atac_panel}"
RUN_ID="${RUN_ID:-$(date -u +%Y%m%dT%H%M%SZ)}"
OUTDIR="${OUTDIR:-${ARTIFACT_ROOT}/${RUN_ID}/100k_smoke}"

R1="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R1_001.fastq.gz"
R2="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R3_001.fastq.gz"
BC="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R2_001.fastq.gz"

mkdir -p \
  "${OUTDIR}/original_chromap/macs3" \
  "${OUTDIR}/original_chromap/compare" \
  "${OUTDIR}/chromap_suite/compare" \
  "${OUTDIR}/logs" \
  "${OUTDIR}/times"

COMMANDS="${OUTDIR}/commands.sh"
SUMMARY="${OUTDIR}/summary.tsv"
VERSIONS="${OUTDIR}/versions.txt"
LOG="${OUTDIR}/logs/run.log"

{
  printf '#!/usr/bin/env bash\n'
  printf 'set -euo pipefail\n\n'
  printf '# Re-run from repo root: %s\n' "${REPO_ROOT}"
  printf '# Output directory from this run: %s\n' "${OUTDIR}"
} >"${COMMANDS}"
chmod +x "${COMMANDS}"

for path in "${INDEX}" "${REF}" "${WHITELIST}"; do
  bench_require_file "${path}" "benchmark input"
done
for lane in L001 L002 L003 L004; do
  bench_require_file "${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_${lane}_R1_001.fastq.gz" "R1 ${lane}"
  bench_require_file "${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_${lane}_R2_001.fastq.gz" "barcode ${lane}"
  bench_require_file "${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_${lane}_R3_001.fastq.gz" "R2 ${lane}"
done
if [[ ! -x "${CHROMAP_ORIGINAL}" && "${CHROMAP_ORIGINAL}" == "${DEFAULT_UNFIXED_CHROMAP}" ]]; then
  echo "[100k-smoke] preparing unfixed upstream Chromap baseline" | tee -a "${LOG}"
  CHROMAP_ORIGINAL="$("${SCRIPT_DIR}/prepare_unfixed_chromap.sh")"
fi
bench_require_exe "${CHROMAP_ORIGINAL}" "unfixed upstream Chromap"
bench_require_exe "${MACS3_BIN}" "MACS3"
bench_require_exe "${SAMTOOLS_BIN}" "samtools"

echo "[100k-smoke] building Chromap-suite binaries" | tee -a "${LOG}"
bench_run_timed "build Chromap-suite chromap + chromap_callpeaks" \
  "${OUTDIR}/times/build.time.txt" "${LOG}" "${COMMANDS}" \
  make -C "${REPO_ROOT}" chromap chromap_callpeaks
bench_require_exe "${CHROMAP_SUITE}" "Chromap-suite"
bench_require_exe "${REPO_ROOT}/chromap_callpeaks" "chromap_callpeaks"

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
  printf 'macs3_pvalue\t%s\n' "${MACS3_PVALUE}"
  printf 'macs3_min_length\t%s\n' "${MACS3_MIN_LENGTH}"
  printf 'macs3_max_gap\t%s\n' "${MACS3_MAX_GAP}"
  printf 'mapping_mode\t%s\n' "explicit_atac_flags_normal_memory_no_preset"
  printf 'chromap_original_unfixed\t%s\n' "${CHROMAP_ORIGINAL}"
  printf 'chromap_suite\t%s\n' "${CHROMAP_SUITE}"
  printf 'macs3\t%s\n' "$(command -v "${MACS3_BIN}")"
  printf 'samtools\t%s\n' "$(command -v "${SAMTOOLS_BIN}")"
} >"${OUTDIR}/manifest.tsv"

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

chromap_common_args=(
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

echo "[100k-smoke] original unfixed Chromap baseline: SAM -> samtools sorted BAM" | tee -a "${LOG}"
original_sam_cmd="$(printf '%q ' "${CHROMAP_ORIGINAL}" -t "${THREADS}" "${chromap_common_args[@]}" --SAM --summary "${OUTDIR}/original_chromap/summary.sam.tsv" -o /dev/stdout)"
original_sam_cmd+=" | $(printf '%q ' "${SAMTOOLS_BIN}" view -@ "${SAMTOOLS_THREADS}" -b -)"
original_sam_cmd+=" | $(printf '%q ' "${SAMTOOLS_BIN}" sort -@ "${SAMTOOLS_THREADS}" -o "${OUTDIR}/original_chromap/possorted.bam" -)"
bench_run_timed_shell "original unfixed Chromap SAM piped to samtools sorted BAM" \
  "${OUTDIR}/times/original_sam_to_bam.time.txt" "${LOG}" "${COMMANDS}" "${original_sam_cmd}"
bench_run_timed "index original baseline BAM" \
  "${OUTDIR}/times/original_bam_index.time.txt" "${LOG}" "${COMMANDS}" \
  "${SAMTOOLS_BIN}" index "${OUTDIR}/original_chromap/possorted.bam"

echo "[100k-smoke] original unfixed Chromap baseline: BED/FRAG rows" | tee -a "${LOG}"
bench_run_timed "original unfixed Chromap BED/FRAG rows" \
  "${OUTDIR}/times/original_fragments.time.txt" "${LOG}" "${COMMANDS}" \
  "${CHROMAP_ORIGINAL}" -t "${THREADS}" "${chromap_common_args[@]}" \
  --BED --summary "${OUTDIR}/original_chromap/summary.fragments.tsv" \
  -o "${OUTDIR}/original_chromap/fragments.tsv"
bench_sort_fragments_5col "${OUTDIR}/original_chromap/fragments.tsv" \
  "${OUTDIR}/original_chromap/compare/fragments.5col.sorted.tsv"
bench_bed3_from_fragments "${OUTDIR}/original_chromap/fragments.tsv" \
  "${OUTDIR}/original_chromap/fragments.bed3"

echo "[100k-smoke] original unfixed Chromap baseline: MACS3 FRAG peaks" | tee -a "${LOG}"
bench_run_timed "MACS3 callpeak on original unfixed Chromap fragments" \
  "${OUTDIR}/times/original_macs3.time.txt" "${LOG}" "${COMMANDS}" \
  "${MACS3_BIN}" callpeak \
  -t "${OUTDIR}/original_chromap/fragments.tsv" \
  -f FRAG -g "${GENOME_SIZE}" \
  -n original_chromap_macs3 \
  -p "${MACS3_PVALUE}" \
  --min-length "${MACS3_MIN_LENGTH}" \
  --max-gap "${MACS3_MAX_GAP}" \
  --outdir "${OUTDIR}/original_chromap/macs3"
bench_bed3_from_narrowpeak "${OUTDIR}/original_chromap/macs3/original_chromap_macs3_peaks.narrowPeak" \
  "${OUTDIR}/original_chromap/macs3/original_chromap_macs3_peaks.bed3"

echo "[100k-smoke] Chromap-suite integrated command" | tee -a "${LOG}"
bench_run_timed "Chromap-suite integrated BAM + fragments + libMACS3 peaks" \
  "${OUTDIR}/times/chromap_suite_integrated.time.txt" "${LOG}" "${COMMANDS}" \
  "${CHROMAP_SUITE}" -t "${THREADS}" "${chromap_common_args[@]}" \
  --BAM --sort-bam --write-index \
  --atac-fragments "${OUTDIR}/chromap_suite/fragments.tsv.gz" \
  --summary "${OUTDIR}/chromap_suite/summary.tsv" \
  --call-macs3-frag-peaks \
  --macs3-frag-peaks-output "${OUTDIR}/chromap_suite/chromap_suite_libmacs3_peaks.narrowPeak" \
  --macs3-frag-summits-output "${OUTDIR}/chromap_suite/chromap_suite_libmacs3_summits.bed" \
  -o "${OUTDIR}/chromap_suite/possorted.bam"
bench_sort_fragments_5col "${OUTDIR}/chromap_suite/fragments.tsv.gz" \
  "${OUTDIR}/chromap_suite/compare/fragments.5col.sorted.tsv"
bench_bed3_from_fragments "${OUTDIR}/chromap_suite/fragments.tsv.gz" \
  "${OUTDIR}/chromap_suite/fragments.bed3"
bench_bed3_from_narrowpeak "${OUTDIR}/chromap_suite/chromap_suite_libmacs3_peaks.narrowPeak" \
  "${OUTDIR}/chromap_suite/chromap_suite_libmacs3_peaks.bed3"

echo "[100k-smoke] standalone chromap_callpeaks reference from original unfixed fragments" | tee -a "${LOG}"
bench_run_timed "chromap_callpeaks on original unfixed Chromap fragments" \
  "${OUTDIR}/times/chromap_callpeaks_original_fragments.time.txt" "${LOG}" "${COMMANDS}" \
  "${REPO_ROOT}/chromap_callpeaks" \
  -i "${OUTDIR}/original_chromap/fragments.tsv" \
  --frag-pileup-macs3-uint8-counts \
  --macs3-frag-narrowpeak "${OUTDIR}/original_chromap/chromap_callpeaks_from_original.narrowPeak" \
  --macs3-frag-summits "${OUTDIR}/original_chromap/chromap_callpeaks_from_original_summits.bed" \
  --bdgpeakcall-cutoff 5 \
  --bdgpeakcall-min-len "${MACS3_MIN_LENGTH}" \
  --bdgpeakcall-max-gap "${MACS3_MAX_GAP}" \
  --frag-score-pseudocount 0
bench_bed3_from_narrowpeak "${OUTDIR}/original_chromap/chromap_callpeaks_from_original.narrowPeak" \
  "${OUTDIR}/original_chromap/chromap_callpeaks_from_original.bed3"

original_frag_md5="$(bench_md5_file "${OUTDIR}/original_chromap/compare/fragments.5col.sorted.tsv")"
suite_frag_md5="$(bench_md5_file "${OUTDIR}/chromap_suite/compare/fragments.5col.sorted.tsv")"
original_bam_body_md5="$(bench_md5_bam_body_sorted "${SAMTOOLS_BIN}" "${OUTDIR}/original_chromap/possorted.bam")"
suite_bam_body_md5="$(bench_md5_bam_body_sorted "${SAMTOOLS_BIN}" "${OUTDIR}/chromap_suite/possorted.bam")"
original_bam_records="$(bench_bam_total_records "${SAMTOOLS_BIN}" "${OUTDIR}/original_chromap/possorted.bam")"
suite_bam_records="$(bench_bam_total_records "${SAMTOOLS_BIN}" "${OUTDIR}/chromap_suite/possorted.bam")"

suite_original_fragments_identical=false
suite_original_bam_body_identical=false
suite_macs3_bed3_identical=false
suite_callpeaks_narrowpeak_identical=false
suite_callpeaks_summits_identical=false

cmp -s "${OUTDIR}/original_chromap/compare/fragments.5col.sorted.tsv" \
       "${OUTDIR}/chromap_suite/compare/fragments.5col.sorted.tsv" &&
  suite_original_fragments_identical=true
[[ "${original_bam_body_md5}" == "${suite_bam_body_md5}" ]] &&
  suite_original_bam_body_identical=true
cmp -s "${OUTDIR}/original_chromap/macs3/original_chromap_macs3_peaks.bed3" \
       "${OUTDIR}/chromap_suite/chromap_suite_libmacs3_peaks.bed3" &&
  suite_macs3_bed3_identical=true
cmp -s "${OUTDIR}/original_chromap/chromap_callpeaks_from_original.narrowPeak" \
       "${OUTDIR}/chromap_suite/chromap_suite_libmacs3_peaks.narrowPeak" &&
  suite_callpeaks_narrowpeak_identical=true
cmp -s "${OUTDIR}/original_chromap/chromap_callpeaks_from_original_summits.bed" \
       "${OUTDIR}/chromap_suite/chromap_suite_libmacs3_summits.bed" &&
  suite_callpeaks_summits_identical=true

macs3_peak_count="$(bench_count_noncomment_lines "${OUTDIR}/original_chromap/macs3/original_chromap_macs3_peaks.narrowPeak")"
suite_peak_count="$(bench_count_noncomment_lines "${OUTDIR}/chromap_suite/chromap_suite_libmacs3_peaks.narrowPeak")"
callpeaks_peak_count="$(bench_count_noncomment_lines "${OUTDIR}/original_chromap/chromap_callpeaks_from_original.narrowPeak")"
fragments_count="$(wc -l <"${OUTDIR}/original_chromap/compare/fragments.5col.sorted.tsv")"
peak_jaccard="$(bench_jaccard_or_na \
  "${OUTDIR}/original_chromap/macs3/original_chromap_macs3_peaks.bed3" \
  "${OUTDIR}/chromap_suite/chromap_suite_libmacs3_peaks.bed3")"
read -r original_only suite_only common_fragments < <(bench_line_set_counts \
  "${OUTDIR}/original_chromap/compare/fragments.5col.sorted.tsv" \
  "${OUTDIR}/chromap_suite/compare/fragments.5col.sorted.tsv" \
  "${OUTDIR}/chromap_suite/compare/original_vs_suite_fragments")

status="PASS"
if [[ "${suite_original_fragments_identical}" != true ||
      "${suite_macs3_bed3_identical}" != true ||
      "${suite_callpeaks_narrowpeak_identical}" != true ||
      "${suite_callpeaks_summits_identical}" != true ]]; then
  status="FAIL"
fi

{
  printf 'metric\tvalue\n'
  printf 'outdir\t%s\n' "${OUTDIR}"
  printf 'status\t%s\n' "${status}"
  printf 'fragments_count\t%s\n' "${fragments_count}"
  printf 'original_fragments_sorted_md5\t%s\n' "${original_frag_md5}"
  printf 'suite_fragments_sorted_md5\t%s\n' "${suite_frag_md5}"
  printf 'suite_vs_original_fragments_identical\t%s\n' "${suite_original_fragments_identical}"
  printf 'suite_vs_original_fragments_original_only\t%s\n' "${original_only}"
  printf 'suite_vs_original_fragments_suite_only\t%s\n' "${suite_only}"
  printf 'suite_vs_original_fragments_common\t%s\n' "${common_fragments}"
  printf 'original_bam_body_sorted_md5\t%s\n' "${original_bam_body_md5}"
  printf 'suite_bam_body_sorted_md5\t%s\n' "${suite_bam_body_md5}"
  printf 'original_bam_records\t%s\n' "${original_bam_records}"
  printf 'suite_bam_records\t%s\n' "${suite_bam_records}"
  printf 'suite_vs_original_bam_body_identical\t%s\n' "${suite_original_bam_body_identical}"
  printf 'macs3_peak_count\t%s\n' "${macs3_peak_count}"
  printf 'chromap_suite_peak_count\t%s\n' "${suite_peak_count}"
  printf 'chromap_callpeaks_peak_count\t%s\n' "${callpeaks_peak_count}"
  printf 'suite_vs_macs3_bed3_identical\t%s\n' "${suite_macs3_bed3_identical}"
  printf 'suite_vs_macs3_bed3_jaccard\t%s\n' "${peak_jaccard}"
  printf 'suite_vs_chromap_callpeaks_narrowpeak_identical\t%s\n' "${suite_callpeaks_narrowpeak_identical}"
  printf 'suite_vs_chromap_callpeaks_summits_identical\t%s\n' "${suite_callpeaks_summits_identical}"
  printf 'commands\t%s\n' "${COMMANDS}"
  printf 'versions\t%s\n' "${VERSIONS}"
  printf 'manifest\t%s\n' "${OUTDIR}/manifest.tsv"
  printf 'log\t%s\n' "${LOG}"
  printf 'note_bam_compare\t%s\n' "standalone SAM/samtools BAM and integrated dual-fragment BAM can vary in non-biological read-level fields, especially read names/tie-selected duplicate representatives; fragment coordinates, barcode, duplicate count, and peak geometry are the mapping/peak parity gates"
} >"${SUMMARY}"

if [[ "${suite_original_fragments_identical}" != true ]]; then
  echo "FAIL: Chromap-suite fragments differ from original Chromap baseline" >&2
  exit 1
fi
if [[ "${suite_macs3_bed3_identical}" != true ]]; then
  echo "FAIL: Chromap-suite peak BED3 differs from MACS3 baseline" >&2
  exit 1
fi
if [[ "${suite_callpeaks_narrowpeak_identical}" != true ]]; then
  echo "FAIL: Chromap-suite narrowPeak differs from chromap_callpeaks reference" >&2
  exit 1
fi
if [[ "${suite_callpeaks_summits_identical}" != true ]]; then
  echo "FAIL: Chromap-suite summits differ from chromap_callpeaks reference" >&2
  exit 1
fi

echo "[100k-smoke] PASS: ${SUMMARY}" | tee -a "${LOG}"
cat "${SUMMARY}"
