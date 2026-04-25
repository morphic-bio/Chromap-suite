#!/usr/bin/env bash
# Full-dataset (unsampled) benchmark: file vs --macs3-frag-peaks-source memory, with
# timing and RSS. Capped --sort-bam-ram is recommended so the BAM coordinate sorter does
# not absorb most of the address space, making the in-memory compact-store cost easier
# to see in "Maximum resident set size" from GNU /usr/bin/time -v.
#
# Defaults target pbmc_unsorted_3k extracted ATAC (same 4 lanes as the 100K harness, full
# reads). Override FIXTURE_ATAC to point at fixture/atac for a faster smoke run.
#
# Environment (optional):
#   CHROMAP, OUTDIR, THREADS, SORT_BAM_RAM (default 512M), FIXTURE_ATAC, INDEX, REF, WHITELIST
#   VERIFY_PARITY=1  — cmp narrowPeak + summits + sorted fragments (default 1)
#   VERIFY_STANDALONE_CPP=0 — skip chromap_callpeaks reference (default 0; can be very slow)
#   SKIP_MAKE=1      — do not `make chromap` before running
#   RUN_SET=file|memory|both  — run only one mode for a quick partial benchmark (default both)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

CHROMAP="${CHROMAP:-${REPO_ROOT}/chromap}"
CALLPEAKS="${CALLPEAKS:-${REPO_ROOT}/chromap_callpeaks}"

# Unsampled ATAC under the benchmark tree (4 lanes × R1/R2/R3; same names as 100K harness)
DEFAULT_BENCH_ROOT="/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k"
FIXTURE_ATAC="${FIXTURE_ATAC:-${DEFAULT_BENCH_ROOT}/extracted/pbmc_unsorted_3k/atac}"
INDEX="${INDEX:-${DEFAULT_BENCH_ROOT}/chromap_index/genome.index}"
REF="${REF:-/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa}"
WHITELIST="${WHITELIST:-${DEFAULT_BENCH_ROOT}/chromap_index/737K-arc-v1_atac.txt}"
OUTDIR="${OUTDIR:-}"
THREADS="${THREADS:-8}"
SORT_BAM_RAM="${SORT_BAM_RAM:-512M}"
MINLEN="${MINLEN:-200}"
MAXGAP="${MAXGAP:-30}"
MACS3_PVAL="${MACS3_PVAL:-1e-5}"
VERIFY_PARITY="${VERIFY_PARITY:-1}"
VERIFY_STANDALONE_CPP="${VERIFY_STANDALONE_CPP:-0}"
SKIP_MAKE="${SKIP_MAKE:-0}"
RUN_SET="${RUN_SET:-both}" # file | memory | both

R1="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R1_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R1_001.fastq.gz"
R2="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R3_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R3_001.fastq.gz"
BC="${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L001_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L002_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L003_R2_001.fastq.gz,${FIXTURE_ATAC}/pbmc_unsorted_3k_S3_L004_R2_001.fastq.gz"

# shellcheck disable=SC2206
chromap_base=(
  -t "${THREADS}" -x "${INDEX}" -r "${REF}" -1 "${R1}" -2 "${R2}" -b "${BC}"
  --barcode-whitelist "${WHITELIST}" -l 2000 --trim-adapters --remove-pcr-duplicates
  --remove-pcr-duplicates-at-cell-level --Tn5-shift --BAM --sort-bam
  --sort-bam-ram "${SORT_BAM_RAM}"
)

# One-line GNU time format (portable; no -v parsing). %e wall sec; %M max RSS KB.
TIME_FMT="wall_sec %e user_sec %U sys_sec %S max_rss_kbytes %M"

parse_time_f() {
  # $1 = path to one-line `time -f` output: wall_sec %e user_sec %U sys_sec %S max_rss_kbytes %M
  if [[ ! -s "$1" ]]; then
    echo "0 0 0 0"
    return
  fi
  read -r _ w _ u _ s _ r _ < <(tr -s ' ' <"$1" | head -1)
  echo "${w:-0} ${u:-0} ${s:-0} ${r:-0}"
}

log_msg() { echo "$@" | tee -a "${LOG}"; }

main() {
  if ! command -v /usr/bin/time >/dev/null; then
    echo "This benchmark expects GNU time at /usr/bin/time" >&2
    exit 2
  fi
  if [[ ! -f "${CHROMAP}" ]]; then
    echo "Build chromap first or set CHROMAP" >&2
    exit 2
  fi
  if [[ ! -d "${FIXTURE_ATAC}" || ! -f "${INDEX}" || ! -f "${REF}" || ! -f "${WHITELIST}" ]]; then
    echo "Missing inputs. Set FIXTURE_ATAC, INDEX, REF, WHITELIST. For a quick path use:" >&2
    echo "  FIXTURE_ATAC=${DEFAULT_BENCH_ROOT}/fixture/atac" >&2
    exit 2
  fi
  for key in S3_L001_R1_001 S3_L004_R3_001; do
    if [[ ! -f "${FIXTURE_ATAC}/pbmc_unsorted_3k_${key}.fastq.gz" ]]; then
      echo "Missing ${FIXTURE_ATAC}/.../${key} — check FIXTURE_ATAC" >&2
      exit 2
    fi
  done
  if [[ -z "${OUTDIR}" ]]; then
    OUTDIR="$(mktemp -d /tmp/chromap_peak_memory_fullset_bench.XXXXXX)"
  fi
  mkdir -p "${OUTDIR}/file" "${OUTDIR}/memory" "${OUTDIR}/logs"
  LOG="${OUTDIR}/benchmark_fullset.log"
  local summary="${OUTDIR}/benchmark_fullset.tsv"

  log_msg "=== Chromap memory-source full-set benchmark ==="
  log_msg "OUTDIR=${OUTDIR}"
  log_msg "FIXTURE_ATAC=${FIXTURE_ATAC}"
  log_msg "THREADS=${THREADS} SORT_BAM_RAM=${SORT_BAM_RAM} RUN_SET=${RUN_SET}"
  log_msg "VERIFY_PARITY=${VERIFY_PARITY} VERIFY_STANDALONE_CPP=${VERIFY_STANDALONE_CPP}"

  if [[ "${SKIP_MAKE}" == "0" ]]; then
    (cd "${REPO_ROOT}" && make chromap chromap_callpeaks) >>"${LOG}" 2>&1
  fi

  run_file() {
    /usr/bin/time -f "${TIME_FMT}" -o "${OUTDIR}/logs/time_v_file.txt" \
      "${CHROMAP}" "${chromap_base[@]}" \
      --atac-fragments "${OUTDIR}/file/fragments.tsv.gz" \
      --summary "${OUTDIR}/file/summary.tsv" \
      -o "${OUTDIR}/file/possorted_bam.bam" \
      --call-macs3-frag-peaks \
      --macs3-frag-peaks-source file \
      --macs3-frag-peaks-output "${OUTDIR}/file/chromap_macs3_frag.narrowPeak" \
      --macs3-frag-summits-output "${OUTDIR}/file/chromap_macs3_frag_summits.bed" \
      --macs3-frag-pvalue "${MACS3_PVAL}" \
      --macs3-frag-min-length "${MINLEN}" \
      --macs3-frag-max-gap "${MAXGAP}" \
      >>"${LOG}" 2>&1
  }

  run_mem() {
    /usr/bin/time -f "${TIME_FMT}" -o "${OUTDIR}/logs/time_v_memory.txt" \
      "${CHROMAP}" "${chromap_base[@]}" \
      --atac-fragments "${OUTDIR}/memory/fragments.tsv.gz" \
      --summary "${OUTDIR}/memory/summary.tsv" \
      -o "${OUTDIR}/memory/possorted_bam.bam" \
      --call-macs3-frag-peaks \
      --macs3-frag-peaks-source memory \
      --macs3-frag-peaks-output "${OUTDIR}/memory/chromap_macs3_frag.narrowPeak" \
      --macs3-frag-summits-output "${OUTDIR}/memory/chromap_macs3_frag_summits.bed" \
      --macs3-frag-pvalue "${MACS3_PVAL}" \
      --macs3-frag-min-length "${MINLEN}" \
      --macs3-frag-max-gap "${MAXGAP}" \
      >>"${LOG}" 2>&1
  }

  if [[ "${RUN_SET}" == "file" || "${RUN_SET}" == "both" ]]; then
    log_msg "[run] file-source peaks + BAM (GNU time -f -> logs/time_v_file.txt)"
    run_file
  fi
  if [[ "${RUN_SET}" == "memory" || "${RUN_SET}" == "both" ]]; then
    log_msg "[run] memory-source peaks + BAM (GNU time -f -> logs/time_v_memory.txt)"
    run_mem
  fi

  local frags_file frags_mem storage_mode
  if [[ -f "${OUTDIR}/file/fragments.tsv.gz" ]]; then
    frags_file="$(zcat "${OUTDIR}/file/fragments.tsv.gz" | wc -l)"
  else
    frags_file="na"
  fi
  if [[ -f "${OUTDIR}/memory/fragments.tsv.gz" ]]; then
    frags_mem="$(zcat "${OUTDIR}/memory/fragments.tsv.gz" | wc -l)"
  else
    frags_mem="na"
  fi
  storage_mode="$(awk '/MACS3 FRAG memory storage mode:/{print $NF; exit}' "${LOG}" 2>/dev/null || true)"
  if [[ -z "${storage_mode}" ]]; then
    storage_mode="na"
  fi

  if [[ "${VERIFY_PARITY}" == "1" && "${RUN_SET}" == "both" ]]; then
    log_msg "[check] sorted fragment lines + narrowPeak + summits"
    if ! cmp -s <(zcat "${OUTDIR}/file/fragments.tsv.gz" | LC_ALL=C sort) \
                <(zcat "${OUTDIR}/memory/fragments.tsv.gz" | LC_ALL=C sort); then
      echo "FAIL: fragments differ" >&2
      exit 1
    fi
    if ! cmp -s "${OUTDIR}/file/chromap_macs3_frag.narrowPeak" \
                "${OUTDIR}/memory/chromap_macs3_frag.narrowPeak"; then
      echo "FAIL: narrowPeak differs" >&2
      exit 1
    fi
    if ! cmp -s "${OUTDIR}/file/chromap_macs3_frag_summits.bed" \
                "${OUTDIR}/memory/chromap_macs3_frag_summits.bed"; then
      echo "FAIL: summits differ" >&2
      exit 1
    fi
    log_msg "parity: OK (fragments, narrowPeak, summits)"
  fi

  if [[ "${VERIFY_STANDALONE_CPP}" == "1" && -f "${CALLPEAKS}" && -f "${OUTDIR}/file/fragments.tsv.gz" ]]; then
    log_msg "[optional] chromap_callpeaks on file fragments (slow)"
    "${CALLPEAKS}" -i "${OUTDIR}/file/fragments.tsv.gz" \
      --frag-pileup-macs3-uint8-counts \
      --macs3-frag-narrowpeak "${OUTDIR}/ref_cpp_macs3_frag.narrowPeak" \
      --macs3-frag-summits "${OUTDIR}/ref_cpp_macs3_frag_summits.bed" \
      --bdgpeakcall-cutoff 5 --bdgpeakcall-min-len "${MINLEN}" --bdgpeakcall-max-gap "${MAXGAP}" \
      --frag-score-pseudocount 0 \
      >>"${LOG}" 2>&1
    cmp -s "${OUTDIR}/file/chromap_macs3_frag.narrowPeak" "${OUTDIR}/ref_cpp_macs3_frag.narrowPeak"
  fi

  local f_wall f_user f_sys f_rss m_wall m_user m_sys m_rss
  f_wall=-1; f_user=-1; f_sys=-1; f_rss=-1
  m_wall=-1; m_user=-1; m_sys=-1; m_rss=-1
  if [[ -f "${OUTDIR}/logs/time_v_file.txt" ]]; then
    read -r f_wall f_user f_sys f_rss < <(parse_time_f "${OUTDIR}/logs/time_v_file.txt")
  fi
  if [[ -f "${OUTDIR}/logs/time_v_memory.txt" ]]; then
    read -r m_wall m_user m_sys m_rss < <(parse_time_f "${OUTDIR}/logs/time_v_memory.txt")
  fi

  local delta_rss="na"
  if [[ "${f_rss}" -ge 0 && "${m_rss}" -ge 0 ]]; then
    delta_rss=$((m_rss - f_rss))
  fi

  {
    echo "metric	file	memory	notes"
    echo "wall_sec	${f_wall}	${m_wall}	GNU time %e (seconds)"
    echo "user_sec	${f_user}	${m_user}	"
    echo "sys_sec	${f_sys}	${m_sys}	"
    echo "max_rss_kbytes	${f_rss}	${m_rss}	GNU time %M (kilobytes)"
    echo "delta_max_rss_kbytes_memory_minus_file	${delta_rss}		+ => memory run higher RSS"
    echo "sort_bam_ram_cap	${SORT_BAM_RAM}		--sort-bam-ram (lower => more sorter spill; can expose compact-store delta)"
    echo "fragment_lines_file	${frags_file}		wc -l of fragments TSV"
    echo "fragment_lines_memory	${frags_mem}		"
    echo "memory_storage_mode		${storage_mode}	from stderr log (memory run)"
    echo "threads	${THREADS}		"
    echo "fixture_atac	${FIXTURE_ATAC}		"
  } | tee "${summary}"

  log_msg "Wrote ${summary} and per-run ${OUTDIR}/logs/time_v_*.txt"
  log_msg "RSS: use max_rss_kbytes (GNU time %M). Tight SORT_BAM_RAM caps sorter RAM so the in-memory fragment store delta is less likely to be hidden by a large default sort buffer."
}

main "$@"
