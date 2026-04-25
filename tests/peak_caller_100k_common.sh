# Shared paired fragments/BAM discovery and BED metrics for 100K peak-caller
# scripts. Sourced by run_peak_caller_100k.sh and run_peak_caller_calibration_100k.sh
# shellcheck shell=bash

: "${BENCH_ROOT:=/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k}"
: "${CHROMAP_PEAK_RUN_ROOT:=}"
: "${FRAGMENTS_TSV_GZ:=}"
: "${ATAC_BAM:=}"
: "${ARC_PEAKS_BED:=}"
: "${RUN_MACS3:=1}"

# Find a BAM in the run tree paired with dual/fragments (same-run harness).
# Prefer possorted_bam under dual/, then ATAC-style names used by STAR/ARC.
bam_in_run_tree() {
  local r=$1
  local p
  for p in \
    "${r}/dual/possorted_bam.bam" \
    "${r}/dual/atac_possorted_bam.bam" \
    "${r}/atac_possorted_bam.bam" \
    "${r}/out/atac_possorted_bam.bam" \
    "${r}/bam_only/possorted_bam.bam"; do
    if [[ -f "$p" ]]; then
      printf '%s' "$p"
      return 0
    fi
  done
  find "$r" -maxdepth 5 \
    \( -name 'possorted_bam.bam' -o -name 'atac_possorted_bam.bam' \) -type f 2>/dev/null |
    LC_ALL=C sort | head -1
}

resolve_fragments_in_run() {
  local r=$1
  local fr
  for fr in "${r}/dual/fragments.tsv.gz" "${r}/dual/fragments.tsv"; do
    if [[ -f "$fr" ]]; then
      printf '%s' "$fr"
      return 0
    fi
  done
  return 1
}

# Return first pair (fragments path, bam) under BENCH; BAM under same run root.
discover_paired_bench() {
  [[ -d "${BENCH_ROOT}" ]] || return 1
  local f r b
  while IFS= read -r f; do
    r="$(cd "$(dirname "$f")/.." && pwd)"
    b="$(bam_in_run_tree "$r")"
    if [[ -n "$b" && -f "$b" ]]; then
      echo "$f"
      echo "$b"
      return 0
    fi
  done < <(find "${BENCH_ROOT}" \( -path '*/dual/fragments.tsv.gz' -o -path '*/dual/fragments.tsv' \) -type f 2>/dev/null | LC_ALL=C sort)
  return 1
}

discover_frags_only_bench() {
  [[ -d "${BENCH_ROOT}" ]] || return 1
  local f
  f="$(find "${BENCH_ROOT}" \( -path '*/dual/fragments.tsv.gz' -o -path '*/dual/fragments.tsv' \) -type f 2>/dev/null | LC_ALL=C sort | head -1 || true)"
  if [[ -n "$f" && -f "$f" ]]; then
    echo "$f"
    return 0
  fi
  return 1
}

discover_arc_peaks() {
  if [[ -n "${ARC_PEAKS_BED}" && -f "${ARC_PEAKS_BED}" ]]; then
    printf '%s' "${ARC_PEAKS_BED}"
    return 0
  fi
  local p
  p="${BENCH_ROOT}/pbmc_unsorted_3k_100k_arc/outs/atac_peaks.bed"
  if [[ -f "${p}" ]]; then
    printf '%s' "${p}"
    return 0
  fi
  return 1
}

peak_100k_run_root_from_frag() {
  local fr=$1
  if [[ "${fr}" == */dual/fragments.tsv* ]]; then
    (cd "$(dirname "${fr}")/.." && pwd)
    return 0
  fi
  dirname "${fr}"
}

peak_100k_assert_bam_under_run() {
  local run_root=$1
  local bam=$2
  [[ -n "${run_root}" && -n "${bam}" && -f "${bam}" ]] || return 0
  local rr br
  rr="$(realpath "${run_root}")"
  br="$(realpath "${bam}")"
  if [[ "${br}" != "${rr}" && "${br}" != "${rr}/"* ]]; then
    echo "PEAK_100K_BAM is not under PEAK_100K_RUN_ROOT=${rr} (bam=${br})" >&2
    return 1
  fi
  return 0
}

# Build MACS3-compatible tagAlign.gz from a paired-end ATAC BAM (same as
# coordination run_macs3_peak_call_from_bam.sh BAM->tagAlign step).
peak_100k_tagalign_from_bam() {
  local bam=$1
  local out_gz=$2
  local threads=${3:-1}
  local samtools_bin="${SAMTOOLS_BIN:-samtools}"
  local bedtools_bin="${BEDTOOLS_BIN:-bedtools}"
  for t in "${samtools_bin}" "${bedtools_bin}"; do
    if ! command -v "${t}" &>/dev/null; then
      echo "peak_100k_tagalign_from_bam: missing ${t}" >&2
      return 1
    fi
  done
  [[ -f "${bam}" ]] || {
    echo "peak_100k_tagalign_from_bam: missing BAM ${bam}" >&2
    return 1
  }
  local work
  work="$(mktemp -d "${TMPDIR:-/tmp}/peak100k_tagalign.XXXXXX")"
  cleanup() {
    rm -rf "${work}"
  }
  trap cleanup EXIT
  local nmsort="${work}/nmsort.bam"
  local bedpe="${work}/bedpe.gz"
  echo "[tagalign] name-sort BAM -> ${out_gz}" >&2
  "${samtools_bin}" sort -n -@ "${threads}" -o "${nmsort}" "${bam}"
  LC_ALL=C "${bedtools_bin}" bamtobed -bedpe -mate1 -i "${nmsort}" | gzip -nc >"${bedpe}"
  zcat -f "${bedpe}" | awk 'BEGIN{OFS="\t"}
    ($1 != "." && $2 >= 0 && $3 > $2 && ($9 == "+" || $9 == "-")) {
      print $1, $2, $3, "N", 1000, $9
    }
    ($4 != "." && $5 >= 0 && $6 > $5 && ($10 == "+" || $10 == "-")) {
      print $4, $5, $6, "N", 1000, $10
    }
  ' | gzip -nc >"${out_gz}"
  trap - EXIT
  cleanup
  return 0
}

# Sets PEAK_100K_RUN_ROOT, PEAK_100K_FRAG and PEAK_100K_BAM (may be empty). Exit 0 on success.
peak_100k_resolve_inputs() {
  PEAK_100K_RUN_ROOT=""
  PEAK_100K_FRAG=""
  PEAK_100K_BAM=""
  if [[ -n "${CHROMAP_PEAK_RUN_ROOT}" ]]; then
    local r
    r="$(cd "${CHROMAP_PEAK_RUN_ROOT}" && pwd)"
    PEAK_100K_RUN_ROOT="${r}"
    if ! PEAK_100K_FRAG="$(resolve_fragments_in_run "$r")"; then
      echo "No dual/fragments.tsv(.gz) under CHROMAP_PEAK_RUN_ROOT=${r}" >&2
      return 1
    fi
    PEAK_100K_BAM="$(bam_in_run_tree "$r")"
    if [[ -n "$PEAK_100K_BAM" && ! -f "$PEAK_100K_BAM" ]]; then
      PEAK_100K_BAM=""
    fi
  elif [[ -n "${FRAGMENTS_TSV_GZ}" && -n "${ATAC_BAM}" ]]; then
    PEAK_100K_FRAG="${FRAGMENTS_TSV_GZ}"
    PEAK_100K_BAM="${ATAC_BAM}"
    [[ -f "$PEAK_100K_FRAG" ]] || {
      echo "Missing FRAGMENTS_TSV_GZ: ${PEAK_100K_FRAG}" >&2
      return 1
    }
    [[ -f "$PEAK_100K_BAM" ]] || {
      echo "Missing ATAC_BAM: ${PEAK_100K_BAM}" >&2
      return 1
    }
    PEAK_100K_RUN_ROOT="$(peak_100k_run_root_from_frag "${PEAK_100K_FRAG}")"
  elif [[ -n "${FRAGMENTS_TSV_GZ}" ]]; then
    PEAK_100K_FRAG="${FRAGMENTS_TSV_GZ}"
    [[ -f "$PEAK_100K_FRAG" ]] || {
      echo "Missing FRAGMENTS_TSV_GZ: ${PEAK_100K_FRAG}" >&2
      return 1
    }
    PEAK_100K_BAM=""
    PEAK_100K_RUN_ROOT="$(peak_100k_run_root_from_frag "${PEAK_100K_FRAG}")"
  else
    if pout="$(discover_paired_bench)"; then
      PEAK_100K_FRAG="$(echo "${pout}" | sed -n '1p')"
      PEAK_100K_BAM="$(echo "${pout}" | sed -n '2p')"
      PEAK_100K_RUN_ROOT="$(peak_100k_run_root_from_frag "${PEAK_100K_FRAG}")"
    elif [[ "${RUN_MACS3}" == "0" ]] && PEAK_100K_FRAG="$(discover_frags_only_bench)"; then
      PEAK_100K_BAM=""
      PEAK_100K_RUN_ROOT="$(peak_100k_run_root_from_frag "${PEAK_100K_FRAG}")"
    else
      cat >&2 <<EOF
No fragment input. Set one of:
  CHROMAP_PEAK_RUN_ROOT=/path/to/run  (directory containing dual/fragments.tsv[.gz]; BAM under same tree for MACS3)
  FRAGMENTS_TSV_GZ=/path and ATAC_BAM=/path  (same run; required for MACS3)
  FRAGMENTS_TSV_GZ=/path only  with RUN_MACS3=0
Or set CHROMAP_100K_BENCH: auto uses first same-run pair (dual/fragments + atac_possorted_bam in that tree), or with RUN_MACS3=0 first dual/fragments only.
BENCH_ROOT=${BENCH_ROOT}
EOF
      return 1
    fi
  fi
  if [[ "${RUN_MACS3}" != "0" && -z "${PEAK_100K_BAM}" ]]; then
    if [[ ! -f "${PEAK_100K_FRAG:-}" ]]; then
      echo "MACS3 enabled but no fragment file and no BAM. Set FRAGMENTS_TSV_GZ or CHROMAP_PEAK_RUN_ROOT." >&2
      return 1
    fi
    # macs3 callpeak -f FRAG uses the fragments file only (BAM/tagAlign not required).
  fi
  if [[ -n "${PEAK_100K_BAM}" ]]; then
    peak_100k_assert_bam_under_run "${PEAK_100K_RUN_ROOT}" "${PEAK_100K_BAM}" || return 1
  fi
  return 0
}

validate_bed3_sorted() {
  local f=$1
  local label=$2
  [[ -f "${f}" ]] || { echo "Missing ${label}: ${f}" >&2; return 1; }
  if [[ ! -s "${f}" ]]; then
    echo "Empty BED: ${f}" >&2
    return 1
  fi
  awk 'BEGIN{FS=OFS="\t"} !/^#/ && !/^track/ { if(NF<3) exit 2; s=$2+0; e=$3+0; if(s>=e) exit 3; print $1, s, e }' "${f}" > "${f}.3clean" || {
    echo "Invalid intervals in ${label} (${f})" >&2
    return 1
  }
  LC_ALL=C sort -c -k1,1 -k2,2n -k3,3n "${f}.3clean" 2>/dev/null || {
    echo "Not sorted: ${f}" >&2
    return 1
  }
  return 0
}

jaccard_bed3() {
  local a=$1
  local b=$2
  if ! command -v bedtools &>/dev/null; then
    echo "NA"
    return 0
  fi
  bedtools sort -i "${a}" > "${a}.s" 2>/dev/null
  bedtools sort -i "${b}" > "${b}.s" 2>/dev/null
  bedtools jaccard -a "${a}.s" -b "${b}.s" 2>/dev/null | awk 'NR==2 {print $3}'
}

median_width() {
  local f=$1
  awk '{print $3-$2}' "${f}" | sort -n | awk '{
    a[NR]=$0
  } END {
    if(NR==0) { print 0; exit }
    if(NR%2) print a[(NR+1)/2]
    else print (a[NR/2]+a[NR/2+1])/2
  }'
}
