#!/usr/bin/env bash
# Shared helpers for Chromap + libMACS3 publication benchmarks.
# shellcheck shell=bash

set -euo pipefail

bench_repo_root() {
  local script_dir
  script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  (cd "${script_dir}/../.." && pwd)
}

bench_require_file() {
  local path=$1
  local label=$2
  if [[ ! -f "${path}" ]]; then
    echo "ERROR: missing ${label}: ${path}" >&2
    exit 2
  fi
}

bench_require_exe() {
  local path=$1
  local label=$2
  if [[ "${path}" == */* ]]; then
    if [[ ! -x "${path}" ]]; then
      echo "ERROR: missing executable ${label}: ${path}" >&2
      exit 2
    fi
  elif ! command -v "${path}" >/dev/null 2>&1; then
    echo "ERROR: missing executable ${label}: ${path}" >&2
    exit 2
  fi
}

bench_quote_cmd() {
  printf '%q ' "$@"
  printf '\n'
}

bench_log_cmd() {
  local commands_file=$1
  shift
  bench_quote_cmd "$@" >>"${commands_file}"
}

bench_run_timed() {
  local label=$1
  local time_file=$2
  local log_file=$3
  local commands_file=$4
  shift 4
  {
    printf '\n# [%s]\n' "${label}"
    bench_quote_cmd "$@"
  } >>"${commands_file}"
  /usr/bin/time -v -o "${time_file}" "$@" >>"${log_file}" 2>&1
}

bench_run_timed_shell() {
  local label=$1
  local time_file=$2
  local log_file=$3
  local commands_file=$4
  local command_text=$5
  {
    printf '\n# [%s]\n' "${label}"
    printf '%s\n' "${command_text}"
  } >>"${commands_file}"
  /usr/bin/time -v -o "${time_file}" bash -o pipefail -c "${command_text}" >>"${log_file}" 2>&1
}

bench_first_existing() {
  local path
  for path in "$@"; do
    if [[ -e "${path}" ]]; then
      printf '%s' "${path}"
      return 0
    fi
  done
  return 1
}

bench_sort_fragments_5col() {
  local input=$1
  local output=$2
  zcat -f "${input}" |
    awk 'BEGIN{FS=OFS="\t"} NF>=5 {print $1,$2+0,$3+0,$4,$5+0}' |
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5n >"${output}"
}

bench_bed3_from_fragments() {
  local input=$1
  local output=$2
  zcat -f "${input}" |
    awk 'BEGIN{FS=OFS="\t"} NF>=3 {print $1,$2+0,$3+0}' |
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n >"${output}"
}

bench_bed3_from_narrowpeak() {
  local input=$1
  local output=$2
  awk 'BEGIN{FS=OFS="\t"} !/^#/ && !/^track/ && NF>=3 {print $1,$2+0,$3+0}' "${input}" |
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n >"${output}"
}

bench_count_noncomment_lines() {
  local input=$1
  awk 'BEGIN{n=0} !/^#/ && !/^track/ && NF {n++} END{print n+0}' "${input}"
}

bench_md5_file() {
  local input=$1
  md5sum "${input}" | awk '{print $1}'
}

bench_md5_bam_body_sorted() {
  local samtools_bin=$1
  local bam=$2
  "${samtools_bin}" view "${bam}" | LC_ALL=C sort | md5sum | awk '{print $1}'
}

bench_jaccard_or_na() {
  local a=$1
  local b=$2
  if [[ -s "${a}" && -s "${b}" ]] && command -v bedtools >/dev/null 2>&1; then
    bedtools jaccard -a "${a}" -b "${b}" | awk 'NR==2{print $3}'
  else
    printf 'na'
  fi
}

bench_line_set_counts() {
  local a=$1
  local b=$2
  local out_prefix=$3
  LC_ALL=C sort "${a}" >"${out_prefix}.a.linesort"
  LC_ALL=C sort "${b}" >"${out_prefix}.b.linesort"
  local a_only b_only common
  a_only="$(comm -23 "${out_prefix}.a.linesort" "${out_prefix}.b.linesort" | wc -l | awk '{print $1}')"
  b_only="$(comm -13 "${out_prefix}.a.linesort" "${out_prefix}.b.linesort" | wc -l | awk '{print $1}')"
  common="$(comm -12 "${out_prefix}.a.linesort" "${out_prefix}.b.linesort" | wc -l | awk '{print $1}')"
  printf '%s\t%s\t%s\n' "${a_only}" "${b_only}" "${common}"
}

bench_bam_total_records() {
  local samtools_bin=$1
  local bam=$2
  "${samtools_bin}" view -c "${bam}"
}

bench_write_git_state() {
  local out=$1
  local label=$2
  local repo=$3
  {
    printf '[%s]\n' "${label}"
    if git -C "${repo}" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
      printf 'repo\t%s\n' "${repo}"
      printf 'commit\t%s\n' "$(git -C "${repo}" rev-parse HEAD 2>/dev/null || true)"
      printf 'status_short_begin\n'
      git -C "${repo}" status --short 2>/dev/null || true
      printf 'status_short_end\n'
    else
      printf 'repo\t%s\n' "${repo}"
      printf 'commit\tna\n'
      printf 'status_short_begin\nstatus_short_end\n'
    fi
    printf '\n'
  } >>"${out}"
}
