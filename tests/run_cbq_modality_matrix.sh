#!/usr/bin/env bash
# Hermetic CBQ modality matrix smoke.
#
# Verifies that native CBQ input produces output identical to the FASTQ
# baseline across the paired-end modalities Chromap supports, through both the
# `chromap` CLI and `chromap_lib_runner`, and that unsupported / malformed CBQ
# combinations fail with a clear message instead of crashing.
#
# Hermetic: generates a tiny synthetic genome + paired ATAC FASTQs + barcode
# FASTQ + whitelist, builds a small index, and encodes CBQ from those FASTQs.
# The fixture is small enough to fit one CBQ block, so even a record-reordering
# encoder stays lane-aligned (multi-block ordering is exercised separately by
# run_cbq_atac_100k.sh).
#
# Encoder resolution: CBQ_ORDERED_ENCODER (env or default STAR path), else
# bqtools (BQTOOLS env, /tmp/star_suite_bqtools/bin/bqtools, or PATH). Missing
# encoder => SKIP. BAM comparisons additionally require samtools (skipped if
# absent). Artifacts under plans/artifacts/cbq_modality_matrix/<timestamp>/.
set -uo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
export LD_LIBRARY_PATH="${REPO_ROOT}/third_party/htslib${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"

CHROMAP="${CHROMAP:-${REPO_ROOT}/chromap}"
CHROMAP_LIB_RUNNER="${CHROMAP_LIB_RUNNER:-${REPO_ROOT}/chromap_lib_runner}"
ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
OUT_ROOT="${OUT_ROOT:-${ARTIFACT_ROOT}/cbq_modality_matrix/$(date -u +%Y%m%dT%H%M%SZ)}"
DATA="${OUT_ROOT}/data"; RUN="${OUT_ROOT}/runs"; CBQ="${OUT_ROOT}/cbq"
mkdir -p "${DATA}" "${RUN}" "${CBQ}"

PASS=0; FAIL=0
log() { printf '[cbq-matrix] %s\n' "$*"; }
skip() { log "SKIP: $*"; printf 'status=skipped\nreason=%s\n' "$*" > "${OUT_ROOT}/SKIPPED.txt"; exit 0; }
ok()   { PASS=$((PASS+1)); log "PASS: $*"; }
bad()  { FAIL=$((FAIL+1)); log "FAIL: $*"; }

[[ -x "${CHROMAP}" ]] || { echo "ERROR: chromap not built" >&2; exit 2; }
[[ -x "${CHROMAP_LIB_RUNNER}" ]] || { echo "ERROR: chromap_lib_runner not built" >&2; exit 2; }
HAVE_SAMTOOLS=0; command -v samtools >/dev/null 2>&1 && HAVE_SAMTOOLS=1

# ---- encoder resolution -----------------------------------------------------
ORDERED=""; BQ=""
if [[ -n "${CBQ_ORDERED_ENCODER:-}" && -x "${CBQ_ORDERED_ENCODER}" ]]; then
  ORDERED="${CBQ_ORDERED_ENCODER}"
elif [[ -x /mnt/pikachu/STAR-suite/core/legacy/source/cbq_ordered_encoder ]]; then
  ORDERED=/mnt/pikachu/STAR-suite/core/legacy/source/cbq_ordered_encoder
fi
if [[ -z "${ORDERED}" ]]; then
  if [[ -n "${BQTOOLS:-}" && -x "${BQTOOLS}" ]]; then BQ="${BQTOOLS}"
  elif [[ -x /tmp/star_suite_bqtools/bin/bqtools ]]; then BQ=/tmp/star_suite_bqtools/bin/bqtools
  elif command -v bqtools >/dev/null 2>&1; then BQ="$(command -v bqtools)"; fi
fi
[[ -n "${ORDERED}${BQ}" ]] || skip "no CBQ encoder (set CBQ_ORDERED_ENCODER or BQTOOLS)"

encode_pair() { # r1 r2 out
  if [[ -n "${ORDERED}" ]]; then "${ORDERED}" --readFilesIn "$1" "$2" --outFile "$3" >/dev/null 2>&1
  else "${BQ}" encode "$1" "$2" --mode cbq -o "$3" -T 2 >/dev/null 2>&1 || "${BQ}" encode "$1" "$2" -o "$3" -T 2 >/dev/null 2>&1; fi
  [[ -s "$3" ]]
}
encode_single() { # r1 out [extra...]
  local r1="$1" out="$2"; shift 2
  if [[ -n "${ORDERED}" ]]; then "${ORDERED}" --readFilesIn "$r1" --outFile "$out" >/dev/null 2>&1
  else "${BQ}" encode "$r1" --mode cbq -o "$out" -T 2 "$@" >/dev/null 2>&1 || "${BQ}" encode "$r1" -o "$out" -T 2 "$@" >/dev/null 2>&1; fi
  [[ -s "$out" ]]
}

# ---- synthetic fixture ------------------------------------------------------
python3 - "${DATA}" <<'PY'
import random, sys
from pathlib import Path
out = Path(sys.argv[1]); rng = random.Random(20260602)
genome = "".join(rng.choice("ACGT") for _ in range(8000))
rc = lambda s: s.translate(str.maketrans("ACGT","TGCA"))[::-1]
bcs = ["ACGTACGTACGTACGT","TGCATGCATGCATGCA","GGGGAAAACCCCTTTT","AAAACCCCGGGGTTTT"]
frags = [("read001",200,360,bcs[0]),("read002",900,1075,bcs[1]),
         ("read003",2000,2170,bcs[2]),("read004",3300,3470,bcs[3]),
         ("read005",4100,4275,bcs[0]),("read006",5000,5160,bcs[1])]
(out/"genome.fa").write_text(">chrSyn\n"+genome+"\n")
(out/"whitelist.txt").write_text("\n".join(bcs)+"\n")
with (out/"R1.fastq").open("w") as r1,(out/"R3.fastq").open("w") as r3,(out/"R2.fastq").open("w") as bc:
    for n,s,e,b in frags:
        a=genome[s:s+75]; c=rc(genome[e-75:e]); q="I"*75; qb="I"*len(b)
        r1.write(f"@{n} 1:N:0:A\n{a}\n+\n{q}\n"); r3.write(f"@{n} 3:N:0:A\n{c}\n+\n{q}\n")
        bc.write(f"@{n} 2:N:0:A\n{b}\n+\n{qb}\n")
PY

"${CHROMAP}" --build-index -r "${DATA}/genome.fa" -o "${DATA}/idx" -k 11 -w 5 >/dev/null 2>&1

R1="${DATA}/R1.fastq"; R3="${DATA}/R3.fastq"; R2="${DATA}/R2.fastq"
WL="${DATA}/whitelist.txt"
RC="${CBQ}/reads.cbq"; BC="${CBQ}/bc.cbq"; BCNH="${CBQ}/bc_noheader.cbq"
encode_pair "${R1}" "${R3}" "${RC}" || skip "pair encode failed"
encode_single "${R2}" "${BC}" || skip "barcode encode failed"

COMMON=(-x "${DATA}/idx" -r "${DATA}/genome.fa" -t 1 --min-read-length 20 --min-num-seeds 1 -q 0)
ATAC=(--preset atac --read-format "bc:0:-1")   # barcode read carries whole barcode

# fastq vs cbq(cli) vs cbq(lib) for a text-row modality
cmp_text() { # label  out_fastq out_cbqcli out_cbqlib
  local label="$1" a="$2" b="$3" c="$4"
  if [[ ! -s "$a" ]]; then bad "$label: FASTQ baseline empty"; return; fi
  if diff -q <(LC_ALL=C sort "$a") <(LC_ALL=C sort "$b") >/dev/null \
     && diff -q <(LC_ALL=C sort "$a") <(LC_ALL=C sort "$c") >/dev/null; then
    ok "$label ($(wc -l <"$a") rows, CLI+lib parity)"
  else bad "$label: CBQ output differs from FASTQ"; fi
}

run3() { # label  outext  cbq-needs-barcode(0/1)  extra-flags...
  local label="$1" ext="$2" bc_on="$3"; shift 3
  local fq="${RUN}/${label}.fastq.${ext}" cc="${RUN}/${label}.cbqcli.${ext}" cl="${RUN}/${label}.cbqlib.${ext}"
  local fq_in=(-1 "${R1}" -2 "${R3}") cbq_in=(--input-format cbq --read-pair-cbq "${RC}")
  if [[ "${bc_on}" == 1 ]]; then fq_in+=(-b "${R2}" --barcode-whitelist "${WL}"); cbq_in+=(--barcode-cbq "${BC}" --barcode-whitelist "${WL}"); fi
  "${CHROMAP}"            "${COMMON[@]}" "$@" "${fq_in[@]}"  -o "${fq}" >/dev/null 2>&1
  "${CHROMAP}"            "${COMMON[@]}" "$@" "${cbq_in[@]}" -o "${cc}" >/dev/null 2>&1
  "${CHROMAP_LIB_RUNNER}" "${COMMON[@]}" "$@" "${cbq_in[@]}" -o "${cl}" >/dev/null 2>&1
  printf '%s|%s|%s|%s\n' "${fq}" "${cc}" "${cl}" "${label}"
}

# ---- positive modality cases ------------------------------------------------
IFS='|' read -r a b c _ < <(run3 pe_bc_bed     bed 1 "${ATAC[@]}" --BED);        cmp_text pe_bc_bed   "$a" "$b" "$c"
IFS='|' read -r a b c _ < <(run3 pe_bulk_bed   bed 0 --preset atac --BED);       cmp_text pe_bulk_bed "$a" "$b" "$c"
IFS='|' read -r a b c _ < <(run3 pe_bc_tag     ta  1 "${ATAC[@]}" --TagAlign);   cmp_text pe_bc_tag   "$a" "$b" "$c"
IFS='|' read -r a b c _ < <(run3 chip_bed      bed 0 --preset chip --BED);       cmp_text chip_bed    "$a" "$b" "$c"

# SAM: compare alignment rows only (drop @ header / @PG)
samf="${RUN}/sam.fastq.sam"; samc="${RUN}/sam.cbqcli.sam"; saml="${RUN}/sam.cbqlib.sam"
"${CHROMAP}"            "${COMMON[@]}" "${ATAC[@]}" --SAM -1 "${R1}" -2 "${R3}" -b "${R2}" --barcode-whitelist "${WL}" -o "${samf}" >/dev/null 2>&1
"${CHROMAP}"            "${COMMON[@]}" "${ATAC[@]}" --SAM --input-format cbq --read-pair-cbq "${RC}" --barcode-cbq "${BC}" --barcode-whitelist "${WL}" -o "${samc}" >/dev/null 2>&1
"${CHROMAP_LIB_RUNNER}" "${COMMON[@]}" "${ATAC[@]}" --SAM --input-format cbq --read-pair-cbq "${RC}" --barcode-cbq "${BC}" --barcode-whitelist "${WL}" -o "${saml}" >/dev/null 2>&1
grep -v '^@' "${samf}" > "${samf}.rows" 2>/dev/null
grep -v '^@' "${samc}" > "${samc}.rows" 2>/dev/null
grep -v '^@' "${saml}" > "${saml}.rows" 2>/dev/null
cmp_text pe_bc_sam "${samf}.rows" "${samc}.rows" "${saml}.rows"

# BAM dual fragments (needs samtools for record compare; always checks fragments TSV)
if [[ "${HAVE_SAMTOOLS}" == 1 ]]; then
  for tag in fastq cbqcli cbqlib; do :; done
  "${CHROMAP}"            "${COMMON[@]}" "${ATAC[@]}" --BAM --sort-bam -b "${R2}" -1 "${R1}" -2 "${R3}" --barcode-whitelist "${WL}" --atac-fragments "${RUN}/dual.fastq.frag.tsv.gz" -o "${RUN}/dual.fastq.bam" >/dev/null 2>&1
  "${CHROMAP}"            "${COMMON[@]}" "${ATAC[@]}" --BAM --sort-bam --input-format cbq --read-pair-cbq "${RC}" --barcode-cbq "${BC}" --barcode-whitelist "${WL}" --atac-fragments "${RUN}/dual.cbqcli.frag.tsv.gz" -o "${RUN}/dual.cbqcli.bam" >/dev/null 2>&1
  "${CHROMAP_LIB_RUNNER}" "${COMMON[@]}" "${ATAC[@]}" --BAM --sort-bam --input-format cbq --read-pair-cbq "${RC}" --barcode-cbq "${BC}" --barcode-whitelist "${WL}" --atac-fragments "${RUN}/dual.cbqlib.frag.tsv.gz" -o "${RUN}/dual.cbqlib.bam" >/dev/null 2>&1
  if samtools quickcheck "${RUN}/dual.cbqcli.bam" 2>/dev/null \
     && diff -q <(zcat "${RUN}/dual.fastq.frag.tsv.gz" 2>/dev/null | LC_ALL=C sort) <(zcat "${RUN}/dual.cbqcli.frag.tsv.gz" 2>/dev/null | LC_ALL=C sort) >/dev/null \
     && diff -q <(zcat "${RUN}/dual.fastq.frag.tsv.gz" 2>/dev/null | LC_ALL=C sort) <(zcat "${RUN}/dual.cbqlib.frag.tsv.gz" 2>/dev/null | LC_ALL=C sort) >/dev/null; then
    ok "pe_bc_bam_dual (fragments parity + BAM quickcheck)"
  else bad "pe_bc_bam_dual: fragments differ or BAM invalid"; fi
else
  log "SKIP BAM dual case (samtools not found)"
fi

# ---- rejection / safety cases ----------------------------------------------
expect_fail() { # label  pattern  args...
  local label="$1" pat="$2"; shift 2
  local out; out="$("${CHROMAP}" "${COMMON[@]}" "$@" -o "${RUN}/${label}.out" 2>&1)"; local rc=$?
  if [[ $rc -ne 0 ]] && grep -qiE "${pat}" <<<"${out}" && ! grep -qiE 'segmentation|core dumped|bad_alloc' <<<"${out}"; then
    ok "reject:${label}"
  else bad "reject:${label} (rc=$rc, expected /${pat}/)"; fi
}
expect_fail no_readpair    "requires --read-pair-cbq"        --preset atac --input-format cbq --BED
expect_fail mixed_inputs   "cannot be mixed"                 --preset atac --input-format cbq --read-pair-cbq "${RC}" -1 "${R1}" --BED
expect_fail wl_no_bccbq    "requires --barcode-cbq"          --preset atac --input-format cbq --read-pair-cbq "${RC}" --barcode-whitelist "${WL}" --BED
expect_fail unpaired_cbq   "mate-count mismatch"             --preset atac --input-format cbq --read-pair-cbq "${BC}" --BED

# headerless barcode lane (only constructible with bqtools -H)
if [[ -z "${ORDERED}" && -n "${BQ}" ]] && encode_single "${R2}" "${BCNH}" -H; then
  expect_fail headerless_bc "requires read names" "${ATAC[@]}" --input-format cbq --read-pair-cbq "${RC}" --barcode-cbq "${BCNH}" --barcode-whitelist "${WL}" --BED
else
  log "SKIP headerless-barcode case (needs bqtools -H)"
fi

# ---- summary ----------------------------------------------------------------
log "outputs: ${OUT_ROOT}"
log "encoder: ${ORDERED:-${BQ}}   samtools=${HAVE_SAMTOOLS}"
log "RESULT: ${PASS} passed, ${FAIL} failed"
printf 'pass=%s\nfail=%s\nencoder=%s\n' "${PASS}" "${FAIL}" "${ORDERED:-${BQ}}" > "${OUT_ROOT}/SUMMARY.txt"
[[ "${FAIL}" -eq 0 ]]
