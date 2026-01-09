#!/bin/bash
# End-to-end tests for low-memory BAM output with Y/noY split + sorting

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Find chromap binary (assume it's in repo root or in PATH)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
CHROMAP=""
if [ -f "$REPO_ROOT/chromap" ]; then
  CHROMAP="$REPO_ROOT/chromap"
elif command -v chromap &> /dev/null; then
  CHROMAP="chromap"
else
  echo -e "${RED}ERROR: chromap binary not found${NC}"
  echo "  Looked in: $REPO_ROOT/chromap and PATH"
  exit 1
fi

if ! command -v python3 &> /dev/null; then
  echo -e "${RED}ERROR: python3 not found (required for test data generation)${NC}"
  exit 1
fi
if ! command -v samtools &> /dev/null; then
  echo -e "${RED}ERROR: samtools not found (required for BAM validation)${NC}"
  exit 1
fi

# Test directory
TEST_DIR=$(mktemp -d)
DATA_DIR="$TEST_DIR/data"
mkdir -p "$DATA_DIR"
cd "$TEST_DIR"

echo "Test directory: $TEST_DIR"
echo "Using chromap: $CHROMAP"
echo "Repo root: $REPO_ROOT"

# Cleanup function
cleanup() {
  cd /
  rm -rf "$TEST_DIR"
}
trap cleanup EXIT

# Generate deterministic test data (chr1 + chrY, paired-end)
python3 "$REPO_ROOT/tests/data/generate_test_data_y.py" -o "$DATA_DIR" > /dev/null

echo "Building index..."
$CHROMAP --build-index -r "$DATA_DIR/test_ref_y.fa" -o "$DATA_DIR/test_ref_y.idx" > /dev/null 2>&1 || {
  echo -e "${RED}FAIL: Index building failed${NC}"
  exit 1
}

count_bam_reads() {
  samtools view -c "$1" 2>/dev/null || echo "0"
}

count_mismatched_contigs() {
  local bam="$1"
  local expected="$2"
  samtools view "$bam" | awk -v expected_contig="$expected" '$3 != expected_contig {bad++} END {print bad+0}'
}

check_sorted_bam() {
  local bam="$1"
  local order_file="$TEST_DIR/contig_order.txt"
  samtools view -H "$bam" | awk '$1=="@SQ"{for(i=1;i<=NF;i++){if($i ~ /^SN:/){sub(/^SN:/,"",$i); print $i}}}' > "$order_file"
  samtools view "$bam" | awk -v order_file="$order_file" '
    BEGIN {
      while ((getline < order_file) > 0) { order[$1]=++idx }
    }
    {
      rname=$3; pos=$4;
      if (rname == "*") { rorder=1000000000; pos=1000000000; }
      else { rorder=order[rname]; }
      if (prev_set && (rorder < prev_order || (rorder == prev_order && pos < prev_pos))) {
        bad=1; exit 1;
      }
      prev_order=rorder; prev_pos=pos; prev_set=1;
    }
    END { exit bad?1:0 }
  '
}

echo ""
echo -e "${GREEN}Baseline: Standard BAM mapping${NC}"

$CHROMAP --BAM --num-threads 2 \
  -r "$DATA_DIR/test_ref_y.fa" -x "$DATA_DIR/test_ref_y.idx" \
  -1 "$DATA_DIR/test_pe_y_R1.fq" -2 "$DATA_DIR/test_pe_y_R2.fq" \
  -o baseline.bam \
  --min-num-seeds 1 --error-threshold 10 --max-insert-size 1000 \
  > /dev/null 2>&1 || {
  echo -e "${RED}FAIL: Baseline BAM mapping failed${NC}"
  exit 1
}

if [ ! -f "baseline.bam" ]; then
  echo -e "${RED}FAIL: baseline.bam not created${NC}"
  exit 1
fi

TOTAL_BASE=$(count_bam_reads "baseline.bam")
if [ "$TOTAL_BASE" -le 0 ]; then
  echo -e "${RED}FAIL: Baseline BAM has no mapped reads${NC}"
  exit 1
fi

echo ""
echo -e "${GREEN}Test 1: Low-memory BAM Y/noY split${NC}"

$CHROMAP --BAM --low-mem --num-threads 2 --emit-noY-bam --emit-Y-bam \
  -r "$DATA_DIR/test_ref_y.fa" -x "$DATA_DIR/test_ref_y.idx" \
  -1 "$DATA_DIR/test_pe_y_R1.fq" -2 "$DATA_DIR/test_pe_y_R2.fq" \
  -o lowmem.bam \
  --min-num-seeds 1 --error-threshold 10 --max-insert-size 1000 \
  > /dev/null 2>&1 || {
  echo -e "${RED}FAIL: Low-memory BAM mapping failed${NC}"
  exit 1
}

if [ ! -f "lowmem.bam" ] || [ ! -f "lowmem.noY.bam" ] || [ ! -f "lowmem.Y.bam" ]; then
  echo -e "${RED}FAIL: Low-memory BAM outputs missing${NC}"
  exit 1
fi

TOTAL_LM=$(count_bam_reads "lowmem.bam")
NOY_LM=$(count_bam_reads "lowmem.noY.bam")
Y_LM=$(count_bam_reads "lowmem.Y.bam")

if [ "$TOTAL_LM" -ne "$TOTAL_BASE" ]; then
  echo -e "${RED}FAIL: Low-memory BAM total mismatch vs baseline: $TOTAL_LM != $TOTAL_BASE${NC}"
  exit 1
fi
if [ "$TOTAL_LM" -ne "$((NOY_LM + Y_LM))" ]; then
  echo -e "${RED}FAIL: Low-memory BAM count mismatch: $TOTAL_LM != $NOY_LM + $Y_LM${NC}"
  exit 1
fi

BAD_Y=$(count_mismatched_contigs "lowmem.Y.bam" "chrY")
BAD_NOY=$(count_mismatched_contigs "lowmem.noY.bam" "chr1")
if [ "$BAD_Y" -ne 0 ] || [ "$BAD_NOY" -ne 0 ]; then
  echo -e "${RED}FAIL: Contig mismatch in Y/noY BAM outputs (Y_bad=$BAD_Y, noY_bad=$BAD_NOY)${NC}"
  exit 1
fi

echo -e "${GREEN}Test 1 PASSED${NC}"

echo ""
echo -e "${GREEN}Test 2: Low-memory BAM + sorting${NC}"

$CHROMAP --BAM --low-mem --num-threads 2 --sort-bam --emit-noY-bam --emit-Y-bam \
  -r "$DATA_DIR/test_ref_y.fa" -x "$DATA_DIR/test_ref_y.idx" \
  -1 "$DATA_DIR/test_pe_y_R1.fq" -2 "$DATA_DIR/test_pe_y_R2.fq" \
  -o lowmem_sorted.bam \
  --min-num-seeds 1 --error-threshold 10 --max-insert-size 1000 \
  > /dev/null 2>&1 || {
  echo -e "${RED}FAIL: Low-memory sorted BAM mapping failed${NC}"
  exit 1
}

if [ ! -f "lowmem_sorted.bam" ] || [ ! -f "lowmem_sorted.noY.bam" ] || [ ! -f "lowmem_sorted.Y.bam" ]; then
  echo -e "${RED}FAIL: Low-memory sorted BAM outputs missing${NC}"
  exit 1
fi

TOTAL_SORTED=$(count_bam_reads "lowmem_sorted.bam")
NOY_SORTED=$(count_bam_reads "lowmem_sorted.noY.bam")
Y_SORTED=$(count_bam_reads "lowmem_sorted.Y.bam")

if [ "$TOTAL_SORTED" -ne "$TOTAL_BASE" ]; then
  echo -e "${RED}FAIL: Sorted BAM total mismatch vs baseline: $TOTAL_SORTED != $TOTAL_BASE${NC}"
  exit 1
fi
if [ "$TOTAL_SORTED" -ne "$((NOY_SORTED + Y_SORTED))" ]; then
  echo -e "${RED}FAIL: Sorted BAM count mismatch: $TOTAL_SORTED != $NOY_SORTED + $Y_SORTED${NC}"
  exit 1
fi

BAD_Y_SORTED=$(count_mismatched_contigs "lowmem_sorted.Y.bam" "chrY")
BAD_NOY_SORTED=$(count_mismatched_contigs "lowmem_sorted.noY.bam" "chr1")
if [ "$BAD_Y_SORTED" -ne 0 ] || [ "$BAD_NOY_SORTED" -ne 0 ]; then
  echo -e "${RED}FAIL: Contig mismatch in sorted Y/noY BAM outputs (Y_bad=$BAD_Y_SORTED, noY_bad=$BAD_NOY_SORTED)${NC}"
  exit 1
fi

check_sorted_bam "lowmem_sorted.bam" || {
  echo -e "${RED}FAIL: lowmem_sorted.bam is not coordinate-sorted${NC}"
  exit 1
}
check_sorted_bam "lowmem_sorted.noY.bam" || {
  echo -e "${RED}FAIL: lowmem_sorted.noY.bam is not coordinate-sorted${NC}"
  exit 1
}
check_sorted_bam "lowmem_sorted.Y.bam" || {
  echo -e "${RED}FAIL: lowmem_sorted.Y.bam is not coordinate-sorted${NC}"
  exit 1
}

echo -e "${GREEN}Test 2 PASSED${NC}"

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}Low-memory BAM Y/noY tests PASSED!${NC}"
echo -e "${GREEN}========================================${NC}"
