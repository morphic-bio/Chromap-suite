#!/bin/bash
# End-to-end tests for Y read names + Y/noY FASTQ emission

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

# Generate deterministic test data
python3 "$REPO_ROOT/tests/data/generate_test_data.py" -o "$DATA_DIR" > /dev/null

echo "Building index..."
$CHROMAP --build-index -r "$DATA_DIR/test_ref.fa" -o "$DATA_DIR/test_ref.idx" > /dev/null 2>&1 || {
  echo -e "${RED}FAIL: Index building failed${NC}"
  exit 1
}

count_fastq_reads() {
  if [ ! -f "$1" ]; then
    echo 0
    return
  fi
  grep -c '^@' "$1" || true
}

SE_Y=15
SE_NOY=35
PE_Y=10
PE_NOY=33

echo ""
echo -e "${GREEN}Test 1: Single-end Y read names + FASTQ split${NC}"

cp "$DATA_DIR/test_se.fq" single_R1.fq

$CHROMAP --SAM --emit-Y-read-names --emit-Y-noY-fastq \
  --emit-Y-noY-fastq-compression none \
  -r "$DATA_DIR/test_ref.fa" -x "$DATA_DIR/test_ref.idx" \
  -1 single_R1.fq -o single.sam \
  --min-num-seeds 1 --error-threshold 10 -t 1 \
  > /dev/null 2>&1 || {
  echo -e "${RED}FAIL: Mapping failed in Test 1${NC}"
  exit 1
}

if [ ! -f "single.Y.names.txt" ]; then
  echo -e "${RED}FAIL: Y read names output missing${NC}"
  exit 1
fi
if [ ! -f "y_separated/single_Y_R1.fq" ] || [ ! -f "y_separated/single_noY_R1.fq" ]; then
  echo -e "${RED}FAIL: Y/noY FASTQ outputs missing${NC}"
  exit 1
fi

Y_COUNT=$(count_fastq_reads "y_separated/single_Y_R1.fq")
NOY_COUNT=$(count_fastq_reads "y_separated/single_noY_R1.fq")
if [ "$Y_COUNT" -ne "$SE_Y" ] || [ "$NOY_COUNT" -ne "$SE_NOY" ]; then
  echo -e "${RED}FAIL: FASTQ counts mismatch in Test 1 (Y=$Y_COUNT, noY=$NOY_COUNT)${NC}"
  exit 1
fi

if [ "$(wc -l < single.Y.names.txt)" -ne "$SE_Y" ]; then
  echo -e "${RED}FAIL: Unexpected Y read names count in Test 1${NC}"
  exit 1
fi

echo -e "${GREEN}Test 1 PASSED${NC}"

echo ""
echo -e "${GREEN}Test 2: Paired-end Y/noY FASTQ split${NC}"

cp "$DATA_DIR/test_pe_R1.fq" paired_R1.fq
cp "$DATA_DIR/test_pe_R2.fq" paired_R2.fq

$CHROMAP --SAM --emit-Y-read-names --emit-Y-noY-fastq \
  --emit-Y-noY-fastq-compression none \
  -r "$DATA_DIR/test_ref.fa" -x "$DATA_DIR/test_ref.idx" \
  -1 paired_R1.fq -2 paired_R2.fq -o paired.sam \
  --min-num-seeds 1 --error-threshold 10 --max-insert-size 1000 -t 1 \
  > /dev/null 2>&1 || {
  echo -e "${RED}FAIL: Mapping failed in Test 2${NC}"
  exit 1
}

if [ ! -f "y_separated/paired_Y_R1.fq" ] || [ ! -f "y_separated/paired_Y_R2.fq" ] || \
   [ ! -f "y_separated/paired_noY_R1.fq" ] || [ ! -f "y_separated/paired_noY_R2.fq" ]; then
  echo -e "${RED}FAIL: Paired-end FASTQ outputs missing${NC}"
  exit 1
fi

Y_R1=$(count_fastq_reads "y_separated/paired_Y_R1.fq")
Y_R2=$(count_fastq_reads "y_separated/paired_Y_R2.fq")
NOY_R1=$(count_fastq_reads "y_separated/paired_noY_R1.fq")
NOY_R2=$(count_fastq_reads "y_separated/paired_noY_R2.fq")

if [ "$Y_R1" -ne "$PE_Y" ] || [ "$Y_R2" -ne "$PE_Y" ] || \
   [ "$NOY_R1" -ne "$PE_NOY" ] || [ "$NOY_R2" -ne "$PE_NOY" ]; then
  echo -e "${RED}FAIL: Paired-end FASTQ counts mismatch (Y_R1=$Y_R1, Y_R2=$Y_R2, noY_R1=$NOY_R1, noY_R2=$NOY_R2)${NC}"
  exit 1
fi

if [ "$(wc -l < paired.Y.names.txt)" -ne "$PE_Y" ]; then
  echo -e "${RED}FAIL: Unexpected Y read names count in Test 2${NC}"
  exit 1
fi

if grep -q '/1' paired.Y.names.txt; then
  echo -e "${RED}FAIL: Unnormalized read names found in Test 2${NC}"
  exit 1
fi

echo -e "${GREEN}Test 2 PASSED${NC}"

echo ""
echo -e "${GREEN}Test 3: Multiple input files naming (.fN)${NC}"

cp "$DATA_DIR/test_se.fq" multiA_R1.fq
cp "$DATA_DIR/test_se.fq" multiB_R1.fq

$CHROMAP --SAM --emit-Y-noY-fastq \
  --emit-Y-noY-fastq-compression none \
  -r "$DATA_DIR/test_ref.fa" -x "$DATA_DIR/test_ref.idx" \
  -1 multiA_R1.fq -1 multiB_R1.fq -o multi.sam \
  --min-num-seeds 1 --error-threshold 10 -t 1 \
  > /dev/null 2>&1 || {
  echo -e "${RED}FAIL: Mapping failed in Test 3${NC}"
  exit 1
}

if [ ! -f "y_separated/multiA_Y_R1.f1.fq" ] || [ ! -f "y_separated/multiA_noY_R1.f1.fq" ] || \
   [ ! -f "y_separated/multiB_Y_R1.f2.fq" ] || [ ! -f "y_separated/multiB_noY_R1.f2.fq" ]; then
  echo -e "${RED}FAIL: Multi-input FASTQ outputs missing${NC}"
  exit 1
fi

Y_F1=$(count_fastq_reads "y_separated/multiA_Y_R1.f1.fq")
NOY_F1=$(count_fastq_reads "y_separated/multiA_noY_R1.f1.fq")
Y_F2=$(count_fastq_reads "y_separated/multiB_Y_R1.f2.fq")
NOY_F2=$(count_fastq_reads "y_separated/multiB_noY_R1.f2.fq")

if [ "$Y_F1" -ne "$SE_Y" ] || [ "$NOY_F1" -ne "$SE_NOY" ] || \
   [ "$Y_F2" -ne "$SE_Y" ] || [ "$NOY_F2" -ne "$SE_NOY" ]; then
  echo -e "${RED}FAIL: Multi-input FASTQ counts mismatch (f1: Y=$Y_F1 noY=$NOY_F1, f2: Y=$Y_F2 noY=$NOY_F2)${NC}"
  exit 1
fi

echo -e "${GREEN}Test 3 PASSED${NC}"

echo ""
echo -e "${GREEN}Test 4: Explicit Y/noY FASTQ output directory${NC}"

$CHROMAP --SAM --emit-Y-noY-fastq \
  --emit-Y-noY-fastq-compression none \
  --Y-noY-fastq-output-dir custom_split \
  -r "$DATA_DIR/test_ref.fa" -x "$DATA_DIR/test_ref.idx" \
  -1 single_R1.fq -o custom_dir.sam \
  --min-num-seeds 1 --error-threshold 10 -t 1 \
  > /dev/null 2>&1 || {
  echo -e "${RED}FAIL: Mapping failed in Test 4${NC}"
  exit 1
}

if [ ! -f "custom_split/single_Y_R1.fq" ] || [ ! -f "custom_split/single_noY_R1.fq" ]; then
  echo -e "${RED}FAIL: Explicit output directory was not used${NC}"
  exit 1
fi

echo -e "${GREEN}Test 4 PASSED${NC}"

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}All fastq_noy end-to-end tests PASSED!${NC}"
echo -e "${GREEN}========================================${NC}"
