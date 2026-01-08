#!/bin/bash
# End-to-end test for Y-chromosome filtering feature
# Tests three-stream SAM output: all, Y-only, and noY

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Find chromap binary (assume it's in repo root or in PATH)
# Do this BEFORE changing directories
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

# Test directory
TEST_DIR=$(mktemp -d)
cd "$TEST_DIR"

echo "Test directory: $TEST_DIR"
if [ -f "$REPO_ROOT/chromap" ]; then
    CHROMAP="$REPO_ROOT/chromap"
elif command -v chromap &> /dev/null; then
    CHROMAP="chromap"
else
    echo -e "${RED}ERROR: chromap binary not found${NC}"
    echo "  Looked in: $REPO_ROOT/chromap and PATH"
    exit 1
fi

echo "Using chromap: $CHROMAP"
echo "Repo root: $REPO_ROOT"

# Check for samtools
if ! command -v samtools &> /dev/null; then
    echo -e "${YELLOW}WARNING: samtools not found, some tests will be skipped${NC}"
    HAS_SAMTOOLS=false
else
    HAS_SAMTOOLS=true
fi

# Cleanup function
cleanup() {
    cd /
    rm -rf "$TEST_DIR"
}
trap cleanup EXIT

# Test 1: Basic three-stream output verification
echo ""
echo -e "${GREEN}Test 1: Basic three-stream output${NC}"

# Create minimal reference with chr1 and chrY
cat > ref.fa <<EOF
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chrY
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF

# Create reads: some mapping to chr1, some to chrY
cat > reads.fq <<EOF
@read1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read3
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
EOF

# Build index
echo "Building index..."
$CHROMAP --build-index -r ref.fa -o index.idx > /dev/null 2>&1 || {
    echo -e "${RED}FAIL: Index building failed${NC}"
    exit 1
}

# Run mapping with Y-filtering
echo "Running mapping with Y-filtering..."
$CHROMAP --SAM --emit-noY-bam --emit-Y-bam \
    -r ref.fa -x index.idx \
    -1 reads.fq -o output.sam \
    --min-num-seeds 1 --error-threshold 20 \
    > /dev/null 2>&1 || {
    echo -e "${RED}FAIL: Mapping failed${NC}"
    exit 1
}

# Verify output files exist
if [ ! -f "output.sam" ]; then
    echo -e "${RED}FAIL: output.sam not created${NC}"
    exit 1
fi
if [ ! -f "output.noY.sam" ]; then
    echo -e "${RED}FAIL: output.noY.sam not created${NC}"
    exit 1
fi
if [ ! -f "output.Y.sam" ]; then
    echo -e "${RED}FAIL: output.Y.sam not created${NC}"
    exit 1
fi

echo "  ✓ All three output files created"

# Verify files are parseable (if samtools available)
if [ "$HAS_SAMTOOLS" = true ]; then
    samtools view -S output.sam > /dev/null 2>&1 || {
        echo -e "${RED}FAIL: output.sam is not valid SAM${NC}"
        exit 1
    }
    samtools view -S output.noY.sam > /dev/null 2>&1 || {
        echo -e "${RED}FAIL: output.noY.sam is not valid SAM${NC}"
        exit 1
    }
    samtools view -S output.Y.sam > /dev/null 2>&1 || {
        echo -e "${RED}FAIL: output.Y.sam is not valid SAM${NC}"
        exit 1
    }
    echo "  ✓ All SAM files are valid"
    
    # Verify all = Y + noY (read count check)
    TOTAL=$(samtools view -c output.sam 2>/dev/null || echo "0")
    NOY=$(samtools view -c output.noY.sam 2>/dev/null || echo "0")
    Y=$(samtools view -c output.Y.sam 2>/dev/null || echo "0")
    
    if [ "$TOTAL" -ne "$((NOY + Y))" ]; then
        echo -e "${RED}FAIL: Read count mismatch: $TOTAL != $NOY + $Y${NC}"
        exit 1
    fi
    echo "  ✓ Read counts match: $TOTAL = $NOY + $Y"
    
    # Verify headers are identical (@SQ lines)
    HEADER_ALL=$(samtools view -H output.sam | grep '^@SQ' | sort)
    HEADER_NOY=$(samtools view -H output.noY.sam | grep '^@SQ' | sort)
    HEADER_Y=$(samtools view -H output.Y.sam | grep '^@SQ' | sort)
    
    if [ "$HEADER_ALL" != "$HEADER_NOY" ] || [ "$HEADER_ALL" != "$HEADER_Y" ]; then
        echo -e "${RED}FAIL: Headers do not match${NC}"
        echo "All headers:"
        echo "$HEADER_ALL"
        echo "noY headers:"
        echo "$HEADER_NOY"
        echo "Y headers:"
        echo "$HEADER_Y"
        exit 1
    fi
    echo "  ✓ Headers are identical across all streams"
fi

echo -e "${GREEN}Test 1 PASSED${NC}"

# Test 2: Path derivation
echo ""
echo -e "${GREEN}Test 2: Path derivation${NC}"

# Test .sam.gz extension
$CHROMAP --SAM --emit-noY-bam --emit-Y-bam \
    -r ref.fa -x index.idx \
    -1 reads.fq -o output.sam.gz \
    --min-num-seeds 1 --error-threshold 20 \
    > /dev/null 2>&1 || {
    echo -e "${RED}FAIL: Mapping with .sam.gz failed${NC}"
    exit 1
}

if [ ! -f "output.noY.sam.gz" ] || [ ! -f "output.Y.sam.gz" ]; then
    echo -e "${RED}FAIL: Derived paths for .sam.gz not created${NC}"
    exit 1
fi
echo "  ✓ .sam.gz path derivation works"

# Clean up for next test
rm -f output.*.sam*

# Test explicit paths
$CHROMAP --SAM --emit-noY-bam --emit-Y-bam \
    --noY-output custom_noY.sam \
    --Y-output custom_Y.sam \
    -r ref.fa -x index.idx \
    -1 reads.fq -o output.sam \
    --min-num-seeds 1 --error-threshold 20 \
    > /dev/null 2>&1 || {
    echo -e "${RED}FAIL: Mapping with explicit paths failed${NC}"
    exit 1
}

if [ ! -f "custom_noY.sam" ] || [ ! -f "custom_Y.sam" ]; then
    echo -e "${RED}FAIL: Explicit paths not used${NC}"
    exit 1
fi
echo "  ✓ Explicit paths work correctly"

echo -e "${GREEN}Test 2 PASSED${NC}"

# Test 3: No Y contigs in reference
echo ""
echo -e "${GREEN}Test 3: Reference with no Y contigs${NC}"

# Create reference without Y
cat > ref_noY.fa <<EOF
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
>chr2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF

# Build index
$CHROMAP --build-index -r ref_noY.fa -o index_noY.idx > /dev/null 2>&1

# Run mapping (should warn but still create files)
$CHROMAP --SAM --emit-noY-bam --emit-Y-bam \
    -r ref_noY.fa -x index_noY.idx \
    -1 reads.fq -o output_noY_ref.sam \
    --min-num-seeds 1 --error-threshold 20 \
    2>&1 | grep -q "WARNING.*No Y chromosome" || {
    echo -e "${YELLOW}WARNING: Expected warning about no Y contigs not found${NC}"
}

if [ ! -f "output_noY_ref.noY.sam" ] || [ ! -f "output_noY_ref.Y.sam" ]; then
    echo -e "${RED}FAIL: Output files not created when no Y contigs present${NC}"
    exit 1
fi

# Verify Y file is empty (or has only headers)
if [ "$HAS_SAMTOOLS" = true ]; then
    Y_COUNT=$(samtools view -c output_noY_ref.Y.sam 2>/dev/null || echo "0")
    if [ "$Y_COUNT" -ne "0" ]; then
        echo -e "${YELLOW}WARNING: Y file should be empty but has $Y_COUNT reads${NC}"
    fi
fi

echo "  ✓ Files created even when no Y contigs present"
echo -e "${GREEN}Test 3 PASSED${NC}"

# Test 4: Low-memory mode
echo ""
echo -e "${GREEN}Test 4: Low-memory mode${NC}"

$CHROMAP --SAM --emit-noY-bam --emit-Y-bam \
    --low-mem --num-threads 2 \
    -r ref.fa -x index.idx \
    -1 reads.fq -o output_lowmem.sam \
    --min-num-seeds 1 --error-threshold 20 \
    > /dev/null 2>&1 || {
    echo -e "${RED}FAIL: Low-memory mode mapping failed${NC}"
    exit 1
}

if [ ! -f "output_lowmem.sam" ] || [ ! -f "output_lowmem.noY.sam" ] || [ ! -f "output_lowmem.Y.sam" ]; then
    echo -e "${RED}FAIL: Low-memory mode output files not created${NC}"
    exit 1
fi

if [ "$HAS_SAMTOOLS" = true ]; then
    TOTAL_LM=$(samtools view -c output_lowmem.sam 2>/dev/null || echo "0")
    NOY_LM=$(samtools view -c output_lowmem.noY.sam 2>/dev/null || echo "0")
    Y_LM=$(samtools view -c output_lowmem.Y.sam 2>/dev/null || echo "0")
    
    if [ "$TOTAL_LM" -ne "$((NOY_LM + Y_LM))" ]; then
        echo -e "${RED}FAIL: Low-memory read count mismatch: $TOTAL_LM != $NOY_LM + $Y_LM${NC}"
        exit 1
    fi
    echo "  ✓ Low-memory mode read counts match: $TOTAL_LM = $NOY_LM + $Y_LM"
fi

echo -e "${GREEN}Test 4 PASSED${NC}"

echo ""
echo -e "${GREEN}Test 5: Y read names + Y/noY FASTQ emission${NC}"
if command -v python3 &> /dev/null; then
    "$SCRIPT_DIR/test_fastq_noy.sh" || {
        echo -e "${RED}FAIL: fastq_noy end-to-end tests failed${NC}"
        exit 1
    }
    echo -e "${GREEN}Test 5 PASSED${NC}"
else
    echo -e "${YELLOW}WARNING: python3 not found, skipping fastq_noy tests${NC}"
fi

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}All end-to-end tests PASSED!${NC}"
echo -e "${GREEN}========================================${NC}"
