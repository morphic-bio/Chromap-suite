# Tests for Y-Chromosome Filtering Feature

This directory contains tests for the three-stream SAM Y-filtering feature.

## Unit Tests

**File**: `test_y_filter.cc`

Tests the Y contig detection logic (`BuildYContigRidMask`) with various contig naming conventions.

### Running Unit Tests

```bash
make test-unit
```

Or manually:
```bash
g++ -std=c++11 -I../src tests/test_y_filter.cc src/sequence_batch.cc -o tests/test_y_filter -lm -lz
./tests/test_y_filter
```

### Test Coverage

- Y contig detection with various naming conventions (chrY, Y, CHRY, etc.)
- Case-insensitive matching
- Exact match requirement (chrY_random should NOT match)
- Multiple contigs including Y
- Reference with no Y contigs

## End-to-End Tests

**File**: `e2e/test_y_streams.sh`

Tests the complete Y-filtering workflow including:
- Three-stream output (all, Y-only, noY)
- Read count verification (`all = Y + noY`)
- Header identity verification
- Path derivation (.sam, .sam.gz, explicit paths)
- No Y contigs handling
- Low-memory mode compatibility

### Running End-to-End Tests

```bash
cd tests/e2e
./test_y_streams.sh
```

### Requirements

- `chromap` binary (built from source or in PATH)
- `samtools` (required for `test_lowmem_bam_y.sh`, optional for SAM validation/read counting in `test_y_streams.sh`)
- `python3` (required for `test_fastq_noy.sh` and `test_lowmem_bam_y.sh`, optional when invoked by `test_y_streams.sh`)

### Test Output

The script creates a temporary directory, runs chromap with various configurations, and verifies:
1. All three output files are created
2. SAM files are valid (if samtools available)
3. Read counts match: `total = Y + noY`
4. Headers are identical across all streams
5. Path derivation works correctly
6. Low-memory mode works correctly

**File**: `e2e/test_fastq_noy.sh`

Tests the Y read names list and Y/noY FASTQ splitting:
- Single-end Y/noY FASTQ emission and name normalization
- Paired-end routing when either mate hits Y
- Multiple input files naming with `.fN` suffixes

### Running End-to-End FASTQ Tests

```bash
cd tests/e2e
./test_fastq_noy.sh
```

**File**: `e2e/test_lowmem_bam_y.sh`

Tests low-memory BAM output with Y/noY split and sorting:
- Low-memory BAM output routing (Y-only vs noY)
- Count conservation (`all = Y + noY`)
- Coordinate-sorted output validation for `--sort-bam`

### Running Low-Memory BAM Tests

```bash
cd tests/e2e
./test_lowmem_bam_y.sh
```

## Test Results

All tests should pass. If any test fails:
1. Check that chromap is built and accessible
2. Verify samtools is installed (for full validation)
3. Check that the test directory has write permissions
4. Review error messages for specific failure points
