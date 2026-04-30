# Tests

## Core libchromap Smoke

`run_libchromap_core_smoke.sh` is the default small regression gate for the
Chromap CLI and `libchromap` runner boundary. It generates synthetic FASTA and
FASTQ inputs under `plans/artifacts/chromap_core_smoke/<timestamp>/`, builds a
temporary index, runs CLI and `chromap_lib_runner` cases, and compares BED,
TagAlign, pairs, BAM, ATAC fragments, barcode summary, and Y/noY outputs.

Run it with:

```bash
make test-libchromap-core-smoke
```

## ENCODE Downsample Smoke

`prepare_encode_downsample_fixtures.sh` creates ignored real-data fixtures from
the candidate ENCODE pairs listed in `encode_downsample_manifest.tsv`. It
downloads full FASTQs only when `ENCODE_ALLOW_DOWNLOAD=1`, then writes
pair-preserving first-N downsampled FASTQs plus a generated manifest under
`plans/artifacts/encode_fixture_cache/`.

`run_encode_downsample_smoke.sh` consumes that generated manifest and runs
Chromap CLI vs `chromap_lib_runner` on real ChIP, ATAC, and Hi-C read pairs.
It requires a matching full-genome reference and Chromap index:

```bash
CHROMAP_GRCH38_REF=/path/to/genome.fa \
CHROMAP_GRCH38_INDEX=/path/to/genome.index \
ENCODE_ALLOW_DOWNLOAD=1 \
make test-encode-downsample-smoke
```

Useful knobs:

- `ENCODE_ASSAYS=chip,atac,hic` selects a subset; default is `all`.
- `ENCODE_DOWNSAMPLE_READS=10000` overrides the manifest default.
- `ENCODE_FIXTURE_CACHE=/path/to/cache` relocates the ignored cache.
- `ENCODE_SKIP_PREPARE=1` reuses an existing generated manifest/cache.

Hi-C coverage stops at Chromap `.pairs` output. Downstream `.hic`, `.cool`, and
TAD callers are outside this smoke gate.

## MCP Recipe Registry Tests

The MCP/Launchpad hardening work adds recipe metadata under
`mcp_server/recipes/registry.yaml`. These tests check that public workflows have
recipe entries, expected outputs, preflight rule ids, smoke coverage, and docs:

```bash
python3 -m pytest mcp_server/tests/test_recipe_registry.py -q
```

Recipe-driven preflight tests cover missing references/indexes, mismatched
paired FASTQ lanes, barcode lane mismatches, output collisions, untrusted output
locations, and the Hi-C `.pairs` boundary:

```bash
python3 -m pytest mcp_server/tests/test_chromap_preflight.py -q
```

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
