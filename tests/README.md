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

## Input Format Smoke

`run_input_format_smoke.sh` verifies the read-input surface independently of
production fixtures. It generates synthetic FASTA/FASTQ inputs, builds a small
Chromap index, and checks:

- paired and single-end FASTQ plain vs `.gz` input parity;
- Y/noY FASTQ sidecar output with `--emit-Y-noY-fastq-compression none` and
  `gz`;
- optional BINSEQ CBQ compatibility through `bqtools decode`, for both default
  compressed and uncompressed (`-l 0`) CBQ files.

Run it with:

```bash
make test-input-format-smoke
```

Set `BQTOOLS=/path/to/bqtools` to exercise the BINSEQ cases when `bqtools` is
not on `PATH`. BINSEQ is currently validated as a decode-to-FASTQ compatibility
path there; native ATAC CBQ ingestion is covered by the smoke below.

## CBQ ATAC Smoke

`run_cbq_atac_smoke.sh` verifies native ATAC CBQ ingestion. It generates a small
paired-end scATAC fixture, encodes read-pair and barcode CBQ files with
`bqtools`, compares CLI FASTQ vs CLI CBQ fragments, and compares CLI FASTQ vs
`chromap_lib_runner` CBQ fragments.

Run it with:

```bash
make test-cbq-atac-smoke
```

Set `BQTOOLS=/path/to/bqtools` when `bqtools` is not on `PATH`.

## CBQ ATAC 100K Parity Gate

`run_cbq_atac_100k.sh` is the pre-merge scale gate. It encodes the 100K PBMC
ATAC fixture lanes to CBQ and checks that barcoded paired-end fragments are
byte-identical (under `LC_ALL=C` sort) across FASTQ baseline, native CBQ via
`chromap`, and native CBQ via `chromap_lib_runner`. It records per-run
wall/user/sys/max-RSS and a manifest under
`plans/artifacts/cbq_atac_100k/<timestamp>/`.

Run it with:

```bash
make test-cbq-atac-100k
```

This gate uses the order-preserving `cbq_ordered_encoder` rather than `bqtools`:
the barcoded path requires the read-pair and barcode CBQ lanes to stay
record-aligned, and stock `bqtools` reorders records across blocks at scale.
Set `CBQ_ORDERED_ENCODER=/path/to/cbq_ordered_encoder` if it is not at the
default location. `LANES="1"` runs a faster single-lane variant. Missing the
encoder or any fixture input is reported as a skip.

The current CBQ reader decodes sequence bases directly into Chromap
`SequenceBatch` buffers, matching the FASTQ path's buffer shape and avoiding a
per-record temporary sequence string. A full 4-lane PBMC 3K ATAC timing run
with warmed cache, 8 threads, and output directed to `/dev/null` is recorded at
`plans/artifacts/cbq_atac_full_timing/20260530T070035Z/`: FASTQ.gz took
`3:08.74`, uncompressed CBQ took `2:56.36`, and both produced `53,969,811`
output mappings with matching core summary totals.

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

## ENCODE Cross-Assay Smoke

`prepare_encode_cross_assay_fixtures.sh` and
`run_encode_cross_assay_smoke.sh` are the S1 real-data regression suite for
Chromap development across ChIP-seq, bulk ATAC-seq, scATAC/snATAC, and Hi-C.
The source rows live in `encode_cross_assay_manifest.tsv`; full ENCODE downloads
and deterministic first-N downsampled FASTQs stay under
`plans/artifacts/encode_cross_assay_cache/`.

Run all cases with:

```bash
CHROMAP_GRCH38_REF=/path/to/genome.fa \
CHROMAP_GRCH38_INDEX=/path/to/genome.index \
ENCODE_ALLOW_DOWNLOAD=1 \
make test-encode-cross-assay-smoke
```

Useful knobs:

- `ENCODE_ASSAYS=chip,atac,scatac,hic` selects a subset; default is `all`.
- `ENCODE_DOWNSAMPLE_READS=10000` overrides the manifest read count.
- `ENCODE_CROSS_ASSAY_CACHE=/path/to/cache` relocates ignored downloads and
  downsampled FASTQs.
- `ENCODE_SKIP_PREPARE=1` reuses an existing generated manifest/cache.
- `SCATAC_WHITELIST=/path/to/whitelist.txt` overrides the manifest whitelist.

The runner compares canonical non-comment, non-empty text rows between
`chromap` and `chromap_lib_runner`. Hi-C stops at `.pairs`; scATAC uses the
barcode/index FASTQ and also checks summary files contain barcode totals. A
whitelist is used only when the manifest provides one; the current ENCODE
snATAC row has no public matching whitelist. See
`plans/2026-05-29-encode-cross-assay-smoke-runbook.md` for the implementation
plan and source-accession rationale.

`run_encode_cbq_cross_assay_smoke.sh` layers CBQ parity on top of the same
generated ENCODE fixtures. It encodes each downsampled paired read set with the
order-preserving `cbq_ordered_encoder`, then compares canonical FASTQ,
`chromap` CBQ, and `chromap_lib_runner` CBQ rows for ChIP, bulk ATAC, scATAC,
and Hi-C `.pairs` output:

```bash
CHROMAP_GRCH38_REF=/path/to/genome.fa \
CHROMAP_GRCH38_INDEX=/path/to/genome.index \
ENCODE_SKIP_PREPARE=1 \
make test-encode-cbq-cross-assay-smoke
```

Use `ENCODE_CBQ_ASSAYS=hic` to run only Hi-C CBQ parity, and
`ENCODE_CBQ_COMPRESSION_LEVEL=0` to store uncompressed CBQs.

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

Recipe run-manifest tests cover dry-run manifest creation, MCP tool exposure,
serial benchmark policy recording, and Launchpad launch log/manifest plumbing:

```bash
python3 -m pytest mcp_server/tests/test_run_manifest.py -q
```

Launchpad tests cover registry-driven recipe listing, metadata-only gating,
recipe form schema generation, registry preflight, render, and dry-run manifest
creation:

```bash
python3 -m pytest mcp_server/tests/test_launchpad_chromap.py -q
```

S0 is the default local gate: MCP tests plus `make test-libchromap-core-smoke`.
S1 ENCODE tests are opt-in because they may download real data and require a
full reference/index; set `ENCODE_ALLOW_DOWNLOAD=1` only when you intend that.
S2 covers longer integration and benchmark runs. Benchmark timings should be
serial unless a runbook explicitly says otherwise.

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
