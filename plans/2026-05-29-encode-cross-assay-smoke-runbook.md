# ENCODE Cross-Assay Smoke Runbook

Date: 2026-05-29

## Goal

Build an opt-in real-data smoke suite for Chromap-suite development across the
assay types Chromap advertises:

- ChIP-seq,
- bulk ATAC-seq,
- single-cell / single-nucleus ATAC-seq,
- Hi-C pairs output.

The suite should catch regressions in shared CLI, `MappingParameters`,
`SequenceBatch`, `chromap_lib_runner`, and output-format plumbing as Chromap
Suite expands beyond ATAC-specific work. It is not a benchmark suite and not a
downstream analysis pipeline.

## Scope

This runbook covers Tier S1 real-data smokes. It should stay opt-in because it
may download ENCODE FASTQs and requires a host with a matching full-genome
reference and Chromap index.

In scope:

- ENCODE source manifest with pinned experiments/files.
- Deterministic first-N downsampling.
- Pair/triplet read-name validation before running Chromap.
- CLI vs `chromap_lib_runner` parity where supported.
- Small nonempty-output and shape checks for each assay.
- Artifacts and generated manifests under `plans/artifacts/`.

Out of scope:

- Committing FASTQs, indexes, alignments, or downloaded ENCODE data.
- Running downstream Hi-C `.cool` / `.hic` / TAD tools.
- STAR-suite multiome orchestration.
- Making this part of the default pre-commit gate.
- CBQ parity for non-ATAC assays. That can be layered on later after FASTQ
  cross-assay coverage is stable.

## Existing Starting Point

The repo already has the first version of this idea:

- `tests/encode_downsample_manifest.tsv`
- `tests/prepare_encode_downsample_fixtures.sh`
- `tests/run_encode_downsample_smoke.sh`
- `make test-encode-downsample-smoke`

Current coverage is paired FASTQ only:

| Assay | Case | Experiment | R1 | R2 | Current args |
| --- | --- | --- | --- | --- | --- |
| ChIP-seq | `C03_chip_tagalign` | `ENCSR089DTY` | `ENCFF851YSO` | `ENCFF451LIG` | `--preset chip --TagAlign` |
| bulk ATAC-seq | `C04_atac_bed` | `ENCSR803FKU` | `ENCFF013DXQ` | `ENCFF473SZG` | `--preset atac --BED` |
| Hi-C | `C10_hic_pairs` | `ENCSR173DMV` | `ENCFF396ZOI` | `ENCFF995XRN` | `--split-alignment --pairs --MAPQ-threshold 1` |

Keep this path working while extending it for scATAC barcode/index FASTQs.

## Source Candidates

All source data should be pulled from ENCODE by accession and recorded in a
generated manifest. Prefer released FASTQ files with `output_type=reads`, and
use ENCODE `paired_with` metadata to confirm mate relationships.

### ChIP-seq

Use the existing manifest row.

- Experiment: `ENCSR089DTY`
- Assay: TF ChIP-seq
- Candidate files:
  - R1: `ENCFF851YSO`
  - R2: `ENCFF451LIG`
- URLs:
  - `https://www.encodeproject.org/files/ENCFF851YSO/@@download/ENCFF851YSO.fastq.gz`
  - `https://www.encodeproject.org/files/ENCFF451LIG/@@download/ENCFF451LIG.fastq.gz`
- Chromap smoke:
  - `--preset chip --TagAlign`
- Checks:
  - CLI output is nonempty.
  - Runner output is nonempty.
  - Sorted non-header rows are byte-identical between CLI and runner.

### Bulk ATAC-seq

Use the existing manifest row.

- Experiment: `ENCSR803FKU`
- Assay: ATAC-seq
- Candidate files:
  - R1: `ENCFF013DXQ`
  - R2: `ENCFF473SZG`
- URLs:
  - `https://www.encodeproject.org/files/ENCFF013DXQ/@@download/ENCFF013DXQ.fastq.gz`
  - `https://www.encodeproject.org/files/ENCFF473SZG/@@download/ENCFF473SZG.fastq.gz`
- Chromap smoke:
  - `--preset atac --BED`
- Checks:
  - CLI output is nonempty.
  - Runner output is nonempty.
  - Sorted non-header rows are byte-identical between CLI and runner.

### scATAC / snATAC

Add a new manifest case. This is the main extension beyond the current paired
FASTQ harness.

- Experiment: `ENCSR308ZGJ`
- Assay: snATAC-seq
- Biosample: K562 nuclear fraction
- Layout: 10x-style ATAC, with genomic read mates in R1/R3 and barcode/index in
  R2.

Start with one lane for a fast smoke:

| Lane | R1 genomic | R3 genomic mate | R2 barcode/index |
| --- | --- | --- | --- |
| L001 | `ENCFF007XKT` | `ENCFF973TAU` | `ENCFF421HKC` |

Optional full-lane extension:

| Lane | R1 genomic | R3 genomic mate | R2 barcode/index |
| --- | --- | --- | --- |
| L001 | `ENCFF007XKT` | `ENCFF973TAU` | `ENCFF421HKC` |
| L002 | `ENCFF797DEE` | `ENCFF687WWF` | `ENCFF208SDB` |
| L003 | `ENCFF538LWE` | `ENCFF754IYJ` | `ENCFF361EJT` |
| L004 | `ENCFF357YOD` | `ENCFF207XDD` | `ENCFF036OLH` |

Download URLs are the standard ENCODE file URLs:

```text
https://www.encodeproject.org/files/<ACCESSION>/@@download/<ACCESSION>.fastq.gz
```

Whitelist source:

```text
https://www.encodeproject.org/files/737K-arc-v1(ATAC)/@@download/737K-arc-v1(ATAC).txt.gz
```

Chromap smoke:

```bash
--preset atac --BED --read-format bc:0:-1 --barcode-whitelist <whitelist>
```

Inputs:

```bash
-1 <R1 fastq>
-2 <R3 fastq>
-b <R2 barcode fastq>
```

Checks:

- R1, R3, and barcode FASTQs have equal downsampled record counts.
- Read-name stems are aligned across all three files after removing `/1`, `/2`,
  `/3`, whitespace suffixes, and lane/read suffixes where present.
- CLI output is nonempty.
- Runner output is nonempty.
- Sorted non-header fragment rows are byte-identical between CLI and runner.
- Summary files, when emitted, report at least one total barcode category.

### Hi-C

Use the existing manifest row.

- Experiment: `ENCSR173DMV`
- Assay: Hi-C
- Candidate files:
  - R1: `ENCFF396ZOI`
  - R2: `ENCFF995XRN`
- URLs:
  - `https://www.encodeproject.org/files/ENCFF396ZOI/@@download/ENCFF396ZOI.fastq.gz`
  - `https://www.encodeproject.org/files/ENCFF995XRN/@@download/ENCFF995XRN.fastq.gz`
- Chromap smoke:
  - `--split-alignment --pairs --MAPQ-threshold 1`
  - or `--preset hic --pairs` once that produces equivalent runner coverage.
- Checks:
  - CLI `.pairs` output is nonempty after removing header/comment lines.
  - Runner `.pairs` output is nonempty.
  - Sorted non-header rows are byte-identical between CLI and runner.
  - Do not run `.cool`, `.hic`, Juicer, Arrowhead, or TAD callers in this suite.

## Manifest Design

Do not overload the existing paired-only TSV by appending columns without updating
the parser. Bash `read` will fold extra fields into the last variable and corrupt
the description.

Recommended implementation: introduce a v2 manifest, for example:

```text
tests/encode_cross_assay_manifest.tsv
```

Suggested columns:

```text
assay
case_id
experiment
layout
r1_accession
r2_accession
barcode_accession
whitelist_id
r1_url
r2_url
barcode_url
whitelist_url
default_read_count
chromap_args
description
r1_accession
r2_accession
barcode_accession
whitelist_id
r1_url
r2_url
barcode_url
whitelist_url
```

Rules:

- `layout=paired` for ChIP, bulk ATAC, and Hi-C.
- `layout=scatac_10x_atac` for R1/R3 genomic plus R2 barcode/index input.
- `barcode_accession`, `barcode_url`, `whitelist_id`, and `whitelist_url` are
  empty for paired-only assays.
- Keep accessions and URLs explicit so the suite remains reproducible even if
  ENCODE search result ordering changes.

The generated manifest should live under:

```text
plans/artifacts/encode_cross_assay_cache/manifest.generated.tsv
```

Generated manifest columns should include source metadata:

```text
assay
case_id
experiment
layout
read_count
r1_fastq
r2_fastq
barcode_fastq
whitelist
r1_md5
r2_md5
barcode_md5
whitelist_md5
chromap_args
source_status
source_file_formats
source_output_types
source_read_lengths
description
```

## Fixture Preparation

Create a new preparation script or carefully refactor the existing one:

```text
tests/prepare_encode_cross_assay_fixtures.sh
```

Behavior:

1. Read `tests/encode_cross_assay_manifest.tsv`.
2. Select rows by `ENCODE_ASSAYS=chip,atac,scatac,hic` or `all`.
3. Download only when `ENCODE_ALLOW_DOWNLOAD=1`.
4. Cache full downloads under:

   ```text
   plans/artifacts/encode_cross_assay_cache/downloads/
   ```

5. Downsample deterministically by first-N FASTQ records.
6. Write downsampled FASTQs under:

   ```text
   plans/artifacts/encode_cross_assay_cache/downsampled/<case_id>_<reads>/
   ```

7. Normalize gzip output with `gzip -n` so checksums are stable.
8. Verify paired or triplet read-name alignment.
9. Download/decompress the whitelist for scATAC if not already cached.
10. Write `manifest.generated.tsv`.

Default read counts:

- `10000` paired reads for ChIP, bulk ATAC, and Hi-C.
- `10000` read triplets for scATAC.
- Allow override with `ENCODE_DOWNSAMPLE_READS`.

## Smoke Runner

Create a new runner or extend the existing runner behind the v2 manifest:

```text
tests/run_encode_cross_assay_smoke.sh
```

Behavior:

1. Require:

   ```bash
   CHROMAP_GRCH38_REF=/path/to/genome.fa
   CHROMAP_GRCH38_INDEX=/path/to/genome.index
   ```

2. Build `chromap` and `chromap_lib_runner` unless `BUILD=0`.
3. Prepare fixtures unless `ENCODE_SKIP_PREPARE=1`.
4. Run cases serially by default.
5. Write logs and outputs under:

   ```text
   plans/artifacts/encode_cross_assay_smoke/<timestamp>/
   ```

6. Record:

   - command line,
   - git state including untracked files,
   - source manifest path,
   - generated manifest path,
   - reference/index paths,
   - thread count,
   - output row counts,
   - pass/fail status per case.

7. Compare CLI and runner output through deterministic canonicalization:

   ```bash
   LC_ALL=C grep -v '^#' output | sed '/^[[:space:]]*$/d' | LC_ALL=C sort
   ```

For BAM/CRAM extensions later, do not use byte compare. Use a semantic comparator
based on `samtools view`, selected fields, and stable sorting.

## Make Targets

Add explicit opt-in targets:

```make
prepare-encode-cross-assay-fixtures:
	./tests/prepare_encode_cross_assay_fixtures.sh

test-encode-cross-assay-smoke: chromap chromap_lib_runner
	./tests/run_encode_cross_assay_smoke.sh
```

Do not include these in default `test-smoke`. They are S1 real-data checks.

## Example Commands

Prepare all fixtures, allowing downloads:

```bash
CHROMAP_GRCH38_REF=/path/to/GRCh38/genome.fa \
CHROMAP_GRCH38_INDEX=/path/to/GRCh38/genome.index \
ENCODE_ALLOW_DOWNLOAD=1 \
ENCODE_ASSAYS=chip,atac,scatac,hic \
ENCODE_DOWNSAMPLE_READS=10000 \
make test-encode-cross-assay-smoke
```

Run only Hi-C from an existing cache:

```bash
CHROMAP_GRCH38_REF=/path/to/GRCh38/genome.fa \
CHROMAP_GRCH38_INDEX=/path/to/GRCh38/genome.index \
ENCODE_ALLOW_DOWNLOAD=0 \
ENCODE_SKIP_PREPARE=0 \
ENCODE_ASSAYS=hic \
make test-encode-cross-assay-smoke
```

Run a faster scATAC-only check:

```bash
CHROMAP_GRCH38_REF=/path/to/GRCh38/genome.fa \
CHROMAP_GRCH38_INDEX=/path/to/GRCh38/genome.index \
ENCODE_ALLOW_DOWNLOAD=1 \
ENCODE_ASSAYS=scatac \
ENCODE_DOWNSAMPLE_READS=1000 \
make test-encode-cross-assay-smoke
```

## Acceptance Criteria

Implementation is complete when:

- `tests/encode_cross_assay_manifest.tsv` exists and includes ChIP, bulk ATAC,
  scATAC, and Hi-C rows.
- Fixture preparation downloads only with `ENCODE_ALLOW_DOWNLOAD=1`.
- Full and downsampled FASTQs stay under `plans/artifacts/`.
- Pair/triplet read-name alignment is verified before Chromap runs.
- Generated manifests and run manifests together record source accessions, URLs,
  file checksums, read counts, reference/index paths, and command lines.
- `make test-encode-cross-assay-smoke` passes on a host with GRCh38 reference and
  index for:
  - `ENCODE_ASSAYS=chip`,
  - `ENCODE_ASSAYS=atac`,
  - `ENCODE_ASSAYS=scatac`,
  - `ENCODE_ASSAYS=hic`,
  - `ENCODE_ASSAYS=all`.
- CLI and `chromap_lib_runner` produce identical canonical text rows for each
  supported text-output case.
- `tests/README.md` documents the new S1 suite and environment variables.
- `docs/` or `mcp_server/README.md` is updated if the suite becomes part of a
  durable public testing workflow.

## Follow-Ups

- Add a compressed CBQ parity case for the scATAC ENCODE downsample after CBQ
  ATAC validation is hardened.
- Add a negative assertion that CBQ plus Hi-C/PAIRS is rejected until CBQ Hi-C is
  intentionally implemented.
- Add optional BAM/CRAM semantic comparison for bulk ATAC or scATAC dual-output
  cases once runner support and comparison tooling are stable.
- Consider a separate STAR-suite runbook for downstream multiome validation
  after Chromap-side S1 coverage is reliable.
