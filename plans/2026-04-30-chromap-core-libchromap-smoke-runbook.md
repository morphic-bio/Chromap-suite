# Chromap Core + libchromap Small-Set Regression Runbook

Date: 2026-04-30
Status: S0 implemented; S1 ENCODE fixture preparation and smoke harness implemented
Scope: Stage 1.5 hardening after MCP/Launchpad Stage 1

## Goal

Add a compact regression matrix proving that Chromap's main user-visible
surfaces still work after MCP/Launchpad/libchromap work. The matrix must stay
small enough for routine local use, but broad enough to catch accidental
breakage in:

- the Chromap CLI command surface,
- the callable `libchromap` boundary,
- assay presets and output formats,
- fork-specific BAM/Y/ATAC-fragment paths.

This is not a benchmark suite and not a downstream analysis pipeline. Hi-C
support stops at Chromap `.pairs` output, analogous to ATAC stopping at
fragments.

## Current Artifact Coverage

Chromap documents or tests name the expected output artifacts for essentially
all core areas, but it does not commit tiny ready-to-run FASTA/FASTQ fixtures
for the full matrix.

Summary count across the proposed 11 smoke areas:

- Named output artifact in README/docs/tests: 11 / 11
- Committed tiny input fixture artifact covering the area: 0 / 11
- Synthetic generator already present for part of the area: 3 / 11
- Existing external/paper fixture path named in tests: 4 / 11
- Needs new fixture sourcing/downsample decision: ChIP real-data, Hi-C real-data,
  optional bulk ATAC real-data

The existing synthetic generators are:

- `tests/data/generate_test_data.py`
- `tests/data/generate_test_data_y.py`

They are useful for hermetic mapping/BAM/Y smoke tests, but they are not
assay-realistic ChIP/ATAC/Hi-C fixtures.

## Artifact Sandbox

All generated fixtures and outputs for this runbook must stay under an ignored
artifact root. Default:

```text
plans/artifacts/
```

Harnesses should accept:

```bash
CHROMAP_ARTIFACT_ROOT="${CHROMAP_ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts}"
```

Recommended subdirectories:

```text
plans/artifacts/chromap_core_smoke/<timestamp>/
plans/artifacts/encode_fixture_cache/
plans/artifacts/libchromap_parity/<timestamp>/
```

Do not write generated FASTQs, Chromap indexes, BAM/CRAM/BED/pairs/fragments,
or downloaded ENCODE files into tracked source directories. Commit only
reproducible scripts, manifests, checksums, and documentation.

## Test Tiers

### Tier S0: Hermetic Synthetic Smoke

Runs with only this repository, `chromap`, and `chromap_lib_runner`.

Fixture strategy:

- Generate synthetic reference and paired FASTQs at runtime.
- Write generated inputs, index, outputs, and logs under
  `CHROMAP_ARTIFACT_ROOT/chromap_core_smoke/<timestamp>/`.
- Build a temporary Chromap index.
- Exercise output format paths and compare CLI vs libchromap where possible.
- Do not fetch external data.

This tier should be the default pre-commit/local smoke.

### Tier S1: ENCODE Downsample Smoke

Runs only when a host has:

- network or pre-downloaded fixture cache,
- a compatible reference FASTA and Chromap index, usually GRCh38,
- enough disk for small downsampled FASTQs and outputs.

Fixture strategy:

- Download selected ENCODE paired FASTQs.
- Downsample pair-preserving to a small read count.
- Store downloads/downsamples under
  `CHROMAP_ARTIFACT_ROOT/encode_fixture_cache/`.
- Store a manifest with ENCODE accessions, URLs, original metadata, read count,
  and downsample checksums. The manifest can be committed; the FASTQs must not.
- Run against a configured full-genome Chromap index.

This tier gives real assay coverage without moving downstream.

### Tier S2: Existing 100K / Paper Fixtures

Runs current larger ATAC/MACS3/paper-oriented tests. This remains explicit and
opt-in; it is not part of the small smoke matrix.

Examples:

- `make test-lowmem-bed-100k`
- `make test-peak-integration-matrix-100k`
- `make test-peak-memory-source-100k`

## Proposed Smoke Matrix

| ID | Area | Chromap names output artifact? | Existing fixture status | S0 synthetic | S1 ENCODE | libchromap parity |
| --- | --- | --- | --- | --- | --- | --- |
| C01 | Index build | Yes: `ref.index` / `genome.index` | No committed fixture | CLI only | optional full-genome index precondition | no current libchromap index API |
| C02 | Basic paired BED | Yes: `test.bed`, `aln.bed` | generator can cover | yes | optional | exact or sorted-text compare |
| C03 | ChIP preset BED/TagAlign | Yes: `aln.bed`, `aln.tagAlign` | no committed ChIP fixture | synthetic command-path only | yes, CTCF ChIP | exact/sorted text compare |
| C04 | ATAC preset BED + Tn5 | Yes: `aln.bed` | synthetic can exercise flags; 100K fixture named | yes | optional ATAC | sorted text compare |
| C05 | scATAC barcode + summary | Yes: whitelist, barcode output, `summary` | 100K fixture paths named; no tiny fixture | generate barcode FASTQ/whitelist | optional | text/summary invariants |
| C06 | Sorted BAM/index | Yes: `output.bam`, `.bam.bai` | external BAM writer paths named | yes, if `samtools` present for validation | optional | semantic BAM compare |
| C07 | Low-mem BED parity | Yes: `test-lowmem-bed-100k` output semantics | 100K fixture named | yes | optional | CLI vs runner plus low-mem parity |
| C08 | ATAC BAM + fragments | Yes: `possorted_bam.bam`, `fragments.tsv.gz`, `summary.tsv` | 100K fixture named | needs runner option support | optional ATAC | after runner supports `--atac-fragments` |
| C09 | MACS3 FRAG peaks | Yes: `*.narrowPeak`, `summits.bed` | 100K/paper fixtures named | no, unless synthetic has enough signal | optional | after runner supports MACS3 flags |
| C10 | Hi-C pairs | Yes: `aln.pairs` | no committed Hi-C fixture | synthetic pairs shape smoke | yes, Hi-C FASTQ | text shape + CLI/runner compare |
| C11 | Y/noY split | Yes: `output.noY.*`, `output.Y.*`, names list/FASTQ prefixes | generator exists | yes | not needed | after runner supports Y flags |

Interpretation:

- "Names output artifact" means docs/tests identify the artifact Chromap should
  write.
- "Existing fixture status" means whether the repository already has runnable
  tiny input files. Today, it generally does not.
- S0 is mandatory for small regression.
- S1 is optional real-data confidence.
- S2 remains for heavier paper gates.

## ENCODE Downsample Candidates

These are proposed S1 fixture sources. Do not download them into the repository.
Create a cache outside git, then commit only a manifest and the downsample
recipe if needed.

### ATAC-seq

Experiment: `ENCSR803FKU`

- Assay: ATAC-seq
- Biosample: human T-helper 17 cell
- Candidate pair:
  - R1: `ENCFF013DXQ`
  - R2: `ENCFF473SZG`
- Download URLs:
  - `https://www.encodeproject.org/files/ENCFF013DXQ/@@download/ENCFF013DXQ.fastq.gz`
  - `https://www.encodeproject.org/files/ENCFF473SZG/@@download/ENCFF473SZG.fastq.gz`

Use only if we want real bulk ATAC coverage beyond synthetic and existing PBMC
100K fixtures.

### ChIP-seq

Experiment: `ENCSR089DTY`

- Assay: TF ChIP-seq
- Target: CTCF
- Biosample: human omental fat pad
- Candidate pair:
  - R1: `ENCFF851YSO`
  - R2: `ENCFF451LIG`
- Download URLs:
  - `https://www.encodeproject.org/files/ENCFF851YSO/@@download/ENCFF851YSO.fastq.gz`
  - `https://www.encodeproject.org/files/ENCFF451LIG/@@download/ENCFF451LIG.fastq.gz`

This is the main gap for real assay coverage because Chromap documents ChIP
usage but the repo does not ship a ChIP fixture.

### Hi-C

Experiment: `ENCSR173DMV`

- Assay: Hi-C
- Biosample: HCT116
- Candidate pair:
  - R1: `ENCFF396ZOI`
  - R2: `ENCFF995XRN`
- Download URLs:
  - `https://www.encodeproject.org/files/ENCFF396ZOI/@@download/ENCFF396ZOI.fastq.gz`
  - `https://www.encodeproject.org/files/ENCFF995XRN/@@download/ENCFF995XRN.fastq.gz`

Hi-C smoke stops at `.pairs`. Matrix generation, `.cool`/`.hic`, Juicer,
Arrowhead, and TAD calls are explicitly downstream and out of scope.

## Downsample Rules

For paired FASTQs, downsample pair-preserving. Prefer first-N records for a
deterministic smoke fixture:

```bash
N_READS=10000
zcat R1.fastq.gz | head -n $((4 * N_READS)) | gzip -nc > R1.${N_READS}.fastq.gz
zcat R2.fastq.gz | head -n $((4 * N_READS)) | gzip -nc > R2.${N_READS}.fastq.gz
```

If using random downsampling, use the same deterministic seed and verify read
names remain paired:

```bash
seqtk sample -s100 R1.fastq.gz "${N_READS}" | gzip -nc > R1.${N_READS}.fastq.gz
seqtk sample -s100 R2.fastq.gz "${N_READS}" | gzip -nc > R2.${N_READS}.fastq.gz
```

Every cached downsample must get a manifest:

```yaml
source: ENCODE
experiment: ENCSR...
r1_accession: ENCFF...
r2_accession: ENCFF...
downloaded_at: YYYY-MM-DD
downsample_method: first_n_fastq_records
read_count: 10000
r1_md5: ...
r2_md5: ...
reference: GRCh38
chromap_index: /path/to/index
```

## libchromap Harness Plan

The repository already has a callable boundary:

- `chromap::RunMapping(const MappingParameters&)`
- `chromap::RunAtacMapping(const MappingParameters&)`

It also has `chromap_lib_runner`, which is a binary harness that calls
`libchromap` while accepting Chromap-like CLI flags. Use that as the first
regression harness instead of inventing a second one.

### Harness Target

Add a future script:

```text
tests/run_libchromap_core_smoke.sh
```

Responsibilities:

1. Build `chromap` and `chromap_lib_runner`.
2. Generate or locate the small fixture.
3. Build a temporary index with the CLI.
4. For each matrix case, run:
   - `chromap ... -o cli.out`
   - `chromap_lib_runner ... -o lib.out`
5. Compare outputs:
   - BED/TagAlign/pairs: exact compare when deterministic, otherwise sorted
     non-header rows plus header sanity.
   - SAM: normalize headers where needed, compare alignment rows.
   - BAM/CRAM: use `samtools view` and optional index validation.
   - summaries: compare required columns and count conservation, not necessarily
     byte-for-byte if ordering differs.
6. Emit one TSV summary with command, status, output counts, md5/semantic hash,
   and failure reason.

### Runner Gaps To Close

`chromap_lib_runner` currently covers many core mapping flags and output
formats, but it should be extended before the full matrix is enforced:

- `--atac-fragments`
- `--call-macs3-frag-peaks`
- `--macs3-frag-*` peak flags
- `--emit-noY-bam`, `--emit-Y-bam`
- `--emit-Y-read-names`
- `--emit-Y-noY-fastq`
- explicit Y/noY output path flags
- optional binary sidecar flags if those remain in scope

Index build has no libchromap API today, so C01 remains CLI-only unless we add
a library index entrypoint later.

## Suggested Implementation Order

1. Documentation-only pass:
   - add this runbook,
   - add recipe reference docs listing purpose, inputs, outputs, and smoke ID.
2. S0 harness:
   - generate synthetic fixture,
   - add CLI smoke for C01, C02, C04, C06, C07, C10, C11,
   - add libchromap parity for cases the runner already supports.
3. Extend `chromap_lib_runner` for missing fork-specific flags.
4. Add C05 and C08 libchromap parity once barcode/fragments are covered.
5. Add optional S1 ENCODE fixture fetch/downsample script and manifest.
6. Keep S2 paper/100K gates separate.

## Acceptance Gates

Small default gate:

```bash
make chromap chromap_lib_runner
tests/run_libchromap_core_smoke.sh
make test-libchromap-core-smoke
python3 -m pytest mcp_server/tests -q
```

Optional real-data gate:

```bash
CHROMAP_GRCH38_REF=/path/to/genome.fa \
CHROMAP_GRCH38_INDEX=/path/to/genome.index \
ENCODE_FIXTURE_CACHE=/path/to/cache \
tests/run_encode_downsample_smoke.sh
```

To create the ignored ENCODE downsample cache explicitly:

```bash
ENCODE_ALLOW_DOWNLOAD=1 \
ENCODE_ASSAYS=chip,atac,hic \
ENCODE_DOWNSAMPLE_READS=10000 \
tests/prepare_encode_downsample_fixtures.sh
```

The committed source manifest is `tests/encode_downsample_manifest.tsv`. The
generated manifest and FASTQs stay under
`plans/artifacts/encode_fixture_cache/` by default.

No test in this plan should call downstream Hi-C TAD tools. The Hi-C handoff
artifact is `.pairs`.

## Validation Record

2026-04-30 S1 implementation check:

- Materialized ignored ENCODE downsample cache at
  `plans/artifacts/encode_fixture_cache/`.
- Read count: 1,000 paired reads per assay.
- Reference:
  `/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa`
- Chromap index:
  `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index`
- Smoke output:
  `plans/artifacts/encode_downsample_smoke/20260430T151521Z/summary.tsv`

Observed real-data smoke result:

| Case | Assay | Status | Non-header rows | Args |
| --- | --- | --- | ---: | --- |
| C03_chip_tagalign | ChIP | PASS | 1516 | `--preset chip --TagAlign` |
| C04_atac_bed | ATAC | PASS | 676 | `--preset atac --BED` |
| C10_hic_pairs | Hi-C | PASS | 763 | `--split-alignment --pairs --MAPQ-threshold 1` |
