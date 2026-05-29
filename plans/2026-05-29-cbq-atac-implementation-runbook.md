# CBQ Paired-Read Implementation Runbook

Date: 2026-05-29

## Scope

This runbook tracks the first Chromap-suite implementation of native CBQ input
for paired-read Chromap workflows. The target branch is:

```bash
feature/cbq-atac-input
```

The implementation started as ATAC/scATAC-only and now covers the shared
paired-read path:

- paired-end read CBQ input,
- optional single barcode CBQ input for scATAC,
- ChIP TagAlign output,
- Hi-C pairs output,
- `chromap` CLI support,
- `chromap_lib_runner` / `libchromap` parity,
- no intermediate FASTQ materialization in the production path.

Do not expand this first pass to single-end CBQ or STAR-suite orchestration
defaults.

## User Interface

FASTQ remains the default input mode.

Native CBQ mode is selected explicitly:

```bash
chromap --preset atac \
  --input-format cbq \
  --read-pair-cbq lane1.reads.cbq,lane2.reads.cbq \
  --barcode-cbq lane1.barcode.cbq,lane2.barcode.cbq \
  --barcode-whitelist whitelist.txt \
  -r ref.fa \
  -x ref.index \
  -o fragments.bed \
  --BED
```

The same options are exposed through `chromap_lib_runner`:

```bash
chromap_lib_runner --preset atac \
  --input-format cbq \
  --read-pair-cbq reads.cbq \
  --barcode-cbq barcodes.cbq \
  --barcode-whitelist whitelist.txt \
  -r ref.fa \
  -x ref.index \
  -o fragments.bed \
  --BED
```

Validation rules:

- `--input-format fastq` is the default.
- `--input-format cbq` requires `--read-pair-cbq`.
- `--barcode-cbq` count must match `--read-pair-cbq` count.
- `--barcode-whitelist` with CBQ requires `--barcode-cbq`.
- FASTQ inputs (`-1`, `-2`, `-b`) cannot be mixed with CBQ inputs.
- CBQ currently rejects `--emit-Y-noY-fastq`.

## Implementation Shape

Core files:

- `src/cbq_reader.h`
- `src/cbq_reader.cc`
- `src/sequence_batch.h`
- `src/sequence_batch.cc`
- `src/mapping_parameters.h`
- `src/chromap.cc`
- `src/chromap.h`
- `src/chromap_driver.cc`
- `src/chromap_lib_runner.cc`
- `src/libchromap.cc`
- `src/mapping_writer.cc`

The production path is:

1. The CBQ reader opens one paired-read CBQ lane and, for scATAC, one barcode
   CBQ lane.
2. Each decoded CBQ record is assigned into the existing `SequenceBatch`
   storage through `AssignLoadedSequence()`.
3. The paired-end mapping loop consumes the filled `SequenceBatch` objects
   exactly like FASTQ-loaded batches.
4. Barcode correction, duplicate handling, BED/TagAlign/pairs/BAM writing,
   low-memory spill handling, summary metadata, and peak sidecars remain
   downstream of the input mode branch.

The CBQ reader decodes CBQ v1 block structure and zstd-compressed columns. It
uses `dlopen()` for `libzstd.so.1` / `libzstd.so`, so the top-level link adds
`-ldl`. This avoids adding a hard link dependency to the existing build while
still requiring zstd at runtime for compressed CBQ files.

## Known Difficulties

The implementation had more integration risk than a normal FASTQ reader swap:

- Chromap previously inferred paired-end and barcode behavior from FASTQ CLI
  fields (`-2`, `-b`). CBQ mode needed semantic helpers such as
  `HasPairedEndInput()`, `HasBarcodeInput()`, and `NumInputLanes()`.
- Barcode pre-sampling and barcode abundance were separate FASTQ-only reads.
  Both now have CBQ branches.
- CBQ is block-oriented. The reader must preserve EOF state across repeated
  fixed-size batch loads.
- BAM read-group auto naming previously assumed `read_file1_paths[0]`; CBQ now
  routes this through `ReadGroupSourcePath()`.

## Current Limitations

These are intentional for the first milestone:

- Paired-end CBQ only.
- No CBQ support for single-end mapping.
- No `--emit-Y-noY-fastq` support in CBQ mode.
- The first reader materializes ASCII sequences into `SequenceBatch`; packed
  CBQ-to-minimizer optimization is deferred.
- Test CBQ generation uses external `bqtools`; released Chromap-suite does not
  depend on STAR-suite or a STAR adapter.

## Build

```bash
make chromap chromap_lib_runner
```

Expected result: both binaries build successfully.

## Synthetic Smoke

Run:

```bash
make test-cbq-atac-smoke
```

The smoke:

1. Generates a small synthetic reference, paired ATAC FASTQs, barcode FASTQ, and
   whitelist.
2. Encodes paired reads and barcodes to CBQ with `bqtools`.
3. Runs FASTQ baseline through `chromap`.
4. Runs native CBQ through `chromap`.
5. Runs native CBQ through `chromap_lib_runner`.
6. Compares sorted fragment rows.

Set `BQTOOLS=/path/to/bqtools` when `bqtools` is not on `PATH`. Artifacts are
written under:

```text
plans/artifacts/cbq_atac_smoke/<timestamp>/
```

Known passing run on this branch:

```text
plans/artifacts/cbq_atac_smoke/20260529T181128Z/
```

## Existing Regression Smoke

Run:

```bash
make test-libchromap-core-smoke
```

This verifies that existing FASTQ CLI vs libchromap behavior still matches
after the input-mode changes.

Known passing run on this branch:

```text
plans/artifacts/chromap_core_smoke/20260529T181111Z/
```

## 100K Fixture Gate

Run:

```bash
make test-cbq-atac-100k
```

This compares barcoded paired-end ATAC fragments over the 100K PBMC fixture
(all four lanes by default; override `LANES="1"` for a faster single-lane run):

```text
/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac
```

- FASTQ baseline through `chromap`,
- CBQ mode through `chromap`,
- CBQ mode through `chromap_lib_runner`,
- canonical comparison: `LC_ALL=C` sorted fragment (BED) rows.

### Encoder ordering requirement

The barcoded ATAC CBQ path consumes a paired-read CBQ lane plus a separate
barcode CBQ lane and requires record `i` of the read lane to correspond to
record `i` of the barcode lane. chromap enforces this with a read/barcode
name-match guard.

Stock `bqtools encode` does **not** preserve input record order at multi-block
scale: the read lane and the barcode lane are independently reordered, so they
no longer align and chromap rejects them (the 4-record synthetic smoke fits one
block and survives; 100K does not). The 100K gate therefore uses the
order-preserving `cbq_ordered_encoder`, which emits records in input FASTQ
order. Point the script at it with `CBQ_ORDERED_ENCODER=/path/to/cbq_ordered_encoder`;
it defaults to:

```text
/mnt/pikachu/STAR-suite/core/legacy/source/cbq_ordered_encoder
```

This is an optional, test-only external dependency. The zstd decompression path
is covered separately by the synthetic smoke, which feeds chromap a compressed
`bqtools` CBQ. Released Chromap-suite does not depend on STAR-suite.

Artifacts (command lines, git state, fixture/index/reference paths, thread
count, output directory, wall/user/sys/max-RSS per chromap run, and the
canonicalization used) are written under:

```text
plans/artifacts/cbq_atac_100k/<timestamp>/MANIFEST.txt
```

Benchmarks run serially.

Known passing run on this branch (4 lanes, 8 threads, 320017 sorted fragment
rows, FASTQ == CLI CBQ == libchromap CBQ):

```text
plans/artifacts/cbq_atac_100k/20260529T202040Z/
```

## ENCODE Cross-Assay CBQ Gate

Run:

```bash
CHROMAP_GRCH38_REF=/path/to/genome.fa \
CHROMAP_GRCH38_INDEX=/path/to/genome.index \
make test-encode-cbq-cross-assay-smoke
```

This gate layers CBQ parity on the ENCODE cross-assay FASTQ fixtures:

- ChIP TagAlign,
- bulk ATAC BED,
- scATAC BED plus barcode summary,
- Hi-C pairs.

It uses `cbq_ordered_encoder` with zstd compression by default and compares
canonical FASTQ, CBQ CLI, and CBQ lib-runner rows.

## Merge Checklist

- `make chromap chromap_lib_runner` passes.
- `make test-cbq-atac-smoke` passes.
- `make test-libchromap-core-smoke` passes.
- `make test-cbq-atac-100k` passes (requires `cbq_ordered_encoder`; see above).
- `make test-encode-cbq-cross-assay-smoke` passes for `chip`, `atac`,
  `scatac`, and `hic` when the ENCODE fixtures and GRCh38 index are available.
- No production CBQ path writes temporary FASTQs.
- Unsupported CBQ combinations fail during argument validation.
- Branch commits exclude unrelated pre-existing dirty files.
