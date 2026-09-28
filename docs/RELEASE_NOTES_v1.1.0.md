# Chromap Suite v1.1.0 Release Notes

Date: draft of 2026-09-28; set when the release is tagged.

Chromap Suite v1.1.0 makes three changes to the ATAC production path. Barcode
learning is bounded again by default. A run can write the AEV1 fragment
sidecar alone, without BAM or fragment text. FASTQ intake reads BGZF with the
reader STAR Suite uses, and reads ordinary gzip files on separate threads in
every batch.

The new FASTQ intake changes no output. The bounded barcode learning does: on
barcoded runs it can change barcode-correction priors, and therefore
corrected barcodes, relative to v1.0.1, which learned from the complete
histogram. `--barcode-sample-limit 0` restores v1.0.1's behavior.

`chromap --version` reports `1.1.0`. The chromap engine version reported by
`chromap --upstream-version` is unchanged.

## Bounded ATAC barcode learning

- Ordinary runs again learn barcode priors from a bounded sample:
  `--barcode-sample-limit N` (default `20000000`) exact whitelist
  observations, stopping at the end of the current input batch, as historical
  Chromap did. `0` scans all barcode inputs. The limit never reduces the reads
  mapped.
- Mergeable-spill workers still collect their complete local histogram while
  mapping. An ordinary single-process run used as an exact parity control for
  a scatter/gather run needs `--deterministic-mapping --barcode-sample-limit 0`.
- Library field: `MappingParameters::barcode_sample_limit`.
- In the full CATATAC run through STAR Suite, the bounded default learned
  from 20,044,513 exact observations in 23.5 s. The complete histogram
  (245,336,345 observations) had taken 254.4 s.

## Sidecar-only ATAC output

- `--atac-sidecar-only` with `--atac-fragment-binary-output PATH` writes only
  the AEV1 sidecar and `PATH.chroms.tsv` (plus `--summary` when requested):
  no `-o` output, no BAM/CRAM and no fragment text rows.
- For the same input and options, the sidecar, its `.chroms.tsv` and the
  summary are byte-identical to those of `--BAM --atac-fragments
  --atac-fragment-binary-output`.
- Invalid combinations (BAM/CRAM or `-o` output, fragment text, sorting or
  indexing, Y/noY BAM streams, mergeable spills, file-source MACS3 peaks,
  single-end input) are rejected by the CLI, `chromap_lib_runner` and
  `RunMapping` alike.
- Library field: `MappingParameters::atac_sidecar_only`.

## FASTQ intake

- **BGZF.** A paired-end lane whose files are all BGZF FASTQ is read with STAR
  Suite's BGZF reader, mirrored under `src/star_input/` (`BgzfBlockReader`,
  `BgzfRangeReader`; see `src/star_input/MIRROR.md` for the source commit,
  checksums and edits). Each file is inflated by its own pool of workers and
  consumed in order on its own thread. Records are paired as STAR Suite pairs
  synchronized BGZF streams: the same ordinal and read-name stem in every
  file, and no file ending before the others. The barcode abundance pass reads
  a BGZF barcode file the same way.
- **Plain gzip.** Read 1, read 2 and the barcode file are decompressed on
  their own threads for every batch from three threads up. Previously only
  the first batch used threads; later batches used OpenMP tasks inside the
  mapping region, and only from twelve threads up. FIFO inputs fed by a
  producer that interleaves records across the files keep working.
- **Options.** `--input-bgzf-mode auto|on|off` (default `auto`) and
  `--input-bgzf-reader-threads N` follow STAR Suite's `--readFilesBgzfMode`
  and `--bgzfReaderThreads`. `auto` uses the BGZF reader when every file of
  the lane is a regular BGZF FASTQ file whose first records pair, and zlib
  otherwise. `on` requires the BGZF reader. `off` always uses zlib.
  `N = 0` derives the inflate workers from `--num-threads`. Library hosts set
  `MappingParameters::input_bgzf_mode` and `input_bgzf_reader_threads`.
- **Limits.** The BGZF reader keeps STAR Suite's fixed capacities: four-line
  FASTQ records, header lines up to 512 characters and sequences up to 650
  bases. Single-end mapping reads with zlib in every mode.

Throughput, from standalone runs on 50,000,000 DOGMA-plex lane-1 ATAC read
pairs at 32 threads: pending (the timed campaign is running; see
docs/handoffs/HANDOFF_STAR_MIRRORED_FASTQ_READER_20260928.md).

## Validation

- Records: Chromap receives the same (name, comment, sequence, quality)
  records in the same order through the old kseq path and the new paths. The
  checks covered DOGMA-plex lane-1 subsets of 4,000,000 and 50,000,000 read
  pairs, as plain gzip and as `bgzip` copies. They also covered the four lanes
  of the PBMC 3k 100k fixture (gzip and `bgzip` copies) and the four
  full-depth PBMC 3k ATAC lanes, 82,781,574 read pairs published as BGZF.
- Standalone output, old (`09b3164`, the commit before the reader change) vs
  new, 16 threads, on the PBMC 3k
  100k fixture, the DOGMA-plex 4M subset and full-depth PBMC 3k, in gzip and
  BGZF: fragments and AEV1 sidecars (with `.chroms.tsv`) are byte-identical
  in every mode (`auto` on gzip and on BGZF, and `off`). `--summary` differs
  only in the candidate-cache statistics (`cachehit`, `fric`, `estfrip`,
  `numcacheslots`), which differ in the same way between two runs of the old
  binary on the same reads. With `--deterministic-mapping` the summaries are
  byte-identical too.
- Embedded: STAR Suite `6e83853`, built unchanged against this branch, ran
  DOGMA-plex lane 1 (2,000,000 reads per arm, sidecar output). Its
  `atac_fragments.bin`, `.chroms.tsv` and every `atac/` output (peaks,
  summits, peak MEX, metrics) are byte-identical to those of the current
  no-BAM composition.
- `make test-fastq-intake-smoke` (new, 27 cases) and the existing smoke
  tiers (`test-smoke`, input-format, libchromap core, sidecar-only, CBQ ATAC,
  CBQ modality matrix, CBQ range reader, runtime spill schema) pass.

## Versioning and provenance

- `chromap --version` -> `1.1.0` (suite); `chromap --upstream-version`
  unchanged.
- `third_party/rapidmacs` is unchanged at `34df448` (RapidMACS v1.0.1 plus
  one documentation commit).
- `src/star_input/` is a copy of STAR Suite `6e83853`
  (`core/legacy/source/input/`, identical at STAR Suite `v1.9.5.a`) with
  namespace, include-path and include-guard edits only, under STAR Suite's
  MIT licences (`src/star_input/LICENSE`).
