# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/), and this
project adheres to [Semantic Versioning](https://semver.org/). Per-release notes
are in [`docs/RELEASE_NOTES_vX.Y.Z.md`](docs/).

## [Unreleased]

### Added

- `--low-mem-finalize-threads N` (`MappingParameters::low_mem_finalize_threads`).
  With `--low-mem`, paired-end barcoded ATAC runs with sidecar-only, BED or
  TagAlign output merge and deduplicate their spill files as one task per
  reference, on N threads, and append the results in reference order.
  - The output is byte-identical to the serial merge: the sidecar, the text,
    the summary CSV, the in-memory peaks and summits, and the stderr counters.
  - 0 (the command-line default) uses `--num-threads`; 1 keeps the serial
    merge. The library default is 1, so embedding hosts opt in.
  - With permit hooks, every merge task holds a host permit; the serial merge
    holds one for its whole duration.
  - Dual BAM/CRAM output and bulk data keep the serial merge.
- `make test-lowmem-parallel-finalize`, added to the release gate (now 16
  targets).

### Changed

- `chromap` and `chromap_lib_runner` raise the soft open-file limit to the hard
  limit before mapping, because the low-memory merge keeps every spill file of
  a reference open. The library never changes process limits.

## [1.1.1] - 2026-09-29

Two fixes to the low-memory (`--low-mem`) merge of spilled mappings. See
[`docs/RELEASE_NOTES_v1.1.1.md`](docs/RELEASE_NOTES_v1.1.1.md).

### Fixed

- A low-memory spill file that cannot be opened now stops the run with an
  error. Previously the merge skipped the file and wrote output without its
  records. When the open-file limit is the cause, the message gives the limit
  and suggests `ulimit -n` or a larger `--low-mem-ram`.
- Bulk-level duplicate removal on barcoded data no longer reads past the
  whitelist table. A barcode outside the table (with
  `--output-mappings-not-in-whitelist`, or when no `--barcode-whitelist` is
  given) has abundance 0 when the kept barcode is chosen. Without a whitelist
  this previously could crash.

### Added

- `make test-lowmem-overflow-edge-cases`, a hermetic regression test for both
  fixes, added to the release gate (now 15 targets).

## [1.1.0] - 2026-09-28

Bounded ATAC barcode learning by default, sidecar-only ATAC output, and FASTQ
intake shared with STAR Suite. See
[`docs/RELEASE_NOTES_v1.1.0.md`](docs/RELEASE_NOTES_v1.1.0.md).

### Added

- `--atac-sidecar-only` (`MappingParameters::atac_sidecar_only`): paired-end
  ATAC runs write only the AEV1 sidecar and its `.chroms.tsv`, byte-identical
  to the sidecar of the BAM dual mode.
- BGZF FASTQ intake with STAR Suite's reader, mirrored under
  `src/star_input/`: each file of a paired-end lane is inflated by parallel
  workers and consumed on its own thread, and records are paired by ordinal
  and read-name stem. `--input-bgzf-mode auto|on|off` (default `auto`) and
  `--input-bgzf-reader-threads N` (`MappingParameters::input_bgzf_mode`,
  `input_bgzf_reader_threads`); the barcode abundance pass reads BGZF too.
- `tests/fastq_intake_harness` and `make test-fastq-intake-smoke`.
- A shared test gate (`make test-release`) and automated Ubuntu 22.04/24.04
  tarball and Debian packaging. Publication requires package runtime checks,
  Debian install/removal checks, a tested source-package rebuild, and a Docker
  smoke check. Manual workflow runs validate without publishing.

### Changed

- Ordinary runs learn barcode priors from `--barcode-sample-limit`
  (default 20,000,000) exact whitelist observations, ending at a batch
  boundary; `0` restores v1.0.1's complete-histogram learning. Barcoded
  outputs can therefore differ from v1.0.1.
- Plain-gzip FASTQ: read 1, read 2 and the barcode file are decompressed on
  their own threads in every batch from three threads up, not only in the
  first batch.
- The stderr log names the FASTQ reader of each lane and its loader time.

## [1.0.1] - 2026-08-23

Chromap Suite patch release embedding RapidMACS v1.0.1 and preserving exact
MACS3 v3.0.3 narrow-peak behavior at score-boundary and fold-enrichment edge
cases. See
[`docs/RELEASE_NOTES_v1.0.1.md`](docs/RELEASE_NOTES_v1.0.1.md).

### Fixed

- Embedded peak calling now uses MACS3 `callpeak`'s strict `score > cutoff`
  predicate in p-value and q-value modes, so a score exactly equal to the
  requested cutoff is excluded from the peak boundary.
- Embedded peak calling now applies MACS3's default fold-enrichment cutoff of
  1.0 after peak construction, excluding statistically significant depletion
  regions from narrowPeak output.
- Direct and materialized ATAC paths now preserve deterministic fragment
  reduction and exact peak-caller input equivalence.

### Changed

- Replaced the former `libMACS3` dependency naming with the public RapidMACS /
  `librapidmacs` project and pinned `third_party/rapidmacs` to RapidMACS
  v1.0.1 (`c22ff439cc66337e1b7fc9b87b547478eb5c0d3b`).
- ATAC mapping releases shared concurrency permits at task boundaries, avoiding
  retained-capacity stalls without changing output.

### Added

- Opt-in mergeable ATAC spill records, deterministic materialization, and
  fixed-width materialized blocks for bounded-memory downstream processing.
- Materialized-reference sidecars and direct index-loading support.
- A ChIP TagAlign workflow schema in the MCP/Launchpad registry.

## [1.0.0] - 2026-06-18

First production Chromap Suite release. `chromap --version` now reports the suite
version (`1.0.0`); the chromap engine version (`0.3.3-r519`) is available via
`chromap --upstream-version`. See [`docs/RELEASE_NOTES_v1.0.0.md`](docs/RELEASE_NOTES_v1.0.0.md).

### Fixed
- SAM writer: emit each record atomically and use length-safe `fwrite` to prevent interleaving and stray bytes.
- Low-memory spill: serialize flush section with `omp critical` and protect memory counter with `omp atomic` to avoid races and temp file name collisions.
- Temp mapping reload: correctly read `size_t` counts into `uint32_t` fields via local `size_t` then cast, preventing in-memory corruption during merges.
- New overflow system: eliminated remaining thread-safety issues in temp file handling by using thread-local writers with coordinated cleanup, preventing malformed SAM lines from unclosed/unflushed temp files.
- Serialization size calculation: fixed `SerializedSize()` methods to account for all variable-length fields, preventing malloc corruption during buffer operations.
- K-way merge: implemented proper k-way merge algorithm for overflow files to ensure correct sorted, deduplicated output with full behavior parity (PCR dedup, MAPQ filtering, Tn5 shift, metadata updates).

### Added
- New CLI flag: `--Tn5-shift-mode STR` to select the Tn5 cut-site offset convention applied when `--Tn5-shift` is active.
  - `classical` = `+4 / -5` (Buenrostro et al. 2013; Cell Ranger ARC / Cell Ranger ATAC default). This is the existing behavior and remains the default when only `--Tn5-shift` is passed.
  - `symmetric` = `+4 / -4` (ChromBPNet convention; symmetric around the 9-bp Tn5 footprint).
  - Passing `--Tn5-shift-mode` implies `--Tn5-shift`.
  - The active offsets are now echoed at startup, e.g. `Perform Tn5 shift (offsets: +4 / -5).`
  - Internal refactor: `Mapping::Tn5Shift()` now takes `(int forward_shift, int reverse_shift)` and the offsets are fields on `MappingParameters` (`Tn5_forward_shift`, `Tn5_reverse_shift`) instead of literals. SAM/BAM and pairs paths remain intentional no-ops (shifting those formats would require coordinated edits to `POS`, `MPOS/PNEXT`, `TLEN`, `CIGAR`, `NM`, `MD`).
- New CLI flag: `--temp-dir DIR` to specify directory for temporary files.
- **New overflow system is now the default**: Thread-safe temp file handling with k-way merge, no compile flags needed.
- `RotateThreadOverflowWriter()`: Per-flush file rotation to ensure one sorted run per overflow file (required for correct k-way merge).
- Enhanced mapping record serialization with precise size calculation and single-write operations.
- Validation scripts under `scripts/`:
  - `validate_sam_fix.sh` for SAM path sanity checks.
  - `validate_low_mem_fix.sh` for low-memory path summary validation.
  - `test_overflow_basic.sh` for new overflow system integration testing.
- Docker support with multi-stage builds and temp directory improvements.

### Changed
- **Default overflow system**: New overflow system with k-way merge is now the default build. Legacy temp file system available via `LEGACY_OVERFLOW=1` compile flag (single-threaded only).
- Temp file operations now use single atomic writes instead of multiple small writes with seeks.
- Thread cleanup is now explicit and coordinated rather than relying on automatic thread-local storage destruction.
- Overflow files are now processed in rid-ascending order to preserve coordinate-sorted output.

### Notes
- Existing CLI flags and behavior remain backward compatible.
- Legacy overflow path (`LEGACY_OVERFLOW=1`) is not thread-safe; use `-t 1` when testing.
- K-way merge implementation includes full deduplication and filtering logic for behavior parity with legacy path.
