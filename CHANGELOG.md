# Changelog

All notable changes to this project will be documented in this file.

## [Unreleased]

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
