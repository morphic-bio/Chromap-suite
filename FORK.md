# Fork Notes: Chromap

This document records detailed changes in this fork, the rationale behind them, and notes for future work.

## Why this fork

- Address rare but impactful output corruption in SAM mode under high concurrency.
- Stabilize and harden the low-memory spill/merge path used for large datasets.

## Major fixes and how they were implemented

### 1) SAM writer: prevent line tearing and unsafe string writes
- Problem: SAM records were emitted via multiple `fprintf()` calls using `%s` with `std::string::data()`. Under concurrency, writes could interleave and `%s` relied on a terminating NUL (not guaranteed pre-C++17), risking over-read.
- Changes:
  - `src/mapping_writer.h`: `AppendMappingOutput` now uses `fwrite(line.data(), 1, line.size(), ...)` (length-safe).
  - `src/mapping_writer.cc`: `MappingWriter<SAMMapping>::AppendMapping` builds the entire record into a single `std::string` and writes it once (atomic at the `FILE*` level).
- Impact: Eliminates interleaved/partial lines and spurious bytes.

### 2) Low-memory spill path: serialize flush and protect counters
- Problem: Output tasks updated `num_mappings_in_mem`, sorted/cleared shared containers, and appended to `temp_mapping_file_handles` concurrently.
- Changes in `src/chromap.h` (both single-end and paired-end paths):
  - Use `#pragma omp atomic` for `num_mappings_in_mem += added_mappings`.
  - Wrap the flush sequence (threshold check → sort → spill → bookkeeping) in `#pragma omp critical(output_flush)`.
  - Temp file naming remains based on vector size but is now inside the critical region, ensuring uniqueness.
- Impact: Prevents data races, container corruption, and temp file name collisions.

### 3) Temp mapping reload: fix size_t → uint32_t read mismatch
- Problem: Spill files store counts as `size_t`. Reload path read directly into a `uint32_t` field, writing 8 bytes into a 4-byte slot on LP64 systems, corrupting adjacent state and causing huge sizes (e.g., 141012992) and truncations during merge.
- Changes in `src/temp_mapping.h`:
  - Read into a local `size_t`, then cast-assign to the `uint32_t` field.
- Impact: Correct in-memory state, preventing downstream corruption.
- This still has problems so we rewrote this and simplified the serialization and spillover

### 4) New overflow system: thread-local writers with k-way merge (now default)
- Problem: The original temp file system had subtle race conditions during thread-local storage cleanup, causing some temp files to remain unclosed/unflushed when the main thread tried to read them, leading to remaining malformed SAM lines.
- Changes (now the **default** build; use `LEGACY_OVERFLOW=1` to revert):
  - `src/overflow_writer.h/cc`: Thread-safe overflow writer using length-prefixed binary format (rid + payload_size + serialized_mapping)
  - `src/overflow_reader.h/cc`: Reader for overflow files that yields (rid, payload) pairs
  - `src/mapping_writer.h/cc`: 
    - Each thread owns its own `OverflowWriter` instance during spill operations
    - `RotateThreadOverflowWriter()`: Closes current writer after each flush to ensure one sorted run per file
    - `ProcessAndOutputMappingsInLowMemoryFromOverflow()`: Full k-way merge implementation
  - `src/chromap.h`: 
    - Calls `RotateThreadOverflowWriter()` immediately after each `OutputTempMappingsToOverflow()` call
    - Added explicit thread cleanup phase where each thread closes its writer and contributes file paths to shared collection
  - All mapping record types: Added precise `SerializedSize()` methods and single-write `WriteToFile()` implementations
- Impact: Eliminates remaining thread-safety issues, ensures all temp files are properly closed before reading, removes complex seek operations for simpler and more robust I/O.

#### K-way merge implementation details

The new overflow system uses a k-way merge algorithm to correctly reconstruct sorted, deduplicated output from multiple overflow files:

1. **Per-flush file rotation**: After each memory threshold spill, `RotateThreadOverflowWriter()` closes the current writer and collects file paths. This ensures each overflow file contains exactly **one sorted run**, which is essential for correct k-way merge.

2. **Rid-group ordering**: Files are scanned to determine which reference sequence IDs (rids) they contain. Rids are then processed in **ascending order** to preserve coordinate-sorted output.

3. **Priority queue merge**: For each rid, a min-heap (priority queue) merges records from all files containing that rid. Records are compared by mapping position to maintain sort order.

4. **Full behavior parity**: The merge includes:
   - PCR duplicate removal (cell-level and bulk-level)
   - MAPQ filtering
   - Tn5 shift (if enabled)
   - Summary metadata updates (dup counts, low-MAPQ counts, mapped counts)

5. **File descriptor management**: Files are opened per-rid and closed after processing to limit FD usage.

6. **Cleanup**: All overflow files are deleted after successful merge.

**Correctness dependency**: If overflow files contained concatenations of multiple sorted runs (rather than one run per file), a simple file-level k-way merge would produce incorrect results. The per-flush rotation ensures this invariant holds.

## Tests and validation

- Scripts (run from repo root):
  - `./scripts/validate_sam_fix.sh`: sanity checks for SAM emission and basic execution.
  - `./scripts/validate_low_mem_fix.sh`: summarizes low-mem concurrency hardening.
  - `./scripts/test_overflow_basic.sh`: validates new overflow system integration.
- Recommended runtime checks:
  - SAM integrity: `awk 'BEGIN{FS="\t"} !/^@/ && NF<11{bad++} END{exit bad!=0}' output.sam`
  - SAM parsing: `samtools view -S output.sam > /dev/null`
- Stress testing: enable low-memory mode with a small threshold and high thread count.

## User-facing behavior and flags

- New flag: `--Tn5-shift-mode {classical|symmetric}` to pick the Tn5 cut-site offset convention. `classical` = `+4 / -5` (Buenrostro 2013; Cell Ranger ARC default; this is the legacy `--Tn5-shift` behavior). `symmetric` = `+4 / -4` (ChromBPNet convention). Implies `--Tn5-shift`. Active offsets are now echoed at startup.
- New flag: `--temp-dir DIR` to specify directory for temporary files (useful for Docker environments and custom temp locations).
- **New overflow system is now the default**: No compile flags needed for the improved overflow system with k-way merge.
- Compile flag: `LEGACY_OVERFLOW=1` to use the legacy temp file system (⚠️ **single-threaded only** - use `-t 1`).
- Existing presets and options work as before.
- Performance: neutral or slightly improved due to fewer stdio calls per record and simplified I/O patterns.

### 5) Configurable Tn5 shift offsets (classical vs symmetric)

- Motivation: the classical ATAC convention (`+4 / -5`, Buenrostro 2013, also used by Cell Ranger ARC/ATAC) is asymmetric around the 9-bp Tn5 duplication. Some downstream tools — notably ChromBPNet — require a symmetric `+4 / -4` cut-site shift. Previously this was always handled post-hoc (e.g. a `+1` adjustment on the minus-strand 5' end in the bigwig pipeline), which adds a silent, pipeline-specific step that is easy to forget or double-apply.
- Changes:
  - `src/mapping_parameters.h`: added `Tn5_forward_shift` (default `4`) and `Tn5_reverse_shift` (default `-5`, signed) alongside the existing `Tn5_shift` boolean.
  - `src/mapping.h`, `src/bed_mapping.h`, `src/paf_mapping.h`, `src/sam_mapping.h`, `src/pairs_mapping.h`: generalized all `Tn5Shift()` overrides to take `(int forward_shift, int reverse_shift)`. The paired-end math now reads `fragment_length -= (forward_shift - reverse_shift)` and `negative_alignment_length += reverse_shift`, which reduces to the previous `-9` / `-5` when `(fs, rs) = (4, -5)`. SAM and pairs keep their intentional no-op semantics (coordinate-shifting those formats would require coordinated edits to `POS`, `MPOS/PNEXT`, `TLEN`, `CIGAR`, `NM`, `MD`).
  - `src/mapping_processor.h`: `ApplyTn5ShiftOnMappings` takes the offsets and forwards them.
  - `src/chromap.h`, `src/mapping_writer.cc`, `src/mapping_writer.h`: the three call sites pass the offsets from `MappingParameters`.
  - `src/chromap_driver.cc`: new `--Tn5-shift-mode` option (`classical` or `symmetric`), implying `--Tn5-shift`; startup echo now reports the active offsets.
- Behavior:
  - `--Tn5-shift` alone → classical `+4 / -5` (legacy behavior, byte-identical to prior builds).
  - `--Tn5-shift-mode classical` → explicit `+4 / -5`.
  - `--Tn5-shift-mode symmetric` → `+4 / -4`. On paired-end BED output, every fragment end shifts by exactly `+1` vs classical, with identical start coordinates.
- Impact: no effect on SAM/BAM output (still intentionally unshifted); BED/BEDPE/PAF output now user-selectable. Eliminates the need for post-hoc `+1` corrections in downstream ChromBPNet pipelines while keeping Cell Ranger ARC parity as the default.

## Compatibility

- Output formats unchanged. SAM/BED/PAF/PAIRS semantics preserved.
- Binary compatible with previous indices and configs.

## Future work

- ~~Make `NEW_OVERFLOW=1` the default build configuration after validation.~~ ✅ **Done** - New overflow system with k-way merge is now the default.
- Consider a dedicated output thread to fully decouple compute from flush operations.
- Add CI checks for SAM field counts and `samtools view` parsing on sample outputs.
- Optional: add `setvbuf` to increase IO buffer size for throughput on very large runs.
- Docker integration improvements for better container-friendly temp directory handling.
