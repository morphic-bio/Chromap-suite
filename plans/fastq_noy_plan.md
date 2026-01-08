# Plan: Y Read Names + Y/noY FASTQ Emission

## Reference Behavior (STAR-Flex)

Sources reviewed:
- ` /mnt/pikachu/STAR-Flex/docs/Y_CHROMOSOME_BAM_SPLIT.md `
- ` /mnt/pikachu/STAR-Flex/docs/flex_parameters.md `
- ` /mnt/pikachu/STAR-Flex/README_flex.md `
- ` /mnt/pikachu/STAR-Flex/.specstory/history/2026-01-07_20-11Z-emityreadnames-feature-testing.md `

Key behaviors to mirror:
- CLI flags:
  - `--emitYReadNames yes|no` + optional `--YReadNamesOutput <path>`
  - `--emitYNoYFastq yes|no`
  - `--emitYNoYFastqCompression gz|none`
  - `--YFastqOutputPrefix <prefix>` / `--noYFastqOutputPrefix <prefix>`
- Y read names:
  - One name per line, normalized (strip leading `@`, stop at whitespace, strip trailing `/1` or `/2`).
  - Names-only mode supported (can run without Y/noY BAM output).
- Y/noY FASTQ emission:
  - Directly emitted during alignment (no second-pass tool required).
- Output naming: insert `_Y` / `_noY` before the last `_<R#>` token actually present
  (e.g., `_R1`, `_R2`, `_R3`, `_R10`).
- If no `_R[0-9]+` token is found, fall back to `Y_reads.mateN.fastq(.gz)` and `noY_reads.mateN.fastq(.gz)`.
  - Optional compression: `gz` (default) or `none`.
  - Preserves header tokens/comments; outputs FASTA if no qualities.
  - Bulk paired-end: if either mate has a Y alignment, both mates route to Y.
  - Single-read mode: route per read.
  - Multiple input files: warn; names derived from first file for each mate.
- FASTQ-only mode in STAR-Flex: `--emitYNoYFastq yes` + `--outSAMtype None`.

## Chromap Current State

- Already supports Y/noY BAM streams via `--emit-noY-bam` / `--emit-Y-bam`.
- Y-hit detection is read-ID based and already pair-aware in `src/mapping_generator.cc`.
- Y contig detection is intentionally strict (exact `Y` after stripping `chr` prefix).
- No existing Y-read-name list output.
- No existing Y/noY FASTQ emission.

## Gaps / Decisions to Resolve

- **Y contig detection parity**: STAR-Flex matches `chrY`, `Y`, `chrY_*`, `Y_*`; Chromap currently matches only exact `Y`.
- **CLI naming**: use Chromap-style flags only (no STAR-Flex aliases).
- **Output naming defaults**:
  - Should Y read names default to `<out>Aligned.out_Y.names.txt` (STAR-Flex) or `<output>.Y.names.txt` (Chromap style)?
  - For FASTQ, should outputs live under the `-o` output directory or alongside input FASTQs?
- **Ordering guarantees**:
  - No ordering guarantees for FASTQ or name-list output (avoid expensive coordination).

## Proposed Implementation Plan

### Phase 1: Parameters + CLI Plumbing

- Add new fields to `src/mapping_parameters.h`:
  - `bool emit_y_read_names`
  - `std::string y_read_names_output_path`
  - `bool emit_y_noy_fastq`
  - `std::string y_noy_fastq_compression` (default `gz`)
  - `std::string y_fastq_output_prefix`
  - `std::string noy_fastq_output_prefix`
  - `std::array<std::string, 2> y_fastq_output_paths` / `noy_fastq_output_paths` (derived)
- Extend `src/chromap_driver.cc` options:
  - `--emit-Y-read-names`
  - `--Y-read-names-output`
  - `--emit-Y-noY-fastq`
  - `--emit-Y-noY-fastq-compression`
  - `--Y-fastq-output-prefix` / `--noY-fastq-output-prefix`
- Validate compression values (`gz` / `none`).
- Validate stdout behavior: if primary output is stdout, require explicit Y-read-names output and FASTQ prefixes (avoid ambiguous defaults).

### Phase 2: Output Naming Helpers

Add helpers (likely in `src/chromap_driver.cc` or new `src/fastq_namer.h`):
- Parse input filenames and insert `_Y` / `_noY` before last `_R[0-9]+` token.
- Handle `.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`, `.gz`.
- `--Y-fastq-output-prefix` / `--noY-fastq-output-prefix` override:
  - Append `mate1` / `mate2` + extension.
- Warn when multiple input files exist; derive names from first file for each mate (STAR-Flex behavior).

### Phase 3: Y-Hit Tracking Expansion

Currently Y-hit collection only runs for Y/noY BAM streams. Expand gating to:
- `emit_y_read_names` OR `emit_y_noy_fastq` OR existing Y/noY BAM flags.
- Update `src/chromap.h` to build Y contig mask and enable tracking when any of these flags are set.

### Phase 4: Y Read Names Emission

Implement a normalized-name writer:
- Add a small helper to normalize names (strip `@`, stop at whitespace, strip `/1` or `/2`).
- Write one line per read with Y alignment.
- Output logic:
  - Per-thread chunk files merged at end (no ordering guarantees).
- Ensure no duplicates (read-level, not alignment-level).

### Phase 5: Y/noY FASTQ Emission

Implement `FastqSplitWriter` (new helper class):
- Supports `gz` or `none`.
- Writes FASTQ if qualities exist, FASTA otherwise (prefix `@` vs `>`).
- Preserve `name + comment` from kseq (add `SequenceBatch::GetSequenceCommentAt`).
- Route logic:
  - Single-end: per read.
  - Paired-end: if pair has Y hit, route both mates to Y.
- Integrate into mapping loop:
  - Create a per-batch `std::vector<uint8_t> y_hit_flags` indexed by `read_index`.
  - After `#pragma omp taskwait` for mapping tasks, write FASTQ and Y-read names for the batch.
  - Keep writers open across all input files; close at end.

### Phase 6: Tests

Add coverage in `tests/`:
- Unit test: name normalization edge cases (`@`, whitespace, `/1`/`/2`).
- Unit test: output naming derivation for `.fastq.gz`, `.fq.gz`, `_R1/_R2` cases.
- Integration test:
  - Small synthetic reference with `chrY` and `chr1`.
  - Reads with some Y alignments, some non-Y.
  - Validate:
    - Y + noY FASTQ counts == total reads.
    - Y read names list contains expected normalized names.
    - Paired-end routing: if either mate hits Y, both mates are in Y FASTQ.

### Phase 7: Documentation

Update:
- `README.md`: new flags, usage examples, output naming rules.
- `chromap.1`: add CLI options + descriptions.
- Any relevant docs in `docs/` (if needed).

## Notes / Open Questions

- Should Chromap broaden Y contig detection to include `chrY_*` / `Y_*` like STAR-Flex?
- What should the default Y-read-names output path be (Chromap-style vs STAR-Flex-style)?
- For FASTQ, should outputs live under the `-o` output directory or alongside input FASTQs?
- Should FASTQ naming include read-file index when multiple input files are provided (avoid collisions)?
