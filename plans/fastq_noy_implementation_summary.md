# Implementation Summary: Y Read Names + Y/noY FASTQ Emission

## Overview

Successfully implemented Y chromosome read filtering features that emit:
1. **Y read names list**: Normalized list of read names with Y-chromosome alignments
2. **Y/noY FASTQ files**: Split FASTQ outputs based on Y-chromosome alignment status

These features work independently or together, and can be used with or without the existing Y/noY BAM stream functionality (`--emit-Y-bam` / `--emit-noY-bam`).

---

## Files Modified/Created

### New Files Created

1. **`src/y_read_names_writer.h`**
   - Helper class for writing normalized Y read names
   - Handles name normalization (strip `@`, stop at whitespace, strip `/1`/`/2`)
   - Deduplicates at read level
   - Thread-safe writing

2. **`src/fastq_split_writer.h`**
   - Helper class for writing Y/noY FASTQ files
   - Supports `gz` or `none` compression
   - Writes FASTQ (with qualities) or FASTA (without qualities)
   - Preserves read name + comment from kseq
   - Routes reads to appropriate streams based on Y-hit status

### Files Modified

1. **`src/mapping_parameters.h`**
   - Added fields to `MappingParameters` struct:
     - `bool emit_y_read_names = false`
     - `std::string y_read_names_output_path`
     - `bool emit_y_noy_fastq = false`
     - `std::string y_noy_fastq_compression = "gz"`  // "gz" or "none"
     - `std::string y_fastq_output_prefix`
     - `std::string noy_fastq_output_prefix`
     - `std::vector<std::vector<std::string>> y_fastq_output_paths_per_file`  // [file][mate]
     - `std::vector<std::vector<std::string>> noy_fastq_output_paths_per_file`  // [file][mate]

2. **`src/chromap_driver.cc`**
   - Added CLI options:
     - `--emit-Y-read-names`: Enable Y read names emission
     - `--Y-read-names-output`: Explicit path for Y read names (default: `<output>.Y.names.txt`)
     - `--emit-Y-noY-fastq`: Enable Y/noY FASTQ emission
     - `--emit-Y-noY-fastq-compression`: Compression type (`gz` or `none`, default: `gz`)
     - `--Y-fastq-output-prefix`: Prefix for Y FASTQ files (overrides auto-naming)
     - `--noY-fastq-output-prefix`: Prefix for noY FASTQ files (overrides auto-naming)
   - Added helper functions:
     - `GetDirectoryFromPath()`: Extract directory from file path
   - `DeriveFastqOutputPath()`: Derive FASTQ output paths from input filenames
      - Inserts `_Y` or `_noY` before last `_R[0-9]+` token (preserves suffixes like `_001`)
      - Preserves base read extension (`.fastq`/`.fq`/`.fasta`/`.fa`/`.fna`) and applies compression (`.gz` when requested)
      - Supports 1-based file index for multiple input files (e.g., `.f1`, `.f2`)
   - Added validation:
     - Compression must be `gz` or `none`
     - Requires explicit paths/prefixes when primary output is stdout
   - Added path derivation logic for FASTQ outputs

3. **`src/sequence_batch.h`**
   - Added new inline methods:
     - `GetSequenceCommentAt()`: Access comment field from kseq_t
     - `GetSequenceCommentLengthAt()`: Get comment length
     - `GetSequenceQualLengthAt()`: Get quality length

4. **`src/chromap.h`**
   - Added includes:
     - `#include "y_read_names_writer.h"`
     - `#include "fastq_split_writer.h"`
   - Expanded Y-hit tracking gates:
     - Updated all checks from `emit_noY_stream || emit_Y_stream` to include `emit_y_read_names || emit_y_noy_fastq`
     - Affects: Y contig mask building, thread-local tracking setup, Y-hit merging
   - Added writer initialization (in both `MapSingleEndReads` and `MapPairedEndReads`):
     - `YReadNamesWriter` creation if `emit_y_read_names` is enabled
     - `FastqSplitWriter` creation if `emit_y_noy_fastq` is enabled
   - Added batch-level output (after `#pragma omp taskwait`):
     - Y read names writing: iterates through batch, writes names for Y-hit reads
     - FASTQ writing: iterates through batch, routes reads to Y/noY streams
     - Both happen before batch swap to ensure read data is still available

---

## Feature Details

### Y Read Names Emission

**Purpose**: Output a normalized list of read names that have Y-chromosome alignments.

**Behavior**:
- One normalized name per line
- Normalization rules:
  - Strip leading `@` if present
  - Stop at first whitespace character
  - Strip trailing `/1` or `/2` if present
- Read-level deduplication (if a read has multiple alignments, name appears once)
- Works independently of Y/noY BAM streams
- Can run without any BAM output (names-only mode)

**Output Format**:
- Plain text file, one name per line
- Default path: `<output>.Y.names.txt` (Chromap-style)
- Example:
  ```
  read123
  read456
  read789
  ```

**CLI Usage**:
```bash
chromap --emit-Y-read-names -o output.sam ...
# Creates output.Y.names.txt

chromap --emit-Y-read-names --Y-read-names-output custom_names.txt -o output.sam ...
# Creates custom_names.txt
```

### Y/noY FASTQ Emission

**Purpose**: Split input FASTQ files into Y and noY streams based on alignment status.

**Behavior**:
- Directly emitted during alignment (no second-pass tool required)
- Output naming:
  - Inserts `_Y` or `_noY` before last `_R[0-9]+` token in input filename
  - Example: `sample_R1.fastq.gz` → `sample_Y_R1.fastq.gz` / `sample_noY_R1.fastq.gz`
  - Fallback: `Y_reads.mateN.<base-ext>(.gz)` / `noY_reads.mateN.<base-ext>(.gz)` if no `_R[0-9]+` token found
  - Base extension is preserved from input (`.fastq`/`.fq`/`.fasta`/`.fa`/`.fna`)
- Compression: `gz` (default) or `none`
- Format:
  - FASTQ if qualities exist (uses `@` prefix)
  - FASTA if no qualities (uses `>` prefix)
- Preserves header tokens/comments from original FASTQ
- Paired-end routing:
  - If either mate has a Y alignment, both mates route to Y FASTQ
  - Single-read mode: route per read
- Multiple input files:
  - Warns user
  - Names derived per input file (base name comes from each file)
  - Includes 1-based file index (e.g., `.f1`, `.f2`) to avoid collisions
- Output location: Under the `-o` output directory (same dir as primary output file)

**CLI Usage**:
```bash
# Basic usage
chromap --emit-Y-noY-fastq -1 read1.fastq.gz -2 read2.fastq.gz -o output.sam ...
# Creates: output_dir/read1_Y_R1.fastq.gz, read1_noY_R1.fastq.gz, etc.

# Custom compression
chromap --emit-Y-noY-fastq --emit-Y-noY-fastq-compression none ...

# Explicit prefixes
chromap --emit-Y-noY-fastq \
  --Y-fastq-output-prefix /path/to/Y_reads \
  --noY-fastq-output-prefix /path/to/noY_reads ...
```

---

## Design Decisions

### 1. Y Contig Detection
**Decision**: Keep strict exact-match behavior (no `chrY_*` / `Y_*` matching)
**Rationale**: Avoids changing existing Y/noY BAM semantics. If STAR-Flex parity is needed, can add opt-in mode later.

### 2. Default Y-read-names Path
**Decision**: Chromap-style `<output>.Y.names.txt` (derived from `-o`)
**Rationale**: Consistent with Chromap naming conventions. Requires explicit path if `-o` is stdout.

### 3. FASTQ Output Location
**Decision**: Under the `-o` output directory (same dir as primary output file)
**Rationale**: Keeps all outputs together. Requires explicit prefixes if `-o` is stdout.

### 4. Multiple Input Files
**Decision**: Include 1-based file index in derived FASTQ names (e.g., `sample_Y_R2.f1.fastq.gz`)
**Rationale**: Avoids collisions when multiple input files are provided.

### 5. Ordering Guarantees
**Decision**: No ordering guarantees for FASTQ or name-list output
**Rationale**: Avoids expensive coordination. Output is functional but order may vary.

### 6. Integration Point
**Decision**: Write FASTQ/Y-read-names after `#pragma omp taskwait` for each batch
**Rationale**: 
- Ensures all mapping tasks for the batch are complete
- Read data is still available before batch swap
- Maintains performance (no blocking on I/O during mapping)

### 7. Thread Safety
**Decision**: 
- Y-read-names: Single writer with read-level deduplication
- FASTQ: Single writer per output stream
- Y-hit tracking: Thread-local pointers to per-thread ID buffers (no shared mutation)
**Rationale**: Simple and safe. Batch-level writing avoids complex thread coordination.

---

## Implementation Flow

### Single-End Reads

1. **Initialization**:
   - Build Y contig mask if any Y-related feature enabled
   - Create `YReadNamesWriter` if `emit_y_read_names` enabled
   - Create `FastqSplitWriter` if `emit_y_noy_fastq` enabled
   - Initialize thread-local Y-hit tracking vectors

2. **Per-Batch Processing**:
   - Load batch of reads
   - Map reads in parallel (Y-hits tracked in thread-local vectors)
   - `#pragma omp taskwait` - wait for all mapping tasks
   - Merge thread-local Y-hits into batch-level set
   - Write Y read names for batch (if enabled)
   - Write FASTQ for batch (if enabled)
   - Swap batches and continue

3. **After All Batches**:
   - Merge all batch Y-hits into global set
   - Set Y-hit filter on writer (for BAM streams)
   - Continue with normal output phase

### Paired-End Reads

Same flow as single-end, with these differences:
- Two read batches (mate1 and mate2)
- Y-hit detection checks both mates (if either hits Y, pair is marked)
- FASTQ writing routes both mates to same stream based on pair-level Y-hit status

---

## Code Structure

### Y-Hit Tracking Expansion

All Y-hit tracking now gates on:
```cpp
if (mapping_parameters_.emit_noY_stream || mapping_parameters_.emit_Y_stream ||
    mapping_parameters_.emit_y_read_names || mapping_parameters_.emit_y_noy_fastq)
```

This ensures Y contig detection and tracking happens when any Y-related feature is enabled.

### Writer Classes

**YReadNamesWriter**:
- Normalizes names using static method `NormalizeReadName()`
- Maintains `written_read_ids_` set for deduplication
- Thread-safe: single writer instance, called from single thread after taskwait

**FastqSplitWriter**:
- Manages multiple output streams (Y and noY, potentially multiple mates)
- Handles both gzip and uncompressed output
- Detects FASTQ vs FASTA based on quality presence
- Thread-safe: single writer instance, called from single thread after taskwait

---

## Testing Recommendations

### Unit Tests Needed

1. **Name Normalization** (`YReadNamesWriter::NormalizeReadName`):
   - `@read123` → `read123`
   - `read123 /1` → `read123`
   - `read123/1` → `read123`
   - `read123/2` → `read123`
   - `@read123 comment` → `read123`
   - Edge cases: empty strings, special characters

2. **FASTQ Path Derivation** (`DeriveFastqOutputPath`):
   - `sample_R1.fastq.gz` → `sample_Y_R1.fastq.gz`
   - `sample_R2.fq` → `sample_Y_R2.fq.gz` (with gz compression)
   - `sample.fastq.gz` → `Y_reads.mate1.fastq.gz` (fallback)
   - Multiple files: `sample_R1.fastq.gz` → `sample_Y_R1.f1.fastq.gz`

3. **FASTQ Writing** (`FastqSplitWriter`):
   - FASTQ format (with qualities)
   - FASTA format (without qualities)
   - Comment preservation
   - Compression (gz vs none)

### Integration Tests Needed

1. **End-to-End Y Read Names**:
   - Small synthetic reference with `chrY` and `chr1`
   - Reads with some Y alignments, some non-Y
   - Validate: Y read names list contains expected normalized names
   - Validate: No duplicates

2. **End-to-End Y/noY FASTQ**:
   - Small synthetic reference with `chrY` and `chr1`
   - Reads with some Y alignments, some non-Y
   - Validate: Y + noY FASTQ counts == total reads
   - Validate: Paired-end routing (if either mate hits Y, both mates in Y FASTQ)
   - Validate: Output format (FASTQ vs FASTA based on input)

3. **Multiple Input Files**:
   - Two input files
   - Validate: File index included in output names
   - Validate: Names derived per input file, no collisions

---

## Known Limitations

1. **No ordering guarantees**: FASTQ and name-list outputs may not preserve input order
2. **Memory**: Y-hit tracking accumulates across all input files (may be large for very large datasets)
3. **Y contig detection**: Strict exact-match only (no `chrY_*` / `Y_*` matching)
4. **Multiple input files**: Outputs use per-file base names with file index suffixes

---

## Future Enhancements

1. **Opt-in Y contig detection**: Add flag for STAR-Flex parity (`chrY_*` / `Y_*` matching)
2. **Ordering option**: Optional flag to preserve input order (with performance cost)
3. **Streaming Y-hit tracking**: For very large datasets, could stream Y-hits to disk instead of keeping in memory
4. **Per-file naming**: Optional alternative naming policies if users want merged outputs

---

## Usage Examples

### Example 1: Y Read Names Only
```bash
chromap \
  -x index.idx \
  -r reference.fa \
  -1 reads_R1.fastq.gz \
  -2 reads_R2.fastq.gz \
  -o output.sam \
  --SAM \
  --emit-Y-read-names
```
Output: `output.Y.names.txt` with normalized read names

### Example 2: Y/noY FASTQ Only
```bash
chromap \
  -x index.idx \
  -r reference.fa \
  -1 reads_R1.fastq.gz \
  -2 reads_R2.fastq.gz \
  -o output.sam \
  --SAM \
  --emit-Y-noY-fastq \
  --emit-Y-noY-fastq-compression gz
```
Output: 
- `output_dir/reads_Y_R1.fastq.gz`
- `output_dir/reads_noY_R1.fastq.gz`
- `output_dir/reads_Y_R2.fastq.gz`
- `output_dir/reads_noY_R2.fastq.gz`

### Example 3: All Features Combined
```bash
chromap \
  -x index.idx \
  -r reference.fa \
  -1 reads_R1.fastq.gz \
  -2 reads_R2.fastq.gz \
  -o output.sam \
  --SAM \
  --emit-Y-bam \
  --emit-noY-bam \
  --emit-Y-read-names \
  --emit-Y-noY-fastq
```
Output: Y/noY BAM streams, Y read names list, and Y/noY FASTQ files

### Example 4: Custom Paths
```bash
chromap \
  -x index.idx \
  -r reference.fa \
  -1 reads_R1.fastq.gz \
  -2 reads_R2.fastq.gz \
  -o output.sam \
  --SAM \
  --emit-Y-read-names \
  --Y-read-names-output /path/to/y_names.txt \
  --emit-Y-noY-fastq \
  --Y-fastq-output-prefix /path/to/Y_reads \
  --noY-fastq-output-prefix /path/to/noY_reads
```

---

## Validation and Error Handling

### Input Validation
- Compression must be `gz` or `none` (exits with error if invalid)
- Requires explicit paths when primary output is stdout
- Validates FASTQ prefixes are provided when stdout is used

### Runtime Behavior
- Warns if no Y contigs found in reference
- Warns if multiple input files (names derived from first file)
- Handles empty Y-hit sets gracefully (creates empty Y files, full noY files)

### Error Handling
- File open failures: Exits with descriptive error message
- Write failures: Exits with descriptive error message
- All file handles properly closed in destructors

---

## Performance Considerations

1. **Y-hit tracking**: Minimal overhead (O(1) per mapping, thread-local vectors)
2. **Y-read-names writing**: Sequential I/O after batch completion (non-blocking)
3. **FASTQ writing**: Sequential I/O after batch completion (non-blocking)
4. **Memory**: Y-hit sets accumulate across batches (may be large for very large datasets)
5. **Compression**: Gzip compression adds CPU overhead but reduces disk I/O

---

## Compatibility

- **Backward compatible**: All new features are opt-in via flags
- **Works with existing features**: Can be combined with `--emit-Y-bam` / `--emit-noY-bam`
- **Independent operation**: Y-read-names and FASTQ can run without BAM output
- **Format support**: Works with all output formats (BED, SAM, BAM, CRAM, etc.)

---

## Documentation Updates Needed

1. **README.md**: Add new flags, usage examples, output naming rules
2. **chromap.1**: Add man page entries for new CLI options
3. **docs/**: Any relevant documentation updates (if needed)

---

## Conclusion

The implementation successfully adds Y read names and Y/noY FASTQ emission features to Chromap, following the design plan and maintaining compatibility with existing functionality. The code is structured, well-commented, and ready for testing and documentation updates.
