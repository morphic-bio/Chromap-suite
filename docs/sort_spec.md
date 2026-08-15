# chromap Coordinate Sort Specification

## Overview

When `--sort-bam` is enabled, chromap outputs BAM/CRAM records in coordinate-sorted order using an internal spill-to-disk sorter. This document describes the exact sorting specification.

## Sort Key (STAR-Flex spec)

BAM/CRAM records are sorted using the following **primary key** (numeric fields, ascending):

1. **tid** (reference ID) - unmapped reads (tid=-1) normalized to INT32_MAX to sort last
2. **pos** (position) - unmapped reads (pos=-1) normalized to INT32_MAX to sort last
3. **flag** (SAM flag)
4. **mtid** (mate reference ID)
5. **mpos** (mate position)
6. **isize** (insert size)

**Tie-break**: `read_id` (uint64, deterministic within run)

The `read_id` is an internal identifier assigned during mapping that is deterministic for identical input across runs (with fixed threading and compression settings).

## Header

When `--sort-bam` is enabled, the output header includes:

```
@HD VN:1.6 SO:coordinate
```

When `--sort-bam` is disabled, the header includes:

```
@HD VN:1.6 SO:unknown
```

## Compatibility

**NOT bit-for-bit identical to samtools sort**:
- samtools uses QNAME (read name) for tie-breaking
- chromap uses internal `read_id` (faster, deterministic within same input order)
- Records at the same coordinate may appear in different order than samtools

**Guarantees**:
- **Deterministic**: same input produces identical output across runs (with `--hts-threads 1` and fixed compression level)
- **Coordinate-sorted**: valid for indexing with `samtools index` or `--write-index`
- **Unmapped reads sort last**: all unmapped reads appear after mapped reads

## Indexing

Use `--sort-bam --write-index` to generate `.bam.bai` or `.cram.crai` indexes.

**`--write-index` requires `--sort-bam`**: htslib's `sam_idx_init()`/`sam_idx_save()` only produce valid indexes when records are written in coordinate order. Without sorting, multithreaded mapping writes records in non-deterministic order, resulting in invalid or empty index files. The `--sort-bam` flag ensures records are buffered, sorted, and written in coordinate order before index finalization.

## Memory Management

The sorter uses spill-to-disk when memory usage exceeds `--sort-bam-ram` (default: 8GB). Sorted chunks are written to temporary files and merged using a k-way merge algorithm.

## Low-Memory Mode Interaction

`--low-mem` and `--sort-bam` are **composable**:
- `--low-mem`: handles mapping overflow during alignment phase
- `--sort-bam`: handles BAM output sorting during output phase

When both are enabled, mappings are processed through low-memory overflow, then converted to BAM and sorted.

## Example

```bash
# Generate coordinate-sorted BAM with index
chromap --BAM --sort-bam --write-index -x index -r ref.fa -1 reads.fq -o output.bam

# Force deterministic output (single-threaded compression)
chromap --BAM --sort-bam -t 1 --hts-threads 1 -x index -r ref.fa -1 reads.fq -o output.bam

# Limit sorting memory to 512MB
chromap --BAM --sort-bam --sort-bam-ram 512M -x index -r ref.fa -1 reads.fq -o output.bam
```

## Alternative: samtools Pipe

For users who prefer exact samtools ordering (QNAME tie-break), use a pipe:

```bash
chromap --SAM -x index -r ref.fa -1 reads.fq -o - | samtools sort -o output.bam
```

Note: This bypasses internal sorting but loses Y/noY stream capability during sort.
