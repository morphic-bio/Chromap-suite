# Mergeable ATAC spill and materializer

This is an opt-in post-alignment scatter/gather path. It is inactive unless a
worker is launched with `--create-mergeable-spill-record`. A worker maps its
own contiguous synchronized read/barcode range, writes one sorted,
pre-deduplication `AtacSpillRecord` file, publishes it atomically, and exits
without writing ordinary mapping outputs. Barcode correction is deliberately
deferred: workers map raw barcodes and carry the raw packed barcode, N mask,
qualities, and a complete local exact-whitelist histogram. The standalone
materializer consumes regular files and has no dependency on a particular
sharder or scheduler.

The worker still loads the reference FASTA and Chromap index because alignment
happens there. BED and BAM materialization load neither: reference names and
lengths, fragment fields, and retained BAM mate payloads are carried by the
spill files. CRAM output alone requires the reference FASTA for reference-based
encoding; it never loads the Chromap index.

## Worker contract

Each worker supplies four identity fields:

- `--mergeable-spill-sample-id`: stable sample identity.
- `--mergeable-spill-input-id`: stable synchronized-input-set identity.
- `--mergeable-spill-shard-ordinal`: zero-based ordinal.
- `--mergeable-spill-shard-count`: exact number of expected shards.

The default native in-place workflow omits
`--mergeable-spill-first-global-read` and
`--mergeable-spill-input-record-count`. Each worker records the synchronized
count it observes at EOF, and the materializer derives global prefixes by
ordinal from the complete ordered spill set. An application that already has a
certified range table may provide both options together; supplying only one is
rejected.

For barcoded ATAC, every worker needs only its own synchronized read/barcode
range and the common whitelist. While mapping that range, it creates a complete
local exact-whitelist histogram from the same barcode batches written to the
spill; mergeable-spill mode does not sample, pre-scan, or reopen the barcode
input. Library callers may instead install a `PairedEndReadProvider`. The
provider fills Chromap's normal paired-end `SequenceBatch` objects, so a caller
can decode its assigned compressed ranges directly from a shared filesystem;
no FIFO, relay process, or decoded shard file is required. The provider API is
generic and contains no sharder or scheduler types, and is currently accepted
only with mergeable-spill output. The barcode length is taken from the common
whitelist and checked against every input barcode. A whitelist-key fingerprint
and the correction policy are recorded in the header. Workers do not correct
or reject barcodes before spill, so no globally resolvable record is lost.

FASTQ example for shard 3 of 8:

```sh
chromap --preset atac \
  -r ref.fa -x ref.index \
  -1 shard3.R1.fastq.gz -2 shard3.R2.fastq.gz \
  -b shard3.barcode.fastq.gz \
  --barcode-whitelist whitelist.txt \
  --create-mergeable-spill-record shard3.atacms \
  --mergeable-spill-sample-id sample-A \
  --mergeable-spill-input-id atac-lanes-A \
  --mergeable-spill-shard-ordinal 3 \
  --mergeable-spill-shard-count 8
```

`--read-format` and all mapping/correction policy options must match the
single-process baseline. Spill creation automatically enables low-memory mode.
The ordinary `-o` output is suppressed.

The materializer's ordinal prefix is not used to order fragments. It rebases
the shard-local read ids embedded in the fragment and both BAM mates, making
them globally unique and restoring the identifiers from the unsharded run.

## Materialization

Inputs may be listed in any order. The materializer validates that ordinals
`0..shard_count-1` occur exactly once, derives late-bound read prefixes (or
validates explicitly declared ranges as contiguous), and requires all policy,
reference, schema, sample/input identity, whitelist, and correction fields to
agree. It sums the shard histograms, applies barcode
correction or rejection to each raw record, restores corrected sort order
within each genomic-coordinate group, and then performs global PCR
deduplication, optional multimapping allocation, MAPQ filtering and Tn5 shift
under the recorded policy. For barcoded bulk-level deduplication, ties by
barcode duplicate count use the gathered global abundance, matching the
ordinary Chromap policy.

```sh
chromap_atac_spill_materializer \
  --spill shard7.atacms --spill shard0.atacms --spill shard1.atacms \
  --spill shard2.atacms --spill shard3.atacms --spill shard4.atacms \
  --spill shard5.atacms --spill shard6.atacms \
  --output-bam atac.bam \
  --fragments fragments.tsv.gz \
  --evidence fragments.aev1 \
  --summary alignment_summary.csv \
  --output-noY atac.noY.bam --output-Y atac.Y.bam \
  --sort-bam --write-index --threads 8
```

Use `--output-bed fragments.bed` when BAM/AEV1 output is not required. BAM is
coordinate-sorted by default; `--no-sort-bam` selects the canonical
fragment-merge order. Use `--output-cram atac.cram --reference ref.fa` for
CRAM.

## Format and failure behavior

The `ATACMS3` envelope versions the shard contract separately from the
`AtacSpillRecord` payload codec. It contains shard identity/range metadata,
reference dictionary, mapping and correction policy, schema mask, record
count, whitelist fingerprint, sorted whitelist keys, and the shard-local
abundance counts. Each mapped payload carries the raw barcode N mask and
quality string. The envelope also carries one fixed input-summary observation
per synchronized pair (20 bytes plus barcode qualities): raw barcode evidence
and the mapped cache-slot hits needed to synthesize total, cache, cardinality,
FRIC, and estimated-FRIP fields after global correction. This input stream is
why unmapped and ultimately rejected reads remain representable at gather.
Writers use a temporary file, flush and `fsync` it, then publish by rename.
Readers reject unsupported versions, truncation, trailing bytes, unsorted
records, schema mismatch, bad local read ids, incomplete summary evidence, and
incomplete or inconsistent shard sets.

Read ids are unsigned 64-bit ordinals throughout worker spill records,
materialization, duplicate tie-breaking, and BAM sorting. The codec is uniform;
there is no conditional 32/64-bit mode. Summary counters are also unsigned
64-bit values, so total, mapped, duplicate, low-MAPQ, cache-hit, and cardinality
fields do not wrap above `UINT32_MAX`. Multimapping allocation intentionally
buffers the gathered, deduplicated candidates because allocation weights
depend on unique mappings across every shard. The normal non-allocation path
remains a bounded k-way streaming merge.

Cache-hit and cache-cardinality fields describe the actual worker caches. They
are merged exactly from shard evidence, but a sharded execution can legitimately
have different cache warm-up behavior than a single-process execution. Mapping,
duplicate, low-MAPQ, total, barcode, fragment, and BAM correctness do not depend
on that diagnostic cache topology.

Barcode correction uses the complete input histogram in both ordinary and
mergeable-spill runs. The historical 20-million exact-whitelist-observation
cutoff is not part of this contract. Because shard ranges are validated as
complete and non-overlapping, summing their histograms reconstructs the same
complete-input model as a single-process run.

Run the hermetic gate with:

```sh
make test-atac-mergeable-spill-materializer
```

The gate covers the v3 header and evidence stream, shuffled ordinal input,
late-bound and explicit read ranges, global barcode correction, summary
synthesis, cell- and bulk-level
deduplication, gathered multimapping allocation, unsigned-64-bit ordinals above
`UINT32_MAX`, summary counts above `UINT32_MAX`, BED, sorted BAM, AEV1, CRAM,
and Y/no-Y BAM routing.
