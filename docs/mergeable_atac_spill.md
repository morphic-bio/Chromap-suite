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
The ordinary `-o` output is suppressed. Spill creation also forces read-local
deterministic mapping: the mutable cross-read candidate cache is disabled and
paired-end multimapping reservoir selection is seeded from the stable read
name. An ordinary one-process run used as an exact parity control must add
`--deterministic-mapping`. Historical ordinary runs retain the candidate cache
unless that flag is supplied.

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

BED is a terminal representation, not materializer working state. The global
merge, barcode correction, and deduplication first write versioned `ATMBLK1`
binary blocks. Only after that binary is complete does the terminal exporter
format BED rows. `--threads` parallelizes this final block export while
preserving byte order with disjoint positioned writes. The binary is normally
temporary; add `--materialized-binary fragments.atmb1` to preserve it for a
later BED export or diagnostics.

Every new worker spill also publishes an `ATACHOT1` companion at
`<spill>.hot`. It is produced during the worker's existing sorted overflow
merge, so records are not sorted a second time. The companion is partitioned
by reference and contains only the fixed-width evidence needed for BED
correction and deduplication: the local 64-bit read id, raw packed barcode and
N mask, coordinates, alignment lengths, MAPQ/direction/uniqueness, Y state,
and fixed-width raw barcode qualities. A bulk record is 32 bytes; a 16-base
scATAC record is 48 bytes. It never carries BAM strings or CIGAR payloads.

For the common BED path, when gathered multimapping allocation and summary
output are not requested, the materializer opens independent positioned-read
cursors for each reference and performs the shard k-way merge, global barcode
correction, and deduplication in parallel.
Each reference task emits a compact binary post-dedup partition; the final
assembly follows reference order and the terminal exporter alone constructs
BED text. For barcodes up to 16 bases, those temporary rows are the same 16-byte
shape as `ATMBLK1`. Empty reference partitions do not create temporary files.
The in-memory decoder mirrors the compact hot record and does not construct the
two `SAMMapping` objects owned by the full spill record.

For 17--32-base barcodes, only the ephemeral per-reference partition widens to
20 bytes so it can retain the corrected 64-bit packed barcode. Ordered assembly
then assigns deterministic 32-bit dictionary ids and writes the same 16-byte
`ATMBLK1` records and 64-bit dictionary footer used by the established long-
barcode path. Neither durable format changes, and the extra four bytes exist
only until the reference partitions are assembled.

Barcode translation remains on this parallel path. Reference tasks retain the
raw packed barcode. During ordered binary assembly, `ATMBLK1` assigns dense
dictionary ids to the raw keys; the terminal exporter translates each distinct
dictionary key once and reuses the resulting text for every fragment. No BED
string or translated-barcode string is created in a worker partition.

The complete `ATACMS3` parent remains authoritative for BAM/CRAM, summary
synthesis, and gathered multimapping allocation.
Those modes use a bounds-checked in-memory decoder for the existing payload
codec. The reader reuses its payload buffer and decodes fields directly,
without constructing a synthetic `FILE` or calling `fmemopen`/`fclose` for
every mapped record. This changes neither the `ATACMS3` envelope nor its record
layout.
The parent advertises the companion only after the companion has been flushed,
`fsync`ed, and atomically renamed; an advertised missing or inconsistent hot
file is a hard error rather than a silent fallback.

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

`ATACHOT1` has its own magic, version, endian marker, parent identity fields,
record width, and per-reference offset/count/start-bound directory. Readers
validate the directory, exact file size, parent record count, shard identity,
whitelist fingerprint, coordinate bounds, barcode bounds, local read ids, and
per-partition sort order. It intentionally has no record checksum pass.

Read ids are unsigned 64-bit ordinals throughout worker spill records,
materialization, duplicate tie-breaking, and BAM sorting. The codec is uniform;
there is no conditional 32/64-bit mode. Summary counters are also unsigned
64-bit values, so total, mapped, duplicate, low-MAPQ, cache-hit, and cardinality
fields do not wrap above `UINT32_MAX`. Multimapping allocation intentionally
buffers the gathered, deduplicated candidates because allocation weights
depend on unique mappings across every shard. The normal non-allocation path
remains a bounded k-way streaming merge.

Mergeable workers disable the mutable candidate cache so mapping decisions are
independent of worker boundaries and scheduling. Their cache-hit/cardinality
diagnostics therefore describe the disabled-cache policy and match an ordinary
`--deterministic-mapping` control. Exact scatter/gather parity must not compare
against a historical-cache run, whose cross-read history is inherently
partition-dependent.

Barcode correction uses the complete input histogram in both ordinary and
mergeable-spill runs. The historical 20-million exact-whitelist-observation
cutoff is not part of this contract. Because shard ranges are validated as
complete and non-overlapping, summing their histograms reconstructs the same
complete-input model as a single-process run.

The post-dedup `ATMBLK1` hot record is 16 bytes: 16-bit reference id, 32-bit
start, 16-bit fragment length, a 32-bit barcode value, and one byte each for
duplicate count, MAPQ, flags, and reserve. The terminal end coordinate is
reconstructed as `start + fragment_length`. The writer rejects zero/out-of-
reference spans and reference dictionaries above 65,536 entries. The worker
also rejects fragment lengths outside 1--65,535 before narrowing to the spill
field, so an oversize insert cannot wrap. Chromap's
packed barcode sequence remains 64-bit so barcodes up to 32 bases are
representable. Barcodes up to 16 bases are stored directly in the 32-bit value,
avoiding a hash lookup in the merge. For 17--32-base barcodes (or translated
barcode output), each distinct materialized key is written once in the file
dictionary and the record value is its dense id. The block directory stores
offset/count and coordinate bounds. It intentionally has no per-block checksum
pass: file structure and bounds are validated, while transfer integrity remains
the responsibility of the workflow's normal artifact digest.

Run the hermetic gate with:

```sh
make test-atac-mergeable-spill-materializer
```

The gate covers the v3 header and evidence stream, shuffled ordinal input,
late-bound and explicit read ranges, global barcode correction, summary
synthesis, cell- and bulk-level
deduplication, gathered multimapping allocation, unsigned-64-bit ordinals above
`UINT32_MAX`, summary counts above `UINT32_MAX`, BED, sorted BAM, AEV1, CRAM,
Y/no-Y BAM routing, compact binary header/dictionary round-trip, and serial
versus parallel multi-block BED byte parity. The hot/full parity gate includes
17- and 32-base raw dictionary barcodes, retained unresolved evidence, and
translated 17-base barcodes. It also rejects a corrupted full BAM payload
through the direct bounded decoder.
