# CBQ Indexed Range Reader Optimization Runbook

Date: 2026-05-31

## Goal

Replace the current serial CBQ batch-loading bottleneck with an indexed,
range-pull reader that lets worker threads decode independent CBQ record ranges.
The intended first target is barcoded paired-end ATAC, because that is the
current production-shaped CBQ path and the easiest place to measure parity
against FASTQ.gz.

The target architecture is:

```text
CBQ block index
  -> record-range scheduler
  -> per-worker CBQ range readers
  -> private SequenceBatch buffers
  -> existing Chromap mapping/output path
```

This is preferred over deepening the single producer-consumer path. A single
producer can hide some read latency, but it cannot scale with mapper-side demand
unless it becomes a multi-producer design. Once there are multiple producers,
range-based ownership is simpler and avoids lock contention around one mutable
file cursor.

## STAR-suite Notes

The latest STAR-suite range-reader work is the most relevant implementation
reference:

- `/tmp/star-suite-cbq-libchromap-atac/docs/RUNBOOK_PF_CBQ_RANGE_COUNTING_20260531.md`
- `/tmp/star-suite-cbq-libchromap-atac/core/legacy/source/input/CbqInputModule.h`
- `/tmp/star-suite-cbq-libchromap-atac/core/legacy/source/input/CbqInputModule.cpp`
- `/tmp/star-suite-cbq-libchromap-atac/core/legacy/source/input/cbq_pf_adapter_harness.cpp`

That branch already has the prerequisite indexed CBQ range reader:

- `CbqInputModule::open_range(lane, first_record, record_count)`;
- tail `CBQINDEX` parsing;
- `current_lane_record_count()`;
- trimmed `CbqReadBatchView` slices for partial first/last blocks;
- `cbq-direct-decode` harness mode proving independent range readers process
  logical ranges exactly once.

The full-lane direct-reader sweep in that runbook is strong evidence that the
range reader itself can scale:

| Mode | 1 thread | 4 threads | 8 threads | 16 threads | 32 threads |
| --- | ---: | ---: | ---: | ---: | ---: |
| count-only | 5.72 s | 1.39 s | 0.86 s | 0.52 s | 0.36 s |
| sequence materialized | 6.56 s | 1.76 s | 0.93 s | 0.64 s | 0.41 s |

Useful lessons to carry over:

- The borrowed `CbqReadBatchView` model is right: decoded records point into
  block-owned storage, and queued consumers retain a shared backing object.
- Production adapters should fill native consumer buffers directly instead of
  materializing FASTQ text.
- The adapter/harness split is valuable: reader correctness is tested before
  full mapper integration.
- Record order is a capability that must be explicit. Chromap ATAC needs a
  stronger read/barcode alignment contract than STAR's generic input contract.
- The next bottleneck after range decode is often the queue/stream handoff.
  Chromap should avoid replacing one serial reader with a shared queue if a
  direct range-to-worker path is feasible.

Important Chromap-specific difference:

- STAR's PF range plan targets `process_features`, where each worker owns a
  PF counting state and finalization merges per-thread hashes. Chromap's mapper
  and writers have different state and determinism constraints, especially
  `read_id` ordering, BAM sorting tie-breaks, Y/noY sidecars, and low-memory
  spill records.

## Current Chromap Baseline

The existing Chromap CBQ path already has:

- native CBQ CLI and libchromap input options,
- paired-read CBQ plus optional barcode CBQ,
- `CbqReadBatchView` borrowed block views,
- direct fill into `SequenceBatch` sequence buffers,
- forward 8-bit and 16-bit CBQ-to-ASCII lookup materialization,
- indexed `OpenRange()` readers with tail `CBQINDEX` parsing,
- an ordered paired-end range batch producer for indexed CBQ mapping,
- optional reverse-complement prefill plumbing, currently disabled because a
  single producer moved parallel mapper work onto the producer critical path,
- parity gates for ATAC, ENCODE cross-assay smoke, BAM/Y-noY paths, and
  libchromap.

The relevant Chromap implementation files are:

- `src/cbq_reader.h`
- `src/cbq_reader.cc`
- `src/cbq_batch_producer.h`
- `src/cbq_batch_producer.cc`
- `src/sequence_batch.h`
- `src/sequence_batch.cc`
- `src/chromap.cc`
- `src/chromap.h`

The single producer-consumer path is now a compatibility fallback for legacy
CBQs without a `CBQINDEX` footer. Do not build more complex scheduling around a
shared `CbqLaneReader`. Indexed CBQs should use worker-owned range readers and
the ordered range batch producer.

## Design Contracts

### Record Range

A record range is half-open:

```text
[global_record_begin, global_record_end)
```

Within one lane, `global_record_begin=0` means the first CBQ record after the
file header. For paired-end CBQ, one record contains both read mates. For ATAC
with barcodes, the paired-read CBQ and barcode CBQ must be read with the same
range boundaries.

### Read/Barcode Alignment

Barcoded ATAC requires:

```text
paired_read_cbq record i == barcode_cbq record i
```

The indexed reader must preserve that contract even when the two files have
different CBQ block boundaries. Workers are assigned record ordinal ranges, not
byte ranges. Each worker resolves those record ranges into the needed blocks for
each file independently.

Headers remain required for barcoded production CBQ unless we add a stronger
out-of-band provenance check. The name-match guard should remain active.

### Sequence IDs

Chromap uses `SequenceBatch::GetSequenceIdAt()` as `read_id` for mapping,
deduplication, BAM sorting tie-breaks, Y/noY sidecars, and low-memory spill
records. Parallel range loading must not let private `SequenceBatch` instances
assign local IDs starting at zero.

Chromap now has a supported way to set sequence IDs from the global record
ordinal:

```text
SequenceBatch::CommitLoadedSequenceBufferWithId(..., uint32_t read_id)
```

The first implementation should fail loudly if a global record ordinal exceeds
`uint32_t`, because Chromap mapping structs currently store `read_id` as
`uint32_t`.

### Output Determinism

Mapping output does not need input batches to be decoded serially, but existing
sort and duplicate paths rely on stable read IDs for deterministic tie-breaks.
The indexed reader must therefore preserve read IDs even if decoded ranges
finish out of order.

For unsorted streaming outputs, keep the implementation conservative: deliver
filled batches to the existing mapper loop in range order. The current range
producer allows worker batches to finish out of order but only releases them to
the mapper by increasing range index.

## Phase 1: Parse CBQ Tail Index Metadata

Goal: read the CBQ tail `CBQINDEX` once and know how to seek to any block.

Add a lightweight index structure in `src/cbq_reader.h`:

```text
CbqBlockIndexEntry
  uint64_t offset
  uint64_t cumulative_records

CbqLaneIndex
  FileHeaderFields file_header
  vector<CbqBlockIndexEntry> blocks
  uint64_t total_records
```

Implementation notes:

- Follow the STAR-suite `/tmp` branch: read the 16-byte footer, verify trailing
  `CBQINDEX`, find the 24-byte index header, decompress the index payload, and
  parse 16-byte `(offset, cumulative_records)` entries.
- Range mode should require the index in the first implementation. A scanner
  fallback can be added later for diagnostics, but it should not be the main
  range path.
- Validate index/header/footer compressed-size consistency.
- Validate offsets are after the 64-byte file header and before the index
  header.
- Validate cumulative records are monotonic.
- Preserve existing file-header validation and mate-count checks.
- Expose `current_lane_record_count()` or the Chromap equivalent so schedulers
  can partition without decoding.

Acceptance:

- Add a small harness or unit test that builds an index for synthetic paired and
  barcode CBQs and verifies total record counts and block offsets.
- Validate compressed and level-0 CBQ files.
- Malformed/truncated fixtures must fail cleanly.

## Phase 2: Add Record-Range Reader

Goal: read a requested record range, even when it spans partial blocks.

Add:

```text
CbqLaneReader::OpenRange(lane_index, first_record, record_count)
CbqLaneReader::CurrentLaneRecordCount()
CbqLaneReader::NextBatch(&batch)
```

The range reader should:

- binary-search the block index for `first_record`;
- seek the worker-owned file handle to the indexed block offset;
- decode each overlapping block;
- expose only the subrange requested from the first and last blocks;
- keep backing storage alive for all returned subranges;
- return fewer records only at EOF;
- report a clear error if `first_record` is past `total_records`.

For the first Chromap integration, returning one or more block-backed subranges
is better than copying records into a synthetic contiguous `CbqReadBatchView`.
The `SequenceBatch` fill layer can iterate subranges and write directly into
its destination buffers.

Implementation notes:

- Mirror STAR's `open_range()` behavior closely: each worker owns a reader
  instance, range mode trims `next_batch()` output, and sequential `Open()` /
  `NextBatch(max_records)` remains available as the oracle path.
- Preferred first pass: each worker owns an `std::ifstream`. Later pass: switch
  to POSIX `pread()` only if profiling shows seek/file-handle overhead matters.
- Do not use one shared `CbqLaneReader` with a mutex. The existing reader owns
  one `ifstream`, one decoded block, and one mutable `record_index`; locking it
  would serialize the hot path.
- Keep the current zstd output-size validation and bounds checks.

Acceptance:

- Range reader dumps for random ranges match the same slice from the sequential
  reader.
- Test ranges that start/end inside a block, exactly on block boundaries, and
  span multiple blocks.
- Add a direct-decode harness mode analogous to STAR's `cbq-direct-decode`.
- Implemented as `make test-cbq-range-reader`.

## Phase 3: ATAC Range Scheduler

Goal: assign independent record ranges to workers while preserving Chromap's
current mapping assumptions.

Initial scheduler:

```text
range_size = read_batch_size_
next_range_begin = atomic counter or serial planner
workers = min(num_threads_for_cbq, number_of_ranges)
```

For each range:

1. Read paired-read CBQ range.
2. Read barcode CBQ range if present.
3. Validate equal record counts.
4. Validate header/name alignment when barcodes are present.
5. Fill private `SequenceBatch` buffers for R1, R2, and barcode.
6. Publish the completed batch with its `range_begin`.

Implementation status:

- `CbqPairedEndRangeBatchProducer` schedules whole `read_batch_size_` record
  ranges.
- Each producer worker opens its own paired-read and barcode `OpenRange()`
  readers.
- Ready batches are keyed by range index and popped in source order.
- Read IDs are assigned from `global_lane_offset + lane_read_ordinal - 1` via
  `SequenceBatch::CommitLoadedSequenceBufferWithId()`.
- Barcoded CBQs still require headers and equal read/barcode record counts.
- CBQs without a `CBQINDEX` footer fall back to the previous sequential
  producer.
- `CHROMAP_REQUIRE_CBQ_INDEX=1` disables that fallback and makes Chromap fail
  loudly if an indexed range producer cannot be initialized.
- The read-pair and barcode lane `CBQINDEX` metadata is parsed once during lane
  setup and passed as immutable cached metadata to worker-owned range readers.
  This avoids rereading and decompressing the tail index for every
  mapper-sized range batch.

Verification run on 2026-05-31:

- `make test-cbq-range-reader`
- `BQTOOLS=/tmp/star_suite_bqtools/bin/bqtools make test-cbq-atac-smoke`
- `make test-libchromap-core-smoke`
- `CBQ_ORDERED_ENCODER=/tmp/cbq_ordered_encoder make test-cbq-atac-100k`
- `CHROMAP_GRCH38_REF=/mnt/pikachu/Homo_sapiens.GRCh38.dna.primary_assembly.fa
  CHROMAP_GRCH38_INDEX=/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index
  CBQ_ORDERED_ENCODER=/tmp/cbq_ordered_encoder ENCODE_SKIP_PREPARE=1
  THREADS=4 make test-encode-cbq-cross-assay-smoke`

The cross-assay smoke covered ChIP, bulk ATAC, scATAC, and Hi-C indexed CBQ
inputs through both `chromap` and `chromap_lib_runner`.

Full 3K PBMC ATAC timing on the SSD fixture with warmed cache, 8 total threads,
and BED output directed to `/dev/null` is recorded under
`plans/artifacts/cbq_atac_full_timing/20260531T081906Z/`:

| Input | Wall | User | Sys | Max RSS | Output mappings |
| --- | ---: | ---: | ---: | ---: | ---: |
| FASTQ.gz | 3:04.47 | 1162.27 s | 10.73 s | 24,699,564 KB | 53,969,811 |
| CBQ level 0 | 2:52.27 | 1133.20 s | 11.26 s | 25,344,496 KB | 53,969,811 |

The CBQ timing run used `CHROMAP_REQUIRE_CBQ_INDEX=1` and logged two indexed
range workers per lane.

The PF range-counting design is the right mental model:

```text
CBQ logical record range -> worker-owned consumer state
```

For Chromap, the first worker-owned state should be the reader plus private
`SequenceBatch` buffers. Deeper mapper-state ownership can come later if the
ordered batch handoff remains a bottleneck.

Keep mapper consumption ordered for the first implementation:

```text
next_range_to_map = 0
completed_ranges[range_begin] -> mapper loop
```

This keeps output behavior close to the existing loop while allowing read,
decompression, and materialization work to happen in parallel ahead of mapping.

Thread budget:

- Start with `num_cbq_workers = max(1, min(4, num_threads / 2))`.
- Reserve enough OpenMP mapper threads that mapping throughput does not drop.
- Record the split in benchmark manifests.

Acceptance:

- Existing FASTQ path is untouched.
- CBQ ATAC smoke passes.
- CBQ 100K parity passes.
- Y/noY FASTQ sidecars and BAM/Y-noY sidecars still pass.

## Phase 4: Benchmark And Revisit Reverse-Complement Prefill

Goal: measure range-pull independently, then decide whether direct
reverse-complement prefill should be enabled.

Benchmark order:

1. Current sequential/batch CBQ reader.
2. Indexed range reader with one worker.
3. Indexed range reader with two workers.
4. Indexed range reader with four workers.
5. Optional: enable direct negative-strand prefill and repeat.

Use the existing under-minute 2M ATAC fixture first:

```text
plans/artifacts/cbq_producer_downsample_2m_20260530T160612Z/
```

Then run the full PBMC 3K ATAC set on SSD once the small fixture shows a clear
signal.

Record:

- command line,
- git commit or dirty-worktree note,
- fixture, reference, and index paths,
- CBQ worker count and mapper thread count,
- output directory under `plans/artifacts/`,
- wall/user/sys/max RSS,
- read/mapping/output counts,
- parity result against FASTQ.gz.

Expected outcome:

- One range worker should be within noise of the current reader.
- Two or more range workers should improve CBQ wall time if CBQ decoding and
  materialization are still visible on the critical path.
- A direct-decode-only benchmark should scale much more strongly than full
  Chromap mapping. If it does not, debug the reader before touching mapper
  scheduling.
- Reverse-complement prefill should only be enabled if range workers keep that
  work off the mapper critical path and improve end-to-end wall time.

## Phase 5: Cleanup And API Decision

Once range-pull passes parity and timing gates:

- remove or retire the experimental `CbqPairedEndBatchProducer` path if it is no
  longer the best implementation;
- keep the sequential `CbqLaneReader` as a simple compatibility reader and test
  oracle;
- document the indexed reader in the CBQ implementation runbook;
- ensure `libchromap.a` includes the same source of truth as the CLI;
- update STAR-suite integration notes if the libchromap ATAC entrypoint changes.

## Future Optimization Notes

If CBQ adoption becomes broad enough to justify deeper Chromap internals work,
the obvious but higher-blast-radius target is a packed/base-view path for the
mapper. CBQ already carries packed sequence words and N positions; consumers
such as seed generation and N checks could use those directly instead of always
reading ASCII. Do not treat this as a CBQ integration cleanup task: ASCII
materialization into `SequenceBatch` is currently cheap, and changing the mapper
base representation would be a Chromap optimization project touching FASTQ,
CBQ, trimming, barcode, mapping, and writer assumptions.

The first simple CBQ-specific optimization has been implemented: Chromap now
parses and validates each lane's `CBQINDEX` once, then passes the cached
metadata to worker-owned readers through `OpenRangeWithIndex()`. Two possible
follow-ups remain if timing shows index/open setup is still visible:

- keep one read-pair and barcode range reader open per producer worker and reset
  it across adjacent assigned ranges;
- schedule coarser worker-owned lane ranges and split them into mapper-sized
  `SequenceBatch` batches inside the worker, matching STAR-suite's direct PF
  range-counting shape more closely.

Acceptance for either should be the existing CBQ range reader smoke, ATAC
smoke, 100K parity, ENCODE CBQ cross-assay smoke, libchromap core smoke, and a
small timing fixture showing the setup cost is actually visible.

## Test Matrix

Required before merge:

```bash
make chromap chromap_lib_runner -j4
BQTOOLS=/tmp/star_suite_bqtools/bin/bqtools make test-cbq-atac-smoke
CBQ_ORDERED_ENCODER=/tmp/cbq_ordered_encoder make test-cbq-atac-100k
CHROMAP_GRCH38_REF=/path/to/genome.fa \
CHROMAP_GRCH38_INDEX=/path/to/genome.index \
CBQ_ORDERED_ENCODER=/tmp/cbq_ordered_encoder \
make test-encode-cbq-cross-assay-smoke
```

Also run the BAM/Y-noY gates used for the previous CBQ merge if they are not
already included in the smoke target for the branch.

## Open Questions

- Does the current encoder write enough small blocks to expose useful range
  parallelism, or do we need to tune CBQ encoder block size for Chromap-sized
  batches?
- Should range workers fill `SequenceBatch` objects directly, or should they
  first publish `CbqReadBatchView` subranges and let the mapper thread fill
  `SequenceBatch`?
- Should range scheduling be lane-local or global across lanes for multi-lane
  ATAC?
- Should indexed range-pull eventually be exposed for ChIP and Hi-C CBQ smoke
  paths, or stay ATAC-only until performance justifies broader integration?

## Stop Conditions

Pause and reassess if:

- read IDs or output ordering change FASTQ-vs-CBQ parity;
- barcode range alignment requires headers to be absent or unverifiable;
- memory grows with worker count because decoded block backings are retained too
  long;
- range-pull improves CBQ decoding time but harms total wall time by starving
  mapper threads;
- block granularity is too coarse for the available CBQ fixtures.
