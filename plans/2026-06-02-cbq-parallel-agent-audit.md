# CBQ Parallel-Agent Commit Audit

Date: 2026-06-02
Branch: `feature/cbq-native-pr-coverage`

## Scope

Adversarial review of the CBQ changes that landed after the minimal ATAC CBQ
milestone:

- `00d0910` - direct CBQ decode into `SequenceBatch` buffers.
- `40ca05a` - indexed CBQ range producer/consumer.
- `3aedd17` - cached CBQ lane indexes for range reads.
- `65f1533` - CBQ Hi-C `.pairs` support.

Assumption for the review: CBQ bytes and `CBQINDEX` metadata are untrusted.

## Findings Fixed

1. `SequenceBatch::AssignKstring()` allocated `len + 1` without first rejecting
   `SIZE_MAX`. A hostile decoded CBQ field length at the platform maximum could
   wrap the allocation size and then write the terminator out of bounds. Fixed
   by rejecting `len == std::numeric_limits<size_t>::max()` before the addition.

2. `CbqPairedEndRangeBatchProducer` computed
   `(total_records + batch_size - 1) / batch_size` without guarding the addition.
   A hostile cached index claiming a near-`UINT64_MAX` record count could wrap
   `total_batches_`. Fixed with an overflow check before the ceil-div.

3. The indexed CBQ path's global read-id guard used
   `UINT32_MAX - (record_count - 1)`, which underflowed if hostile metadata
   supplied a record count above the `uint32_t` id space. Fixed by checking the
   lane offset/count against `UINT32_MAX + 1` records using subtraction only
   after proving the offset is in range.

4. The headerless barcoded-CBQ rejection text still said to re-encode
   "without --skip-headers". Fixed to say "re-encode with headers".

## Reviewed Without Further Findings

- Direct-buffer decode (`00d0910`): decoded CBQ spans are copied into kseq-owned
  `SequenceBatch` buffers before the reader advances; no borrowed CBQ block
  storage is retained by mapping code. `ResetLoadedSequences()` clears loaded
  counts and negative-sequence state before reuse.

- Indexed range producer (`40ca05a`): worker slots are owned by the producer,
  popped by the consumer, swapped into mapping batches, and returned only after
  the swap. Producer shutdown joins workers. The range loader creates local
  readers per range, so reader state is not shared across worker threads.

- Cached lane indexes (`3aedd17`): cached metadata is copied into
  `shared_ptr`-owned lane indexes captured by producer workers. `OpenRangeWithIndex()`
  re-checks header presence against the file, rejects missing blocks, invalid
  offsets, and non-monotonic cumulative counts. The read/barcode metadata path
  still rejects count mismatches and headerless barcoded input.

- Hi-C pairs (`65f1533`): the previous CBQ rejection was removed from both
  front-ends. The expanded matrix now compares FASTQ, CLI CBQ, and
  `chromap_lib_runner` CBQ `.pairs` rows.

- Existing `src/cbq_reader.cc` hardening remains present: block/header size
  caps, checked offset math, remaining-file-byte validation, bounded vector
  reserve, BitVector word-count checks, and zstd decompressed-size validation.

## Verification

- `make test-cbq-atac-smoke` passed.
- `make test-cbq-modality-matrix` passed 25/25, including Hi-C pairs,
  Y/noY FASTQ, bulk BAM, ChIP BAM, CRAM, `--read-group auto`, headerless CBQ,
  and read/barcode count-mismatch rejection.
- `make test-libchromap-core-smoke` passed.
- `make test-cbq-atac-100k` passed, with `CHROMAP_REQUIRE_CBQ_INDEX=1` on the
  CBQ CLI and libchromap runs, exercising the indexed range producer.

Artifacts are under `plans/artifacts/`:

- `cbq_atac_smoke/20260602T174708Z`
- `cbq_modality_matrix/20260602T174715Z`
- `chromap_core_smoke/20260602T174752Z`
- `cbq_atac_100k/20260602T174807Z`
