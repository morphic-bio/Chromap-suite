# fqgzip paired FASTQ input adapter

Status: opt-in adapter implemented and synthetic-parity tested on 2026-07-30.
No non-synthetic or 100K execution was launched for this implementation.

## Design

Chromap Suite can consume ordinary paired `fastq.gz` files through fqgzip's
progressive record-aligned shard API. R1 and R2 use independent decode/index
readers. They coordinate only at canonical read-name synchronization and each
subsequent pair is validated before it enters Chromap.

The adapter runs as one bounded background producer. It converts validated
fqgzip pair views into Chromap's existing `SequenceBatch` objects; alignment,
insert-size estimation, mapping policy, ordering, output, and peak calling stay
inside one unchanged Chromap process. It does not run one mapper per shard.

The build pins fqgzip revision
`bcdc92d8d6e7f1258724195c2789c677eb33a1a6` and rejects a different or
tracked-dirty fqgzip checkout.

## Build and use

Ordinary builds remain independent of fqgzip:

```bash
make clean
make -j8 WITH_FQGZIP=0 chromap
```

Build the opt-in adapter:

```bash
make clean
make -j8 WITH_FQGZIP=1 FQGZIP_DIR=/mnt/pikachu/fgqzip chromap
```

Use it with direct bulk paired gzip FASTQ inputs:

```bash
./chromap --fqgzip --fqgzip-shards 8 --fqgzip-threads 8 \
  -r reference.fa -x reference.index \
  -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz \
  --SAM -o alignments.sam -t 24
```

When `--fqgzip-threads` is omitted, R1 and R2 receive one independent worker
each. Larger explicit budgets are divided across the mates. Reader workers are
additional to the `-t` mapping workers, so benchmarks must record both and keep
their combined host allocation comparable. When `--fqgzip-shards` is omitted,
the requested shard count follows the reader-worker budget.

## Current scope

- Direct paired ordinary-gzip FASTQ input is supported, including multiple
  paired lanes.
- The first adapter is bulk paired-end only. Separate barcode-file/scATAC
  three-stream input, single-end input, stdin, and plain FASTQ stay on Chromap's
  existing reader.
- `--fqgzip` cannot be combined with `--input-format cbq`.
- fqgzip continuously rejects missing, reordered, duplicated, or mismatched
  mates rather than silently handing incomplete pairs to Chromap.
- The current fqgzip implementation materializes a bounded shard before its
  record iterator. Lazy batch delivery is a later performance refinement.

## Verification

Both the fqgzip-enabled and normal builds succeed. The opt-in synthetic input
smoke compares normalized SAM output from native gzip and fqgzip paths and
rejects a deliberately dropped R2 record:

```bash
make clean
make -j8 WITH_FQGZIP=1 FQGZIP_DIR=/mnt/pikachu/fgqzip chromap
TEST_FQGZIP_INPUT=1 \
  CHROMAP_ARTIFACT_ROOT=/tmp/chromap_fqgzip_smoke_20260730 \
  bash tests/run_input_format_smoke.sh
```

Generated FASTQs, indexes, alignments, and logs remain untracked under the
artifact root. A 100K or larger validation requires a frozen command and
explicit authorization under the repository repeat-run policy.
