# ATAC Runtime Spill Schema Runbook

This runbook captures the design and test plan for replacing the current
output-specific ATAC low-memory spill records with one runtime-configured spill
schema. The immediate trigger is STAR-suite multiome production: BAM plus binary
sidecar ATAC output caused much higher RSS than the fragment-only low-memory
path because the current dual record carries both fragment evidence and two
SAM/BAM-side payloads before spilling.

## Decision

Implement and validate this in Chromap-suite first. Do not restart STAR-suite
multiome production on the BAM plus sidecar path until the Chromap-side harness
passes and STAR-suite has been rebuilt against the tested `libchromap.a`.

The Y/noY path is not the root cause of this issue. Y routing can add output
streams, but the observed memory growth came from the ATAC low-memory mapping
buffer holding large dual records before spill.

## Goals

- Keep one internal source of truth for ATAC mappings regardless of requested
  outputs.
- Decide at runtime which optional payload sections are present in the spill
  record.
- Sort and deduplicate using a compact, mandatory fragment-evidence prefix.
- Materialize requested outputs from that one record:
  - fragments / binary sidecar / peak evidence from the prefix,
  - BAM/CRAM from an optional BAM-pair payload.
- Make low-memory spilling byte-based so enabling BAM cannot silently multiply
  the real memory footprint.
- Test this in Chromap-suite independently before STAR-suite consumes it.

## Non-Goals

- Do not derive full BAM from the binary sidecar. The sidecar intentionally lacks
  CIGAR, sequence/quality, NM/MD, flags, and mate fields.
- Do not keep separate fragment-only and BAM-plus-fragment spill paths once the
  runtime schema lands.
- Do not change STAR-suite multiome production behavior before Chromap-side
  validation is complete.

## Current Failure Mode

Current dual output mode uses `PairedEndAtacDualMapping` to emit both:

- primary BAM/CRAM output, and
- secondary ATAC fragment output or binary sidecar.

That object contains the fragment evidence plus two `SAMMapping` payloads. The
low-memory threshold was calibrated around smaller fragment-style records. With
the dual object, a nominal low-memory run can hold many more bytes before
spilling.

Small one-batch smoke tests can also skip the mid-run spill branch and rely on
the final drain. The harness must explicitly test both:

- forced mid-run overflow spills, and
- final-drain-only runs.

## Proposed Record Model

Introduce one ATAC spill record with a fixed mandatory prefix and optional
sections selected by a per-run schema mask.

### Mandatory Prefix

The prefix must be sufficient for sort, dedup, fragment/sidecar output, and
peak-MEX construction:

```text
schema_version
rid
fragment_start
fragment_end or fragment_length
barcode_key
read_id
mapq
is_unique
direction_or_strand
num_dups
positive_alignment_length
negative_alignment_length
flags
```

Recommended flags:

```text
HAS_BAM_PAIR
HAS_Y_HIT
IS_BULK
```

The sort comparator and duplicate equivalence must use only the prefix fields
needed by the existing fragment semantics. Optional payloads must not change
fragment deduplication.

### Optional BAM Pair Section

Present only when BAM/CRAM output is requested.

```text
mate1:
  read_name
  sequence
  quality
  rid
  pos
  mate_rid
  mate_pos
  tlen
  flag
  mapq
  NM
  MD
  cigar

mate2:
  same fields
```

The section is variable-size because names, qualities, CIGAR, and MD tags are
variable-size. The spill record should therefore expose:

```text
encoded_size()
write(FILE*)
read(FILE*, schema_mask)
```

## Spill File Contract

Each overflow file should start with a run-level header:

```text
magic
version
schema_mask
record_codec_version
```

Each record should be length-prefixed:

```text
uint32_t rid
uint32_t byte_len
payload bytes
```

The `rid` wrapper keeps compatibility with the existing k-way merge shape. The
payload reader uses the file header's `schema_mask` to decode optional sections.

Prefer a per-run schema mask over arbitrary per-record TLV. It is simpler,
faster, and still supports variable-size fields inside the optional BAM section.

## Memory Accounting

Replace spill decisions based on record counts or `sizeof(MappingRecord)` with
byte accounting.

The buffer should track approximate bytes currently held in memory:

```text
current_buffer_bytes += record.memory_size_or_encoded_size()
flush when current_buffer_bytes >= low_mem_ram_limit
```

For vectors that hold records before encoding, include the real heap payload for
strings, qualities, CIGAR vectors, and MD tags. An encoded-size estimate is
acceptable if it is conservative.

Default thresholds should be expressed in bytes and should apply consistently
across:

- fragment-only output,
- sidecar output,
- BAM-only output,
- BAM plus sidecar output.

## Output Adapters

After sort/dedup, `AppendMapping` should become output-adapter driven:

```text
write_fragment_text(prefix)
write_binary_sidecar(prefix)
write_macs3_memory_record(prefix)
write_bam_pair(prefix, bam_pair_section)
```

Adapters must validate required payloads:

- sidecar/fragments require prefix only,
- BAM/CRAM requires `HAS_BAM_PAIR`,
- Y split for ATAC BAM should use `HAS_Y_HIT` or the existing read-id routing
  state if that remains the chosen implementation.

## Implementation Phases

1. Add `AtacSpillRecord` and schema-mask definitions.
2. Add serializer/deserializer tests for prefix-only and prefix-plus-BAM records.
3. Move dual ATAC low-memory output to `AtacSpillRecord`.
4. Convert fragment-only ATAC low-memory output to the same record.
5. Replace spill thresholds with byte accounting.
6. Preserve current output parity:
   - fragment-only output matches old fragment-only output,
   - dual fragments match fragment-only output,
   - dual BAM record count is two records per emitted fragment,
   - sorted BAM passes `samtools quickcheck` and coordinate-sort validation.
7. Rebuild STAR-suite against the tested Chromap archive and rerun the
   STAR-suite multiome smoke.

## Harness Plan

Add a Chromap-suite harness under `tests/`, with generated artifacts under
`plans/artifacts/` by default. The harness should accept:

```text
CHROMAP_ARTIFACT_ROOT
CHROMAP
THREADS
LOW_MEM_RAM
FIXTURE_ATAC
REF
INDEX
WHITELIST
```

Recommended script name:

```text
tests/run_atac_runtime_spill_schema_harness.sh
```

Required cases:

1. Fragment-only baseline, no low-mem.
2. Fragment-only, low-mem forced spill.
3. BAM plus sidecar, low-mem forced spill.
4. BAM plus sidecar, low-mem final-drain-only.
5. BAM plus sidecar, sorted BAM enabled.
6. If the synthetic Y fixture is practical for ATAC dual mode, BAM plus sidecar
   with Y/noY streams.

Required checks:

```text
fragment tuple multiset equality:
  fragment-only baseline == fragment-only low-mem
  fragment-only baseline == dual sidecar decoded to fragments

BAM checks:
  samtools quickcheck
  coordinate-sort check when --sort-bam is enabled
  BAM record count == 2 * emitted fragment count
  mate-pair validator passes

sidecar checks:
  binary header record_count > 0
  chroms TSV exists and matches reference sequence count/name order
  decoded record count matches emitted fragment count after dedup

memory checks:
  forced-spill stderr includes "Processing N overflow files"
  measured max RSS is recorded for each case
  dual BAM plus sidecar RSS stays within the configured low-memory envelope
```

Use `/usr/bin/time -v` when available and write a summary file:

```text
plans/artifacts/atac_runtime_spill_schema/<timestamp>/HARNESS_SUMMARY.txt
```

## Acceptance Criteria

Chromap-side acceptance:

- `make libchromap.a` succeeds.
- Existing low-memory BED smoke still passes.
- Existing ATAC dual output regression still passes.
- New runtime spill schema harness passes forced-spill and final-drain cases.
- No generated FASTQ, BAM, sidecar, index, or peak artifacts are tracked.

STAR-suite acceptance after Chromap integration:

- STAR-suite rebuilds cleanly with `WITH_CHROMAP=1`.
- The STAR-suite multiome lane smoke produces:
  - non-empty `atac_possorted.bam`,
  - non-empty `atac_fragments.bin`,
  - local ATAC peaks/MEX/metrics,
  - GEX GeneFull and Velocyto MEX,
  - Y-removal artifacts for the JAX KOLF2 production line.
- Only after the smoke passes should production be restarted.

