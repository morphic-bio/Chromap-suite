# Chromap-suite CBQ ATAC Integration Handoff

Date: 2026-05-29

## Goal

Add production CBQ input support to Chromap-suite for ATAC/scATAC without
materializing intermediate FASTQs. The existing FASTQ path should remain
unchanged. CBQ should enter as an in-memory reader/view that fills the same
Chromap mapping surface used downstream today.

Scope for this handoff is ATAC only:

- paired-end ATAC reads,
- optional single barcode read for scATAC,
- `chromap` CLI and `libchromap` runner parity,
- STAR-suite multiome consumption after Chromap-suite support is validated.

Do not expand this to ChIP, Hi-C, or other modalities in the first pass.

## Repo And Branch

Work in Chromap-suite, preferably from a clean worktree because the current
checkout has local dirty state.

```bash
git -C /mnt/pikachu/Chromap-suite status --short --branch
git -C /mnt/pikachu/Chromap-suite worktree add /tmp/chromap_cbq_atac origin/master
cd /tmp/chromap_cbq_atac
git switch -c feature/cbq-atac-input
```

Read `/mnt/pikachu/Chromap-suite/AGENTS.md` before editing. Do not revert
unrelated dirty files in `/mnt/pikachu/Chromap-suite`.

## Existing Reference Work

STAR-suite already has the CBQ pieces needed as a reference:

- `/mnt/pikachu/STAR-suite/core/legacy/source/input/CbqInputModule.h`
- `/mnt/pikachu/STAR-suite/core/legacy/source/input/CbqInputModule.cpp`
- `/mnt/pikachu/STAR-suite/core/legacy/source/input/InputContract.h`
- `/mnt/pikachu/STAR-suite/core/legacy/source/input/CbqChromapAdapter.h`
- `/mnt/pikachu/STAR-suite/core/legacy/source/input/CbqChromapAdapter.cpp`
- `/mnt/pikachu/STAR-suite/tests/run_cbq_chromap_adapter_smoke.sh`

The existing `CbqChromapAdapter` is correctness proof only: it writes temporary
FASTQs and then calls Chromap. Do not use that as the production integration
surface. Production means CBQ reader to in-memory Chromap consumer, with no
intermediate FASTQ files.

## Chromap Integration Points

Primary files in Chromap-suite:

- `src/mapping_parameters.h`: add explicit CBQ input fields.
- `src/chromap_driver.cc`: add CLI parsing and validation.
- `src/chromap_lib_runner.cc`: expose the same options for libchromap parity.
- `src/sequence_batch.h` and `src/sequence_batch.cc`: add a controlled way to
  fill `SequenceBatch` records from a non-kseq source.
- `src/chromap.h`: branch the ATAC read-loading loop when CBQ input is active.
- `src/chromap.cc`: barcode sampling/correction helpers may need a CBQ branch.
- `Makefile`: compile the CBQ reader and add smoke targets.
- `tests/`: add CBQ input parity smokes.

Current Chromap loads reads through `SequenceBatch` and `kseq`:

- `SequenceBatch::InitializeLoading()`
- `SequenceBatch::LoadBatch()`
- `Chromap::MapPairedEndReads()`
- `Chromap::LoadPairedEndReadsWithBarcodes()`

The low-risk approach is to keep all downstream mapping, barcode correction,
fragment writing, BAM writing, duplicate logic, and peak sidecars unchanged.
Only replace the file loading step.

## Proposed API Shape

Avoid overloading the existing FASTQ vectors with ambiguous semantics. Add an
explicit input mode and explicit CBQ path lists.

```cpp
enum class ReadInputFormat { kFastq, kCbq };

struct MappingParameters {
  ReadInputFormat read_input_format = ReadInputFormat::kFastq;
  std::vector<std::string> read_pair_cbq_paths;
  std::vector<std::string> barcode_cbq_paths;
};
```

Suggested CLI:

```bash
chromap --preset atac \
  --input-format cbq \
  --read-pair-cbq lane1.reads.cbq,lane2.reads.cbq \
  --barcode-cbq lane1.barcode.cbq,lane2.barcode.cbq \
  --barcode-whitelist whitelist.txt \
  ...
```

Keep `-1/-2/-b` as FASTQ fields. If a later UI wants aliases, add them only
after the explicit surface is working and tested.

Validation rules:

- `--input-format fastq` is the default and requires the existing `-1/-2`
  path behavior.
- `--input-format cbq` requires `--read-pair-cbq`.
- `--barcode-cbq` count must match `--read-pair-cbq` count when barcode input
  is used.
- `--barcode-whitelist` with CBQ still requires barcode CBQ input.
- Reject mixed FASTQ and CBQ inputs in one command.
- Initially reject `--emit-Y-noY-fastq` in CBQ mode, or require explicit output
  prefixes and test it separately. The first production milestone should not
  depend on Y/noY FASTQ sidecar behavior.

## Reader Strategy

Port a narrow C++ CBQ reader into Chromap-suite. Keep it independent from
STAR-specific namespaces and build rules. It should expose one lane reader for:

- paired read CBQ: two mates, emitted as read1/read2,
- barcode CBQ: one mate, emitted as barcode.

For phase 1, decode into `SequenceBatch` slots. This is still production
because it is in-memory and goes straight into the Chromap consumer. It avoids
touching minimizer, barcode correction, writer, and duplicate code.

Add methods to `SequenceBatch` instead of mutating its internals from outside,
for example:

```cpp
void ResetLoadedSequences();
void AssignLoadedSequence(uint32_t index,
                          const char* name,
                          size_t name_len,
                          const char* comment,
                          size_t comment_len,
                          const char* seq,
                          size_t seq_len,
                          const char* qual,
                          size_t qual_len,
                          uint32_t source_id);
```

`AssignLoadedSequence()` should apply the existing effective range behavior so
CBQ and FASTQ reach the same downstream state.

Later optimization can avoid ASCII sequence materialization in minimizer
generation, but do not do that in the first pass. The first target is exact
parity with a small code surface.

## Mapping Loop Shape

In `Chromap::MapPairedEndReads()`:

1. Keep reference and index loading exactly as-is.
2. Keep barcode whitelist loading and all downstream objects exactly as-is.
3. At the per-file loop, branch only at the loader:
   - FASTQ mode: current `InitializeLoading()` and `LoadPairedEndReadsWithBarcodes()`.
   - CBQ mode: open paired read CBQ plus optional barcode CBQ for the lane and
     fill `read_batch1_for_loading`, `read_batch2_for_loading`, and
     `barcode_batch_for_loading`.
4. After batches are filled, keep the same `SwapSequenceBatch()` and mapping
   path.

Barcode pre-sampling is a separate required check. `SampleInputBarcodesAndExamineLength()`
currently reads barcode FASTQ paths; add a CBQ-aware branch or helper so scATAC
CBQ mode computes barcode length and abundance from the barcode CBQ.

## Test Plan

Add a small synthetic smoke first:

- Generate synthetic FASTA, paired ATAC FASTQs, barcode FASTQ, and whitelist.
- Encode `reads_R1 + reads_R2` to one paired CBQ.
- Encode barcode FASTQ to one single-end CBQ.
- Run FASTQ baseline through `chromap`.
- Run CBQ mode through `chromap`.
- Compare sorted fragment/BED rows and summary metrics.
- Run the same comparison through `chromap_lib_runner`.

Then add a 100K fixture smoke using the existing PBMC ATAC fixture:

- fixture root: `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac`
- current scripts to mirror:
  - `tests/run_chromap_lowmem_bed_smoke_100k.sh`
  - `tests/run_atac_dual_output_100k.sh`
  - `tests/run_atac_runtime_spill_schema_harness.sh`

The CBQ smoke can use the STAR-suite ordered encoder during development:

```bash
CBQ_ORDERED_ENCODER_BIN=/mnt/pikachu/STAR-suite/core/legacy/source/cbq_ordered_encoder
```

Do not make the released Chromap-suite depend on STAR-suite. Either vendor a
small test encoder, import the needed reader/encoder code, or make the external
encoder an optional test-only dependency.

## Acceptance Criteria

- `make chromap chromap_lib_runner` succeeds.
- Existing FASTQ tests still pass.
- Synthetic FASTQ vs CBQ ATAC smoke passes for CLI and libchromap.
- 100K PBMC ATAC FASTQ vs CBQ smoke passes for at least barcoded paired-end
  ATAC fragments.
- No production CBQ path writes temporary FASTQs.
- Threaded CBQ output is deterministic relative to FASTQ for the comparison
  surface, or the runbook documents the exact canonicalization used.
- STAR-suite can call the updated `libchromap.a` for ATAC after the Chromap
  smoke passes.

## Suggested Agent Prompt

Implement ATAC-only CBQ input in Chromap-suite. Keep FASTQ paths unchanged.
Use a clean worktree and do not revert unrelated dirty files. Add explicit CBQ
input fields to `MappingParameters`, CLI/lib runner options, a narrow CBQ reader
that fills `SequenceBatch` in memory, and a loader branch in paired-end ATAC
mapping. Do not materialize FASTQ files in the production path. Add synthetic
and 100K fixture parity smokes comparing FASTQ vs CBQ outputs. Start from the
STAR-suite CBQ reader and `tests/run_cbq_chromap_adapter_smoke.sh` as references,
but replace the temporary-FASTQ adapter with an in-memory reader.
