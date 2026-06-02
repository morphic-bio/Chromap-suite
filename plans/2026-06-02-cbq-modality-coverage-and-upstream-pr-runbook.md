# CBQ Modality Coverage + Upstream-Chromap PR Runbook

Date: 2026-06-02

## Purpose

Two related goals:

1. **Modality coverage** — native CBQ input now works far beyond the original
   ATAC-only milestone, but the committed smoke tests only exercise ATAC BED.
   Define and add a full-modality CBQ parity smoke so every supported output is
   gated, and every unsupported combination has a verified rejection.
2. **Upstream-scoped PR** — produce a version of the CBQ feature that depends
   only on surfaces present in upstream Chromap (`haowenz/chromap`) or the
   `morphic-bio/chromap` fork, so it can be offered as a clean PR without
   dragging in Chromap-suite-only machinery (libchromap, libMACS3, the AEV1
   sidecar, Y/noY, the low-mem spill rewrite).

## Current State (verified 2026-06-02 on `master` @ `3aedd17`)

`feature/cbq-atac-input` (the ATAC milestone + the security hardening) was
merged into `master` (`57e6f9a`). Parallel CBQ work then landed on top:

- `65f1533` Add CBQ Hi-C pairs parity support (removed the old PAIRS rejection)
- `00d0910` Decode CBQ directly into SequenceBatch buffers
- `40ca05a` Add indexed CBQ range producer
- `3aedd17` Cache CBQ lane indexes for range reads

Net effect: CBQ is now effectively a **drop-in replacement for FASTQ across all
paired-end modalities**, not just ATAC. The earlier "ATAC-only" framing and the
PAIRS / `--emit-Y-noY-fastq` rejections are obsolete.

### What was probed (single lane L001 of the 100K PBMC ATAC fixture, current master binary)

| Modality / output | CBQ result vs FASTQ baseline |
|---|---|
| scATAC BED + barcode | sorted rows byte-identical |
| scATAC BAM + `--sort-bam` + `--atac-fragments` (dual) | fragments byte-identical |
| scATAC TagAlign + barcode | 155,640 rows, byte-identical |
| scATAC SAM + barcode | alignment rows byte-identical |
| bulk ATAC BED (no barcode) | 79,030 rows, byte-identical |
| ChIP (`--preset chip`) paired BED | 77,868 rows, byte-identical |
| Hi-C (`--split-alignment --pairs`) | produces a valid `.pairs` file |
| `--emit-Y-noY-fastq` (with CBQ) | emits Y / noY mate FASTQs |
| unpaired CBQ as `--read-pair-cbq` | rejected at open: "CBQ mate-count mismatch … file paired=false but requested mate_count=2" |

These were ad-hoc probes, not committed tests. The only committed CBQ gates are
`make test-cbq-atac-smoke` (synthetic ATAC BED, CLI + libchromap) and
`make test-cbq-atac-100k` (100K ATAC BED, CLI + libchromap).

### Coverage gap

No committed smoke covers, per front-end (CLI `chromap` **and** `chromap_lib_runner`):

- BAM / sorted-BAM / dual fragments
- TagAlign, SAM
- bulk (no-barcode) paired-end
- ChIP preset
- Hi-C / PAIRS
- `--emit-Y-noY-fastq` under CBQ
- the **rejection** paths (unpaired CBQ; barcoded CBQ with headers stripped;
  read/barcode count mismatch; FASTQ+CBQ mixed inputs).

## Part A — Full-Modality CBQ Smoke Matrix

### Design principles

- **Hermetic.** Reuse the synthetic-fixture pattern from
  `tests/run_cbq_atac_smoke.sh`: a tiny generated genome + paired FASTQs +
  barcode FASTQ + whitelist, a small `--build-index`, and CBQ files encoded
  from those FASTQs by the vendored C++ writer under `tests/`. No STAR-suite,
  no `bqtools`, no 12 GB index, and no external fixture.
- **Canonical comparison.** `LC_ALL=C sort` over text rows for BED / BEDPE /
  TagAlign / pairs / fragments; for BAM, compare `samtools view` records (or
  fragment TSV) and `samtools quickcheck`; for SAM, compare non-`@` rows.
- **Both front-ends.** Every positive case runs through `chromap` and
  `chromap_lib_runner` and is compared to the FASTQ baseline produced by
  `chromap`.
- **Rejections are tests too.** Each unsupported/guarded combination asserts a
  non-zero exit and an expected message, and (for the oversized/garbage cases)
  that the process does not crash or OOM.

### Positive cases (CBQ == FASTQ, both front-ends)

| Case | Key flags |
|---|---|
| `pe_bc_bed` | `--preset atac --BED -b/--barcode-cbq --barcode-whitelist` |
| `pe_bulk_bed` | paired, no barcode, `--BED` |
| `pe_bc_tagalign` | `--TagAlign` + barcode |
| `pe_bc_sam` | `--SAM` + barcode (compare non-`@` rows) |
| `pe_bc_bam` | `--BAM --sort-bam` + barcode (samtools view rows + quickcheck) |
| `pe_bc_dual` | `--BAM --sort-bam --atac-fragments frags.tsv.gz` (compare fragments) |
| `pe_bulk_bam` | paired, no barcode, `--BAM --sort-bam` (samtools view rows + quickcheck) |
| `chip_bam` | `--preset chip --BAM --sort-bam`, no barcode (samtools view rows + quickcheck) |
| `pe_bulk_cram` | paired, no barcode, `--CRAM --sort-bam` (samtools view rows + quickcheck) |
| `read_group_auto` | `--BAM --read-group auto`; normalize RG tags for row parity and assert source-derived IDs |
| `chip_bed` | `--preset chip --BED`, no barcode |
| `hic_pairs` | `--split-alignment --pairs` (compare sorted pairs rows) |
| `emit_y_noy_fastq` | `--emit-Y-noY-fastq` (compare SAM rows and noY FASTQ payloads) |

### Rejection / safety cases (assert non-zero exit + message, no crash)

| Case | Expectation |
|---|---|
| `--input-format cbq` without `--read-pair-cbq` | "requires --read-pair-cbq" |
| `--barcode-cbq` count != `--read-pair-cbq` count | "count must match" |
| `--barcode-whitelist` without `--barcode-cbq` | "requires --barcode-cbq" |
| FASTQ (`-1/-2/-b`) mixed with CBQ | "cannot be mixed" |
| unpaired (1-mate) CBQ as `--read-pair-cbq` | "mate-count mismatch" |
| barcoded CBQ with headers stripped (`tests/cbq_ordered_encoder --strip-headers`) | "requires read names (headers)" |
| read lane longer than barcode lane and vice-versa | "Numbers of reads and barcodes don't match!" |
| garbage / truncated / oversized-dimension CBQ | clean error, no crash, RSS bounded |

### Deliverable

`tests/run_cbq_modality_matrix.sh` + `make test-cbq-modality-matrix`. Builds and
uses `tests/cbq_ordered_encoder` by default; `CBQ_ORDERED_ENCODER` is only for
intentional local overrides. Artifacts under
`plans/artifacts/cbq_modality_matrix/<timestamp>/`.

Wire it into the merge checklist next to `test-cbq-atac-smoke`.

## Part B — Upstream-Chromap-Only PR

### Goal

A PR that adds native CBQ input touching **only files and symbols that exist in
upstream Chromap**, so it applies cleanly to `morphic-bio/chromap` (and largely
to `haowenz/chromap`) without any Chromap-suite-only dependency.

### Fork facts (verified)

- `/mnt/pikachu/chromap-orig` = `haowenz/chromap` (true upstream; last suite sync
  was `73ec0c2`, 2021-10-30 per `HISTORY.md`).
- `/mnt/pikachu/chromap` = `morphic-bio/chromap` fork.
- **Neither** has `src/libchromap.cc`, `src/chromap_lib_runner.cc`, or
  `third_party/libMACS3`. Those are Chromap-suite-only.
- Both **do** have `mapping_parameters.h`, `chromap_driver.cc`, `chromap.h`,
  `chromap.cc`, `sequence_batch.{h,cc}`, `mapping_writer.cc`, and the core
  functions the CBQ branch hooks (`MapPairedEndReads`,
  `LoadPairedEndReadsWithBarcodes`, `SampleInputBarcodesAndExamineLength`,
  `ComputeBarcodeAbundance`, `ReplaceByEffectiveRange`).
- `AtacDualFragmentAndBam` does **not** exist upstream (suite-only helper).

### File-by-file scope

**Portable — include in the upstream PR:**

- `src/cbq_reader.h`, `src/cbq_reader.cc` — new, self-contained. Includes the
  full security hardening (header/block bounds, overflow-safe offset/ceil-div,
  decompressed- and compressed-size caps, remaining-file-bytes check, bounded
  reserve, BitVector bounds). **Must travel with the PR.**
- `src/mapping_parameters.h` — additive only: `enum class ReadInputFormat`, the
  `read_pair_cbq_paths` / `barcode_cbq_paths` fields, and the
  `UsesCbqInput / HasPairedEndInput / HasBarcodeInput / NumInputLanes /
  ReadGroupSourcePath` helpers. Do **not** port the `AtacDualFragmentAndBam`
  edit (suite-only).
- `src/chromap_driver.cc` — `--input-format`, `--read-pair-cbq`,
  `--barcode-cbq` options + validation, the mapping-branch gate change, and the
  run-summary print. Drop any branch that references suite-only output (e.g.
  the AEV1 sidecar / MACS3 frag-peaks block).
- `src/chromap.h`, `src/chromap.cc` — `LoadPairedEndReadsFromCbq`,
  `LoadBarcodesFromCbq`, the CBQ branch in `MapPairedEndReads`, and the CBQ
  branches in `SampleInputBarcodesAndExamineLength` / `ComputeBarcodeAbundance`.
  The `NumInputLanes()` loop substitution is portable.
- `src/sequence_batch.h`, `src/sequence_batch.cc` — `ResetLoadedSequences`,
  `AssignLoadedSequence`, `AssignKstring`.
- `Makefile` — add `cbq_reader.cc` to the object list and `-ldl` to the link
  line (for the `dlopen` of libzstd).

**Exclude — Chromap-suite-only:**

- `src/libchromap.cc`, `src/chromap_lib_runner.cc` — do not exist upstream. The
  upstream PR is **CLI-only**; libchromap parity is a Chromap-suite concern.
- `src/mapping_writer.cc` `ReadGroupSourcePath()` change — only needed for BAM
  read-group auto-naming, and native BAM is a suite addition (per `HISTORY.md`).
  Port only if the target fork already has BAM RG auto-naming.
- Hi-C/PAIRS CBQ, `--emit-Y-noY-fastq` CBQ, AEV1 sidecar, MACS3 frag-peaks,
  low-mem spill interplay — out of scope for the first upstream PR. Keep the
  upstream PR to paired-end BED/BEDPE/TagAlign/SAM (whatever the target fork
  supports), matching the original minimal milestone.

### Recommended baseline

Base the upstream PR on the **minimal ATAC milestone** (`feature/cbq-atac-input`
content: `bcd9e72` + the two hardening commits `5d65ab1`, `f6dad06`), **not**
current master. Master's CBQ has since grown suite-specific tendrils (Hi-C
parity, Y/noY emission, direct-buffer decode against suite-specific batch
internals, indexed range producer/consumer, index caching) that would bloat and
entangle an upstream PR. Cherry-pick the minimal core, then re-apply hardening.

### Runtime dependency

The reader `dlopen`s `libzstd.so.1` / `libzstd.so` at run time (no link-time
dep beyond `-ldl`). Document that compressed CBQ requires libzstd present at run
time. This keeps the upstream build change to a single `-ldl`.

### Test self-containment (the real blocker)

The original smokes generated CBQ with **external** tools:

- `cbq_ordered_encoder` lived in STAR-suite and was **ephemeral** — it
  disappeared from `core/legacy/source/` mid-session after a STAR rebuild.
- `bqtools` is a separate external binary and **reorders records across blocks
  at scale**, so it only works for tiny single-block fixtures.

An upstream PR cannot depend on STAR-suite or a floating `bqtools`. Vendor the
order-preserving C++ CBQ v1 writer source under `tests/` and build it as part of
the smoke target. This is the only acceptable path for a gated upstream checkout:
the test writer is local source, emits the exact block layout `cbq_reader.cc`
parses, appends the CBQINDEX footer used by the indexed producer, and preserves
read/barcode lane order across blocks.

### Alignment contract (carry into PR docs)

CBQ barcoded input requires the read-pair and barcode lanes to be
**record-aligned** (record *i* = same original read) and to **carry headers**
so the per-record name match can verify alignment. The reader enforces the
header requirement and a per-record name check; document that inputs must be
produced by an order-preserving encoder.

### Suggested PR shape

1. Commit 1: `cbq_reader.{h,cc}` (with hardening) + Makefile `-ldl` + object.
2. Commit 2: `SequenceBatch::{ResetLoadedSequences,AssignLoadedSequence}`.
3. Commit 3: `MappingParameters` CBQ fields/helpers (additive).
4. Commit 4: `chromap.{h,cc}` CBQ loader branch + barcode sampling/abundance.
5. Commit 5: `chromap_driver.cc` CLI + validation.
6. Commit 6: vendored test writer + hermetic modality smoke (paired-end subset
   the target fork supports) + docs (README option block + alignment contract).

Target: `morphic-bio/chromap` (active fork; `haowenz/chromap` is unmaintained
and Chromap-suite no longer tracks it). Keep single-end, Hi-C, BAM/sidecar,
Y/noY out of the first PR.

## Acceptance / Checklist

- [ ] `tests/run_cbq_modality_matrix.sh` added; `make test-cbq-modality-matrix`
      passes (all positive cases CBQ==FASTQ on both front-ends; all rejection
      cases fail with the expected message and no crash).
- [ ] Modality matrix wired into the merge checklist.
- [ ] Upstream branch cut from the minimal CBQ core (no libchromap /
      lib_runner / sidecar / MACS3), builds against `morphic-bio/chromap` with
      only `-ldl` added.
- [ ] Vendored hermetic CBQ test writer (no STAR / bqtools dependency).
- [ ] Security hardening present in the upstream `cbq_reader.cc`.
- [ ] README/docs in the PR cover the CBQ options and the lane-alignment +
      header contract.
