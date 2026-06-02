# CBQ Upstream-PR + Full-Coverage Handoff

Date: 2026-06-02

This is the entry point for finishing the native-CBQ work. Read this, then the
companion runbook
[`plans/2026-06-02-cbq-modality-coverage-and-upstream-pr-runbook.md`](2026-06-02-cbq-modality-coverage-and-upstream-pr-runbook.md)
for the detailed plan, and `AGENTS.md` before editing.

## TL;DR

Native CBQ input is implemented, hardened, and works as a FASTQ drop-in across
paired-end modalities on `master`. Four things remain:

1. **Vendor a self-contained CBQ test writer** (no STAR-suite, no `bqtools`).
2. **Extend the modality smoke** to all modalities (add Hi-C/PAIRS,
   `--emit-Y-noY-fastq`, read/barcode count-mismatch) using that writer.
3. **Adversarially audit the parallel-agent CBQ code** that landed after the
   ATAC milestone (direct-buffer decode, indexed range producer/consumer, lane
   index cache, Hi-C pairs) — it has not been reviewed.
4. **Cut and prepare an upstream-Chromap-only PR** from the minimal CBQ core.

## Repo state (as of 2026-06-02)

- Working repo: `/mnt/pikachu/Chromap-suite`, default branch `master` @ `3aedd17`.
- **This repo is under concurrent multi-agent development.** During the last
  session, branches `feature/cbq-{producer-consumer,index-cache,direct-sequence-buffers}`
  and a `/tmp/chromap-suite-atac-cbq-worker` worktree (`codex/atac-cbq-worker`)
  were committing and fast-forward-merging into `master`. Re-check `git log`/
  `git worktree list` before starting and **work on your own branch**.
- Deliverable branch from the last session: `feature/cbq-modality-matrix`
  (`a1f616c`) — adds the modality smoke + the runbook + Makefile target. Not
  pushed.

### What is already done (in `master`)

- CBQ reader `src/cbq_reader.{h,cc}`; loader branch in `MapPairedEndReads`;
  `SequenceBatch::{ResetLoadedSequences,AssignLoadedSequence}`; CLI
  (`--input-format`, `--read-pair-cbq`, `--barcode-cbq`) + validation in
  `chromap_driver.cc` and `chromap_lib_runner.cc`. (origin: branch
  `feature/cbq-atac-input`, commits `bcd9e72`, `66ed152`.)
- Security hardening of the reader (header/block bounds, overflow-safe
  offsets/ceil-div, decompressed- and compressed-size caps, remaining-file-bytes
  check, bounded reserve, BitVector bounds, headers-required guard for barcoded
  CBQ). Commits `5d65ab1`, `f6dad06`. **Verified still present in master.**
- Gates: `make test-cbq-atac-smoke` (synthetic ATAC BED), `make test-cbq-atac-100k`
  (100K ATAC BED), both CLI + `chromap_lib_runner`. STAR-suite
  `star-libchromap-contract` rebuilds and maps ATAC byte-identically against the
  updated `libchromap.a` (validated, not committed).
- Modality smoke `tests/run_cbq_modality_matrix.sh` + `make test-cbq-modality-matrix`
  (on `feature/cbq-modality-matrix`): CBQ vs FASTQ parity for BED, bulk BED,
  TagAlign, SAM, ChIP, BAM-dual (CLI + lib) + 5 rejection cases. Passes 11/11.
- Parallel work merged to master (NOT authored or reviewed here): Hi-C pairs
  parity (`65f1533`), direct-buffer decode (`00d0910`), indexed range producer
  (`40ca05a`), lane index cache (`3aedd17`).

### Verified current behavior (probed on master)

CBQ is a FASTQ drop-in with byte-identical (sorted) parity for scATAC
BED/BAM-dual/TagAlign/SAM, bulk ATAC, ChIP; Hi-C `--pairs` produces a valid
`.pairs` file; `--emit-Y-noY-fastq` emits Y/noY mate FASTQs; an unpaired
(1-mate) CBQ passed as `--read-pair-cbq` is rejected at open.

## Key facts and gotchas

### Upstream forks (for the PR)

- `/mnt/pikachu/chromap-orig` = `haowenz/chromap` (true upstream; last suite sync
  `73ec0c2`, 2021). Unmaintained; suite no longer tracks it.
- `/mnt/pikachu/chromap` = `morphic-bio/chromap` fork — **recommended PR target**.
- **Neither** has `src/libchromap.cc`, `src/chromap_lib_runner.cc`, or
  `third_party/libMACS3`. Those, plus the AEV1 sidecar, `--emit-Y-noY-fastq`,
  native BAM, and `AtacDualFragmentAndBam`, are **Chromap-suite-only**.
- Both **do** have `mapping_parameters.h`, `chromap_driver.cc`, `chromap.{h,cc}`,
  `sequence_batch.{h,cc}`, `mapping_writer.cc`, and the hook functions
  (`MapPairedEndReads`, `LoadPairedEndReadsWithBarcodes`,
  `SampleInputBarcodesAndExamineLength`, `ComputeBarcodeAbundance`,
  `ReplaceByEffectiveRange`).

### Encoder situation (the crux of testability)

- `cbq_ordered_encoder` lives in STAR-suite (`/mnt/pikachu/STAR-suite/core/legacy/source/`)
  and is **ephemeral** — it disappeared mid-session after a STAR rebuild. Do not
  depend on it.
- `bqtools` (`/tmp/star_suite_bqtools/bin/bqtools`) is external and **reorders
  records across blocks at scale**, breaking read/barcode lane alignment for
  multi-block files. Only safe for tiny single-block fixtures.
- The CBQ v1 block layout the reader expects is fully specified by the parser in
  `src/cbq_reader.cc` (`ParseFileHeader` 64-byte header; `ParseBlockHeader`
  96-byte block header; 7 zstd columns: `z_seq_len, z_header_len, z_npos, z_seq,
  z_flags, z_headers, z_qual`; uncompressed columns are accepted — `Decompress`
  validates size, and the reader handles uncompressed input). A vendored writer
  emitting **uncompressed** columns is sufficient and removes the libzstd
  runtime requirement from the test path.

### Alignment contract

Barcoded CBQ requires the read-pair and barcode lanes to be **record-aligned**
(record *i* = same original read) and to **carry headers** (the reader enforces a
per-record name match and rejects headerless barcoded input).

### Runtime dependency

The reader `dlopen`s `libzstd.so.1`/`libzstd.so`; the only build change needed
is `-ldl`. Uncompressed CBQ needs no libzstd at run time.

### Test fixtures

- 100K PBMC ATAC: `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac`
  (`*_R1`=read1, `*_R3`=read2, `*_R2`=barcode).
- Index `…/pbmc_unsorted_3k_100k/chromap_index/genome.index`; whitelist
  `…/chromap_index/737K-arc-v1_atac.txt`; reference
  `/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa`.
- Keep all generated artifacts under `plans/artifacts/` (git-ignored).

## Open work — acceptance criteria

1. **Vendored CBQ writer** (`tests/` ): emits CBQ v1 (uncompressed columns OK)
   matching `cbq_reader.cc`; round-trips through `chromap --input-format cbq`
   with FASTQ parity; no STAR-suite / `bqtools` dependency.
2. **Modality smoke extended**: `run_cbq_modality_matrix.sh` uses the vendored
   writer (drops the external-encoder skip), and adds positive cases for Hi-C
   `--pairs` and `--emit-Y-noY-fastq`, plus the read/barcode count-mismatch
   rejection. All CBQ==FASTQ on both front-ends; all rejections exit non-zero
   with the documented message and no crash.
3. **Audit of parallel CBQ code**: adversarial correctness + memory-safety review
   of the diffs in `00d0910`, `40ca05a`, `3aedd17`, `65f1533` (direct-buffer
   decode, indexed range producer/consumer, lane index cache, Hi-C pairs).
   Findings fixed or filed; the existing reader hardening must remain intact.
4. **Upstream PR branch**: cut from the minimal CBQ core (NOT current master —
   see runbook "Recommended baseline"), build against `/mnt/pikachu/chromap`
   with only `-ldl` added, CLI-only (no libchromap/lib_runner/sidecar/MACS3),
   hardening preserved, hermetic vendored test included, README/docs cover the
   options + alignment contract. Builds clean; CBQ==FASTQ for the paired-end
   subset the fork supports.

Existing gates (`test-cbq-atac-smoke`, `test-cbq-atac-100k`,
`test-cbq-modality-matrix`, `test-libchromap-core-smoke`) must still pass.

## Constraints

- Do not depend on STAR-suite or a floating `bqtools` in released/PR code or its
  tests. STAR-suite is downstream only.
- Work on a feature branch; do not commit to `master`. No `Co-authored-by` /
  generated-by trailers (see `AGENTS.md`).
- Benchmarks serial; record command line, git state, fixtures, threads, output
  dir, and wall/user/sys/max-RSS for any timing run.

---

## Suggested Agent Prompt

> Finish the native-CBQ work in `/mnt/pikachu/Chromap-suite`. Read
> `plans/2026-06-02-cbq-upstream-pr-and-coverage-handoff.md`, the companion
> runbook `plans/2026-06-02-cbq-modality-coverage-and-upstream-pr-runbook.md`,
> and `AGENTS.md` first. The repo is under concurrent multi-agent development —
> check `git log`/`git worktree list`, branch off `master`, and never commit to
> `master`.
>
> Deliver, each verified and committed on your branch:
>
> 1. **A vendored, self-contained CBQ v1 writer under `tests/`** (Python or a
>    tiny C++ tool) that emits the exact block layout `src/cbq_reader.cc`
>    parses. Uncompressed columns are fine (the reader accepts them; this also
>    avoids the libzstd test dependency). Prove it: encode a small FASTQ pair +
>    barcode, run `chromap --input-format cbq`, and show byte-identical (sorted)
>    fragments vs the FASTQ baseline. It must NOT depend on STAR-suite or
>    `bqtools`.
>
> 2. **Extend `tests/run_cbq_modality_matrix.sh`** to use the vendored writer
>    (remove the external-encoder skip) and add committed positive cases for
>    Hi-C `--pairs` and `--emit-Y-noY-fastq`, plus a read/barcode count-mismatch
>    rejection case. Every positive case compares CBQ vs FASTQ on both `chromap`
>    and `chromap_lib_runner`; every rejection asserts non-zero exit, the
>    documented message, and no crash/OOM. Keep it hermetic.
>
> 3. **Adversarially audit the parallel-agent CBQ code** added after the ATAC
>    milestone — the diffs in commits `00d0910` (direct-buffer decode),
>    `40ca05a` (indexed range producer), `3aedd17` (lane index cache), and
>    `65f1533` (Hi-C pairs). Treat CBQ bytes as untrusted. Focus on bounds,
>    integer overflow, lifetime/aliasing of the direct `SequenceBatch` buffers,
>    range/offset math in the indexed producer/consumer, and cache invalidation.
>    Fix real findings (or file them with repro) and confirm the existing reader
>    hardening (the guards in `src/cbq_reader.cc`) is intact.
>
> 4. **Cut an upstream-Chromap-only PR branch.** Base it on the minimal CBQ
>    core, not current master (see the runbook "Recommended baseline" and
>    "File-by-file scope"): `cbq_reader.{h,cc}` with hardening; additive
>    `MappingParameters` fields/helpers (no `AtacDualFragmentAndBam`); the
>    `chromap.{h,cc}` loader + barcode-sampling/abundance branches;
>    `sequence_batch.{h,cc}` `Reset/AssignLoadedSequence`; `chromap_driver.cc`
>    CLI + validation; Makefile object + `-ldl`. Exclude everything that does
>    not exist upstream (`libchromap.cc`, `chromap_lib_runner.cc`, BAM
>    read-group `ReadGroupSourcePath`, AEV1 sidecar, MACS3, Y/noY, low-mem
>    interplay). Verify it builds against `/mnt/pikachu/chromap`
>    (`morphic-bio/chromap`) with only `-ldl` added and that CBQ==FASTQ for the
>    paired-end output formats that fork supports, using the vendored hermetic
>    test. Add a README option block documenting `--input-format cbq` /
>    `--read-pair-cbq` / `--barcode-cbq`, the lane-alignment + headers contract,
>    and the libzstd-at-runtime note.
>
> Keep `make test-cbq-atac-smoke`, `test-cbq-atac-100k`,
> `test-cbq-modality-matrix`, and `test-libchromap-core-smoke` green. Put
> artifacts under `plans/artifacts/`. Do not depend on STAR-suite. Report what
> you did, what passed, and anything you deferred.
