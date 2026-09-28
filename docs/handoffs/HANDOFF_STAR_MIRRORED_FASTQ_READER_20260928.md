# Handoff: STAR Suite FASTQ reader mirrored into Chromap Suite (2026-09-28)

Runbook: `docs/runbooks/RUNBOOK_STAR_MIRRORED_FASTQ_READER_20260928.md`.
Validation root: `/mnt/pikachu/fastq_mirror_validation_20260928` (below `V`).

## State

Branch `feat/star-mirrored-fastq-reader` in
`/mnt/pikachu/Chromap-suite-fastq-mirror-20260928`, on top of `09b3164`:

| Commit | Content |
|---|---|
| `baeb3ef` | Mirror (`src/star_input/`), glue and BGZF provider (`src/fastq_bgzf_input.{h,cc}`), options, gzip threads per batch, harness, intake smoke |
| `66d4537` | `src/version.h` -> 1.1.0 (code validated at this commit) |
| `b4bcc51` | README, tests/README, CHANGELOG, `docs/RELEASE_NOTES_v1.1.0.md` draft (throughput pending) |
| this commit | runbook and this handoff |

Done:

- Mirror: `BgzfBlockReader.{h,cpp}`, `BgzfRangeReader.{h,cpp}` from STAR
  `6e83853` (identical at STAR `v1.9.5.a`); namespace/include/guard edits only
  (`src/star_input/MIRROR.md`, `src/star_input/LICENSE`).
- Provider: each file inflated by its own workers and consumed on its own
  thread; pairing by ordinal and read-name stem, as STAR Suite does.
  `--input-bgzf-mode auto|on|off`, `--input-bgzf-reader-threads N`;
  `MappingParameters::input_bgzf_mode` / `input_bgzf_reader_threads`.
- Plain gzip: thread per file for every batch (was: first batch only;
  later batches OpenMP tasks from 12 threads).
- Version 1.1.0; release notes and CHANGELOG drafted.
- Tests (all exit 0 at `66d4537`): `test-fastq-intake-smoke` (27 cases),
  `test-libchromap-core-smoke`, `test-atac-sidecar-only-smoke`,
  `test-input-format-smoke`, `test-cbq-atac-smoke`; at the pre-fold commit
  (same zlib code, lockstep provider) also `test-smoke`,
  `test-cbq-modality-matrix` (25/0), `test-cbq-range-reader`,
  `test-atac-runtime-spill-schema-harness`.
- Record identity (`V/records/*.tsv`): one record-dump sha per dataset across
  kseq serial, kseq threads, kseq on BGZF, the BGZF provider (lockstep and
  per-file threads) and `auto`: DOGMA 4M (12,000,000 records,
  `cf975355...`), DOGMA 50M (150,000,000 records incl. the plain extract,
  `82da1df7...`), PBMC 100k four lanes, PBMC full depth four lanes native BGZF
  (248,344,722 records = 82,781,574 pairs).
- Standalone output identity (`V/identity/identity.json`,
  `V/identity_det/identity.json`), old `09b3164` vs new `66d4537`, 16 threads,
  PBMC 100k gz/bgz, DOGMA 4M gz/bgz, PBMC full bgz, BED and sidecar-only:
  fragments and sidecar+chroms identical in every mode; `--summary` differs
  only in cache columns (also old vs old); with `--deterministic-mapping`
  summaries identical. Superseded runs of the pre-fold binary `a021ce7` are in
  `V/identity_superseded_a021ce7/` (same results).
- Embedded (`V/embedded/`): STAR `6e83853` built unchanged against
  `66d4537` (STAR sha `5e65d2e5...`, build record
  `/mnt/pikachu/STAR-suite-readercheck-20260928/multiomics_build.json`,
  logs `V/embedded/build_logs/`). DOGMA-plex lane 1, 2M reads, 16 threads,
  sidecar: `atac_fragments.bin` (1,587,839 records), `.chroms.tsv` and all
  six `atac/` files byte-identical to the no-BAM binary
  (`V/embedded/compare_readercheck_66d4537.json`, `VERDICT: PASS`).

In progress:

- Throughput campaign, PID `1009455`:
  `cd V/throughput && ./run_throughput.sh pass2 harness_bgzf_bgz harness_kseq_threads_gz harness_kseq_serial_gz chromap_new_bgz chromap_old_bgz chromap_new_gz chromap_old_gz`
  (log `V/throughput/campaign.log`, runs in `V/throughput/pass{1,2}/<case>/`).
  It waits for a quiet host before each case. Remaining in pass2:
  `chromap_old_bgz`, `chromap_new_gz`, `chromap_old_gz`.

Clean results so far (`python3 V/throughput/summarize.py` ->
`throughput.tsv`), 50,000,000 lane-1 pairs, 32 threads, inputs evicted:

| Case | Wall (s) | Verdict |
|---|---|---|
| chromap_old_gz (pass1) | 145.2 | clean |
| chromap_new_bgz (pass1, pass2) | 122.8, 122.3 | clean |
| harness_kseq_serial_gz (pass2) | 82.0 | clean |
| harness_bgzf_bgz (pass2) | 12.5 | clean |
| chromap_new_gz, chromap_old_bgz, harness_kseq_threads_gz | pass1/2 contaminated | rerun |

## Next

1. Let pass2 finish (`tail V/throughput/campaign.log`).
2. Rerun every case without a clean verdict under a new pass label, e.g.
   `cd V/throughput && nohup ./run_throughput.sh pass3 chromap_new_gz chromap_old_bgz harness_kseq_threads_gz > campaign3.log 2>&1 &`.
3. `python3 V/throughput/summarize.py`; after the campaign (not during it)
   check the 50M sidecars old vs new are identical
   (`sha256sum V/throughput/pass*/chromap_*/atac_fragments.bin`).
4. Put the numbers in `docs/RELEASE_NOTES_v1.1.0.md` (replace "pending"),
   update this handoff, commit, report.

## Open problems and decisions

- zlib-ng not added: no system package; compat mode would replace zlib for
  htslib and STAR's link line, and vendoring the native API is a cmake build
  with SIMD dispatch. Not small and clean.
- Single-end mapping keeps its loader (zlib; OpenMP tasks); only the barcode
  abundance pass and paired-end lanes use the new intake.
- With BGZF ATAC input inside STAR, inflate workers run outside STAR's permit
  hooks. The mirrored reader accepts `BgzfWorkPermitHooks` with the same
  signature as `MappingParameters`' permit hooks; wiring them is a follow-up.
  Production DOGMA-plex inputs are plain gzip.
- Standalone at 32 threads Chromap is mostly mapping-bound on plain gzip, so
  the gzip threading changes wall time little there; the harness shows the
  loader itself.
- BGZF timings use `bgzip` copies of the lane-1 subset (test inputs; the
  production files are plain gzip).
- Dockerfile `CHROMAP_SUITE_VERSION`/`REVISION` ARGs still pin v1.0.1.
- Left in place for inspection: `/mnt/pikachu/STAR-suite-readercheck-20260928`
  (detached, built, untracked `multiomics_build.json`) and
  `/mnt/pikachu/multiomics-suite-readercheck-20260928` (detached, edited
  `compatibility_manifest.json`, uncommitted by design).
