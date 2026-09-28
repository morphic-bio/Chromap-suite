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
| `a88fa44` | runbook and initial handoff |
| `b30347e` | completed throughput results, 50M output identity and release-note measurements |
| `fc47c2f` | integration merge onto master `a7d2ba1`, preserving the feature history |
| release finalization | dated release notes/changelog, Docker version/revision and release-build provenance |

The user subsequently authorized integration. It was prepared in
`/mnt/pikachu/Chromap-suite-v110-integration-20260928`, branch
`integration/chromap-v1.1.0-20260928`. Merge `fc47c2f` combines master
`a7d2ba1` and feature tip `b30347e` without conflicts. It retains the barcode
learning and sidecar prerequisites and every feature commit, including
validated source revision `66d4537`; the source and Makefile are identical
to that revision. The feature worktree and its validation artifacts remain
at their original revisions.

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

- Throughput: all seven cases have at least one clean result. The campaign
  is no longer running; PID `1009455` is gone. Pass2 completed
  `chromap_old_bgz` and `chromap_new_gz`, then stopped while waiting before
  its redundant `chromap_old_gz` run (that directory contains only `tmp/`;
  use the clean pass1 baseline). Pass3 completed the missing clean
  `harness_kseq_threads_gz` case at 18:13:34 UTC on 2026-09-28. No further
  timing runs are required by this runbook.
- Final summary: `V/throughput/throughput.{tsv,json}` contains 14 completed
  timed runs, eight clean and six contaminated. All completed runs exited
  zero and recorded lock acquisition/release. The three timed binaries'
  SHA256 hashes were checked against the archived values. A tracked copy of
  all measurements is in
  [`docs/benchmarks/fastq_intake_20260928/throughput.tsv`](../benchmarks/fastq_intake_20260928/throughput.tsv).
- Post-campaign output comparison: `nice -n 10 python3
  V/throughput/compare_outputs.py` -> `V/throughput/identity.json`,
  `VERDICT: PASS`. All seven completed 50M-pair mapping runs, including
  contaminated timings, have identical sidecar and `.chroms.tsv` hashes.
  Summary differences are limited to `cachehit`, `fric`, `estfrip` and
  `numcacheslots`. The report with full hashes and binary provenance is
  [tracked here](../benchmarks/fastq_intake_20260928/identity.json).
- Release notes now contain the throughput numbers and 50M output check.
  Source and build files remain unchanged from validated commit `66d4537`.
  The existing reports were also checked: all 26 standalone comparisons,
  six deterministic comparisons and eight required embedded files pass.

Final clean results, 50,000,000 lane-1 pairs, 32-thread setting, inputs
evicted and reference/index warmed, Intel Core i9-13900KF (32 logical CPUs):

| Case | Pass | Wall (s) | Max RSS (GiB) |
|---|---|---:|---:|
| chromap_old_gz | pass1 | 145.2 | 22.9 |
| chromap_new_gz | pass2 | 138.6 | 23.0 |
| chromap_old_bgz | pass2 | 141.7 | 23.1 |
| chromap_new_bgz | pass1, pass2 | 122.8, 122.3 | 22.4, 22.1 |
| harness_kseq_serial_gz | pass2 | 82.0 | 1.2 |
| harness_kseq_threads_gz | pass3 | 30.3 | 1.2 |
| harness_bgzf_bgz | pass2 | 12.5 | 1.0 |

Full mapping wall time decreased 4.5% on gzip and 13.3–13.7% on BGZF.
Harness loading rates are 0.612M, 1.667M and 4.042M pairs/s for serial
gzip, threaded gzip and BGZF, respectively. Harness rates use the internal
loading timer and compare reader strategies; full-run version comparisons
use process wall time. These are individual observations, with two clean
new-BGZF runs and one clean run for every other case. Contaminated runs are
retained as evidence but excluded from these comparisons.

## Integration validation

Artifacts are under
`/mnt/pikachu/Chromap-suite-v110-integration-20260928/plans/artifacts/v1.1.0-integration-20260928/`.
Builds and checks were invoked at `nice -n 10`; smoke targets ran serially.

- Fresh `make -j8 all chromap_lib_runner tests/fastq_intake_harness`: pass.
- `make test-smoke`: pass, including all 12 peak-integration matrix cells
  with identical narrowPeak and summit outputs, and all three low-memory
  BED parity cases.
- `make test-fastq-intake-smoke`: 27/27 pass.
- `make test-atac-sidecar-only-smoke`: 17/17 pass.
- `make test-libchromap-core-smoke`: 10/10 pass.
- Installed-header/static-library SDK link check: pass.
- Release tarball builds and its `chromap --version` is `1.1.0`; the engine
  remains `0.3.3-r519`. The tarball's `licenses/STAR-input-LICENSE` matches
  the mirrored source notices byte for byte.
- Release workflow YAML, Docker build arguments, ancestry and the repository's
  attribution-history check: pass. Historical commits and validation evidence
  are preserved.
- `docker buildx build --builder default --platform linux/amd64 --load
  -t local/chromap-suite:v1.1.0 .`: pass, including the installed-SDK link
  check in the Dockerfile. Local image ID:
  `sha256:28a256637244ae4d7eea205e0add927879fdd22ba4fcc293b8a8b0c88b1b4624`.
- `tests/run_container_smoke.sh local/chromap-suite:v1.1.0`: pass, 16
  identical CLI/library alignment rows, version `1.1.0`, source revision
  `fc47c2f`. The installed STAR notices match `src/star_input/LICENSE`.
  Container validation for this integration covers amd64; the container
  guide separately retains historical v1.0.1 arm64 evidence.

`run_checks.sh`, `checks.tsv`, individual logs, SDK staging, the local tarball
and its checksums are retained in that artifact directory. The container
smoke now derives its expected version/revision from Dockerfile, with
`EXPECTED_VERSION`/`EXPECTED_REVISION` overrides for tagged images.

## Next

Local integration includes merge `fc47c2f` and release-finalization commit
`0a4c71c`. The subsequent user-requested release pipeline adds test-gated
tarballs, Debian binaries and source packages, clean-runtime validation, and
Docker publication dependencies. See the [release procedure](../releasing.md)
and [pipeline validation runbook](../runbooks/RUNBOOK_TEST_GATED_RELEASE_20260928.md).
The unpublished local `v1.1.0` tag is advanced with that integration; the
earlier container and tarball above remain evidence for the initial merge.
The primary checkout's pre-existing untracked files are preserved. Publishing
remains a separate step; nothing has been pushed to origin or an image registry.

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
- Integration updated Dockerfile `CHROMAP_SUITE_VERSION` to `1.1.0` and
  `REVISION` to merge `fc47c2f`. The release workflow passes the version
  derived from the tag and the exact tagged commit as build arguments.
- Left in place for inspection: `/mnt/pikachu/STAR-suite-readercheck-20260928`
  (detached, built, untracked `multiomics_build.json`) and
  `/mnt/pikachu/multiomics-suite-readercheck-20260928` (detached, edited
  `compatibility_manifest.json`, uncommitted by design).
