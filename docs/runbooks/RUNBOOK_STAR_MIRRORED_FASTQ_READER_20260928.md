# Runbook: STAR Suite FASTQ reader mirrored into Chromap Suite (2026-09-28)

## Goal and scope

Give Chromap the FASTQ intake STAR Suite uses, for Chromap Suite v1.1.0:

- mirror STAR Suite's in-tree BGZF FASTQ reader into `src/star_input/` and
  read BGZF paired-end lanes through a `PairedEndReadProvider` over it
  (`--input-bgzf-mode auto|on|off`, `--input-bgzf-reader-threads N`, also
  `MappingParameters` fields for embedding hosts);
- read ordinary gzip with zlib, but with read 1, read 2 and the barcode file on
  their own threads in every batch;
- bump the version to 1.1.0 and draft `docs/RELEASE_NOTES_v1.1.0.md`
  (with the base commits `7265546` bounded barcode learning and `09b3164`
  sidecar-only output);
- validate record identity, output identity (standalone and embedded in STAR)
  and measure throughput.

A later shared open `libfastq` replaces the mirrored copy; only
`src/fastq_bgzf_input.{h,cc}`, `src/star_input/` and one Makefile line see it.

Out of scope: any change to STAR-suite (STAR no longer hosts Chromap;
Multiomics Suite owns the single binary). No push, merge or tag; the
coordinator releases.

## Rules

- Excluded material: follow the exclusions in the maintainers' private notes.
  Mirror only from STAR Suite `core/legacy/source/input/`. If something seems
  to exist only in excluded material, stop and report.
- Clean room: never read 10x Genomics code.
- STAR-suite: no source changes. A detached worktree for the embedded build
  is allowed.
- Git: plain commit messages, no AI attribution, no push/merge/tag.
- Do not touch other agents' worktrees (`multiomics-suite-analysis`,
  `rapidglm`, the e2e bench, nm-refresh, production checkouts).
- Timed-run lock: every timed run takes
  `/mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock` (`flock`) and holds it
  from the page-cache eviction to the end of the run. Untimed builds and tests
  run at `nice -n 10` or lower and need no lock.

## Locations

| What | Path |
|---|---|
| Chromap worktree, branch `feat/star-mirrored-fastq-reader` | `/mnt/pikachu/Chromap-suite-fastq-mirror-20260928` |
| Base (must stay unchanged) | `/mnt/pikachu/Chromap-suite-atac-nobam-20260928` (`09b3164`) |
| STAR source for the mirror (read only) | `/mnt/pikachu/STAR-suite` or `/mnt/pikachu/STAR-suite-atac-nobam-20260928`, commit `6e83853` |
| STAR check worktree (detached `6e83853`) | `/mnt/pikachu/STAR-suite-readercheck-20260928` |
| Recipe copy (detached multiomics `5791d78`, manifest edited, uncommitted) | `/mnt/pikachu/multiomics-suite-readercheck-20260928` |
| Validation root (data, binaries, scripts, results) | `/mnt/pikachu/fastq_mirror_validation_20260928` |
| Archived binaries + `SHA256SUMS` | `.../bin/` (`chromap_old_09b3164`, `chromap_new_66d4537`, `fastq_intake_harness_66d4537`, ...) |
| DOGMA-plex lane-1 subsets | `.../lane1_4M/`, `.../lane1_50M/` (`R{1,2,3}.fastq.gz` pigz, `.fastq.bgz` bgzip copies) |
| PBMC 100k fixture bgzip copies | `.../pbmc_100k_bgzf/` |
| Record identity | `.../records/` |
| Standalone identity | `.../identity/`, `.../identity_det/` (+ `identity.json`) |
| Embedded check | `.../embedded/` |
| Throughput | `.../throughput/` |
| Final tracked throughput and 50M output identity | `docs/benchmarks/fastq_intake_20260928/{throughput.tsv,identity.json}` |

## Steps

1. **Worktree.**
   `git -C /mnt/pikachu/Chromap-suite worktree add -b feat/star-mirrored-fastq-reader /mnt/pikachu/Chromap-suite-fastq-mirror-20260928 09b3164`,
   then initialise RapidMACS from the local checkout without changing the
   shared config:
   `git -c protocol.file.allow=always -c submodule.third_party/rapidmacs.url=/mnt/pikachu/Chromap-suite-atac-nobam-20260928/third_party/rapidmacs submodule update --init third_party/rapidmacs`.
   Check: `git submodule status` shows `34df448`; the shared config still
   has the GitHub URL.
2. **Old binary.** Build at `09b3164` (`make -j16 chromap chromap_lib_runner`)
   and keep it as `bin/chromap_old_09b3164` (sha256 `975c9225...`).
3. **Mirror.** `git show 6e83853:core/legacy/source/input/<f>` for
   `BgzfBlockReader.{h,cpp}` and `BgzfRangeReader.{h,cpp}` into
   `src/star_input/`, then edit namespaces (`star::input` ->
   `chromap::star_input`), include paths (`input/` -> `star_input/`) and
   include guards only. Check: `diff` against the source shows only those
   lines; `src/star_input/MIRROR.md` has the per-file sha256 table; the four
   files are identical at STAR `master` and `v1.9.5.a`.
4. **Glue, provider, options.** `src/fastq_bgzf_input.{h,cc}`
   (`SelectBgzfFastqInput`, `BgzfInflateWorkersPerStream`, `BgzfFastqStream`,
   `BgzfPairedEndReadProvider`), `MapPairedEndReads` lane selection,
   barcode abundance pass, CLI and `chromap_lib_runner` options,
   `MappingParameters` fields. Check: `nice -n 10 make -j8 chromap
   chromap_lib_runner tests/fastq_intake_harness` builds with no warnings from
   the new files; STAR's link contains both `star::input` and
   `chromap::star_input` symbols (`nm -C STAR`).
5. **Plain gzip threads.** `Chromap::LoadPairedEndReadsWithBarcodes` uses a
   `std::thread` per file for every batch; later batches pass `num_threads >= 3`.
   Check: FIFO cases in the smoke test (step 6).
6. **Tests.** `make test-fastq-intake-smoke` (27 cases), `make test-smoke`,
   `test-input-format-smoke`, `test-libchromap-core-smoke`,
   `test-atac-sidecar-only-smoke`, `test-cbq-atac-smoke`,
   `test-cbq-modality-matrix`, `test-cbq-range-reader`,
   `test-atac-runtime-spill-schema-harness`, with
   `CHROMAP_ARTIFACT_ROOT` outside the repository. Check: every target exits 0.
   `make test-smoke` rebuilds the tracked `tests/test_frag_compact_store`;
   restore it with `git checkout -- tests/test_frag_compact_store`.
7. **Subsets.** `make_subsets.sh` (first 50,000,000 records of lane 1) and
   `make_subsets_step2.sh` (4M and 50M, pigz and bgzip). Check:
   `records/dogma50M_records.tsv`: the plain extract, the gzip and the BGZF
   copy give the same record dump (150,000,000 records).
8. **Record identity.** `bin/fastq_intake_harness_66d4537 --reader
   kseq-serial|kseq-threads|bgzf|auto --threads N --dump - R1 R3 R2 | sha256sum`
   (`records/*.tsv`, `records/bgzf_66d4537.sh`). Check: one sha per dataset.
9. **Standalone output identity.**
   `run_identity.sh OLD NEW identity pbmc100k dogma4M pbmcfull`,
   `run_identity_det.sh`, then `python3 compare_identity.py identity` and
   `identity_det`. Check: fragments and sidecars identical; `--summary`
   differs only in `cachehit,fric,estfrip,numcacheslots` (also old vs old);
   identical under `--deterministic-mapping`.
10. **Embedded (no STAR changes).**
    `git -C /mnt/pikachu/STAR-suite worktree add --detach /mnt/pikachu/STAR-suite-readercheck-20260928 6e83853`;
    `git -C /mnt/pikachu/multiomics-suite-atac-nobam-20260928 worktree add --detach /mnt/pikachu/multiomics-suite-readercheck-20260928 5791d78`
    and edit only its `compatibility_manifest.json` (star path, chromap
    path/revision/tree, rapidmacs path, composed_artifact shas; copy in
    `embedded/compatibility_manifest_readercheck.json`). Build:
    `nice -n 10 python3 scripts/build_star_composition.py --log-dir .../embedded/build_logs --jobs 8`
    in the recipe copy (the Chromap worktree must be clean at the pinned
    commit). Run the baseline with the no-BAM recipe and the new build with
    the copy:
    `python3 <recipe>/scripts/run_dogmaplex_lane.py --lane 1 --fastq-root /mnt/pikachu/dogmaplex_gse309834/production_20260926/acquisition/fastq --max-reads 2000000 --threads 16 --skip-active-check --atac-output sidecar --output-root .../embedded/runs --run-id <id>`
    with `MULTIOMICS_STAR_SUITE_DIR` / `MULTIOMICS_STAR_BIN` set for the new
    build. Check: `python3 embedded/compare_embedded.py runs/nobam_baseline runs/readercheck_66d4537 out.json`
    prints `VERDICT: PASS`.
11. **Throughput.** `throughput/run_throughput.sh PASS CASE...` with cases
    `chromap_{old,new}_{gz,bgz}` and `harness_{kseq_serial_gz,kseq_threads_gz,bgzf_bgz}`.
    Each case waits for a quiet host (no STAR process, load1 <= 4) before
    taking the lock, evicts the three inputs (`evict_cache.py`), warms the
    index, and runs under
    `run_timed.py --record-dir DIR --time-output DIR/time.txt`.
    `python3 throughput/summarize.py` writes `throughput.tsv`/`.json`
    (`NA` in the TSV marks fields not reported by that case).
    Check: each case needs at least one `clean` verdict; rerun cases lacking
    a clean result under a new pass label. The recorded clean criteria are
    start load1 <= 4, mean outside CPU fraction <= 0.05, outside md0 I/O
    fraction <= 0.05 and no other STAR process seen.
    After all timing jobs have stopped, run
    `nice -n 10 python3 throughput/compare_outputs.py` from the validation
    root. It requires clean coverage of all seven cases, successful exits
    and lock records for every completed run, and the archived binary
    hashes. It then compares every completed mapping run's sidecar and
    `.chroms.tsv` to `pass1/chromap_old_gz`, allowing only candidate-cache
    column differences in summaries. Check: `throughput/identity.json`
    reports `PASS`. Copy the final `throughput.tsv` and `identity.json` to
    `docs/benchmarks/fastq_intake_20260928/` and replace the release notes'
    pending measurements. The final campaign has clean coverage in
    pass1/pass2/pass3; see the handoff for selected runs.
12. **Version and notes.** `src/version.h` -> `1.1.0`;
    `docs/RELEASE_NOTES_v1.1.0.md`, `CHANGELOG.md`. The initial reader branch
    left Docker version/revision ARGs for the coordinator. The subsequent
    user-authorized integration sets the release date to 2026-09-28,
    `CHROMAP_SUITE_VERSION` to `1.1.0` and the default source revision to
    integration merge `fc47c2f`. Tagged Docker builds override both ARGs with
    the release version and actual tagged commit. See the handoff for the
    integration checks and publication state.
