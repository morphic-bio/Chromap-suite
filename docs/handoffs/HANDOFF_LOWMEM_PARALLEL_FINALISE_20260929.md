# Handoff: parallel low-memory finalisation for paired-end ATAC (2026-09-29)

- Runbook: `docs/runbooks/RUNBOOK_LOWMEM_PARALLEL_FINALISE_20260929.md`. §10
  has the M-1 results; §11 has the 1.2.0 results through M4.
- Design note: `docs/design/LOWMEM_PARALLEL_FINALISE_20260929.md`.

## State: M4 and D17 complete (29 September); stopped for the author

### Branches

| Branch | Worktree | State |
|---|---|---|
| `feat/lowmem-parallel-finalize-20260929` (1.2.0) | `/mnt/pikachu/Chromap-suite-lowmem-feat-20260929` | Local only: nothing pushed, tagged or merged to master. Contains the 1.1.1 fix branch (merged, `9d3c30b`) and `origin/master` `128ead3` (v1.1.1, merged, `600dce4`). Code commits: `f2f1b9f` (per-reference path), `e14c8d3` (lean decode), `ebd474b` (docs, 16-target release gate), `6787063` (D17 summary aggregation). |
| `fix/v1.1.1-lowmem-edge-cases` | `/mnt/pikachu/Chromap-suite-v111-fix-20260929` | Released as v1.1.1 (`128ead3` on `origin/master`). |
| `design/lowmem-parallel-finalise-20260929` | — | No longer updated. |

### Validation root

`/mnt/pikachu/lowmem_parallel_validation_20260929`:

| Path | Contents |
|---|---|
| `bin/` | Baseline `chromap_base_v111` (`d0ae41d2…`); final build `chromap_new_e14c8d3` (`5572046f…`); `SHA256SUMS` |
| `compare/m4_summary.tsv` | 51 of 51 comparisons identical, plus the expected open-file error |
| `unit/` | Unit harness logs, baseline goldens, mutation checks |
| `measure/` | strace check, lane-1 per-reference table, D2 spill listing |
| `logs/` | Release gate and smoke runs |

### Results

- **Identity.** Every case is byte-identical to the 1.1.1 serial merge:
  - the synthetic fixture as 400 lanes;
  - PBMC 100k, and PBMC 100k as 256 lanes;
  - DOGMA 2M, with and without `--output-mappings-not-in-whitelist`;
  - DOGMA 50M;
  - PBMC 3k at full depth;
  - DOGMA lane 1 at full depth. Its sidecar is also byte-identical to the
    v1.1.0 production sidecar that Multiomics 0.9.0 wrote for the same lane.
- **Tests.**
  - `make test-release`: 16 of 16 targets pass.
  - The release-pipeline unittest passes.
  - The fixture smokes pass, including the runs forced to the serial merge.
  - The unit harness passes 143 runs.
- **Informal finalisation times** (not benchmarks; shared host, single runs):

  | Case | Serial | Parallel |
  |---|---|---|
  | F1 (DOGMA lane 1, full depth) | 680 s | 211 s |
  | D2 (DOGMA 50M) | 32 s | 13.5 s |
  | P3 (PBMC 3k, full depth) with `--summary` | 43 s | **76 s** |
  | P3 without `--summary` | — | 7.6 s |

## D17 (29 September): summary contention fixed

- Approved by the author and implemented in `6787063`. The harness was
  broadened in `e07edd7`.
- **Change.** Tasks aggregate summary deltas per barcode and apply them with
  one atomic add per barcode and field when they end.
- **Re-verification.** All of these pass:
  - the unit harness: 143 runs, plus the regenerated baseline goldens;
  - the mutation checks;
  - `make test-release` (16 targets) and the pipeline unittest;
  - all 51 identity comparisons, which include full-depth PBMC 3k, DOGMA 50M
    and the full-depth lane. The lane sidecar also equals the Multiomics 0.9.0
    production sidecar.
- **Timing.** PBMC 3k full with `--summary` now finalises in 2.3 s in parallel,
  against 43 s serial (informal).
- **Final binary.** `V/bin/chromap_new_6787063` (`f4a6d570…`).

## Waiting on the author

1. Releasing 1.2.0:
   - bump `src/version.h` to 1.2.0;
   - date the release notes and the changelog;
   - set the Dockerfile version and source revision;
   - tag and push.
2. Multiomics adoption, in multiomics-suite:
   - `--chromapAtacLowMemFinalizeThreads` (the library default is 1);
   - G-M1 on a build pinned to 1.2.0.

## Excluded material

Follow the exclusions in the maintainers' private notes; never read or copy the
material they name.
