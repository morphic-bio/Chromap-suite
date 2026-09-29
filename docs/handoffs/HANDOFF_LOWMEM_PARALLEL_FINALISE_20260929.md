# Handoff: parallel low-memory finalisation for paired-end ATAC (2026-09-29)

- Runbook: `docs/runbooks/RUNBOOK_LOWMEM_PARALLEL_FINALISE_20260929.md`. §10
  has the M-1 results; §11 has the 1.2.0 results.
- Design note: `docs/design/LOWMEM_PARALLEL_FINALISE_20260929.md`.

## State (update of 29 September, during M4)

### Branches

| Branch | Worktree | Head / commits | Notes |
|---|---|---|---|
| `feat/lowmem-parallel-finalize-20260929` (implementation) | `/mnt/pikachu/Chromap-suite-lowmem-feat-20260929` | from the design branch, with `fix/v1.1.1-lowmem-edge-cases` merged in (`9d3c30b`) | Local commits up to `ebd474b`; the runbook and handoff are kept here from now on |
| `fix/v1.1.1-lowmem-edge-cases` (1.1.1) | `/mnt/pikachu/Chromap-suite-v111-fix-20260929` | `25afe27` | Public `origin/master` `71030e9` merged in; whole-tree check clean. Waiting for the author's release |
| `design/lowmem-parallel-finalise-20260929` | — | `31a10ee` | No longer updated |

### Validation root

`/mnt/pikachu/lowmem_parallel_validation_20260929`:

- `bin/`: baseline 1.1.1 `chromap_base_v111` (`d0ae41d2…`) and final build
  `chromap_new_e14c8d3` (`5572046f…`), with `SHA256SUMS`.
- `runs/` and `compare/`: end-to-end identity runs and their comparisons.
- `unit/`: unit logs and the baseline goldens.
- `measure/`: the measurements.
- `scripts/`: the drivers. These include `run_baseline.sh` and `run_new.sh`,
  which run in the background with the shared lock.

### Milestones

- **M0:** done (runbook §11).
- **M1/M2:** done, commit `f2f1b9f`.
- **M3:** done, commit `e14c8d3`.
  - Unit harness: 143 runs pass, including a check against the baseline
    goldens.
  - Mutation checks: the tests catch every ordering mistake that matters. The
    one mutation they miss is harmless; runbook §11 explains why.
- **M4:** in progress.
  - Documentation and the 16-target release gate are committed (`ebd474b`).
  - `make test-release` and the end-to-end runs are running.

### Decisions applied

- D13-D16 are recorded in runbook §7.
- 1.2.0 builds on 1.1.1 by merge.
- The N = 1 serial merge holds one permit when hooks are present.
- Runs that load the genome index take the shared lock.

## Waiting on the author

1. Releasing 1.1.1: the tag, the push and the GitHub release.
2. Releasing 1.2.0 after M4 reports: the version bump, release notes, tag and
   push.
3. Multiomics adoption. In multiomics-suite:
   - pass `--chromapAtacLowMemFinalizeThreads` (the library default is 1);
   - run G-M1 with a build pinned to the 1.2.0 commit.

## Excluded material

Follow the exclusions in the maintainers' private notes; never read or copy the
material they name.

## If usage runs short

Stop and update this handoff. The background drivers write each run to
`V/runs/<case>/<mode>/<binary>` and log progress to `V/logs/run_new.log`.
Compare finished runs with `V/scripts/compare_case.py BASE NEW MODE OUT.json`.
