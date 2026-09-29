# Handoff: parallel low-memory finalisation for paired-end ATAC (2026-09-29)

Runbook: `docs/runbooks/RUNBOOK_LOWMEM_PARALLEL_FINALISE_20260929.md`.
Design note: `docs/design/LOWMEM_PARALLEL_FINALISE_20260929.md`.

## State

- **Design branch** `design/lowmem-parallel-finalise-20260929`, in worktree
  `/mnt/pikachu/Chromap-suite-lowmem-parallel-20260929`, starts from `master`
  `a47f077` (v1.1.0). It holds the design note, the runbook and this handoff,
  committed on the local branch.
- **Author decisions (29 September)** are in runbook §7; the TODO list is in
  §7.1. The main points:
  - Every per-reference task runs under a host permit (D4).
  - The summary uses atomic adds into existing khash entries and logs only
    new barcodes, in order (D8).
  - TagAlign is in scope; bulk is out (D6).
  - Two v1.1.0 bugs were fixed first, as Chromap Suite 1.1.1 (D11).
- **M-1 (Chromap Suite 1.1.1) is complete.** Results are in runbook §10.
  - **Branch** `fix/v1.1.1-lowmem-edge-cases`, in worktree
    `/mnt/pikachu/Chromap-suite-v111-fix-20260929`. Local commits `777462b`,
    `73cafcd` and `32e8808` sit on `v1.1.0`.
  - **Validation root:** `/mnt/pikachu/chromap_v111_validation_20260929`.
  - **Regression test:** fails 4 of 5 cases on v1.1.0 and passes 5 of 5 on
    1.1.1.
  - **Release gate:** `make test-release` passes 15 of 15.
  - **Identity:** all 12 comparisons against v1.1.0 are identical, plus the
    v1.1.0 control.
  - Nothing has been pushed, tagged, released, merged or rebased between the
    branches.
- **M0-M4 (1.2.0)** have not started. They wait for the coordinator.

## Waiting on the author

1. **Release of 1.1.1:** the tag and the push. The release date in
   `CHANGELOG.md` and the release notes, and the Dockerfile's default source
   revision, are set in that step.
2. **How 1.2.0 picks up 1.1.1:** merge the fix branch into the design branch,
   or rebase the design branch onto it.
3. **The N = 1 permit question** (TODO §7.1): should the N = 1 serial loop
   hold a permit when hooks are present?
4. **Public-repository text.**
   - None of this work's commits contains the excluded names or terms.
     Checked with `git log -p a47f077..<branch>` on both branches.
   - The already-public `origin/master` has two such lines, in
     `docs/runbooks/RUNBOOK_STAR_MIRRORED_FASTQ_READER_20260928.md`, from
     commit `a88fa44`, which is in `v1.1.0`.
   - A whole-tree grep of any branch therefore still finds them.
   - Whether to neutralise that file in a new commit is for the coordinator
     and the author.

## Excluded material

Follow the exclusions in the private coordinator handoff; never read or copy
the directories it names.

## Estimate for the rest

About 4.25 agent-days for M0-M4, plus about 11-12 hours of untimed identity
runs, run one at a time (runbook §8).
