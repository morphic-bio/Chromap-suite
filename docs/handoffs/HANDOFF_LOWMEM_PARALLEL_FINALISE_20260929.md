# Handoff: parallel low-memory finalisation for paired-end ATAC (2026-09-29)

Runbook: `docs/runbooks/RUNBOOK_LOWMEM_PARALLEL_FINALISE_20260929.md`.
Design note: `docs/design/LOWMEM_PARALLEL_FINALISE_20260929.md`.

## State

**Design branch.**

- Branch `design/lowmem-parallel-finalise-20260929` in the worktree
  `/mnt/pikachu/Chromap-suite-lowmem-parallel-20260929`, from `master`
  `a47f077` (v1.1.0).
- It holds the design note, the runbook and this handoff, committed on the
  local branch.

**Author decisions (29 September).** They are in runbook §7; the TODO list is
in §7.1. The main changes from the first draft:

- Every per-reference task runs under a host permit (D4).
- The summary uses atomic adds into existing khash entries, and logs only new
  barcodes, in order (D8).
- TagAlign is in scope and bulk is out (D6).
- Two v1.1.0 bugs are fixed first, as Chromap Suite 1.1.1 (D11).

**Current milestone: M-1, Chromap Suite 1.1.1** (runbook §4.0 and §5.8).

- Branch `fix/v1.1.1-lowmem-edge-cases` in the worktree
  `/mnt/pikachu/Chromap-suite-v111-fix-20260929`, from `v1.1.0`.
- Validation root: `/mnt/pikachu/chromap_v111_validation_20260929`.

**M0-M4 (1.2.0)** have not started. They wait for the coordinator after M-1.

## Waiting on the author

- Approval to release 1.1.1. The tag and push are the author's.
- How 1.2.0 work picks up 1.1.1: merge it into the design branch, or rebase
  onto it.
- The question in TODO §7.1: should the N = 1 serial loop hold a permit when
  hooks are present?

## Estimate

About 4.75 agent-days in total, including M-1 (range 3.5-6). Add about
12-14 hours of untimed identity runs, one at a time (runbook §8).
