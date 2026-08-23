# Chromap Suite v1.0.1 Release Notes

Date: 2026-08-23

Chromap Suite v1.0.1 is a correctness and reliability release. It embeds
RapidMACS v1.0.1 and makes its MACS3 v3.0.3 compatibility corrections available
through both the standalone Chromap Suite peak-calling surface and applications
embedding `libchromap` + `librapidmacs`.

`chromap --version` reports `1.0.1`. The underlying chromap engine version is
unchanged and remains available through `chromap --upstream-version`.

## Peak-calling correctness

The embedded RapidMACS dependency is pinned to public release v1.0.1
(`c22ff439cc66337e1b7fc9b87b547478eb5c0d3b`). That release corrects two
MACS3 `callpeak` edge cases:

- Peak regions use MACS3's strict `score > cutoff` predicate in p-value and
  q-value modes. A score exactly equal to the requested cutoff is no longer
  included in a peak boundary.
- MACS3's default fold-enrichment cutoff of 1.0 is applied after peak
  construction, excluding statistically significant depletion regions from
  narrowPeak output.

The corrections apply to FRAG, BAMPE, BEDPE, and BED callers and to in-process
peak calling from Chromap. The historical custom streaming caller is unchanged,
and the separate `bdgpeakcall` diagnostic retains its inclusive threshold to
match MACS3's `bdgpeakcall` command.

There are no public command-line or library API compatibility breaks.

## ATAC reliability

- Direct and materialized ATAC paths now use the same deterministic fragment
  reduction before peak calling, preserving exact peak-caller input equivalence.
- The opt-in mergeable spill materializer validates shard identity and coverage,
  performs global barcode correction and deduplication, and emits canonical
  fragments plus supported materialized outputs.
- Fixed-width materialized blocks and reference sidecars reduce downstream
  parsing and index-loading overhead while retaining fail-closed format and
  provenance checks.
- ATAC mapping releases shared concurrency permits at task boundaries so other
  work can use idle capacity without changing mapping or peak outputs.

## Naming and workflow updates

- The peak-calling dependency and user-facing companion executable use the
  canonical RapidMACS / `librapidmacs` names. `chromap_callpeaks` remains as a
  compatibility symlink.
- The MCP/Launchpad registry includes a ChIP TagAlign workflow schema.

## Validation

- RapidMACS exact-cutoff unit regressions cover p-value and optimized q-value
  modes.
- Available full-input ATAC, ChIP-seq, and CUT&RUN checks reproduce complete
  narrowPeak and summit outputs byte-for-byte against MACS3 v3.0.3.
- Direct-versus-materialized ATAC parity, mergeable-spill materialization,
  fragment storage, q-value CLI behavior, and the hermetic Chromap smoke tier
  pass on the release candidate.

## Versioning and provenance

- `chromap --version` -> `1.0.1` (suite).
- `chromap --upstream-version` -> the unchanged chromap engine version.
- `third_party/rapidmacs` -> RapidMACS v1.0.1 at
  `c22ff439cc66337e1b7fc9b87b547478eb5c0d3b`.
- The annotated `v1.0.1` tag identifies the release source; release CI verifies
  that the binary version matches the tag before publishing artifacts.
