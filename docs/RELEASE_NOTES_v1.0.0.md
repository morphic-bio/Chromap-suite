# Chromap Suite v1.0.0 Release Notes

Date: 2026-06-18

`v1.0.0` is the first production Chromap Suite release. The release artifact
version is `1.0.0`, and `chromap --version` now reports the suite version
`1.0.0`. The underlying chromap engine version is available through
`chromap --upstream-version` (`0.3.3-r519`) and is written to the BAM/SAM `@PG`
`VN` tag for provenance. Chromap Suite was spun off from
[Chromap](https://github.com/haowenz/chromap) in 2026 and no longer tracks
upstream; see [`HISTORY.md`](../HISTORY.md) for lineage.

## Production Surfaces

- **ATAC-seq alignment with native BAM output** — coordinate-sorted, indexed BAM
  produced directly by the aligner, alongside the existing BED/BEDPE/PAF and
  fragments outputs. Bulk ATAC, scATAC, ChIP-seq, and Hi-C all run through the
  same `chromap` binary.
- **In-process libMACS3 narrow peak calling** (`--call-macs3-frag-peaks`) — a
  single `chromap` invocation produces sorted/indexed BAM, fragments, and
  MACS3-equivalent narrowPeak + summit outputs with no intermediate
  `fragments.tsv.gz` write. Output is **byte-identical to standalone MACS3
  v3.0.3** in the default p-mode path (50,274 peaks, md5 `34f9f991…` on the 3K
  PBMC ATAC channel). `--macs3-frag-qvalue Q` switches to q-value/FDR
  thresholding for `macs3 callpeak -q` compatibility; bulk ATAC via
  `--macs3-frag-peaks-source memory`.
- **`libchromap.a` callable library** — the full Chromap Suite ATAC pipeline
  exposed through the `libchromap` API (`RunAtacMapping()`, `ChromapAtacConfig`,
  `ChromapPermitHooks`). The same library backs the `chromap` CLI and is the
  integration point used by STAR Suite for multiomic processing.
- **Multiomic integration with STAR Suite** — `libchromap.a` runs as a STAR
  Suite worker thread for concurrent ATAC + GEX processing in a single `STAR`
  invocation, with a shared-thread-budget permit allocator. End-to-end on the
  public 3K PBMC Multiome: **18:17 / 64.8 GB peak RSS, 2.19× faster and ~18%
  lower memory than Cell Ranger ARC v2.2.0** (40:04 / 79.1 GB).
- **ATAC fragment sidecar** (AEV1 binary, `--atac-fragment-binary-output`) — a
  compact 32-byte-header + 24-byte-record file emitted alongside BAM/fragments,
  letting STAR Suite's `libscrna` call ATAC cells without re-parsing the gzipped
  fragments TSV (53,969,811 records, md5 `a4251bbc…` on the headline run).
- **Native CBQ input** (`--input-format cbq`) — maps paired-end ATAC/scATAC
  reads straight from BINSEQ CBQ files with no intermediate FASTQ, producing
  fragments byte-identical (under canonical sort) to the FASTQ path. ATAC-only
  for this milestone.

## Reliability

- **Rewritten low-memory spillover** — per-thread overflow writers feeding a
  *k*-way merge on read-back replace the prior shared-buffer + atomic-write
  design, and are now the **default build** (no compile flag). Supports the full
  cross-product of `--low-mem` with `--atac-fragments`, BAM output,
  `--macs3-frag-low-mem`, and Y-filtering at production scale (≳10⁹ reads). A
  pre-existing legacy race that produced silent read drops (~1/10⁴) and
  intermittent hangs at scale is resolved as a side effect. The legacy path
  remains available via `LEGACY_OVERFLOW=1` (single-threaded only).
- **SAM/serialization hardening** — atomic, length-safe record writes;
  corrected `SerializedSize()` accounting; coordinated thread cleanup; and
  rid-ascending overflow processing to preserve coordinate-sorted output.

## Additional Surfaces

- **Y-chromosome filtering** (`--emit-Y-bam`, `--emit-noY-bam`,
  `--emit-Y-noY-fastq`, `--emit-Y-read-names`) for sex-aware analyses; works
  with `--sort-bam`.
- **`--Tn5-shift-mode {classical|symmetric}`** selects the Tn5 cut-site offset
  convention on BED/BEDPE/PAF (`classical` `+4/-5` default; `symmetric` `+4/-4`).
- **`--temp-dir DIR`** for a custom temporary directory.

## Tooling, MCP, and Tests

- **MCP server + Launchpad** (`mcp_server/`) — a schema-driven workflow renderer
  for agents plus a browser recipe UI for humans, both consuming the same
  parameterised YAML recipe registry.
- **Regression suite (C01–C11)** — an 11-area parity matrix (index build, paired
  BED, ChIP/ATAC presets, scATAC barcodes, sorted BAM + index, low-mem BED
  parity, ATAC BAM + fragments, libMACS3 narrow peaks, Hi-C pairs, Y/noY split)
  across three tiers (S0 hermetic synthetic, S1 ENCODE downsample, S2 heavy
  fixtures). S0 is mandatory for pre-commit checks.

## Versioning and Provenance

- `chromap --version` → `1.0.0` (suite); `chromap --upstream-version` →
  `0.3.3-r519` (engine). The version constant lives in
  [`src/version.h`](../src/version.h) (`CHROMAP_SUITE_VERSION`).
- Releases are cut as annotated `vX.Y.Z` tags; each carries a
  `docs/RELEASE_NOTES_vX.Y.Z.md` and a rolled `CHANGELOG.md` entry.
- Release CI (`.github/workflows/release.yml`) builds and validates the release
  artifacts and publishes them on the tagged commit; the final tagged commit
  hash and artifact checksums are recorded in the GitHub release.
