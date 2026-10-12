# Chromap Suite v1.2.1 Release Notes

Date: 2026-10-12 (draft; the date is set when the release is made).

Chromap Suite v1.2.1 moves the embedded RapidMACS from `34df448` (v1.0.1 plus
one documentation commit) to v1.0.2. Chromap's own source is unchanged.

- The peaks and summits that `chromap` writes with `--call-macs3-frag-peaks`
  are byte-identical to those of v1.2.0.
- The packaged `rapidmacs` command (also installed as `chromap_callpeaks`) and
  `librapidmacs` carry the three RapidMACS v1.0.2 fixes described below.

`chromap --version` reports `1.2.1`. The chromap engine version reported by
`chromap --upstream-version` is unchanged.

## RapidMACS v1.0.2

RapidMACS v1.0.2 fixes three cases in which v1.0.1 differed from MACS3 3.0.3
`callpeak`. Its changelog and release notes
(`third_party/rapidmacs/CHANGELOG.md` and
`third_party/rapidmacs/docs/RELEASE_NOTES_v1.0.2.md`) describe them in full,
with their validation against MACS3.

### Single-end BED input: the tag total under a numeric `--keep-dup`

With a numeric `--keep-dup`, MACS3 keeps the tag of a chromosome strand that
holds exactly one tag in the pile-up, and leaves it out of the tag total that
scales the genome-wide background.

**v1.0.1** counted that tag in the total.

**v1.0.2** leaves it out, as MACS3 does. Fold enrichment, P and q scores can
change, and with them the results that depend on a threshold. `--keep-dup all`
still counts every tag.

### Single-end BED input: summit scores under `--pvalue`

**v1.0.1** wrote -log10 q in the fifth column of the summits file under either
threshold.

**v1.0.2** writes -log10 P under `--pvalue`, as MACS3 does. Under `--qvalue`
the column is -log10 q, as before.

### Fragment input: a count cap and the placing of summits

`--frag-pileup-max-count N` (`frag_pileup_max_count` in the library) caps the
count of a fragment row, as MACS3's `--max-count N` does.

**v1.0.1** applied the cap to the pile-up and to the background total. The
fragment ends that place summits were still repeated by the count without the
cap.

**v1.0.2** applies the cap to those ends as well. Summit positions, and the
scores reported at them, can change. RapidMACS's reference outputs without a
cap are unchanged.

## Chromap's own peak calls are unchanged

`chromap --call-macs3-frag-peaks`, and the same call through `libchromap`,
reaches RapidMACS only through the fragment pipeline:

- It does not call the single-end BED path, so neither BED fix applies.
- It never sets a count cap.
- The fragment pipeline writes its summits file as before, under
  `--macs3-frag-pvalue` and under `--macs3-frag-qvalue`.

## Validation

All checks ran on one host. Outputs were compared between clean builds of
v1.2.0 and of the same source with RapidMACS moved to v1.0.2, before the
version was set. The tests ran on the release source.

- **Output identity with v1.2.0.** The peak tests and eight further `chromap`
  runs wrote the same bytes before and after the move. That covers all 286
  files that hold peaks, summits, fragments, bedGraphs or BAM records. The
  runs are:
  - the eleven `test-peak-*-100k` targets;
  - `test-peak-integration-matrix-100k`: 12 combinations of peak source,
    workspace and `--low-mem`, with BAM and with BED output (2,047 peaks in
    every one);
  - `test-lowmem-bed-100k`, `test-macs3-bed-callpeak-cli` and
    `test-macs3-frag-qvalue-cli`;
  - eight runs on the PBMC 3k 100k fixture with `--call-macs3-frag-peaks` and
    `--deterministic-mapping`: P and q thresholds, barcoded and bulk input,
    cell-level and bulk-level duplicate removal, the sweep workspace,
    `--low-mem` and `--macs3-frag-no-uint8-counts`.

  454 files were compared. The others that differ are 41 logs and timing
  files, and 18 mapping summaries that are equal once the four
  candidate-cache columns are removed. Those columns also differ between two
  runs of one build.
- **The release source.** The integration matrix, run again on the release
  source, wrote the same 12 narrowPeak and 12 summits files as v1.2.0.
- **Tests.** The 16-target release gate, `make test-smoke` and
  `make test-unit` pass.
- **RapidMACS regression tests.** RapidMACS's two new tests pass with the
  packaged `rapidmacs` command. They compare complete narrowPeak and summits
  files with RapidMACS's stored MACS3 3.0.3 outputs, byte for byte:
  - single-end BED input with strands that hold one tag, under P and q
    thresholds and `--keep-dup` 1, 2 and `all`;
  - fragment input with count caps 1, 2 and 5, and without a cap.

## Versioning and provenance

- `chromap --version` reports `1.2.1`. `chromap --upstream-version` is
  unchanged.
- `third_party/rapidmacs` moves from `34df448` to
  `1021e0eb0a781dbcccf00ebc69b72caf5d625dc8` (RapidMACS v1.0.2).
- Programs that link `librapidmacs` must be rebuilt and relinked. The shared
  narrowPeak and summits writer has one more optional argument; existing
  source calls compile unchanged.
