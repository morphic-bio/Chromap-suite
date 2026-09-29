# Chromap Suite v1.1.1 Release Notes

Date: 2026-09-29

Chromap Suite v1.1.1 fixes two errors in the low-memory (`--low-mem`) merge of
spilled mappings. When neither error applies, the output is byte-identical to
v1.1.0.

`chromap --version` reports `1.1.1`. The chromap engine version reported by
`chromap --upstream-version` is unchanged.

## Fixes

### A spill file that cannot be opened stops the run

In low-memory mode, Chromap writes sorted runs of mappings to temporary spill
files, one per reference and flush. It then merges each reference's files.

**v1.1.0:** a spill file that could not be opened was skipped, both while the
files were grouped by reference and while a reference was merged. The run then
finished normally, without that file's records.

The likeliest cause is the open-file limit. The merge keeps all spill files of
one reference open at once, and a reference gets one file per flush.

**v1.1.1:** the run stops with an error that names the file and the reason.

- If the open-file limit is the cause, the message gives the number of files
  for the reference and the limit. It suggests raising the limit
  (`ulimit -n`) or using a larger `--low-mem-ram`, so that fewer spill files
  are written.
- The spill files are removed before exiting.

The fix covers every output format that uses the low-memory merge.

### Bulk-level duplicate removal stays inside the whitelist table

With `--remove-pcr-duplicates` on barcoded data, duplicates are removed at bulk
level unless `--remove-pcr-duplicates-at-cell-level` or `--preset atac` is
given. At a duplicated position, the low-memory merge keeps one barcode. It
chooses by duplicate count first, then by the barcode's abundance in the
whitelist table.

**v1.1.0** looked up that abundance without checking that the barcode was in
the table:

- With `--output-mappings-not-in-whitelist`, a barcode outside the whitelist
  read a value past the end of the table's storage, so the choice was
  undefined.
- With barcodes but no `--barcode-whitelist`, the table is empty, and the
  lookup crashed the run.

**v1.1.1:** a barcode that is not in the table has abundance 0. The other
selection rules are unchanged:
- more duplicates wins;
- then higher abundance;
- then the first in sort order.

Runs whose barcodes are all in the table are unaffected. That covers every run
with a whitelist and without `--output-mappings-not-in-whitelist`.

## Validation

- **New regression test.** `make test-lowmem-overflow-edge-cases` is hermetic
  and runs each case in a separate process. The cases are:
  - a spill file removed before the merge;
  - an open-file limit below one reference's spill files;
  - the same spills under the normal limit;
  - bulk-level duplicate removal with an empty whitelist table;
  - bulk-level duplicate removal with barcodes outside the table.

  Linked against v1.1.0, four of the five cases fail:
  - with one of six spill files removed, the run writes 5 of 6 records and
    exits normally;
  - with an open-file limit of 24 and 40 spill files for one reference, the
    run writes 20 of 40 records and exits normally;
  - with an empty whitelist table, the run crashes (segmentation fault);
  - with barcodes outside the table, it keeps a barcode outside the table
    over a whitelisted barcode with abundance 7.

  All five pass on v1.1.1. Under valgrind, v1.1.0 shows invalid reads in the
  whitelist lookup and v1.1.1 shows none. The test is now part of the release
  gate, which has 15 targets and passes.
- **Command-line runs** on the small synthetic scATAC fixture from the
  sidecar-only smoke test:
  - Barcodes, no whitelist, bulk-level duplicate removal and `--low-mem`:
    v1.1.0 crashes; v1.1.1 writes 399 fragments.
  - The same reads given as 100 lanes, with `--low-mem-ram 1K` (100 spill
    files per reference) and `ulimit -n 64`:
    - v1.1.0 exits normally with the same 3,241 fragments as a run without the
      limit, but 3,082 of them carry lower duplicate counts: 303,921 in total,
      instead of 460,645.
    - v1.1.1 stops with the open-file message and leaves no temporary files.
- **Output identity with v1.1.0.** Twelve comparisons are byte-identical: the
  AEV1 sidecar and `.chroms.tsv`, the fragment text, the peaks, the summits,
  the summary and the stderr counters. The cases are:
  - the synthetic fixture as 100 lanes, sidecar output, with cell-level and
    with bulk-level duplicate removal (100 spill flushes each);
  - the PBMC 3k 100k fixture:
    - sidecar and BED output, each at the default `--low-mem-ram` and at `1K`;
    - sidecar output with bulk-level duplicate removal;
    - paired-end BED without barcodes;
    - single-end barcoded BED;
  - DOGMA-plex lane 1 at 2,000,000 read pairs, sidecar and BED output at
    `--low-mem-ram 1K`.

  These eleven runs used `--deterministic-mapping`.

  The twelfth comparison ran DOGMA-plex 2M sidecar output with the candidate
  cache on:
  - Its summary is compared without the four candidate-cache columns. Those
    columns differ between two v1.1.0 runs as well.
  - The MACS3 peak-metrics file records each run's own output paths, so it is
    compared after replacing the run directory; no other byte differs.

## Versioning and provenance

- `chromap --version` reports `1.1.1`. `chromap --upstream-version` is
  unchanged.
- `third_party/rapidmacs` is unchanged at `34df448`.
