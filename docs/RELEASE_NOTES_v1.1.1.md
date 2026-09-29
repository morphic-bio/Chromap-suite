# Chromap Suite v1.1.1 Release Notes

Date: 2026-09-29 (draft; the date is set when the release is made).

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
  - two runs finish without the missing records;
  - one crashes;
  - one keeps a barcode outside the table over a whitelisted one.

  All five pass on v1.1.1. The test is now part of the release gate, which has
  15 targets.
- **Output identity with v1.1.0.** Results are pending; see the handoff.

## Versioning and provenance

- `chromap --version` reports `1.1.1`. `chromap --upstream-version` is
  unchanged.
- `third_party/rapidmacs` is unchanged at `34df448`.
