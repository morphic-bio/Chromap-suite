# Chromap Suite v1.2.0 Release Notes

Date: 2026-09-29

Chromap Suite v1.2.0 runs the low-memory (`--low-mem`) merge of paired-end
barcoded ATAC output on several threads, one reference per task. The output
is byte-identical to the serial merge of v1.1.1:

- the AEV1 sidecar and its `.chroms.tsv`;
- BED and TagAlign rows;
- the summary CSV;
- in-memory peaks and summits;
- the stderr counters.

`chromap --version` reports `1.2.0`. The chromap engine version reported by
`chromap --upstream-version` is unchanged.

## Parallel low-memory ATAC merge

With `--low-mem`, Chromap writes sorted runs of mappings to temporary spill
files, one per reference and flush. After mapping, it merges each
reference's files, removes duplicates and writes the output. v1.1.1 did this
on one thread. On a full-depth DOGMA-plex lane with 32 threads, this step took
about as long as the mapping itself.

v1.2.0 merges and deduplicates each reference in its own task, writes the
task's output to a temporary partition file, and appends the partitions in
reference order.

**Byte identity.** The last group of each reference is emitted exactly as the
serial merge emits it, including the rule the serial merge applies to the
last reference. Summary rows keep their values and their order:
- each task totals its summary counts per barcode;
- when a task ends, it adds the totals to existing rows with one atomic add
  per barcode and field;
- barcodes without a row are inserted afterwards, in the order the serial
  merge inserts them.

**Scope.** Paired-end barcoded ATAC runs with sidecar-only
(`--atac-sidecar-only`), BED or TagAlign output. Dual BAM/CRAM output and bulk
data (no barcodes) keep the serial merge.

**Options.**

| Setting | Meaning |
|---|---|
| `--low-mem-finalize-threads N` (`chromap`, `chromap_lib_runner`) | 0, the default, uses `--num-threads`; 1 keeps the serial merge |
| `MappingParameters::low_mem_finalize_threads` (library) | defaults to 1, so embedding hosts opt in |

**Host permits.** When the host sets permit hooks
(`MappingParameters::permit_acquire_hook` and `permit_release_hook`):
- every merge task runs under a host permit and reports its work (records
  merged, spill bytes, elapsed time) when it releases the permit;
- the serial merge holds one permit for its whole duration;
- concurrency is bounded by the permits the host grants.

**Open files.** A task keeps every spill file of its reference open, so the
number of concurrent tasks is capped by the open-file limit.
- `chromap` and `chromap_lib_runner` raise the soft limit to the hard limit
  before mapping; the library never changes process limits.
- If even one reference needs more files than the limit allows, the run stops
  and suggests `ulimit -n` or a larger `--low-mem-ram`.

**Spill decoding.** Spill records are decoded straight into the fragment
fields, with no per-record payload string or BAM fields.

**Memory.**
- Task partitions go to disk in the spill directory and are removed as they
  are appended.
- A task holds its per-barcode summary totals only while it runs.

## Validation

- **Unit harness.** `make test-lowmem-parallel-finalize` (new, part of the
  release gate) runs synthetic spills through every variant:
  - the serial merge;
  - 2, 7, 32 and 64 threads;
  - the automatic setting;
  - host permit pools of 1 and 3;
  - reduced open-file limits.

  It covers:
  - 400 flushes;
  - duplicates across reference boundaries;
  - identical records in several files;
  - bulk-level duplicate removal;
  - the end-of-stream rule;
  - empty references;
  - more than 255 duplicates;
  - pending summary-table resizes and new summary rows;
  - Tn5 shift off, and no duplicate removal;
  - BED with a translate table, and TagAlign;
  - a failing task.

  Every variant writes the same bytes as the serial merge. The serial outputs
  match outputs saved from the unmodified v1.1.1 library.
- **Output identity with v1.1.1.** 51 comparisons of the parallel merge
  (32 threads, the command-line default, and 1, 2 and 7) against the serial
  v1.1.1 merge, on the same inputs, are byte-identical. All 51 were repeated
  after the per-task summary aggregation was added. They compare the
  sidecar, the fragment text, peaks, summits, the summary and the stderr
  counters. The inputs are:
  - a synthetic scATAC fixture given as 400 lanes (400 spill flushes), with
    sidecar, BED and TagAlign output and with cell- and bulk-level duplicate
    removal;
  - the PBMC 3k 100k fixture, and the same fixture as 256 lanes (25,600 spill
    files);
  - DOGMA-plex lane 1 at 2 million read pairs, including
    `--output-mappings-not-in-whitelist`, and at 50 million read pairs;
  - the full PBMC 3k (82.8 million pairs, 168 flushes);
  - a full DOGMA-plex lane (335 flushes, 40,049 spill files, 316,453,702
    fragments).

  Summaries were compared byte for byte with `--deterministic-mapping`. In
  default mode the candidate-cache columns were left out, since they also
  differ between two runs of v1.1.1. The full-lane sidecar is also
  byte-identical to the v1.1.0 sidecar written by Multiomics Suite 0.9.0 for
  the same lane.
- **Open-file limits.** Under `ulimit -n 1024` the merge ran 2 tasks at a
  time and wrote the same bytes. Under `ulimit -n 256`, below one reference's
  400 spill files, it stopped with the open-file message and left no
  temporary files.
- **Tests.**
  - The 16-target release gate passes.
  - The fixture smokes pass, including the sidecar-only and libchromap smokes
    forced to the serial merge.

## Versioning and provenance

- `chromap --version` reports `1.2.0`. `chromap --upstream-version` is
  unchanged.
- `third_party/rapidmacs` is unchanged at `34df448`.
