# Design: parallel low-memory finalisation for paired-end ATAC (Chromap Suite 1.2.0)

Status: draft for author review, 29 September 2026. Nothing has been implemented.
Branch `design/lowmem-parallel-finalise-20260929` from `master` (`a47f077`, v1.1.0).

## 1. Problem

On a full-depth DOGMA-plex lane 1 on pikachu (32 threads, `--low-mem`, sidecar-only output), Chromap's work after
ATAC mapping runs on one thread. This was measured on 29 September by the Multiomics performance work, with binary
`baf1fef5`, untimed, from the `Log.out` phases:

| Phase | Seconds |
|---|---|
| ATAC mapping | 1,592 |
| Low-memory finalisation: k-way merge, dedup and output, one thread, 335 mid-batch spill flushes | 1,552 |

At 20M read pairs the same step takes 8.8 s, with 8 flushes. The output is 316,453,702 fragments over 159
references. The largest reference is chr1, with 9.6%; the next are chr19 7.3%, chr2 7.0% and chr17 6.2%.

## 2. What already exists, and why it isn't used here

The August mergeable-spill work (`a2f39aa`, `84fd604`; `docs/mergeable_atac_spill.md`) built a parallel finalisation
for the split-and-gather path:
- per-reference tasks run the k-way merge, barcode correction and dedup in parallel;
- they read compact fixed-width records and never construct `SAMMapping` objects;
- each task writes a binary partition, and the partitions are assembled in reference order.

That code can't be called on an ordinary run:
- It reads `ATACMS3`/`ATACMS4` worker spills and `ATACHOT1` companions, and corrects barcodes after the merge.
- Its workers map with the candidate cache disabled and use the complete barcode histogram.
- An ordinary run corrects barcodes during mapping, uses the cache, and learns barcodes from a sample. Routing an
  ordinary run through the materializer would change its output.

The plan reuses the materializer's design (per-reference tasks, partition files, ordered assembly, fixed-width
decoding), not its entry point.

## 3. The current ordinary path

`MappingWriter<MappingRecord>::ProcessAndOutputMappingsInLowMemoryFromOverflow`, `src/mapping_writer.cc:1274`
(called from `src/chromap.h:2182`):

1. Each flush writes one overflow file per reference present in the buffer. The k-way file header carries the
   reference id, so the scan pass is cheap.
2. References are processed one after another in ascending id order. For each reference, a `std::priority_queue`
   merges all of its files. Each step does the following:
   - reads the next record into a `std::string` payload;
   - decodes the payload into a full `AtacSpillRecord`, which holds two `SAMMapping` members;
   - copies the record into and out of the heap.
3. Dedup:
   - Equal records (barcode, start, length) collapse, keeping the last in sort order.
   - With bulk-level dedup for single-cell data, `FindBestMappingIndexFromDuplicates` picks the kept record.
   - After dedup, MAPQ filtering and the Tn5 shift apply.
4. `AppendMapping` in sidecar mode:
   - writes one 24-byte AEV1 record with `fwrite`;
   - pushes a `macs3::FragmentRecord` into `buckets[rid]`.
   - The summary counts (`SummaryMetadata`, a khash keyed by barcode) are updated.

**What makes this parallel.** Duplicates never cross references (the check is `last_rid == min_rid`), so every
reference is independent. AEV1 records are fixed-width, and the header is written at finalize with the record count.
Writing each reference's records to its own partition and concatenating them in reference order therefore gives
the same bytes.

## 4. Design

1. **Per-reference tasks.**
   - Run on a pool of `--num-threads` workers.
   - Schedule longest first, sized by each reference's total spill bytes, which is known from file sizes.
   - Each task runs the current loop unchanged, with its own dedup state.
2. **Task-local outputs:**
   - AEV1 records go to a per-reference partition file in the temp directory.
   - MACS3 fragment records go to `buckets[rid]`. The bucket vector is sized before the parallel region, because
     today's `resize` inside `AppendMapping` would race.
   - Summary updates go to a task-local list in first-seen order.
   - The unique, multi and passing counters are task-local.
3. **Ordered assembly:**
   - Partitions are appended to the sidecar in reference order with `copy_file_range`, then the existing finalize
     writes the header.
   - Summary updates are applied to the global khash reference by reference, in first-seen order.
   - The serial path inserts barcodes in exactly that order, so khash slot order, and with it the summary CSV row
     order, is unchanged.
4. **End-of-stream case.**
   - Today the last record of the last reference goes through a separate block. That block checks MAPQ before
     bulk-level best-duplicate selection; the loop does it the other way round for every other reference.
   - Tasks must reproduce this: loop semantics for every reference except the globally last one.
   - It only matters with bulk-level dedup of single-cell data, but it has to be exact.
5. **Lean decode (same release, second step).**
   - When the spill schema has no BAM pair section (sidecar-only and BED output), decode straight from the reader's
     buffer into `PairedEndMappingWithBarcode` plus the prefix fields.
   - No `std::string` payload is allocated and no `SAMMapping` is constructed. The heap holds small fixed-size
     entries.
   - This is a constant-factor gain on every thread.
6. **File descriptors.**
   - Open files = concurrent tasks × files per reference (about 335 at full depth).
   - Pikachu's limit is effectively unlimited, but other hosts may allow 1,024–4,096.
   - At start, raise the soft limit to the hard limit and cap concurrency to the descriptor budget. If the budget is
     small, fall back to fewer tasks, down to one.
7. **Memory stays bounded:** each task holds one record and one read buffer per file, and partitions go to disk.
8. **Out of scope for 1.2.0:**
   - dual BAM/CRAM output (BAM sort order and writer), which stays serial;
   - bulk and non-ATAC record types;
   - mergeable-spill mode;
   - splitting one reference into ranges. Duplicates share a start, so ranges by start position would be safe, but
     that needs a block index in the spill format. Deferred unless chr1 proves to be the limit.
9. **Option:** `--low-mem-finalize-threads N`, with 1 meaning today's serial path.
   - Default is an author decision; see section 8.
   - Multiomics passes it through as `--chromapAtacLowMemFinalizeThreads`.

## 5. Expected effect

- **This step:** chr1 holds 9.6% of the output, which caps parallelism across references at about 10×. With 32
  threads, the 1,552 s step should drop to roughly 150–250 s, before lean decode. This is an estimate, not a
  measurement.
- **Wall time of a full DOGMA-plex lane on pikachu: small.**
  - Today the ATAC side (mapping plus finalisation) ends at about 3,150 s. STAR's side (RNA mapping, reader-bound,
    then Solo) ends at about 3,270 s.
  - After this change, ATAC would end at about 1,850 s, but the wall stays set by the RNA side. The gain on that run
    is about 2 minutes.
- **Where it helps directly:**
  - Chromap-only runs at depth;
  - hosts where ATAC is the long pole;
  - any later RNA-side improvement, which only shows up in the wall once this is done.

## 6. Gates for Chromap 1.2.0

1. **Unit tests: serial and parallel byte-identical at 1, 2, 7 and 32 threads.** Synthetic spills cover:
   - duplicates on both sides of a reference boundary;
   - equal-key ties across files;
   - bulk-level dedup;
   - a MAPQ-threshold record as the very last record;
   - empty references.
2. **Byte identity against v1.1.0** (same inputs, `--low-mem`) for the AEV1 sidecar, the BED fragments path, the
   summary CSV, embedded peaks and summits, and the stderr counters, on:
   - the test fixtures;
   - PBMC 3k;
   - DOGMA-plex lane-1 2M and 50M subsets;
   - forced-spill runs with a small `--low-mem-ram`, to get hundreds of flushes on small inputs.
3. **One full-depth lane-1 run, untimed**, to confirm identity at 335 flushes.
4. **Multiomics G-M1** with a binary pinned to 1.2.0.
5. **No timing claims** until benchmarking resumes.

## 7. Work plan (one agent, own worktree, runbook and handoff; no push or tag)

0. Measure: pre-dedup record counts and spill files per reference, from a forced-spill run.
1. Move the per-reference loop into a function; serial, byte-identical.
2. Parallel tasks and ordered assembly; unit tests.
3. Lean decode.
4. Gates 1–4.
5. Report. The release (1.2.0 tag, public push) waits for the author.

Estimate: 2–4 days.

## 8. Author decisions

1. Default on (outputs are byte-identical) or opt-in for 1.2.0?
2. Leave dual BAM/CRAM output serial in 1.2.0, as proposed?
3. Start the implementing agent after this note is approved?
