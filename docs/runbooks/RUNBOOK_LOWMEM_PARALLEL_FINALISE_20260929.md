# Runbook: parallel low-memory finalisation for paired-end ATAC (Chromap Suite 1.2.0), 2026-09-29

## Status

- **Author review:** 29 September. The decisions are recorded in §7, and the
  TODO list is in §7.1.
- **Order of work:**
  1. **M-1**, a Chromap Suite 1.1.1 patch for two v1.1.0 bugs, on its own
     branch and worktree (§4.0).
  2. M0-M4 for 1.2.0. M0 starts only when the coordinator says so.
- **Where the documents are:** the design note
  (`docs/design/LOWMEM_PARALLEL_FINALISE_20260929.md`), this runbook and the
  handoff are on branch `design/lowmem-parallel-finalise-20260929`, from
  `master` `a47f077` (v1.1.0).
- **Read the design note first.** This runbook adds the code-level plan, the
  corrections found while checking the note against the code, and the
  validation plan.

## 0. Rules

From the author for this work:

- No push, tag or release. The author approves the 1.1.1 release and the 1.2.0
  release separately.
- Commits on local branches are approved. Use plain messages, with no AI
  attribution and no co-author trailers.
- Update this runbook and the handoff
  (`docs/handoffs/HANDOFF_LOWMEM_PARALLEL_FINALISE_20260929.md`) at each
  milestone.
- Near 90% of the usage limit, stop, write the handoff and report.
- Stop and report after M-1. Do not start M0 until the coordinator says so.

Standing project rules:

- **Clean room.** Never read 10x Genomics code. From
  `/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0`, use only
  `fasta/genome.fa`, as the existing tests do. Use the ATAC barcode whitelist
  as reference data. Do not open other files in that package.
- **Excluded material.** Follow the exclusions in the private coordinator
  handoff; never read or copy the directories it names.
- **Spill reader.** `OverflowReader` (the reader for the low-memory spill
  temp files) may change in M3.
- **No timing.** No benchmarks or timing claims. Runs are for output identity
  only. Do not record, compare or quote wall times. The comparisons ignore
  timing lines in stderr.
- **Worktrees.** Work only in the worktrees this runbook names: this design
  worktree and the M-1 worktree (§4.0). Do not modify other agents'
  worktrees. Reading from them is allowed where this runbook says so.
- **Running jobs:**
  - Builds and tests run at `nice -n 10`.
  - A timed Multiomics speed check (v0.10.0, G-M2) may be running on pikachu.
    The author says contention does not matter for it, so untimed builds and
    identity runs are fine. Never run anything timed.
  - Large identity runs (P3, D2, F1) take
    `flock /mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock`, so they do
    not overlap someone else's timed run. They are still untimed.
  - Run large jobs one at a time (`AGENTS.md`, benchmarking policy).
- **Artifacts.** Keep generated artifacts out of tracked paths. Point
  `CHROMAP_ARTIFACT_ROOT` into a validation root.

## 1. Goal and non-goals

Goal: in a `--low-mem` run, the work after mapping (k-way merge, dedup,
output) runs as one task per reference on several threads, for paired-end
ATAC `AtacSpillRecord`. Outputs stay byte-identical to v1.1.0. In an embedded
run, every task runs under the host's permits.

In scope for 1.2.0:

- `MappingWriter<AtacSpillRecord>::ProcessAndOutputMappingsInLowMemoryFromOverflow`
  for single-cell (barcoded) runs, when the output is one of:
  - sidecar-only (`--atac-sidecar-only`: AEV1 sidecar plus `.chroms.tsv`);
  - barcoded BED fragment text;
  - barcoded TagAlign text (D6).
- The side outputs of those runs:
  - the `--summary` CSV;
  - the in-memory MACS3 fragment buckets (peaks and summits);
  - the stderr counters.

Non-goals for 1.2.0 (unchanged by construction; see the TODO list, §7.1):

- **Dual BAM/CRAM output** (`AtacDualFragmentAndBam()`) stays on the v1.1.0
  serial loop (D2).
- **Bulk data** (no barcodes, `is_bulk_data`) stays on the v1.1.0 loop (D6).
- **Non-ATAC record types.** The generic
  `ProcessAndOutputMappingsInLowMemoryFromOverflow` template keeps its body.
  It serves `SAMMapping`, `PAFMapping`, the BED types, `PairsMapping` and
  `PairedPAFMapping`.
- **Mergeable spill.** `WriteMergeableAtacSpillFromOverflow` and the
  materializer (`src/atac_spill_materializer.cc`) are unchanged.
- **Everything else:** the in-memory (non-`--low-mem`) path, the spill writer
  (`src/overflow_writer.*`), the spill format, and splitting one reference
  into ranges (design §4.8).

## 2. Locations

| What | Path |
|---|---|
| Design worktree (branch `design/lowmem-parallel-finalise-20260929`) | `/mnt/pikachu/Chromap-suite-lowmem-parallel-20260929` (below `W`) |
| M-1 worktree (branch `fix/v1.1.1-lowmem-edge-cases`, from `v1.1.0`) | `/mnt/pikachu/Chromap-suite-v111-fix-20260929` (below `W111`) |
| M-1 validation root | `/mnt/pikachu/chromap_v111_validation_20260929` (below `V111`) |
| 1.2.0 validation root (created in M0) | `/mnt/pikachu/lowmem_parallel_validation_20260929` (below `V`) |
| Main checkout (do not modify) | `/mnt/pikachu/Chromap-suite` |
| RapidMACS clone source, pinned `34df448` (read only) | `/mnt/pikachu/Chromap-suite/third_party/rapidmacs` |
| Index (GRCh38 ARC) | `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/genome.index` |
| Reference FASTA | `/mnt/pikachu/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/fasta/genome.fa` |
| ATAC whitelist | `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/737K-arc-v1_atac.txt` |
| ATAC-to-GEX translate table | `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/chromap_index/atac2gex.tsv` |
| PBMC 3k 100k fixture (4 lanes, gzip) | `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac` |
| PBMC 3k full depth (4 lanes, BGZF, 82,781,574 pairs) | `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/extracted/pbmc_unsorted_3k/atac` |
| DOGMA-plex lane 1, 2M (G-M1 reference inputs; read only) | `/mnt/pikachu/single_binary_reference_20260928/dogmaplex_lane1_2m_sidecar/run1/input_fastqs/atac/DP01_ATAC_S1_L001_R{1,2,3}_001.fastq.gz` |
| DOGMA-plex lane 1, 50M (gzip; read only) | `/mnt/pikachu/fastq_mirror_validation_20260928/lane1_50M/R{1,2,3}.fastq.gz` |
| DOGMA-plex lane 1, full depth (read only) | `/mnt/pikachu/dogmaplex_gse309834/production_20260926/acquisition/fastq/lane_01/atac/DP01_ATAC_S1_L001_R{1,2,3}_001.fastq.gz` |
| Lane-1 full-depth v1.1.0 output (read only; Multiomics 0.9.0 embeds Chromap `a47f077`) | `/mnt/pikachu/perf_lane_throughput_20260929/full_baseline/out/lane01/{atac_fragments.bin,atac_fragments.bin.chroms.tsv,chromap_summary.csv}` and `../../phases.json` |
| Earlier identity scripts (templates; read only) | `/mnt/pikachu/fastq_mirror_validation_20260928/run_identity.sh`, `run_identity_det.sh`, `compare_identity.py` |
| Timed-run lock | `/mnt/pikachu/e2e_bench_20260926/pikachu_timed.lock` |

Host facts, checked read-only on 29 September:

- 32 cores, 125 GiB RAM, glibc 2.35.
- `ulimit -Sn` and `-Hn` are both 1,073,741,816.
- `/mnt/pikachu` is ext4 with a 4 KiB I/O block, and 3.7 TB free.

DOGMA-plex read roles: R1 and R3 are the mates, and R2 carries the barcode
(`--read-format bc:8:23:-`).

## 3. Corrections to the design note

Each claim in the note was checked against the code at `a47f077`. Claims not
listed here hold. The lane-1 numbers in §1 match `phases.json` (1,592.36 s,
1,552.51 s, 335 flushes). The AEV1 header reports 316,453,702 records.

1. **Forced spills (gate 6.2).** Hundreds of flushes on a small input can't
   come from `--low-mem-ram` alone.
   - A batch is 500,000 pairs (`read_batch_size_`, `src/chromap.h:227`).
   - The flush check runs once per batch, in the per-batch output task
     (`src/chromap.h:2021-2053`).
   - Each input lane starts its own batches (`src/chromap.h:1099`).
   - So flushes ≤ batches: at most 4 for 2M, at most 100 for 50M, about 166
     for PBMC 3k at full depth.
   - For hundreds on small inputs, pass the same small lane many times with
     `--low-mem-ram 1K`, so every batch flushes (§5.4).
2. **Summary CSV identity (gate 6.2).** Byte identity holds only with
   `--deterministic-mapping`.
   - Without it, `cachehit`, `fric`, `estfrip` and `numcacheslots` vary
     between two runs of the same binary. The v1.1.0 reader validation showed
     this.
   - In default mode, compare with those four columns removed, but only after
     a v1.1.0-vs-v1.1.0 control confirms that policy.
3. **Pre-sizing the MACS3 buckets (§4.2) is already done.**
   - `OutputHeader` for `AtacSpillRecord` assigns `num_reference_sequences`
     empty buckets in both branches (`src/mapping_writer.cc:2118-2121`,
     `2198-2201`).
   - So the `resize` in `AppendMapping` (`2249-2251`) does not run today.
   - Keep one serial size check before the parallel region. Tasks never
     resize.
4. **Summary replay (§4.2-4.3): replaced by the D8 variant.**
   - `kh_put` checks for a resize before its lookup, even for a key already
     present (`src/khash.h:309-320`). A threshold crossing is acted on at the
     next put, whatever its key.
   - So the table layout depends on two things: the sequence of new-key
     insertions, and a put following each insertion.
   - Existing-key puts only move a resize earlier. That has no effect before
     the next insertion.
   - Per-reference delta lists would have held about 0.4-0.8 GB at lane-1 full
     depth. That bound comes from the 691,757-row summary: at most 11M-25M
     (reference, barcode) entries.
   - D8 replaces them with atomic adds into the existing table and an ordered
     log of new barcodes only (§4.3 item 2).
5. **Appending with `copy_file_range` (§4.3) is incomplete.**
   - The sidecar and text outputs are stdio `FILE*`s. A direct descriptor copy
     needs `fflush` before and `fseeko(fp, 0, SEEK_END)` after.
   - It returns `EXDEV` across filesystems on some kernels, and `EINVAL` for a
     pipe.
   - The plan uses a plain `fread`/`fwrite` copy. `copy_file_range` is a TODO
     (§7.1).
6. **Partitions must carry text as well as AEV1.**
   - The non-sidecar outputs on this path are text: BED rows through
     `mapping_output_file_`, and TagAlign rows
     (`src/mapping_writer.cc:2230-2285`).
   - Barcoded TagAlign with `--low-mem` also uses `AtacSpillRecord`
     (`src/chromap_driver.cc:1715-1728`, `src/libchromap.cc:101-117`). It is
     now in scope (D6).
7. **Header count.**
   - The AEV1 header count comes from `atac_evidence_records_written_`.
     `FinalizeAtacEvidenceBinaryOutput` writes it, called from
     `~MappingWriter` (`src/mapping_writer.h:137-138`, `406-438`).
   - Assembly must add each partition's record count.
8. **The file-descriptor budget (§4.6) is incomplete.**
   - Files per reference = mid-batch flushes plus the final drain, about 336
     at lane-1 full depth.
   - The budget must also cover partitions, the sidecar and an embedding
     host's own files.
   - Only the CLI raises the limit (D7).
   - v1.1.0 silently skips a spill file that fails to open. M-1 fixes that
     (§4.0).
9. **Memory (§4.7) is incomplete.**
   - `OverflowReader` calls `setvbuf(file_, nullptr, _IOFBF, 8 MiB)`
     (`src/overflow_reader.cc:14`). With a NULL buffer, glibc ignores the size
     and uses its default buffer, which is `st_blksize`-based (4 KiB here).
     Confirm with `strace -e trace=read` in M0.
   - Do not raise it to 8 MiB per file in the parallel path: 32 × 336 × 8 MiB
     is about 84 GiB.
   - The MACS3 buckets (16 bytes per output fragment, about 5 GB at lane-1
     full depth) are unchanged.
10. **The scan pass is serial over every spill file (§3.1).** It is cheap per
    file, but there are tens of thousands of files at full depth.
    - A header probe gives the rid, schema and size without stdio: `open`,
      `pread` of the 32-byte `AtacKwaySpillFileHeaderV1`, `fstat`, `close`.
11. **End of stream (§4.4): correct as stated.** In detail:
    - The two blocks are the loop block (`src/mapping_writer.cc:1459-1495`)
      and the end block (`1543-1578`).
    - At end of stream, the MAPQ gate tests the last record in sort order, not
      the bulk-selected best.
    - If the gate fails, LOWMAPQ and MAPPED are credited to that record's
      barcode, and the uni/multi counter uses that record.
    - If it passes, the best record is emitted even if its own MAPQ is below
      the threshold.
    - Bulk-level dedup of single-cell data is the CLI default when
      `--remove-pcr-duplicates` is given without
      `--remove-pcr-duplicates-at-cell-level` and without `--preset atac`
      (`src/mapping_parameters.h:105`, `src/chromap_driver.cc:603-608, 819-827`).
      Multiomics uses cell-level dedup.
    - The "globally last" reference is the highest rid with a spill file.
      Every k-way spill file holds at least one record.
12. **A task can't run "the current loop unchanged" (§4.1).**
    - The loop carries `last_rid`, `last_mapping` and the global
      `is_first_iteration` across references.
    - It emits a reference's final group only when the next reference's first
      record arrives.
    - So each task returns its final group as a tail for the assembler (§4.3).
13. **Worker threads must not exit.**
    - The v1.1.0 loop calls `ExitWithMessage` (`exit(-1)`).
    - Calling `exit` from a worker inside an OpenMP region is unsafe.
    - Tasks return errors instead, as `MaterializeHotPartitions` does
      (`src/atac_spill_materializer.cc:110-480`).
14. **Embedded concurrency was not covered.**
    - Multiomics runs Chromap alongside STAR under a permit allocator
      (`MappingParameters::permit_acquire_hook` and `permit_release_hook`,
      `src/mapping_parameters.h:268-284`).
    - The mapping loop releases all permits before finalisation
      (`src/chromap.h:2067-2074`).
    - The author requires permits for every task (D4, §4.3 item 5).
15. **"Bulk out of scope" needs an explicit rule.**
    - Sidecar-only output for bulk data also uses `AtacSpillRecord`
      (`kAtacSpillSchemaIsBulk`).
    - The eligibility check excludes it (D6).
16. **Standalone peaks.**
    - "Embedded peaks and summits" needs
      `--call-macs3-frag-peaks --macs3-frag-peaks-source memory`.
    - The default source is `file`, and sidecar-only requires `memory`.
17. **A v1.1.0 full-depth output already exists (gate 6.3).**
    - Multiomics 0.9.0 pins Chromap `a47f077`.
    - The existing lane-1 full-depth output is therefore v1.1.0 output: a
      cross-check, not the gate. The embedded configuration is not the
      standalone one.
18. **Not re-verified:**
    - "159 references" with fragments (the AEV1 header lists 194);
    - the per-reference shares;
    - the 20M-pair figure.

    M0 measures the first two. §5 stays an unmeasured estimate.
19. **The materializer (§2) is consistent with `docs/mergeable_atac_spill.md`
    and the code.**
    - Its tasks run in ascending rid order.
    - It assembles after the region with `fread`.
    - It leaves the summary to the `ATACMS3` parent.

## 4. Code-level plan

Names are suggestions; the structure is what matters.

### 4.0 M-1: Chromap Suite 1.1.1 (before M0)

Two v1.1.0 bugs, fixed on their own branch before any 1.2.0 work:

- **Silent skip.**
  - The low-memory merge skips a spill file it cannot open, in both the scan
    and the merge (`src/mapping_writer.cc:1323-1326`, `1386-1389`).
  - A shortage of file descriptors therefore drops records without an error.
  - The code is in the generic template, so every record type is affected.
- **Out-of-bounds read.**
  - `FindBestMappingIndexFromDuplicates` (`src/mapping_writer.h:530-568`)
    reads `kh_value` at the iterator `kh_get` returns.
  - For a barcode not in the whitelist table, that iterator is `kh_end`, one
    element past the value array.
  - Reachable with bulk-level dedup and `--output-mappings-not-in-whitelist`.
  - Also reachable with barcodes and no `--barcode-whitelist`: the table is
    empty, and the read goes through a NULL value array.

Plan:

1. **Worktree.**
   - `git -C /mnt/pikachu/Chromap-suite worktree add -b fix/v1.1.1-lowmem-edge-cases /mnt/pikachu/Chromap-suite-v111-fix-20260929 v1.1.0`.
   - Initialise RapidMACS from the local clone as in §5.2.
   - Build the v1.1.0 baseline there before any edit, and archive it in
     `V111/bin/` with its sha256.
2. **Fix A: unopenable spill files stop the run.**
   - In `ProcessAndOutputMappingsInLowMemoryFromOverflow` (generic template),
     a spill file that cannot be opened stops the run with its path and
     `strerror(errno)`, in both places.
   - For `EMFILE`/`ENFILE`, the message also gives the number of files for
     the reference and suggests `ulimit -n` or a larger `--low-mem-ram`.
   - The spill files are unlinked first.
   - Successful runs are unchanged.
3. **Fix B: defined abundance for absent barcodes.**
   - In `FindBestMappingIndexFromDuplicates`, a barcode missing from the
     whitelist table has abundance 0. Nothing is read at `kh_end`.
   - Selection is otherwise unchanged: more duplicates first, then higher
     abundance, then the first in sort order.
   - When every barcode is present (every run with a whitelist and without
     `--output-mappings-not-in-whitelist`), the result is the same.
4. **Regression test** `tests/test_lowmem_overflow_edge_cases.cc`
   (`make test-lowmem-overflow-edge-cases`).
   - It links `libchromap` and needs no fixture. It drives
     `MappingWriter<AtacSpillRecord>` in sidecar-only mode through its public
     low-memory API, with synthetic spills. Each case runs in a forked child,
     so the exit status and stderr can be checked.
   - Cases:
     - a spill file deleted before finalisation: must exit non-zero with the
       message;
     - `RLIMIT_NOFILE` below the files for one reference: must exit non-zero
       with the message;
     - the same with a normal limit: exits 0, with every expected record;
     - bulk-level dedup with an empty whitelist table: v1.1.0 crashes; the fix
       exits 0 with the defined selection;
     - bulk-level dedup with barcodes absent from a non-empty table: the fix
       gives the defined selection.
   - It must fail when linked against v1.1.0 and pass after the fixes.
   - It joins the hermetic release gate as the 15th target. Update the counts
     in `AGENTS.md`, `README.md` and `docs/releasing.md`, and list the target
     in `tests/README.md`.
5. **Identity against v1.1.0** on normal cases (§5.8). A valgrind check of the
   absent-barcode case is evidence only.
6. **Release files.**
   - Version 1.1.1 in `src/version.h`, the Dockerfile `CHROMAP_SUITE_VERSION`
     `ARG`, `debian/changelog`, and the manual-release default in
     `.github/workflows/release.yml`.
   - A `CHANGELOG.md` `[1.1.1]` entry and `docs/RELEASE_NOTES_v1.1.1.md`.
   - The Dockerfile's default source revision is set when the release is
     finalised, as for 1.1.0.
7. **Commit** on the fix branch. Stop and report. No tag, push or merge.

1.2.0 work then starts from the 1.1.1 commit, merged or rebased as the author
directs, so the silent skip no longer exists anywhere.

### 4.1 Entry point and dispatch

- **Keep the v1.1.0 loop verbatim (D5).** Move the body of the generic template
  `MappingWriter<MappingRecord>::ProcessAndOutputMappingsInLowMemoryFromOverflow`
  (as fixed in 1.1.1), unchanged, into a private member template
  `ProcessLowMemOverflowSerial(...)`. The generic public function calls it, so
  every record type other than `AtacSpillRecord` runs exactly that code.
- **Add an explicit specialisation** for
  `MappingWriter<AtacSpillRecord>::ProcessAndOutputMappingsInLowMemoryFromOverflow`.
  - Declare it in `src/mapping_writer.h` next to the other `AtacSpillRecord`
    specialisations, and define it in `src/mapping_writer.cc`.
  - Remove the explicit instantiation for `AtacSpillRecord` at
    `src/mapping_writer.cc:1614-1615`.
- **Dispatch.** The specialisation calls `ProcessLowMemOverflowSerial` unless
  all of these hold:
  - resolved workers N ≥ 2;
  - `!AtacDualFragmentAndBam()`;
  - `!CreatesMergeableAtacSpill()`;
  - `!is_bulk_data`;
  - output format BED or TagAlign;
  - every spill file has a k-way header.
- **Option (D1, D3):** `MappingParameters::low_mem_finalize_threads`.
  - 0 means auto (`num_threads`); 1 means the v1.1.0 loop; N ≥ 2 means
    per-reference tasks.
  - The library field defaults to 1, so an embedding host opts in.
  - The CLI (`chromap` and `chromap_lib_runner`) defaults to 0, through
    `--low-mem-finalize-threads N`.
- **Unchanged:** the call site `src/chromap.h:2182` and its signature.
- **C++11 build** (`Makefile:18`): no `if constexpr`, no `std::atomic_ref`, no
  `std::make_unique`. Use the GCC `__atomic` builtins.

### 4.2 M1: per-reference function, serial and byte-identical

1. **Probe the spills (correction 10).** New `src/atac_lowmem_finalize.{h,cc}`,
   added to `core_cpp_source` in the `Makefile`.
   - `ProbeAtacKwaySpillFile(path, &rid, &schema_mask, &bytes, &error)`
     validates the same header fields as
     `OverflowReader::ConsumeAtacSpillFilePrefixIfPresent`.
   - Build `std::vector<AtacRefJob>`: rid, file indices ascending by global
     index into `shared_overflow_file_paths_` (this keeps the heap tie-break),
     and total bytes.
   - Enforce one schema mask across all files, with v1.1.0's message.
   - Print the two v1.1.0 lines unchanged and in the same order: "Processing N
     overflow files for k-way merge" and "Low-memory overflow: mid-batch flush
     count: N".
2. **The per-reference function:**
   `bool MergeDedupAtacReference(const AtacRefJob&, const khash_t(k64_seq)* wl, AtacRefSink*, AtacRefResult*)`.
   - Copy the per-rid loop (`src/mapping_writer.cc:1377-1540`), with its own
     dedup state and a per-reference first-record flag.
   - It opens its own readers. A reader that fails to open is an error.
   - It emits every group except the last, and returns the last as
     `AtacRefResult::tail`: `last_mapping`, `num_last_mapping_dups` and the
     bulk-dedup list.
   - It asserts that the reference yielded at least one record.
   - It reports errors through `result->error`, never through
     `ExitWithMessage`.
3. **Group emission:** `EmitAtacGroup(rid, reference, state, bool end_of_stream, sink, counters)`
   reproduces the two v1.1.0 blocks exactly:
   - `false` is the loop block: best-duplicate selection, then the MAPQ gate;
   - `true` is the end block: the MAPQ gate, then best-duplicate selection.
4. **Serial driver (direct sink).** For each rid in ascending order:
   - emit the pending tail of the previous non-empty reference with
     `end_of_stream = false`;
   - merge this reference: groups go to `AppendMapping` and
     `summary_metadata_.UpdateCount` as in v1.1.0;
   - keep this reference's tail pending.

   After the loop:
   - emit the pending tail with `end_of_stream = true`;
   - unlink the spills and clear the path list;
   - print the three trailing lines as v1.1.0 does.

   In M1, any N ≥ 2 selects this path but runs it on one thread.
5. **Check:** the M1 gates (§6).

### 4.3 M2: task-local outputs, ordered assembly, permits, the pool

1. **Extract the non-dual branch** of `AppendMapping` for `AtacSpillRecord`
   (`src/mapping_writer.cc:2209-2287`) into
   `AppendAtacNonDualFragment(const AtacFragmentSink&, rid, reference, const PairedEndMappingWithBarcode&)`.
   - The sink holds the sidecar `FILE*` and its record count, the text
     `FILE*`, and the bucket pointer.
   - `AppendMapping` calls it with the member streams.
   - The task sink never resizes buckets; it asserts `rid < buckets.size()`.
   - These are read-only and safe to call concurrently:
     `FindBestMappingIndexFromDuplicates`, `BarcodeTranslator::Translate`,
     `SequenceBatch::GetSequenceNameAt` and the whitelist khash.
2. **Summary: the exact low-memory variant (D8).**

   New `SummaryMetadata` methods (`src/summary_metadata.h`):
   - `bool ResizePendingOnNextPut() const`: `n_occupied >= upper_bound`.
   - `void ApplyPendingResizeLikePut()`: the same resize step `kh_put` takes
     before its lookup. Copy both branches from `src/khash.h:312-320`.
   - `khiter_t FindExisting(uint64_t)`: `kh_get`.
   - `void AddExistingAtomic(khiter_t, int type, uint64_t)`:
     `__atomic_fetch_add(&kh_value(h, k).counts[type], change, __ATOMIC_RELAXED)`.

   How the counts are handled:
   - **Before the region.** Act only when `--summary` is set and at least one
     spill record exists. If a resize is pending, apply it: that is exactly
     what the first serial `UpdateCount` would do before its own lookup.
   - **During the region.** Nothing is inserted.
     - For a barcode already in the table, the task adds its counts
       atomically.
     - For a barcode not in the table, the task logs it, with its DUP, LOWMAPQ
       and MAPPED deltas, in a task-local list: first-seen order, one entry
       per barcode, aggregated.
     - The log is usually empty. New keys appear only when mapped barcodes
       were never credited TOTAL during mapping, for example with
       `--output-mappings-not-in-whitelist`.
   - **At assembly**, the log is replayed in the ordered sequence (item 7).
     Each barcode gets three `UpdateCount` calls (DUP, LOWMAPQ, MAPPED), so a
     put always follows each insertion.

   Why it is exact:
   - The serial path's layout is set by its sequence of new-key insertions
     (global first-seen order, with the tails in their positions), with a
     resize at the first put after each threshold crossing.
   - This variant inserts the same keys in the same order.
   - It applies the pending start resize at the same point relative to those
     insertions.
   - It keeps a put after every insertion.
   - Counts are integer sums.

   Verify with byte-identical summaries under `--deterministic-mapping` (§5).
   Multiomics plans to make its summary opt-in; that is in its own repository.
3. **Task sink:**
   - a partition file (AEV1 records, or text rows in the same format);
   - a pointer to `buckets[rid]`;
   - its own counters;
   - its new-barcode log.
4. **Temp partitions:**
   - Directory: `dirname(shared_overflow_file_paths_[0])`, the resolved spill
     directory.
   - Name: `mkstemp("<dir>/chromap_lowmem_part_<pid>_<rid>_XXXXXX")`, one per
     reference with output. Each is closed at the end of its task.
   - When a task succeeds, unlink its reference's spill files.
   - At assembly, unlink each partition right after appending it.
   - On any error, unlink every partition and every remaining spill file, then
     `ExitWithMessage` with the first error in rid order.
5. **Permits (D4, applied to every task).** When
   `MappingParameters::PermitHooksEnabled()`:
   - Before a per-reference task opens any file, it calls
     `wait_ns = permit_acquire_hook(permit_hook_ctx)` and records the start
     time.
   - When it finishes, it calls
     `permit_release_hook(ctx, wait_ns, work_units, work_bytes, work_ns)`,
     the shape the PE mapping loop uses (`src/chromap.h:1379-1416`):
     - work units = spill records merged;
     - work bytes = that reference's spill bytes;
     - work ns = the task's elapsed time, as telemetry for the host. It is
       not recorded or reported.
   - The release also happens when the task fails.
   - Concurrency is bounded by three things:
     - the permits the host grants (a worker blocks in `acquire`);
     - the resolved workers N;
     - the descriptor budget.
   - With no hooks (standalone), N workers run. No pool runs outside the
     host's permits.
   - The acquire happens outside any lock.
6. **File-descriptor budget (D7):**
   `ResolveLowMemFinalizeWorkers(requested, max_files_per_ref, &note)`.
   - The CLI (`chromap_driver.cc`, `chromap_lib_runner.cc`) raises the soft
     `RLIMIT_NOFILE` to the hard limit at start. The library never changes
     limits.
   - `open_now` = entries in `/proc/self/fd`, or 64 if that is unreadable.
   - `budget = soft - open_now - max(64, open_now)`.
   - `workers = min(requested, non-empty references, budget / (max_files_per_ref + 2))`.
   - If that is 0, clean up and fail with a message. The message names the
     files needed, the limit and the remedies (`ulimit -n`, a larger
     `--low-mem-ram`).
7. **Scheduling and ordered assembly.**

   Scheduling:
   - Sort jobs by spill bytes descending, then rid ascending.
   - Run them under
     `#pragma omp parallel for schedule(dynamic, 1) num_threads(workers)`,
     with an `std::atomic<bool> failed`.
   - Results go in a vector indexed by rid.

   Assembly (serial, after the region, rid ascending). For each rid with a
   job:
   1. Emit the pending tail (`end_of_stream = false`) through the direct sink.
   2. Append the partition to `atac_evidence_fp_` or `mapping_output_file_`
      with a 4 MiB `fread`/`fwrite` loop. Add its count to
      `atac_evidence_records_written_`, then unlink it.
   3. Replay its new-barcode log (item 2).
   4. Add its counters.
   5. Make its tail the pending tail.

   After the loop:
   - Emit the pending tail with `end_of_stream = true`.
   - Print the v1.1.0 trailing lines unchanged.
   - When N ≥ 2, add one line (D9), for example "Low-memory finalization: R
     reference tasks on T threads".
8. **Option plumbing:**
   - `src/mapping_parameters.h`;
   - `src/chromap_driver.cc`: option list near line 88, parsing near line 899;
     reject N < 0;
   - `src/chromap_lib_runner.cc`: near lines 394 and 661;
   - `docs/chromap.html`, `README.md`, and `CHANGELOG.md` `[Unreleased]`.

   The Multiomics pass-through (`--chromapAtacLowMemFinalizeThreads`) belongs
   to multiomics-suite.
9. **Check:** the M2 gates (§6).

### 4.4 M3: lean decode

1. **Lean decoder.** In `src/atac_kway_spill.{h,cc}`, add
   `DecodeAtacKwaySpillRecordLean(const AtacKwaySpillRecordHeaderV1&, uint16_t schema, PairedEndMappingWithBarcode*, std::string*)`.
   It validates and extracts exactly what the full decoder does
   (`src/atac_kway_spill.cc:616-685`).
2. **Record reads.** In the spill-file reader (`src/overflow_reader.{h,cc}`),
   add `ReadNextAtacRecordHeader(AtacKwaySpillRecordHeaderV1*, std::string* error)`.
   - It applies `ReadNext`'s block framing checks, but returns errors instead
     of exiting.
   - It requires a 48-byte length prefix.
   - It reads straight into the struct.
   - Stdio buffering is unchanged.
3. **Templating.** Template the per-reference function, `EmitAtacGroup` and the
   tail on the record type.
   - The lean heap entry is `{PairedEndMappingWithBarcode, file_index}`, with
     the same comparator.
   - Add `FindBestIndexFromDuplicatesT<R>`, which carries the 1.1.1 fix.
4. **When it applies.** Use the lean decode only for schemas with no BAM pair
   and no raw barcode evidence. That is always true for eligible runs; assert
   it.
5. **Check:** the M3 gates (§6).

### 4.5 Files (1.2.0)

| File | Change |
|---|---|
| `src/mapping_writer.h` | Serial helper declaration; `AtacSpillRecord` specialisation declaration; per-reference members; `AppendAtacNonDualFragment`; `FindBestIndexFromDuplicatesT` |
| `src/mapping_writer.cc` | Loop body moved verbatim; specialisation, driver and assembly; non-dual branch extracted; the `AtacSpillRecord` explicit instantiation of this function removed |
| `src/summary_metadata.h` | Pending-resize, find and atomic-add methods (D8) |
| `src/atac_lowmem_finalize.{h,cc}` (new) | Probe, jobs, partition I/O and cleanup, new-barcode log, descriptor budget, permit wrapper |
| `src/atac_kway_spill.{h,cc}`, `src/overflow_reader.{h,cc}` | M3 only |
| `src/mapping_parameters.h`, `src/chromap_driver.cc`, `src/chromap_lib_runner.cc` | Option; limit raised in the CLI only |
| `Makefile`, `tests/test_lowmem_parallel_finalize.cc` (new), `tests/README.md` | Unit harness; target in the release gate (D10) |
| `scripts/release/run_release_tests.sh`, `docs/releasing.md`, `AGENTS.md`, `README.md` | Release target and counts (D10) |
| `docs/chromap.html`, `CHANGELOG.md` | Documentation |

Not touched:

- `src/chromap.h` and `src/overflow_writer.*`;
- dual BAM/CRAM, the mergeable spill and the materializer;
- `src/version.h`, which is bumped only at release.

## 5. Test and validation plan

### 5.1 Validation root (1.2.0)

```text
V=/mnt/pikachu/lowmem_parallel_validation_20260929
V/bin/        archived binaries, SHA256SUMS, per-build patch (git diff) and git status
V/inputs/     synthetic fixture and index, many-lane lists
V/unit/       unit goldens from the baseline, unit run logs
V/runs/<case>/<mode>/<binary>/   one directory per run
V/compare/    per-case comparison JSON, summary.tsv
V/measure/    spill listings, per-reference tables
V/scripts/    run_case.sh, compare_case.py, capture_spill_listing.sh, measure_sidecar.py
V/artifacts/  CHROMAP_ARTIFACT_ROOT for make targets
```

The 1.2.0 baseline is the 1.1.1 binary. On normal cases 1.1.1 is
byte-identical to v1.1.0 (M-1), so identity against it is identity against
v1.1.0. Keep the v1.1.0 binary from `V111/bin/` as well.

### 5.2 Builds and provenance

```bash
cd <worktree>
git -c protocol.file.allow=always \
  -c submodule.third_party/rapidmacs.url=/mnt/pikachu/Chromap-suite/third_party/rapidmacs \
  submodule update --init third_party/rapidmacs
git submodule status            # expect 34df44818853a1beb59e98b0848813e38a4a8c51
nice -n 10 make -j16 chromap chromap_lib_runner
sha256sum chromap chromap_lib_runner
```

- **Provenance.** Archive every binary with its commit and
  `git status --porcelain`. For an uncommitted build, add `git diff` and the
  patch's sha256.
- **Library path.** Run with `LD_LIBRARY_PATH=<worktree>/third_party/htslib`
  if `ldd` needs it.
- **Stale binary.** The binary in `/mnt/pikachu/Chromap-suite-v110-integration-20260928`
  predates `71b82b7`. Do not use it.
- **Tracked test binary.** `make test-smoke` rebuilds the tracked
  `tests/test_frag_compact_store`. Restore it with
  `git checkout -- tests/test_frag_compact_store`.

### 5.3 Unit tests (`tests/test_lowmem_parallel_finalize.cc`, `make test-lowmem-parallel-finalize`)

The harness drives the public API and needs no FASTQ or index.

Setup and run:

- `MappingWriter<AtacSpillRecord>` with case parameters: output (sidecar, or
  BED/TagAlign text), temp dir, summary path, MACS3 buffer, dedup flags, MAPQ
  threshold, Tn5.
- A `SequenceBatch` built with `AssignLoadedReferenceMetadata`.
- `OutputHeader`.
- For each flush: sort each rid's vector, then `OutputTempMappingsToOverflow`,
  then `RotateThreadOverflowWriter()`.
- Optional pre-seeding with `UpdateSummaryMetadata(bc, SUMMARY_METADATA_TOTAL, n)`.
- Finalisation, `OutputSummaryMetadata`, then destruction of the writer.

Outputs per case:

- the sidecar and `.chroms.tsv`, or the text file;
- the summary CSV;
- a dump of the buckets;
- the stderr counter lines;
- a listing of the temp directory, which must be empty.

Records come from `std::mt19937_64` with fixed seeds, using raw engine output.

| ID | Content |
|---|---|
| U01 | One reference, one file, no duplicates |
| U02 | 3 references × 400 flushes; cell-level dedup; duplicates across files with different MAPQ and read ids around the threshold of 30 |
| U03 | Same (barcode, start, length) at the end of reference r and the start of r+1: must not merge |
| U04 | Identical sort tuples in several files (heap tie-break) |
| U05 | Bulk-level dedup: several barcodes per position; `num_dups` ties broken by abundance; abundance ties; barcodes absent from the table (1.1.1 semantics) |
| U06 | End of stream A: last record MAPQ ≥ 30 and bulk-selected best < 30, at the end of the last non-empty reference |
| U07 | End of stream B: last record < 30, best ≥ 30 |
| U08 | The U06/U07 patterns at the end of a non-last reference, and before trailing empty references |
| U09 | Empty references first, middle and last; all but one empty; no records at all |
| U10 | One group of 300 duplicates (AEV1 count 255; summary DUP 299) |
| U11 | Summary: pre-seeded and new keys; exactly 788 pre-seeded keys (1,024 buckets, upper bound 788), so a resize is pending at the first finalisation put; a new key crossing a threshold as the last insertion; the same new key in several references |
| U12 | Tn5 shift on and off; `remove_pcr_duplicates` off |
| U13 | Outputs: sidecar with buckets; BED text with and without a translate table; TagAlign |
| U14 | Descriptors, in forked children: a limit that reduces workers (identical output); a limit below one reference's files (defined error, empty temp directory) |
| U15 | Cleanup on failure: partition directory unwritable after the spills are written |
| U16 (M3) | Lean decode matches the full decode field for field; both reject the same crafted invalid headers |
| U17 | Permits: fake hooks granting P permits (P = 1, 3). Concurrent tasks never exceed P. Acquires = releases = non-empty references. Work units total the spill records. A failing task still releases. No hooks means N workers. Output identical in every case |

Oracles:

- **Goldens (M0).** Build the harness against the unmodified baseline library
  (1.1.1), adding only the test file and its target. Record each output's
  sha256 in `V/unit/goldens_baseline.tsv`.
- **From M1 on.** Run every case with the v1.1.0 loop (N = 1), the
  per-reference path on one worker, and N = 2, 7, 32 and 64. Each must match
  the goldens and the N = 1 run of the same binary.

### 5.4 Many-lane and forced-spill inputs

With `--low-mem-ram 1K`, every batch flushes. The flush count then equals the
batch count, which is one per lane for lanes of up to 500,000 pairs.

```bash
rep() { local list=$1 n=$2 out="" i; for ((i=0;i<n;i++)); do out+="${out:+,}$list"; done; echo "$out"; }
lanes() { local d=$1 r=$2 out="" l; for l in 1 2 3 4; do out+="${out:+,}$d/pbmc_unsorted_3k_S3_L00${l}_${r}_001.fastq.gz"; done; echo "$out"; }
P100K=/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/fixture/atac
for r in R1 R2 R3; do rep "$(lanes $P100K $r)" 64 > $V/inputs/pbmc100k_x64_$r.txt; done   # 256 lanes
```

Synthetic fixture:

- Run the existing smoke once with the baseline binaries:
  `BUILD=0 CHROMAP_BIN=<baseline chromap> LIBRUNNER_BIN=<baseline chromap_lib_runner> OUTROOT=$V/inputs/synth_smoke bash <worktree>/tests/run_atac_sidecar_only_smoke.sh`.
- It writes `fixture/{ref.fa,read1.fq,read2.fq,barcode.fq,whitelist.txt}`
  and `index/ref.index`.
- Repeat each FASTQ 400 times with `rep`.

Expected mid-batch flush counts, at least: S1 300, P2 200, D2 80, P3 120. If a
count falls short, stop and report; do not change the comparison.

### 5.5 End-to-end cases (baseline against the new build)

Common to every run: `-t 32 --temp-dir $D/tmp`.

Flag sets:

- `PBMC=(-l 2000 --trim-adapters --remove-pcr-duplicates --remove-pcr-duplicates-at-cell-level --Tn5-shift --barcode-whitelist $WL --low-mem)`.
- `PBMC_BULKDEDUP`: `PBMC` without `--remove-pcr-duplicates-at-cell-level`.
- `DOGMA=(--preset atac --read-format bc:8:23:- --barcode-whitelist $WL --barcode-translate $A2G --barcode-translate-from-first)`.
- `SYNTH=(-x $V/inputs/synth_smoke/index/ref.index -r $V/inputs/synth_smoke/fixture/ref.fa --barcode-whitelist $V/inputs/synth_smoke/fixture/whitelist.txt -l 2000 --trim-adapters --remove-pcr-duplicates --Tn5-shift --low-mem)`.
  Add `--remove-pcr-duplicates-at-cell-level` for the cell-level variant.

Outputs:

- `SIDECAR=(--atac-sidecar-only --atac-fragment-binary-output $D/atac_fragments.bin --summary $D/summary.csv --call-macs3-frag-peaks --macs3-frag-peaks-source memory --macs3-frag-peaks-output $D/peaks.narrowPeak --macs3-frag-summits-output $D/summits.bed)`.
  DOGMA runs add `--macs3-frag-low-mem`.
- `BEDTXT=(-o $D/fragments.tsv --summary $D/summary.csv)`.
- `TAGALIGN=(--TagAlign -o $D/tagalign.txt --summary $D/summary.csv)`.

Modes: `default`, and `det` (adds `--deterministic-mapping`).

Inputs:

- Synthetic × 400 and PBMC × 64: the comma lists in `$V/inputs`.
- PBMC: `lanes <dir> R1|R3|R2`, with `-x` and `-r` from §2.
- DOGMA: `-1 R1 -2 R3 -b R2`.

| ID | Input | Flags / output | `--low-mem-ram` | Modes | Expected flushes |
|---|---|---|---|---|---|
| S1c | synthetic × 400 lanes | SYNTH cell-level; SIDECAR, BEDTXT, TAGALIGN | 1K | default, det | ~400 |
| S1b | synthetic × 400 lanes | SYNTH bulk-level; SIDECAR, BEDTXT | 1K | default, det | ~400 |
| P1 | PBMC 100k, 4 lanes | PBMC; SIDECAR and BEDTXT | default and 1K | default, det | ≤1 / 4 |
| P2c | PBMC 100k × 64 (256 lanes) | PBMC; SIDECAR | 1K | default, det | ~256 |
| P2b | PBMC 100k × 64 | PBMC_BULKDEDUP; SIDECAR | 1K | default, det | ~256 |
| D1 | DOGMA 2M | DOGMA; SIDECAR; BEDTXT at 1K | default and 1K | default, det | ≤1 / 4 |
| D1w | DOGMA 2M | DOGMA `--output-mappings-not-in-whitelist`; SIDECAR (new summary keys during finalisation) | 1K | det | 4 |
| D2 | DOGMA 50M | DOGMA; SIDECAR | 1K | default, det | ~100 |
| P3 | PBMC 3k full depth | PBMC; SIDECAR | 1K | default | ~166 |
| F1 | DOGMA lane 1 full depth | DOGMA; SIDECAR | default (1 GiB) | default | ~335 |

Thread settings for the new build:

- Every case runs with `--low-mem-finalize-threads 32` and with the CLI
  default (0).
- S1, P1 and D1 also run with 2 and 7.
- P1 also runs with 1. It must equal the baseline through the v1.1.0 loop.
- Descriptor cases (new build only):
  - S1c under `ulimit -n 1024`: fewer workers, identical output.
  - S1c under `ulimit -n 256`: the defined error, and an empty temp
    directory.

Spill measurement (M0, baseline only, on D2 at 1K and on F1):

- `capture_spill_listing.sh` waits for "overflow files for k-way merge" in
  stderr, then runs
  `find $D/tmp -maxdepth 1 -name 'chromap_*' -printf '%f\t%s\n'`. Spill file
  names are `chromap_<pid>_<thread>_<rid>_<counter>.tmp`.
- `measure_sidecar.py` reads the existing lane-1 full-depth sidecar with numpy
  (read only) for per-reference output records and distinct barcodes.

### 5.6 What is compared (`V/scripts/compare_case.py BASE NEW`)

Byte for byte:

- `atac_fragments.bin` and `atac_fragments.bin.chroms.tsv`;
- `fragments.tsv` and `tagalign.txt`;
- `peaks.narrowPeak` and `summits.bed`;
- in `det` mode, `summary.csv`.

Default-mode summaries:

- `summary.csv` is compared with `cachehit`, `fric`, `estfrip` and
  `numcacheslots` removed by header name. Row order and every other column
  must match.
- This is allowed only after a baseline-vs-baseline control (P1 and D1,
  default mode) differs in those columns alone.
- Any change to this policy needs the author.

Stderr lines that must match exactly:

- `^Number of reads:`, `^Number of reads have multi-mappings:`,
  `^Number of uni-mappings:`;
- `^Processing [0-9]+ overflow files for k-way merge`;
- `^Low-memory overflow: mid-batch flush count:`;
- `^# uni-mappings:`;
- `^Number of output mappings \(passed filters\):`.

The timing line must be present; its value is ignored. The D9 line is ignored.

Also required: the same file set, both exits 0, an empty temp directory after
both runs, and no leftover `atac_fragments.bin.tmp`.

Record each comparison in `V/compare/`. The F1 cross-check against the
embedded v1.1.0 output is informational.

### 5.7 Existing tests, at each code milestone

Under `CHROMAP_ARTIFACT_ROOT=<validation root>/artifacts`:

- `make test-release` (the hermetic release gate);
- `make test-lowmem-bed-100k`;
- `make test-atac-runtime-spill-schema-harness`;
- `make test-lowmem-parallel-finalize`.

The CLI default is on (D1), so the smokes exercise the new path. Also run
`test-atac-sidecar-only-smoke` and `test-libchromap-core-smoke` with wrapper
binaries that append `--low-mem-finalize-threads 1`, so the v1.1.0 loop stays
covered.

### 5.8 M-1 validation (1.1.1 against v1.1.0)

Root: `V111`.

- **Regression test.** Build the new test against the v1.1.0 library, then
  against the fixed one. It must fail before and pass after.
- **Release gate.** `make test-release` (15 targets) on the fixed build.
- **Identity in `det` mode, byte for byte.** Compare the sidecar,
  `.chroms.tsv`, text, peaks, summits, summary and stderr counters on:
  - P1 sidecar and BED, at the default `--low-mem-ram` and at 1K;
  - P1 with `PBMC_BULKDEDUP`;
  - PBMC 100k paired-end without barcodes, BED `--low-mem`, which exercises
    the generic template;
  - PBMC 100k single-end barcoded BED `--low-mem`, also the generic
    template;
  - D1 sidecar and BED at 1K;
  - S1c and S1b sidecar at 1K, with a 100-lane synthetic list.
- **Identity in default mode.** A D1 sidecar pair, plus a v1.1.0-vs-v1.1.0
  control. The summary is compared with the cache columns removed.
- **Evidence, not a gate:**
  - Under valgrind, the regression test's absent-barcode case reports an
    invalid read on v1.1.0 and none on 1.1.1.
  - A synthetic CLI run with barcodes, no whitelist and bulk-level dedup
    crashes on v1.1.0 and completes on 1.1.1.

### 5.9 Multiomics G-M1 (described only; the harness lives in multiomics-suite)

G-M1 compares five fixtures file by file against the 0.9.0 reference outputs
under `/mnt/pikachu/single_binary_reference_20260928`:

- DOGMA-plex lane 1 at 2M reads (sidecar);
- DOGMA-plex five-arm 100k;
- PBMC 3k 100k;
- HIV DOGMA four-arm 100k;
- CAT-ATAC trimodal 100k.

It uses `scripts/compare_composition_outputs.py` and the policies in
`tests/reference/`.

For 1.2.0, a Multiomics agent:

- pins the approved Chromap commit in `compatibility_manifest.json`;
- chooses `--chromapAtacLowMemFinalizeThreads` and passes its permit hooks,
  which already exist;
- runs its G-M1 driver as described in
  `docs/handoffs/HANDOFF_SINGLE_BINARY_M1_M5_20260928.md`.

## 6. Milestones, stop points and stop conditions

Each milestone ends in a stop-and-report: what was done, the evidence paths,
any differences, the commits made (local branches only), and an updated
runbook and handoff.

| Milestone | Work | Gate to pass before reporting |
|---|---|---|
| M-1 Chromap 1.1.1 | §4.0 | The regression test fails on v1.1.0 and passes on 1.1.1; `make test-release` (15 targets) passes; the §5.8 identity cases are byte-identical; release notes, changelog and version drafted. **Stop; wait for the coordinator.** |
| M0 Baseline | Validation root; baseline (1.1.1) build and sha; `strace` check (correction 9); inputs; unit harness and goldens; measurements; all baseline runs in §5.5; control pair | Goldens written; baseline runs exit 0 with the expected flush counts; the control differs only in the cache columns |
| M1 Per-reference, serial | §4.1, §4.2 | Unit: one worker equals the goldens and N = 1; §5.7 tests pass; S1c, S1b, P1 and D1 at N = 2 (serial in M1) identical |
| M2 Parallel | §4.3 | Unit at N = 1, one worker, 2, 7, 32 and 64, plus U14, U15 and U17; §5.7 tests pass; S1c, S1b, P1, D1 and D1w at 2, 7, 32 and 0; the descriptor cases |
| M3 Lean decode | §4.4 | U16 and all unit cases; §5.7 tests pass; S1c, S1b, P1, D1, P2c and P2b at 32 |
| M4 Gates | All §5.5 cases, including D2, P3 and F1; documentation; release-note and changelog drafts | Every comparison passes; G-M1 is ready to hand over |
| Release | 1.1.1 and 1.2.0 tags and pushes | Waits for the author |

Stop and report immediately if any of these happens:

- **An output differs.** Any compared output differs from the baseline, or
  from the N = 1 run in the same binary.
  - Record the case, the file and the first differing byte or line.
  - Never adjust a comparison, a policy, the inputs or the goldens to make a
    case pass.
- **A control differs** outside the cache columns.
- **A run fails:** it exits non-zero, crashes, leaves temp files, or reports
  an unexpected flush count.
- **Out of scope.** The work would need changes outside §4.5's file list (for
  example the spill format, the spill writer, dual BAM/CRAM or the mergeable
  spill), excluded material, or vendor source.
- **Disk.** `/mnt/pikachu` has less than 500 GB free before D2, P3 or F1.
- **Memory.** Peak memory is clearly above the baseline. `/usr/bin/time -v` is
  allowed for max RSS only; never record times.
- **Usage** nears 90%.

## 7. Author decisions (29 September)

| ID | Decision |
|---|---|
| D1 | Accepted. On by default in the CLI (0 = auto). The library field defaults to 1, so hosts opt in. |
| D2 | Dual BAM/CRAM stays serial in 1.2.0. A parallel dual path is a TODO. |
| D3 | Accepted: `--low-mem-finalize-threads N`, with the Multiomics pass-through `--chromapAtacLowMemFinalizeThreads`. |
| D4 | Changed. Permits apply to every task. Each task acquires a permit through the existing hooks before it starts, and releases it reporting work units, bytes and ns. Concurrency is bounded by the granted permits, the workers and the descriptor budget. Without hooks, N workers run. No pool runs outside the host's permits (§4.3 item 5). |
| D5 | Keep the v1.1.0 loop verbatim for N = 1 for now. Removing the duplicate path is a TODO. |
| D6 | TagAlign is included. Bulk is excluded in 1.2.0 and is a TODO. |
| D7 | Accepted. The CLI raises the soft limit; the library does not. When one reference needs more files than the limit allows, fail with a message. |
| D8 | Changed to the exact low-memory variant: atomic adds into existing entries; new barcodes logged and inserted at assembly; the pending resize triggered as the first serial put would; no delta lists (§4.3 item 2). Multiomics will make its summary opt-in, in its own repository. |
| D9 | Accepted. One informational stderr line when N ≥ 2. |
| D10 | Accepted. `test-lowmem-parallel-finalize` joins the release gate. |
| D11 | Fix first. The two bugs ship as 1.1.1 (M-1) before 1.2.0 work. |
| D12 | The documents are committed on the design branch. M-1 now; M0 waits for the coordinator. |

### 7.1 TODO (after 1.2.0)

- **Parallel dual BAM/CRAM finalisation (D2).** The BAM writer and sorter are
  a single stream.
- **Remove the duplicate serial path (D5).** Two copies of the merge loop will
  diverge. Once the per-reference path has shipped, make it the only path.
- **Parallel finalisation for bulk data (D6)**, for the same reason.
- **N = 1 and permits.** Should the N = 1 serial loop (not a pool) hold a
  permit when hooks are present? It would take one permit for its whole
  duration, as a wrapper outside the verbatim loop. Settle this with the
  author.
- **`copy_file_range` for partition assembly**, with the discipline in
  correction 5 and a fallback. Only if it proves useful once benchmarking
  resumes.
- **Splitting one reference into ranges** (design §4.8), if chr1 proves to be
  the limit.
- **Multiomics:** the pass-through option, permit use, and its summary
  becoming opt-in. All three belong to multiomics-suite.

## 8. Effort estimate

Agent-days of work. Machine time is for untimed identity runs, one at a time;
it is a planning figure, not a measurement.

| Milestone | Agent effort | Machine time (untimed runs) |
|---|---|---|
| M-1 | 0.5 day | about 1-2 h |
| M0 | 0.75 day | about 5 h, mostly F1, P3 and D2 |
| M1 | 0.75 day | about 1 h |
| M2 | 1.5 days (permits and the D8 variant add about 0.25) | about 1 h |
| M3 | 0.5 day | about 1 h |
| M4 | 0.75 day | about 3-4 h |
| Total | about 4.75 days (range 3.5-6) | about 12-14 h |

## 9. Pre-existing issues

The two issues found while checking the design note are fixed in M-1 (§4.0).

The `setvbuf(..., nullptr, _IOFBF, size)` calls in `src/overflow_reader.cc:14`
and `src/mapping_writer.h:132` and `:358` do not set the requested size under
glibc. They have no effect on correctness. This note is here so nobody
"fixes" them into large per-file buffers.

## 10. M-1 results

To be filled in when M-1 finishes.
