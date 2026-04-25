# Runbook: MACS3-compatible C++ peak caller speed parity

## Decision

Use external peak calling from the written `fragments.tsv.gz` as the operational path for now.

The in-process / memory-source C++ MACS3-compatible caller stays useful for parity diagnostics and future integration, but it should not be the default end-to-end path until it is within range of external `macs3 callpeak` speed on full data.

## Why this runbook exists

Current full-set timing shows the Chromap mapping phase is not the bottleneck.

Observed full-set run:

- Run output: `/tmp/chromap_peak_memory_fullset_20260425_062130/benchmark_fullset.tsv`
- Fragment rows: `53,969,811`
- Chromap mapping/sort/dedup phase: `507.19s`
- Integrated C++ peak file-source wall: `3548.04s`
- Integrated C++ peak memory-source wall: `3467.33s`
- C++ peak stage estimate: roughly `49-51 min` after subtracting Chromap mapping/sort/dedup.
- External `macs3 callpeak` full-set reference from earlier validation: roughly `12 min`.
- ARC `_perf` reference for this dataset: `DETECT_PEAKS` about `64s`, `GENERATE_PEAK_MATRIX` about `41s`.

Interpretation: the C++ implementation is currently parity-first, not performance-first. It should be instrumented and optimized against external MACS3 before being reintroduced as a production in-process option.

## Current implementation shape

Primary files:

- `src/peak_caller/macs3_frag_peak_pipeline.{h,cc}`
- `src/peak_caller/peak_io.{h,cc}`
- `src/peak_caller/fragment_input.{h,cc}`
- `src/chromap_callpeaks.cc`
- `src/chromap_driver.cc`

Current C++ pipeline stages:

1. Load or materialize fragments into `std::vector<ChromFragments>`.
2. Write MACS3 FRAG treat pileup bedGraph.
3. Write MACS3 no-control lambda bedGraph.
4. Write ppois score bedGraph.
5. Run bdgpeakcall / narrowPeak / summits generation.

Known inefficiency to verify with timers:

- `RunMacs3FragPeakPipelineFromFragments` writes treat and lambda intermediates.
- `WriteFragMacs3BdgcmpScoreBedGraphs` currently regenerates scratch treat/lambda bedGraphs and parses them back before scoring.
- The code serializes large intermediate tracks even when the final output only needs narrowPeak and summits.
- The implementation is mostly serial across chromosomes.

## Target operating model for now

Operational multiome path:

1. Chromap emits `fragments.tsv.gz`.
2. External `macs3 callpeak -f FRAG` runs from that file.
3. The cell-calling evidence layer reads the selected peak set and computes per-barcode ATAC peak-region metrics.
4. The STAR/Chromap multiome pipeline remains file-boundary based until C++ speed parity is reached.

Do not block current cell-calling work on in-process C++ peak calling.

## Phase 0: freeze benchmark surfaces

Create fixed benchmark inputs:

- `100K`: existing 100K fixture for fast correctness / regression.
- `1M` or `5M`: intermediate subset large enough to expose stage costs without hour-long iteration.
- `full`: PBMC 3k multiome full fragments file.

Required metadata for each benchmark:

- Fragment file path.
- Fragment line count.
- MD5 of fragment file.
- Command line.
- Threads.
- Hostname and date.
- Wall/user/sys/max RSS from `/usr/bin/time`.
- Peak count, total peak bp, median peak width.
- BED3 Jaccard and pValue summary vs external MACS3.

Suggested output layout:

```text
out/peak_speed_parity_<date>/
  100k/
  1m/
  5m/
  full/
  summary.tsv
```

## Phase 1: add instrumentation

Add opt-in stage timers to `chromap_callpeaks` and the shared C++ pipeline.

Suggested flag:

```text
--profile-stages <stage_profile.tsv>
```

Minimum profile columns:

```text
stage	wall_sec	user_sec	sys_sec	max_rss_kb	input_rows	output_rows	output_bytes	notes
```

Stages to time:

- `load_fragments`
- `materialize_memory_fragments` if memory source is used
- `sort_fragments_by_chrom`
- `write_treat_pileup`
- `write_control_lambda`
- `write_score_ppois`
- `bdgpeakcall_regions`
- `narrowpeak_summits`
- `cleanup`
- `total`

Implementation notes:

- Use `clock_gettime(CLOCK_MONOTONIC)` or `std::chrono::steady_clock` for wall time.
- Use `getrusage(RUSAGE_SELF)` for user/sys/RSS snapshots.
- Keep profiling disabled by default.
- Emit profile TSV even on failure where possible.
- Include row counts and output byte counts for intermediate bedGraphs.

Acceptance checks:

- Existing 100K parity tests still pass with profiling disabled.
- Profiling enabled does not change narrowPeak or summits bytes on the 100K fixture.
- Profile TSV contains every expected stage exactly once.

## Phase 2: benchmark external MACS3 vs current C++

For each benchmark surface, run:

1. External `macs3 callpeak`.
2. `chromap_callpeaks` current C++ implementation from file fragments.
3. Optional integrated Chromap run only for end-to-end accounting, not for optimization.

Recommended commands:

```bash
/usr/bin/time -v macs3 callpeak \
  -t fragments.tsv.gz \
  -f FRAG \
  -g hs \
  -p 1e-5 \
  --min-length 200 \
  --max-gap 30 \
  -n macs3_frag_p1e5 \
  --outdir macs3
```

```bash
/usr/bin/time -v ./chromap_callpeaks \
  -i fragments.tsv.gz \
  --frag-pileup-macs3-uint8-counts \
  --macs3-frag-narrowpeak cpp_macs3_frag.narrowPeak \
  --macs3-frag-summits cpp_macs3_frag_summits.bed \
  --bdgpeakcall-cutoff 5 \
  --bdgpeakcall-min-len 200 \
  --bdgpeakcall-max-gap 30 \
  --frag-score-pseudocount 0 \
  --profile-stages cpp_profile.tsv
```

Metrics:

- C++ wall / MACS3 wall.
- Stage-level wall percent.
- Peak BED3 identity or Jaccard.
- pValue median/mean/max absolute difference on matched BED3.
- Summit exact rate and distance distribution.
- Max RSS.
- Intermediate bedGraph byte counts.

Gate for moving to optimization:

- Profile identifies the top two wall-time stages on `1M` or larger.
- Full-set baseline is recorded once, but do not repeat full-set runs until a change improves an intermediate benchmark.

## Phase 3: remove duplicated work

First optimization target: avoid regenerating and reparsing treat/lambda tracks.

Current likely hotspot:

- `RunMacs3FragPeakPipelineFromFragments` writes treat and lambda once.
- `WriteFragMacs3BdgcmpScoreBedGraphs` writes scratch treat/lambda again and parses them.

Preferred change:

- Split the score writer so it can consume in-memory per-chromosome treat/lambda run arrays directly.
- Keep the bedGraph-writing mode for diagnostics.
- Add a no-intermediate fast path for narrowPeak/summits.

Acceptance:

- 100K narrowPeak and summits remain byte-identical to the pre-optimization C++ output.
- 100K and intermediate benchmark BED3 remains identical to external MACS3 where it was previously identical.
- Stage profile shows the duplicated treat/lambda work removed or sharply reduced.

## Phase 4: chromosome-level parallelism

Once duplicate work is removed, parallelize independent chromosome work.

Candidate stages:

- Build fragment-span treat pileup by chromosome.
- Build control lambda by chromosome after global totals are known.
- Compute ppois scores by chromosome.
- Call regions / summits by chromosome.

Constraints:

- Global totals for no-control lambda must be computed before per-chromosome scoring.
- Output order must remain deterministic: sort chromosomes lexicographically before final write.
- Do not parallel-write to the same output file; collect per-chromosome buffers or write temp shards then concatenate in deterministic order.

Acceptance:

- `THREADS=1` and `THREADS=8` produce identical output on 100K.
- Intermediate benchmark scales materially with threads.
- Full-set C++ runtime moves toward external MACS3.

## Phase 5: reduce intermediate I/O

After correctness and parallelism are stable, make bedGraph emission optional.

Modes:

- Diagnostic mode: write treat/lambda/ppois bedGraphs for parity debugging.
- Fast mode: keep tracks in memory and emit only narrowPeak and summits.

Suggested flags:

```text
--keep-intermediates <dir>
--no-intermediate-bedgraphs
```

Default for `chromap_callpeaks` can remain diagnostic until speed parity is proven. Default for future pipeline integration should be fast mode.

Acceptance:

- Fast mode and diagnostic mode produce identical narrowPeak and summits on 100K.
- Fast mode reduces full-set disk writes and wall time.

## Phase 6: re-evaluate in-process integration

Only revisit in-process Chromap integration after file-source `chromap_callpeaks` is close to external MACS3 speed.

Re-entry criteria:

- Full-set C++ file-source runtime within `1.5x` external MACS3 wall time, or a clear user-approved reason to accept slower runtime.
- Stage profile shows no obvious duplicated whole-genome pass.
- Memory footprint is documented.
- File-source and memory-source outputs are byte-identical on 100K.

Until then:

- Keep `--macs3-frag-peaks-source file` as the stable path.
- Keep `--macs3-frag-peaks-source memory` as experimental / diagnostic.
- Do not require in-process peak calling for STAR/Chromap end-to-end benchmarks.

## Deliverables for the next coding agent

1. Add `--profile-stages` to `chromap_callpeaks`.
2. Instrument the current stages without changing outputs.
3. Add a benchmark script that runs external MACS3 and C++ on `100K` and an intermediate subset.
4. Produce `summary.tsv` with wall/RSS/parity metrics.
5. Do not optimize until the profile identifies the dominant stage.

## Stop conditions

Stop and report before changing algorithms if:

- C++ output diverges from the previous C++ output on 100K.
- External MACS3 parity regresses.
- Full-set benchmark requires repeated hour-long runs without an intermediate benchmark signal.
- A proposed optimization would make MACS3 parity impossible to verify.
