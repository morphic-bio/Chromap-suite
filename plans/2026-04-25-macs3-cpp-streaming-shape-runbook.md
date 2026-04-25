# Runbook: MACS3-compatible C++ peak caller streaming shape

## Decision

Move the C++ MACS3-compatible FRAG caller away from bedGraph files as internal
state. The first target is not micro-optimization. The first target is the
right execution shape:

1. Read or receive fragments once.
2. Accumulate all fragment-dependent state needed by later phases.
3. Finalize treat, lambda, ppois, regions, narrowPeak, and summits from in-memory
   structures.
4. Keep diagnostic bedGraphs as optional exports of those structures.

Use the existing 100K fixture as the correctness surface at every step. Only
after this shape is working should we benchmark 5M/full and tune memory layout,
parallelism, and streaming dependencies.

## Why this runbook exists

The current implementation mirrors the MACS3 command-line/file workflow:

1. Load fragments into `std::vector<ChromFragments>`.
2. Write treat pileup bedGraph.
3. Write control lambda bedGraph.
4. Read treat/lambda to write ppois.
5. Read ppois/treat/lambda again to call regions, narrowPeak, and summits.

That file-boundary shape was useful for parity diagnostics. It is not the desired
end state for integration with Chromap, where fragments are generated in process.

The desired end state is that fragment generation feeds a peak-calling workspace
directly. Files become debug outputs, not the transport between phases.

## Current five uses of fragment-derived information

These are the logical uses we need to preserve while changing the implementation
shape:

1. **Treat pileup construction**: fragment-span pileup over `[start, end)`.
2. **Lambda construction**: local end-window events plus global totals for MACS3
   no-control lambda scaling.
3. **Ppois scoring**: min-merge treat and lambda breakpoints and compute MACS3
   `bdgcmp -m ppois` scores.
4. **Region calling**: run MACS3-style `bdgpeakcall` on ppois.
5. **Summit and peak annotation**: use treat to choose summits and treat/lambda at
   the summit locus to compute signal, p-value, fold-enrichment, and q-value
   fields.

The current code may read or write bedGraphs to satisfy these uses. The streaming
shape should satisfy them from in-memory state.

## Key constraint

This is not fully online peak emission yet. Some finalization is required after
the fragment stream:

- Input fragments may be unsorted.
- Lambda needs global totals (`total_len`, `tot_c`) before scaled lambda values
  can be finalized.
- Deterministic chromosome ordering must be preserved.

So the immediate target is **single fragment read plus in-memory accumulation and
finalization**, not necessarily emitting final peaks while the first fragment file
read is still active.

## Proposed core abstraction

Introduce a MACS3 FRAG peak workspace, tentatively:

```text
Macs3FragPeakWorkspace
  params
  chrom_order / chrom metadata
  per_chrom fragment-derived event buffers
  global totals for lambda
  finalized treat tracks
  finalized lambda tracks
  finalized ppois tracks
```

Recommended track representation for correctness-first development:

```text
Track
  vector<int32_t> end_positions  // MACS3 p-array semantics
  vector<float> values           // value on previous_end..end_positions[i]
```

Recommended event representation for accumulation:

```text
EventBuffer
  vector<pair<int32_t, double>> treat_span_events
  vector<pair<int32_t, double>> lambda_end_window_events
```

Do not optimize these structures first. Choose clarity and parity. Compact
storage, deduplication, chunking, and parallel per-chrom finalization come later.

## Phase 0: Freeze correctness oracles

Keep the existing file-based path as the oracle while building the streaming
workspace.

Required 100K checks:

- Treat bedGraph byte-identical or accepted exact parity metric vs current path.
- Lambda bedGraph byte-identical or accepted exact parity metric vs current path.
- Ppois bedGraph byte-identical or accepted exact parity metric vs current path.
- NarrowPeak byte-identical vs current C++ path.
- Summits byte-identical vs current C++ path.
- Existing MACS3 BED3 parity target still passes.

Use 100K only in this phase. Do not use 5M/full to drive design decisions until
the phase-level parity checks are stable.

## Phase 1: Build workspace from materialized fragments

Start from `std::vector<ChromFragments>` so the first implementation is isolated
from Chromap mapping internals.

Add functions along these lines:

```text
BuildMacs3FragWorkspaceFromFragments(by_chrom, params, workspace)
FinalizeTreatTracks(workspace)
FinalizeLambdaTracks(workspace)
FinalizePpoisTracks(workspace)
CallPeaksFromWorkspace(workspace, narrowPeak, summits)
```

This phase may still load fragments from TSV once using `LoadFragmentsFromTsv`.
The goal is to remove bedGraph reads/writes between peak-caller phases.

Acceptance:

- Streaming workspace treat export matches current `WriteFragSpanPileupBedGraph`.
- Streaming workspace lambda export matches current
  `WriteFragMacs3NoControlLambdaBedGraph`.
- No change to default CLI behavior.

## Phase 2: Replace ppois file dependency

Compute ppois directly from finalized in-memory treat/lambda tracks.

Keep the existing `WriteFragMacs3BdgcmpScoreBedGraphs` and
`WriteMacs3BdgcmpScoreBedGraphsFromTreatLambdaBdgs` functions for diagnostics and
standalone file workflows.

Add an in-memory scoring function:

```text
BuildPpoisTrackFromTreatLambda(treat_track, lambda_track, pseudocount, ppois_track)
```

It must use the same min-merge and rounding semantics as
`EmitMacs3BdgcmpScoresForChrom`.

Acceptance:

- Exported ppois bedGraph from the in-memory ppois track matches current C++ ppois.
- 100K narrowPeak/summits remain byte-identical after switching the pipeline to
  use in-memory ppois.

## Phase 3: Replace narrowPeak file dependency

Refactor `RunMacs3FragPpoisNarrowPeaks` so the algorithm can consume in-memory
tracks directly:

```text
RunMacs3FragPpoisNarrowPeaksFromTracks(ppois, treat, lambda, params, outputs)
```

The current file-reading function should become a thin adapter:

1. Read bedGraphs into tracks.
2. Call the track-based function.

The track-based function should own the actual algorithm:

- Call merged regions from ppois.
- Use treat to choose summits.
- Use treat/lambda at summit locus for p-value and fold-enrichment fields.
- Run BH q-value calculation.
- Write narrowPeak and summits.

Acceptance:

- 100K narrowPeak and summits byte-identical vs the file-based current C++ output.
- Existing `--bdgpeakcall-*` diagnostic mode still works.

## Phase 4: Use workspace in `RunMacs3FragPeakPipelineFromFragments`

Change the shared pipeline to:

1. Build workspace from `by_chrom`.
2. Finalize treat and lambda tracks.
3. Optionally export treat/lambda bedGraphs if requested or needed for diagnostics.
4. Finalize ppois track in memory.
5. Optionally export ppois bedGraph if requested.
6. Call narrowPeak/summits from in-memory tracks.

At the end of this phase, bedGraph files are no longer required internally for
the normal narrowPeak/summits path.

Acceptance:

- 100K C++ outputs byte-identical to pre-refactor C++ outputs.
- All existing 100K parity targets still pass.
- `--frag-span-pileup-bdg`, `--frag-lambda-bdg`, and `--frag-score-ppois-bdg`
  still emit diagnostics with the same semantics.

## Phase 5: Feed workspace from fragment generation

After the `by_chrom`-based workspace is correct, connect it to Chromap’s fragment
generation path.

Current memory source (`FragPeakMemoryAccumulator`) stores compact fragment rows
and later materializes `ChromFragments`. Replace or augment it with direct
workspace accumulation:

```text
OnFragment(chrom_id, start, end, count):
  workspace.AddFragment(chrom_id, chrom_name, start, end, count)
```

This should accumulate:

- Treat span events.
- Lambda end-window events.
- Lambda global totals.
- Any metadata needed for deterministic finalization.

Acceptance:

- Memory-source and file-source outputs byte-identical on 100K.
- No requirement to reduce RAM yet.
- Finalization still produces deterministic chromosome order.

## Phase 6: Benchmark the shape

Only after phases 1-5 are correct:

- Run the 5M fixture.
- Compare stage profile vs current file-boundary pipeline.
- Record MACS3 wall vs C++ wall.
- Identify whether remaining cost is event finalization, ppois scoring,
  region calling, or output writing.

Acceptance:

- 5M C++ narrowPeak and summits byte-identical to pre-refactor C++ output.
- `stage_profile.tsv` shows bedGraph file-read stages removed from the internal
  path.

## Phase 7: Optimize structures and parallelism

Only after the streaming shape is correct:

- Choose compact event/track representations.
- Avoid storing redundant event and finalized track forms where possible.
- Parallelize per-chrom finalization after global lambda totals are known.
- Consider streaming sorted shards or per-chrom chunking if RAM becomes limiting.

Do not do these before parity is locked.

## Non-goals during streaming-shape work

- Do not change MACS3 thresholds or scoring semantics.
- Do not remove diagnostic bedGraph outputs.
- Do not optimize memory layout first.
- Do not run full datasets repeatedly.
- Do not change default CLI behavior until the new path is proven.

## Suggested development order

1. Add shared in-memory track types and bedGraph export helpers.
2. Generate treat from workspace and export it; compare 100K.
3. Generate lambda from workspace and export it; compare 100K.
4. Generate ppois from in-memory treat/lambda; compare 100K.
5. Refactor narrowPeak/summits to consume tracks; compare 100K.
6. Switch `RunMacs3FragPeakPipelineFromFragments` to workspace path; compare 100K.
7. Connect Chromap fragment generation to direct workspace accumulation; compare
   memory source vs file source on 100K.
8. Benchmark 5M and then plan memory/parallel optimization.

## Stop conditions

Stop and report before continuing if:

- Any 100K C++ output changes unexpectedly.
- Diagnostic bedGraph exports diverge from the current implementation.
- A proposed change requires changing MACS3 semantics.
- The implementation starts optimizing storage before the phase parity checks are
  in place.

## Staged simplification implementation plan

This section tracks the follow-on cleanup after the first workspace-backed path
landed. The priority is still correctness and streaming shape before memory
tuning.

### Stage 1: Make the workspace API truly streaming-shaped

Goal: stop treating `std::vector<ChromFragments>` as the only population path.

Add explicit workspace lifecycle APIs:

```text
InitMacs3FragWorkspace(chrom_names, params, workspace)
AddMacs3FragWorkspaceFragment(workspace, chrom_id, start, end, count)
```

Then make `BuildMacs3FragWorkspaceFromFragments` a compatibility adapter:

```text
InitMacs3FragWorkspace(...)
for chrom in by_chrom:
  for fragment in chrom.frags:
    AddMacs3FragWorkspaceFragment(...)
```

Acceptance:

- `make chromap_callpeaks`
- `make chromap`
- `make test-peak-narrowpeak-100k`
- 100K direct diagnostic export compare: workspace pipeline treat/lambda/ppois
  byte-identical vs standalone diagnostic writers.

### Stage 2: Separate finalization from export

Goal: bedGraphs become optional outputs, not internal pipeline stages.

Refactor `RunMacs3FragPeakPipelineFromFragments` into explicit phases:

```text
Build/initialize workspace
FinalizeTreatTracks
FinalizeLambdaTracks
FinalizePpoisTracks
RunMacs3FragPpoisNarrowPeaksFromTracks
ExportDiagnosticsIfRequested
Cleanup
```

Compatibility rules:

- If the user passed `--frag-span-pileup-bdg`, write treat.
- If the user passed `--frag-lambda-bdg`, write lambda.
- If the user passed `--frag-score-ppois-bdg` or an override path, write ppois.
- If `keep_intermediates_dir` is set, write all intermediates there.
- Do not create a temp workdir unless an export path needs it.

Stage profile should split compute from export:

- `finalize_treat_tracks`
- `finalize_lambda_tracks`
- `finalize_ppois_tracks`
- `export_treat_pileup`
- `export_control_lambda`
- `export_score_ppois`

Acceptance:

- Default narrowPeak/summits behavior unchanged.
- Diagnostic bedGraph flags still emit identical files.
- Profiling on/off narrowPeak/summits byte-identical on 100K.

### Stage 3: Deduplicate narrowPeak logic

Goal: one algorithm implementation for file and in-memory callers.

Make `RunMacs3FragPpoisNarrowPeaksFromTracks` canonical. Change
`RunMacs3FragPpoisNarrowPeaks` to read ppois/treat/lambda into track
collections and delegate to the track-backed function.

Acceptance:

- `make test-peak-bdgpeakcall-100k`
- `make test-peak-narrowpeak-100k`
- Track-backed and file-backed narrowPeak/summits outputs match on 100K.

### Stage 4: Deduplicate MACS3 track helpers

Goal: remove parity-sensitive duplicate helper code between `peak_io.cc` and
`macs3_frag_workspace.cc`.

Extract shared helpers into a module such as:

```text
src/peak_caller/macs3_track_ops.{h,cc}
```

Candidate helpers:

- effective fragment weight/count handling
- event coalescing
- event-to-track conversion
- treat/control breakpoint merge
- bedGraph run merge
- ppois track construction from treat/lambda
- bedGraph track export

Acceptance:

- Treat/lambda/ppois diagnostic exports byte-identical on 100K.
- Score parity target still passes.
- No behavior changes beyond code consolidation.

### Stage 5: Wire direct fragment generation into the workspace

Goal: let Chromap fragment generation feed the workspace directly.

Add a path where each retained ATAC fragment calls:

```text
AddMacs3FragWorkspaceFragment(workspace, rid, start, end, count)
```

Initially this can coexist with `FragPeakMemoryAccumulator`. After parity is
locked, peak finalization should use the workspace directly rather than
materializing `ChromFragments`.

Acceptance:

- 100K memory-source vs file-source narrowPeak/summits byte-identical.
- Existing integrated 100K target passes.
- No RAM optimization required yet.

### Stage 6: Benchmark and optimize

Only after the shape is correct, run the 5M benchmark and compare:

- old file-boundary pipeline
- workspace no-intermediate-export path
- external MACS3

Then optimize based on data:

1. `bdgpeakcall_regions`
2. event finalization
3. memory / duplicate event-track storage
4. chromosome-level parallel finalization

Acceptance:

- 5M narrowPeak/summits byte-identical to pre-simplification C++ output.
- Stage profile shows intermediate bedGraph export removed from the default
  critical path.
- Performance claims are backed by 5M data, not 100K.
