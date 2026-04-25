# Summary: MACS3-compatible C++ peak caller — profiling, 5M benchmark, first optimization

Date: 2026-04-25
Workspace: `/mnt/pikachu/Chromap-suite`

## What shipped

1. **Instrumentation** — `--profile-stages` on `chromap_callpeaks`, stage TSV (`wall_sec`, `user_sec`, `sys_sec`, `max_rss_kb`, row/byte counts, notes). Stage wall time excludes full-file line scans; companion rows `profile_file_metrics_*` record scan/stat cost.
2. **Harness** — `tests/run_macs3_cpp_speed_profile.sh` (MACS3 + C++, parity check profiling on vs off). `RUN_MACS3=0` skips external `macs3` even if on `PATH`.
3. **5M subset** — `tests/build_peak_speed_fragments_5m.sh` → `out/peak_speed_parity_20260425/fragments_5m.tsv.gz` + `fragments_5m.metadata.tsv` (not committed; `out/` in `.gitignore`).
4. **First optimization** — `WriteMacs3BdgcmpScoreBedGraphsFromTreatLambdaBdgs` in `peak_io`: ppois from **existing** treat + lambda bedGraphs (no scratch rebuild from fragments). `RunMacs3FragPeakPipelineFromFragments` uses this after writing treat/lambda. `WriteFragMacs3BdgcmpScoreBedGraphs` kept for fragment-only diagnostic callers.

## Key files touched

- `src/peak_caller/stage_profile.{h,cc}`, `src/chromap_callpeaks.cc`, `src/peak_caller/macs3_frag_peak_pipeline.{h,cc}`, `src/peak_caller/bdgpeakcall.cc`, `src/peak_caller/peak_io.{h,cc}`, `src/chromap_driver.cc`, `Makefile`
- `tests/run_macs3_cpp_speed_profile.sh`, `tests/build_peak_speed_fragments_5m.sh`, `tests/peak_caller_100k_common.sh` (MACS3 + FRAG without BAM)
- `.gitignore` — `out/`

## 5M subset (local build)

| Field | Value |
|-------|--------|
| Path | `out/peak_speed_parity_20260425/fragments_5m.tsv.gz` |
| Rows | 5,000,000 |
| MD5 (gzip) | `d0e787fbdc2d51d10c79ec7d09a7da4f` |
| Source TSV | `.../star_chromap_concurrent_full_20260424_174426/chromap_fragments.tsv` |

Example harness:

```bash
FRAGMENTS_TSV_GZ=/path/to/Chromap-suite/out/peak_speed_parity_20260425/fragments_5m.tsv.gz \
OUTDIR=/path/to/Chromap-suite/out/peak_speed_parity_20260425/5m \
RUN_MACS3=1 bash tests/run_macs3_cpp_speed_profile.sh
```

## Validation

- `make chromap_callpeaks` — OK
- `make test-peak-narrowpeak-100k` — OK (BED3 vs MACS3 unchanged)
- Profiling on vs off — C++ narrowPeak and summits **byte-identical** (100K and 5M)
- 5M C++ outputs **byte-identical** before vs after optimization (`cmp` on off-profile narrowPeak + summits)

## 5M timing (stage `wall_sec`, from `stage_profile.tsv`)

Compared to a temporary baseline where the pipeline still called `WriteFragMacs3BdgcmpScoreBedGraphs` (rebuild + parse scratch treat/lambda):

| Stage | Before | After |
|-------|--------|--------|
| `write_score_ppois` | ~18.36 s | ~13.69 s |
| `total` (profile sum) | ~503.7 s | ~496.5 s |
| GNU time C++ wall (prof off) | ~8:20 | ~8:15 |

Most wall time remains in **`bdgpeakcall_regions`** (~475–478 s): triple bedGraph read + region merge.

On one run with MACS3 enabled on the same 5M input: **~1:07.5** MACS3 wall vs **~8:20** C++ wall (~7.5×).

## Next optimization target

**`bdgpeakcall_regions` / `RunMacs3FragPpoisNarrowPeaks`:** avoid parsing treat, lambda, and ppois from disk again when tracks are already produced in-process; keep deterministic chrom order and current thresholds.

## Related plan

- `plans/2026-04-25-macs3-cpp-speed-parity-runbook.md`
