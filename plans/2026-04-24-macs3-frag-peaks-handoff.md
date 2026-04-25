# Handoff: Chromap MACS3-compatible FRAG peaks + follow-on design discussion

**Audience:** Reviewing agent (implementation, architecture, or validation).
**Repos:** Primary work is **`/mnt/pikachu/Chromap-suite`** (not `multiomic-atac-scrna` unless docs live there).

---

## 1. What was implemented (baseline)

- **Opt-in** MACS3-compatible **FRAG** narrowPeak/summits in the main **`chromap`** binary after mapping.
- **Not default**; requires **`--call-macs3-frag-peaks`** plus **`--macs3-frag-peaks-output`** and **`--macs3-frag-summits-output`**, with **`--atac-fragments`** and the same dual-ATAC constraints as today (paired-end, barcodes, BAM/CRAM; **`--low-mem` still incompatible with `--atac-fragments`**).
- **First integration** reads the **written** fragments file (`.tsv` / `.tsv.gz`) and runs the **shared** C++ pipeline (`RunMacs3FragPeakPipelineFromFragments` in `src/peak_caller/macs3_frag_peak_pipeline.{h,cc}`). **`chromap_callpeaks`** uses the same implementation (refactor from standalone-only duplication).
- **Sidecar:** if **`--summary`** is set, parameters + qValue/signalValue non-parity note go to **`<summary>.macs3_frag_peaks.tsv`**.
- **Harness:** `tests/run_chromap_peak_integration_100k.sh`, **`make test-peak-integration-100k`** (builds `chromap` + `chromap_callpeaks`).
- **Validation context:** Full-set notes from the user (BED3 identity, Jaccard 1.0, ~50k peaks, q/signal not byte parity) vs **100K fixture** subset in CI-style runs.

Key files touched earlier in the thread: `chromap_driver.cc`, `mapping_parameters.h`, `chromap_callpeaks.cc`, `Makefile`, new `macs3_frag_peak_pipeline.*`, integration test script.

---

## 2. Performance (measured on 100K fixture, 8 threads)

- **Integrated** (`chromap` map + peaks): ~**11.7 s** wall.
- **Sequential** (`chromap` map-only ~**9.55 s** + `chromap_callpeaks` ~**2.12 s**): ~**11.7 s** total.
- **Conclusion:** End-to-end wall time is ~**parity**; peak stage is the **same C++ code**; integration mainly avoids a **second process**, not a cheaper algorithm. **C++ path vs real `macs3`** on **full** data remains **slower** (user’s prior note); 100K is not representative of that gap.

---

## 3. “True integration” / peaks before or without re-reading fragments

**Question:** Run MACS3-style peaks **before** fragments hit disk or **without** `LoadFragmentsFromTsv`.

**Facts:**

- The validated pipeline needs a **complete retained fragment set** (or equivalent) because **no-control λ** uses a **dataset-wide** term (e.g. total fragment mass vs effective genome). **Prefixes of a stream are insufficient** for correct final λ/ppois/peaks.
- **Options:**
  - **(A)** Current: write fragments → **read file again** for peaks (extra I/O).
  - **(B)** **In-memory accumulate** during `MappingWriter<PairedEndAtacDualMapping>::AppendMapping` (same rows as TSV), then call **`RunMacs3FragPeakPipelineFromFragments`** at end of mapping → **saves decompress + parse + second read**; **RAM = Θ(N fragments)** in a compact struct form (still linear in row count).
  - **(C)** Buffer **all** fragments, run peaks, **then** write user fragments file → peaks **before** user-visible file; **highest RAM**.

**External `macs3` binary:** Still needs a **materialized FRAG** input (file/pipe/temp); not a free in-mapper stream without buffering or replay.

---

## 4. “Streaming fragments + ingest MACS3 at the same time”

- **Parallel wall-clock** (producer writes / queues fragments, consumer builds tracks) is possible.
- **Finishing** correct MACS3-style peaks still waits on **global** statistics and **full** score/pileup context — consumer cannot produce **final** peaks from only a **prefix** unless the algorithm is changed or **two-pass** / **buffer** is used.
- **Sorted ingest:** **Not required** for correctness of the **multiset** (additive pileup); **sorted / grouped chrom** is used for **simple linear** bedGraph-style code and parity with MACS3’s sequential readers. Unsorted parallel ingest is OK if implementation **normalizes** (per-chrom sort or equivalent) before track build.

---

## 5. RAM, “efficient structure,” `--low-mem`

- In-memory **`ChromFragments`** is **smaller than TSV text** but still **O(N)** fragments.
- **`--low-mem`** is **not** a drop-in fix: driver **forbids** **`--atac-fragments` with `--low-mem`** today. Supporting dual ATAC + peaks + low-mem would need **new** design (e.g. spill binary fragment temp, or extend overflow path).

---

## 6. Suggested next steps (if product wants them)

1. **(B)** Accumulate fragments in **`AppendMapping`** (thread-safe if ever multi-writer), run shared pipeline post-map, **keep** writing the same fragments file → **remove extra read**; document **RAM** tradeoff vs file-based peak path.
2. Optional **flag:** `--macs3-frag-peaks-from-file` (current) vs `--macs3-frag-peaks-in-memory` (new) for operator control.
3. **Low-mem + dual:** large separate project; clarify with stakeholders before scoping.

---

## 7. Open questions for reviewer

- Is **peak-from-memory** acceptable default when `--call-macs3-frag-peaks`, or must **file path** remain the only source for audit/repro?
- Should **integration test** be extended to assert **byte-identical** peaks when comparing in-memory vs file path (once (B) exists)?

---

*End handoff.*
