# Chromap Suite

Chromap Suite is an open-source C++ chromatin-accessibility platform: ATAC-seq alignment with native BAM output and coordinate sorting, an in-process narrow peak caller (`libMACS3`) byte-identical to MACS3 v3.0.3, an embeddable callable-library API (`libchromap`), and an MCP server with a browser Launchpad. Chromap Suite produces the MACS3 narrow peaks the field already uses and composes with [STAR Suite](https://github.com/morphic-bio/STAR-suite) for an end-to-end multiomic (RNA + ATAC) single-binary pipeline. Bulk ATAC, scATAC, ChIP-seq, and Hi-C all run through the same binary; multiomic processing is delivered through `libchromap` embedded in STAR Suite as a worker thread.

On the public 3K PBMC Multiome at 32 threads, the integrated multiomic pipeline (STAR Suite + libchromap + libMACS3) completes in **18:17 / 64.8 GB peak RSS** vs **40:04 / 79.1 GB** for Cell Ranger ARC v2.2.0 (**2.19× faster, ~18% lower memory**), producing **50,274 narrowPeak peaks byte-identical** to MACS3 v3.0.3 (md5 `34f9f991…`). Standalone `libMACS3` runs **7.8× / 10.3× / 11.3× faster** than Cython MACS3 v3.0.3 at 1 / 4 / 24 threads.

Agent quickstart: see [`AGENTS.md`](AGENTS.md) for repo-specific guardrails, tests, and recent changes.

Chromap Suite was spun off from [Chromap](https://github.com/haowenz/chromap) in 2026 and no longer tracks upstream; see [HISTORY.md](HISTORY.md) for lineage.

## Chromap Suite capabilities

The full set of capabilities is organised by scope, mirroring Table 2 of the [Chromap Suite preprint](https://github.com/morphic-bio/chromap_suite_paper):

### Core capabilities

- **Native BAM output + coordinate sort** (`--BAM --sort-bam --write-index`). Replaces the conventional `chromap | samtools sort | samtools index` chain via htslib + a *k*-way disk-merge spillover. Produces `@HD VN:1.6 SO:coordinate` headers and `.bam.bai` indexes in one process. Sort key: `(tid, pos, flag, mtid, mpos, isize)` with `read_id` tie-break (note: differs from `samtools sort` QNAME tie-break). See [`docs/sort_spec.md`](docs/sort_spec.md). Compatible with `--low-mem`.
- **In-process libMACS3 narrow peak calling** (`--call-macs3-frag-peaks`). Fragments are handed to `libMACS3` through the `FragmentIterator` API without an intermediate `fragments.tsv.gz` write, so a single chromap invocation produces sorted-and-indexed BAM, fragments, and MACS3-equivalent narrowPeak and summit outputs. The default threshold remains MACS-style p-mode (`--macs3-frag-pvalue 1e-5`); `--macs3-frag-qvalue Q` switches to q-value/FDR thresholding for `macs3 callpeak -q` compatibility. Bulk ATAC supported via `--macs3-frag-peaks-source memory` (no barcode required); scATAC/multiomic via barcoded fragments. Output is byte-identical to standalone MACS3 v3.0.3 in the default p-mode path (50,274 peaks, md5 `34f9f991…` on the 3K PBMC ATAC channel).
- **Y-chromosome filtering** (`--emit-Y-bam`, `--emit-noY-bam`, `--emit-Y-noY-fastq`, `--emit-Y-read-names`). Three-stream output (all / Y-only / noY) for sex-aware analyses. Works with `--sort-bam`. Detection is case-insensitive and matches `Y`, `chrY`, `CHR_Y`, `chr_y`; decoy/random/alt contigs (`chrY_random`, etc.) are intentionally excluded. See the [Y-chromosome filtering](#y-chromosome-filtering) section below.
- **`libchromap.a` callable library**. The full Chromap Suite ATAC pipeline is exposed through the `libchromap` API (`RunAtacMapping()`, `ChromapAtacConfig`, `ChromapPermitHooks`). The same library backs the `chromap` CLI binary and is the integration point used by STAR Suite for multiomic processing. See [`src/libchromap.h`](src/libchromap.h) and [`src/libchromap.cc`](src/libchromap.cc).
- **`--Tn5-shift-mode {classical|symmetric}`** picks the Tn5 cut-site offset convention on BED/BEDPE/PAF. `classical` (`+4 / -5`; Buenrostro 2013 / Cell Ranger ARC) is the default; `symmetric` (`+4 / -4`; ChromBPNet) is the alternative. Implies `--Tn5-shift`. Active offsets are echoed at startup. SAM/BAM output remains intentionally unshifted (shifting would require coordinated edits to `POS`, `MPOS`, `TLEN`, `CIGAR`, `NM`, `MD`).
- **`--temp-dir DIR`** for custom temporary directory (helpful in Docker/container environments).

### Reliability and tooling

- **Low-memory spillover (rewritten architecture)**. Per-thread overflow writers feeding a *k*-way merge on read-back **replace** the prior shared-buffer + atomic-write design. Supports the full cross-product of `--low-mem` with `--atac-fragments`, BAM output, `--macs3-frag-low-mem`, and Y-filtering modes; production-scale runs (≳10⁹ reads) handled cleanly. A pre-existing race condition in the legacy spillover that produced silent read drops at ~1/10⁴ on typical datasets and intermittent hangs at production scale is resolved as a side effect of the rewrite. The legacy temp-file system is available via `LEGACY_OVERFLOW=1` at compile time (single-threaded only). See [`HISTORY.md`](HISTORY.md) for file references and validation history.
- **Concurrency coordination across new collaborators**. Worker threads coordinate with the native BAM writer, STAR Suite's permit allocator (when embedded), and libMACS3 peak-call paths under the existing OpenMP scheduler. The smoke matrix exercises all combinations of low-mem / BAM / peak-call / Y-filter modes.
- **Regression suite (C01–C11)**. An 11-area parity matrix covering the main user-visible surfaces: index build, paired BED output, ChIP and ATAC presets, scATAC barcode handling, sorted BAM and index, low-memory BED parity, ATAC BAM + fragments, libMACS3 narrow peak calling, Hi-C pairs, and Y/noY split. Three tiers: **S0** hermetic synthetic for pre-commit smoke; **S1** ENCODE downsample for real-assay confidence (paired ENCODE accessions and downsample manifests committed; FASTQs cached out-of-tree); **S2** the 100K and paper-fixture tier reserved for heavier gates. The S0 tier is mandatory for pre-commit checks; S1 and S2 are opt-in.
- **MCP server + Launchpad** (`mcp_server/`). Schema-driven workflow renderer for agents plus a browser recipe UI for humans. Both surfaces consume the same parameterised YAML recipe registry (`mcp_server/recipes/registry.yaml`). The MCP face exposes recipes to agents as schema-validated tool calls; the Launchpad face renders the same recipes as web forms. See the [Chromap Launchpad](#chromap-launchpad-recipe-builder) section below.

### ATAC-multiomic

- **ATAC fragment sidecar** (AEV1 format, `--atac-fragment-binary-output`). Compact binary file emitted alongside BAM/CRAM plus fragments TSV: a 32-byte header followed by 24-byte fragment records keyed by the run's barcode length, with `(size - 32) / 24 == record-count` parity. Chromosome ids are paired with names in `<sidecar>.chroms.tsv`. The runtime spill harness decodes sidecar records back to full `(chrom, start, end, barcode, count)` tuples and compares them with the fragments TSV/BED baseline. Lets STAR Suite's `libscrna` empty-cells function call cells without re-parsing the gzipped fragments TSV. On the 3K PBMC headline run the sidecar contains 53,969,811 records (md5 `a4251bbc…`). See [`src/mapping_writer.h`](src/mapping_writer.h) and [`docs/atac_runtime_spill_schema_runbook.md`](docs/atac_runtime_spill_schema_runbook.md).
- **Multiomic integration with STAR Suite**. `libchromap.a` runs as a STAR Suite worker thread for concurrent ATAC + GEX processing in a single `STAR` invocation. STAR Suite's permit allocator partitions a shared thread budget across GEX mapping, feature processing, and ATAC mapping by observing per-domain drain rates and rebalancing toward simultaneous completion. End-to-end on 3K PBMC: 18:17 / 64.8 GB / 2.19× faster than Cell Ranger ARC v2.2.0. See [STAR Suite](https://github.com/morphic-bio/STAR-suite) for the integration entry point.
- **Native CBQ input** (`--input-format cbq`, `--read-pair-cbq`, `--barcode-cbq`). Maps paired-end ATAC/scATAC reads straight from BINSEQ CBQ files with no intermediate FASTQ, producing fragments byte-identical (under canonical sort) to the FASTQ path. CBQ sequence decode writes directly into Chromap's `SequenceBatch` buffers, including Y/noY FASTQ sidecar emission. FASTQ remains the default. ATAC-only for this milestone; see the [Native CBQ input (ATAC)](#native-cbq-input-atac) sample command below.

## Folder Structure

```
src/                  # Chromap Suite C++ sources and public libchromap header
libchromap.a          # Built static library
bin/                  # CLI binaries (chromap, chromap_callpeaks)
third_party/
  htslib/             # Vendored htslib (BAM read/write/sort/index)
  libMACS3/           # Vendored libMACS3 submodule (narrow peak caller)
mcp_server/           # MCP server + Launchpad UI
  recipes/            # YAML recipe registry consumed by both MCP and Launchpad
  workflows/          # Workflow schemas
  launchpad/          # Browser UI (api.py + static assets)
tests/                # Regression matrix C01–C11, smoke harnesses
plans/                # Design notes and runbooks
docs/                 # User-facing docs (sort_spec, chromap_launchpad, ...)
scripts/              # Helper scripts (validators, launchpad runner, ...)
```

## Modules

- **chromap** (`src/`, `bin/chromap`) — the CLI binary covering bulk ATAC, scATAC, ChIP-seq, and Hi-C. Supports the standard ATAC, ChIP, scATAC, and Hi-C flag surfaces, plus the Chromap Suite additions above (BAM, sort, peak-call, sidecar, Y-filtering). Build: `make`.
- **libchromap** (`libchromap.a`, `src/libchromap.h`) — the same code as a callable C++ library. The `RunAtacMapping()` entry point + `ChromapAtacConfig` + `ChromapPermitHooks` structs let a host process embed Chromap Suite's ATAC pipeline directly with shared thread coordination. Build: produced as a side effect of `make`.
- **libMACS3** (`third_party/libMACS3/lib/libmacs3.a`) — vendored submodule providing a portable C++ implementation of MACS3's narrow peak-calling capability. Standalone CLI `macs3frag` reads a fragments TSV and produces byte-identical narrowPeak. Build: `make libmacs3` (or transitively via `make`). See [the libMACS3 repo](https://github.com/morphic-bio/libMACS3).
- **chromap_callpeaks** (`bin/chromap_callpeaks`) — alias for libMACS3's `macs3frag` binary, kept for harness compatibility. Reads standard fragments files and produces narrowPeak; **7.8× faster** than Cython MACS3 single-threaded, **10.3× faster at 4 threads**.
- **MCP server + Launchpad** (`mcp_server/`) — agent automation service plus browser recipe UI. Initial recipes: `chromap_index`, `chromap_atac_bed`, `chromap_atac_bam_fragments`, `chromap_hic_pairs`. See [`mcp_server/README.md`](mcp_server/README.md).

## Benchmarks

Benchmarks below use public 10x PBMC 3K Multiome data unless noted. The STAR
Suite headline run uses a 32-thread host (i9-13900KF, 126 GB RAM). The Chromap
CLI parity benchmark uses 16 Chromap threads and 16 samtools threads, records
`/usr/bin/time -v` for every step, and is intended to isolate the integrated
Chromap Suite command from the established Chromap + samtools + MACS3 workflow.

| Workflow | Surface | Baseline | Result | Note |
|---|---|---|---|---|
| Multiome (RNA + ATAC) | Single `STAR` invocation: alignment + Solo + GEX `EmptyDrops_CR` + concurrent libchromap + libMACS3 narrow peaks + ATAC `evidence-from-peaks` (T=32, lease=256, pool=48) | Cell Ranger ARC v2.2.0 (`cellranger-arc count --create-bam=true --nosecondary --disable-cell-annotation --localcores=32`); **40:04 / 79.07 GB peak RSS** | **18:17.52 / 64.80 GB; 2.19× faster, ~18% lower memory** | Apples-to-apples scope: alignment + GEX UMI counting + GEX `EmptyDrops_CR` + ATAC mapping + narrow peak calling + ATAC per-barcode cell call. ARC RSS read from cellranger's own `_perf` JSON (ARC forks 252 child processes at peak, so parent-process `/usr/bin/time -v` is misleading). |
| Narrow peaks vs MACS3 | Same multiomic pipeline output | MACS3 v3.0.3 standalone on the same input | **50,274 peaks, byte-identical** (md5 `34f9f991…`); summits also byte-identical (md5 `b54a556…`) | The pipeline delivers the MACS3 narrow peaks downstream ATAC analyses already use; ARC produces 81,157 peaks via its own custom wide caller, which is not a parameterisation of MACS. |
| libMACS3 standalone | Fragments TSV in → narrowPeak out | MACS3 v3.0.3 single-threaded (721 s, 2.50 GB) | **92 s / 2.77 GB at 1 thread (7.8×)**; **70 s at 4 threads (10.3×)**; 64 s at 24 threads (11.3×) | All thread counts produce byte-identical narrowPeak. Memory parity at the practical 4-thread budget; sublinear scaling beyond 4 threads because the per-chromosome serial sweep dominates the remaining wall time. |
| ATAC fragment sidecar | AEV1 binary, 32-byte header + 24-byte records + `.chroms.tsv` metadata | gzipped fragments TSV re-parse | 53,969,811 records, 1,295,275,496 B (md5 `a4251bbc…`) | Consumed by `libscrna::atac::RunAtacEvidenceFromBinary` for ATAC cell calling; integration guard fires when `atacAcquire=0` after `--chromapAtacEnable 1 --dynamicThreadInterface 1`. |

The end-to-end multiomic comparison anchors on a fully public dataset for direct reproducibility.

### Chromap CLI + libMACS3 parity benchmark

This benchmark compares the integrated Chromap Suite command against the
established command-chain baseline:

```sh
# Baseline BAM path
chromap ... --SAM -o /dev/stdout | samtools view -b - | samtools sort -o possorted.bam -
samtools index possorted.bam

# Baseline peak path
chromap ... --BED -o fragments.tsv
macs3 callpeak -t fragments.tsv -f FRAG -g hs -p 1e-5 --min-length 200 --max-gap 30

# Integrated path
chromap ... --BAM --sort-bam --write-index \
  --atac-fragments fragments.tsv.gz \
  --call-macs3-frag-peaks \
  -o possorted.bam
```

Use `--macs3-frag-qvalue 0.05` on the integrated path, or
`MACS3_QVALUE=0.05` in the benchmark scripts, when comparing against workflows
that use MACS `-q` semantics.

The baseline uses unfixed upstream Chromap `0.3.3-r519` and standalone MACS3
v3.0.3. The integrated path uses Chromap Suite `0.3.3-r519`. Both normal and
`--low-mem` modes were run on the full 3K PBMC ATAC fixture with 16 threads.

| Mode | Integrated wall / RSS | Full baseline wall / peak RSS | Integrated speedup | Fragment parity | Peak parity |
|---|---:|---:|---:|---|---|
| normal | **7:34.23 / 101.1 GB** | 23:05.60 / 127.0 GB | **3.05x faster; 67.2% less wall time** | **53,969,811 vs 53,969,811; exact sorted 5-col MD5 match** | **50,274 vs 50,274 peaks; BED3 Jaccard 1.0; all summit distances 0 bp** |
| `--low-mem` | **8:49.98 / 42.6 GB** | 22:39.71 / 24.2 GB | **2.57x faster; 61.0% less wall time** | **53,969,811 vs 53,969,811; exact sorted 5-col MD5 match** | **50,274 vs 50,274 peaks; BED3 Jaccard 1.0; all summit distances 0 bp** |

The full baseline wall time includes Chromap SAM to sorted BAM, BAM indexing,
Chromap FRAG output, and standalone MACS3 FRAG peak calling.

BAM record counts differ between the standalone SAM/samtools path and the
integrated dual-fragment BAM path (`108,066,582` baseline records vs
`107,939,622` integrated records). This is reported as variation in
non-biological read-level fields, especially read names and tie-selected
duplicate representatives. Fragment coordinates, barcode, duplicate count, peak
geometry, and summit positions are the parity gates, and they match exactly.

Reusable scripts and recorded outputs:

- Full benchmark runner: [`scripts/benchmarks/run_full_publication_benchmark.sh`](scripts/benchmarks/run_full_publication_benchmark.sh)
- 100K smoke runner: [`scripts/benchmarks/run_100k_publication_smoke.sh`](scripts/benchmarks/run_100k_publication_smoke.sh)
- Benchmark panel and run note: [`plans/2026-05-06-libmacs3-chromap-atac-benchmark-panel.md`](plans/2026-05-06-libmacs3-chromap-atac-benchmark-panel.md)
- Full-run summary table: `plans/artifacts/libmacs3_chromap_atac_panel/20260506T_full_pbmc3k_16t/full_pbmc3k/full_summary.tsv`

## Building & Installing

### From source

```sh
# Initialise submodules (htslib + libMACS3)
git submodule update --init --recursive

# Build chromap CLI + libchromap.a + chromap_callpeaks (libMACS3 standalone)
make
```

Requirements: GCC ≥ 7.3.0, GNU make, OpenMP, zlib, htslib build dependencies (libcurl, libcrypto, libbz2, liblzma, libdeflate). The vendored htslib submodule is built from source as part of `make` so no system htslib is required.

Compile-time switches:

```sh
# Address sanitizer
make asan=1

# Use the legacy single-threaded overflow path instead of the k-way merge default
make LEGACY_OVERFLOW=1
```

### Validation

```sh
# Hermetic S0 smoke matrix (regression cases C01–C11, synthetic fixtures)
make test-smoke

# Single-component validators (pre-spinoff correctness checks; still useful)
./scripts/validate_sam_fix.sh
./scripts/validate_low_mem_fix.sh
./scripts/test_overflow_basic.sh
```

## Sample Commands

### Index build

```sh
chromap -i -r ref.fa -o ref.index
```

### Bulk ATAC (BED + inline narrow peaks, no barcode)

```sh
chromap --preset atac \
  -x ref.index -r ref.fa \
  -1 read1.fq.gz -2 read2.fq.gz \
  -o aln.bed \
  --call-macs3-frag-peaks \
  --macs3-frag-peaks-source memory
```

`--macs3-frag-peaks-source memory` is required for the bulk path because the file-source loader expects a 5-col fragments TSV (duplicate count in column 5), but bulk BED output puts MAPQ in column 5.

### scATAC (BED + inline narrow peaks, barcoded)

```sh
chromap --preset atac \
  -x ref.index -r ref.fa \
  -1 read1.fq.gz -2 read2.fq.gz \
  -b barcode.fq.gz --barcode-whitelist whitelist.txt \
  -o aln.bed \
  --call-macs3-frag-peaks
```

### scATAC (BAM dual + fragments TSV + sidecar + inline peaks)

```sh
chromap --preset atac \
  -x ref.index -r ref.fa \
  -1 read1.fq.gz -2 read2.fq.gz \
  -b barcode.fq.gz --barcode-whitelist whitelist.txt \
  --BAM --sort-bam --write-index \
  --atac-fragments fragments.tsv.gz \
  --atac-fragment-binary-output fragments.bin \
  --call-macs3-frag-peaks \
  -o aln.bam
```

### Native CBQ input (ATAC)

Maps paired-end ATAC/scATAC directly from BINSEQ CBQ without writing intermediate FASTQ. The same options work with `chromap_lib_runner`.

```sh
chromap --preset atac \
  --input-format cbq \
  --read-pair-cbq lane1.reads.cbq,lane2.reads.cbq \
  --barcode-cbq lane1.barcode.cbq,lane2.barcode.cbq \
  --barcode-whitelist whitelist.txt \
  -x ref.index -r ref.fa \
  -o fragments.bed --BED
```

Rules: `--input-format fastq` is the default; `--input-format cbq` requires `--read-pair-cbq`; `--barcode-cbq` count must match `--read-pair-cbq` count; `--barcode-whitelist` requires `--barcode-cbq`; and FASTQ inputs (`-1/-2/-b`) cannot be mixed with CBQ. `--emit-Y-noY-fastq` is supported for CBQ inputs and writes FASTQ sidecars from the decoded sequence/quality buffers. The hermetic CBQ modality matrix covers ChIP, bulk paired, scATAC, Hi-C `.pairs`, BAM/CRAM writer output, `--read-group auto`, and Y/noY FASTQ parity.

Alignment contract: the read-pair and barcode CBQ lanes must be **record-aligned** — record *i* of each lane must be the same original read, in the same order. To make this verifiable, barcoded CBQ inputs must include read names (headers); both lanes are checked to carry headers at startup and read/barcode names are compared per record. Encode with an order-preserving encoder: encoders that reorder records across blocks under parallelism (e.g. `bqtools` at scale) break this alignment.

Full 3K PBMC ATAC timing on the 4-lane fixture, 8 threads, warmed cache, and
BED output directed to `/dev/null`: FASTQ.gz took `3:04.47`, uncompressed CBQ
took `2:52.27`, with identical read/mapping/output counts (`53,969,811`
output mappings). The manifest is under
`plans/artifacts/cbq_atac_full_timing/20260531T081906Z/`.

### ChIP-seq

```sh
chromap --preset chip \
  -x ref.index -r ref.fa \
  -1 read1.fq.gz -2 read2.fq.gz \
  -o aln.bed
```

### Hi-C (pairs format)

```sh
chromap --preset hic \
  -x ref.index -r ref.fa \
  -1 read1.fa -2 read2.fa \
  -o aln.pairs
```

### Multiomic (RNA + ATAC) via STAR Suite

The multiomic integration is invoked through STAR Suite (libchromap is embedded as a worker thread). See the [STAR Suite README — Multiome ATAC row](https://github.com/morphic-bio/STAR-suite#benchmarks) for the full command.

### Coordinate-sorted BAM standalone

```sh
chromap --BAM --sort-bam --write-index \
  -x ref.index -r ref.fa \
  -1 reads.fq -o output.bam

# Cap sort RAM
chromap --BAM --sort-bam --sort-bam-ram 512M \
  -x ref.index -r ref.fa -1 reads.fq -o output.bam

# Deterministic single-thread output across runs
chromap --BAM --sort-bam -t 1 --hts-threads 1 \
  -x ref.index -r ref.fa -1 reads.fq -o output.bam
```

`--write-index` requires `--sort-bam`. Without coordinate sorting, htslib on-the-fly indexing produces invalid/empty index files because records are written in non-deterministic order from multithreaded mapping. See [`docs/sort_spec.md`](docs/sort_spec.md).

### Y-chromosome filtering

```sh
chromap --BAM --sort-bam --write-index \
  --emit-Y-bam --emit-noY-bam \
  -x ref.index -r ref.fa -1 reads_R1.fq -2 reads_R2.fq -o output.bam
```

Produces:

| File | Contents |
|---|---|
| `output.bam` | All reads, coordinate-sorted |
| `output.bam.bai` | Index for main BAM |
| `output.noY.bam` | Non-Y reads only, coordinate-sorted |
| `output.Y.bam` | Y-chromosome reads only, coordinate-sorted |

Y-chromosome detection is case-insensitive and matches: `Y`, `chrY`, `CHR_Y`, `chr_y`. Decoy/random/alt contigs (`chrY_random`, `Y_alt`, `chrY_hap1`) are intentionally excluded.

To emit a normalised list of read names that hit Y, plus split FASTQ output:

```sh
chromap --SAM --emit-Y-read-names --emit-Y-noY-fastq \
  --emit-Y-noY-fastq-compression none \
  -x ref.index -r ref.fa -1 reads_R1.fq -2 reads_R2.fq -o output.sam
```

Read-name list: `<output>.Y.names.txt` by default, normalised (strip leading `@`, stop at first whitespace, strip trailing `/1` or `/2`); no ordering guarantees.

FASTQ split naming:
- Auto-named FASTQ sidecars are written under `y_separated/` next to the primary output path, matching STAR Suite's Y-removal layout. Override this with `--Y-noY-fastq-output-dir`.
- Inserts `_Y` / `_noY` before the last `_R[0-9]+` token if present (e.g., `sample_R1.fastq.gz` → `sample_Y_R1.fastq.gz`).
- Falls back to `Y_reads.mateN.<ext>(.gz)` / `noY_reads.mateN.<ext>(.gz)` otherwise.
- Compression: `gz` (default) or `none`.
- Paired-end routing: if either mate hits Y, both mates go to Y outputs.

Related flags: `--emit-Y-read-names`, `--Y-read-names-output`, `--emit-Y-noY-fastq`, `--emit-Y-noY-fastq-compression {gz|none}`, `--Y-noY-fastq-output-dir`, `--Y-fastq-output-prefix`, `--noY-fastq-output-prefix`.

## Chromap Launchpad (Recipe Builder)

Browser-based recipe builder served from the Chromap Suite MCP server. Pick a workflow, fill in parameters through a guided form, and get a validated, copy-pasteable shell command. The same recipes are exposed as MCP tool calls for agents.

### Quick start

```sh
# Install dependencies (once)
pip install -r mcp_server/requirements.txt

# Start Launchpad + MCP on one port
bash scripts/launchpad_server.sh up

# Open in your browser
# http://127.0.0.1:8765/launchpad/
```

### MCP endpoints

When the server is running on port `8765`:

- `http://127.0.0.1:8765/launchpad/` — Launchpad UI
- `http://127.0.0.1:8765/` — MCP streamable-HTTP endpoint
- `http://127.0.0.1:8765/sse` — MCP SSE endpoint

### MCP client setup

#### VS Code / GitHub Copilot

`.vscode/mcp.json`:

```json
{
  "servers": {
    "chromapSuite": {
      "type": "http",
      "url": "http://127.0.0.1:8765/"
    }
  }
}
```

#### Cursor

```json
{
  "mcpServers": {
    "chromap-suite": {
      "url": "http://127.0.0.1:8765/sse"
    }
  }
}
```

#### Claude

```sh
claude mcp add --transport http chromap-suite http://127.0.0.1:8765/
```

Initial public recipes: `chromap_index`, `chromap_atac_bed`, `chromap_atac_bam_fragments`, `chromap_hic_pairs`. See [`docs/chromap_launchpad.md`](docs/chromap_launchpad.md) and [`mcp_server/README.md`](mcp_server/README.md).

## More Detail

- Project lineage and pre-spinoff fork notes: [HISTORY.md](HISTORY.md)
- Changes by release: [CHANGELOG.md](CHANGELOG.md)
- BAM sort specification: [docs/sort_spec.md](docs/sort_spec.md)
- Launchpad design: [docs/chromap_launchpad.md](docs/chromap_launchpad.md)
- libMACS3 repository: https://github.com/morphic-bio/libMACS3
- STAR Suite (multiomic integration entry point): https://github.com/morphic-bio/STAR-suite
- Chromap Suite preprint: https://github.com/morphic-bio/chromap_suite_paper

## Citing

If you use Chromap Suite, please cite the preprint above plus the Chromap and MACS3 papers it builds on:

> Zhang, H., Song, L., Wang, X., Cheng, H., Wang, C., Meyer, C. A., …, Liu, X. S., Li, H. (2021). Fast alignment and preprocessing of chromatin profiles with Chromap. *Nature Communications*, 12(1), 1-6. https://doi.org/10.1038/s41467-021-26865-w

> Zhang, Y., Liu, T., Meyer, C. A., Eeckhoute, J., Johnson, D. S., Bernstein, B. E., Nusbaum, C., Myers, R. M., Brown, M., Li, W., Liu, X. S. (2008). Model-based Analysis of ChIP-Seq (MACS). *Genome Biology*, 9, R137. https://doi.org/10.1186/gb-2008-9-9-r137

The Chromap QC summary file format (carried forward unchanged in Chromap Suite) is described in:

> Ahmed, O., Zhang, H., Langmead, B., Song, L. (2025). Quality control of single-cell ATAC-seq data without peak calling using Chromap. *bioRxiv*. https://doi.org/10.1101/2025.07.15.664951

## Authors

Chromap Suite extensions copyright Ling-Hong Hung. Preprint co-authors: Ling-Hong Hung, Ka Yee Yeung. The GitHub Contributors graph reflects inherited git history from the pre-spinoff fork; see [AUTHORS.md](AUTHORS.md) for current authorship and acknowledgments, and [HISTORY.md](HISTORY.md) for project lineage.

## Licence

- `libchromap` and the `chromap` CLI: MIT (see [LICENSE](LICENSE)).
- `libMACS3`: BSD-3 — a single source-inspected adaptation for summit edge cases necessitates BSD-3 distribution; see [the libMACS3 repo](https://github.com/morphic-bio/libMACS3) for the methodology.
- MCP server + Launchpad (`mcp_server/`): MIT.
