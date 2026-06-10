# Agent-era strategy + future directions — handoff

Date: 2026-06-09

Entry point for the strategic direction that emerged while selecting a benchmark
dataset for the suites. Read this, then the durable writeup
[`docs/agent_era_composable_bioinformatics.md`](../docs/agent_era_composable_bioinformatics.md),
and `AGENTS.md` before editing. **Next working focus is the CAT-ATAC analysis
(§5), not the items below** — these are the roadmap to return to.

## TL;DR

We converged on a product/paper strategy:
1. **nf-core participation** — package the transcriptomics/ATAC capabilities as
   **single-binary nf-core pipelines** (one input→output unit; DAG-free for the
   consumer). Reaches the human/workflow community *and* makes the tools agent-shaped.
2. **Agent side** — an **MCP server exposing a small set of major composable
   recipes**, agnostic to workflow engine and scripting language; the LLM handles
   glue/exceptions. MCP is the lead surface (converging cross-vendor protocol).
3. **Canonical base + translation engine** — keep a simple vendor-neutral recipe as
   the source of truth; project it to MCP / nf-core / (later) per-vendor skills. The
   engine also ingests already-curated workflows (nf-core) → recipes, bootstrapping a
   discovery index ("Pedro's page for the agent era") without a curation treadmill.
4. **One opinionated implementation paper** — for visibility, NIH (R01/R21/SBIR,
   ITCR-shaped) and VC. Thesis baked into the Discussion; CAT-ATAC as the demo.

## Status (be honest about what exists)

- **Built:** Chromap-suite (deterministic alignment, byte-identical libMACS3 narrow
  peaks, fragments/matrices, CBQ input, Hi-C `.pairs`, `libchromap`, AEV1 sidecar,
  MCP server, Launchpad recipe registry); STAR-suite (STARsolo `GeneFull`, feature
  barcoding/guides, MuData); Bwb (NCI-funded workflow engine); a scheduler + common
  interface (preprint, per team).
- **Planned / in progress:** single-binary nf-core pipelines; the multi-target
  *recipe → {MCP, nf-core, skills}* projector; the *workflow → recipe* translation
  at scale + fidelity validation; the agent-discovery index as a named resource.
- **Pending empirical:** the CAT-ATAC end-to-end run (data downloaded, not yet run).

## Decisions made this session

- **Dataset:** GSE288996 (CAT-ATAC) is the benchmark — open, complete, cell-line,
  **multiome Perturb-seq** with **MACS2 (≡MACS3) peaks** (Signac `CallPeaks`,
  confirmed in `ucsf-lgr/catatac_public`). It is the open/complete analog of Norman
  (whose ATAC R3 was never deposited). Norman download **deleted**.
- **Format strategy:** vendor-neutral canonical recipe is the source of truth; **MCP
  now** (protocol converges), **nf-core as a cheap projection** (human reach + hedge),
  **vendor "skills" deferred** (packaging diverges across Cursor/Gemini/OpenAI;
  generate only when one consolidates). Skills add to discovery only at large corpus
  scale (progressive disclosure), and that benefit is also available from a deferred
  MCP server; the MCP server is mandatory as the **execution** endpoint regardless.
- **Granularity rule for translation:** one input→output unit per workflow by default
  (DAG below the interface) + exposed seams for cross-workflow composition; inherit
  determinism/provenance from source; validate fidelity. Avoid both the opaque
  monolith and the 40-piece explosion.
- **Hi-C/ChIP are recipe-orchestrated, not new C++.** Chromap maps Hi-C (`.pairs`)
  and ChIP (BAM) via legacy presets; ChIP *peak calling* is full MACS3 (input
  control) or public ENCODE tracks — **not** libMACS3 (which is ATAC-narrow,
  treatment-only). State this honestly; do not claim libMACS3 calls ChIP peaks.

## Future directions (roadmap)

1. **nf-core single-binary pipelines** for transcriptomics + ATAC; treat as a
   generated surface of the canonical recipe.
2. **MCP recipe layer**: curate the *major* composable recipes (ATAC→peaks/fragments,
   RNA→matrix, guide→counts, MuData assembly, Hi-C→pairs→cooler, ChIP→peaks);
   deferred/progressive discovery; recipe = intent + typed I/O + determinism +
   footprint + execution binding.
3. **Translation/projection engine**: *workflow → recipe* (seed from nf-core) and
   *recipe → {MCP, nf-core, skills}*; fidelity validation harness.
4. **Discovery index** ("Pedro's page for agents"): machine-readable, AI/usage-curated,
   pointing at deterministic units; seed narrow.
5. **Paper**: implementation + perspective (see writeup). Suggested path —
   **bioRxiv preprint now** (immediate funding utility) → Patterns / Genome Biology
   Opinion (accessible, read by NIH + industry) → Nature Methods Comment as a
   visibility stretch. Bundle with the suites' methods papers + scheduler preprint as
   a portfolio. Lead with the practitioner credential (Bwb) + evidence (IGVF
   seqspec-not-Kallisto, byte-identical MACS3, CAT-ATAC); include the self-critique
   (provenance drift). Two funder faces from one thesis: reproducibility (NIH) and
   blast-radius adoption / discovery-layer (SBIR/VC).

## CAT-ATAC dataset (next focus)

- Location: `/mnt/pikachu/catatac_gse288996/` (subset, ~102 GB).
  - `fastq/ATAC/` — SRR32265760 (DMSO rep1), SRR32265762 (DASA rep1); each
    `_1`=100bp genomic, `_2`=24bp **barcode**, `_3`=100bp genomic (complete).
  - `fastq/GEX/` — SRR32265752 (DMSO), SRR32265754 (DASA); `_1/_2` = 101/101.
  - `fastq/guide/` — SRR32265756 (DMSO), SRR32265758 (DASA); `_1/_2` = 100/100.
  - `metadata/subset_manifest.tsv`, `metadata/ena_PRJNA1220572.tsv`,
    `download_subset.sh` (resumable; tops up rep2 / iPSC if re-pointed).
- Full set is 5 samples × {RNA, ATAC, guide} (K562 DMSO×2, DASA×2 + iPSC); only the
  DMSO/DASA rep1 subset was pulled. Add reps/iPSC via the script if needed.
- **Downstream the paper authors did** (after our processing layer; informs what the
  MuData handoff must feed): custom ZINB guide→cell assignment, Seurat SCTransform +
  cell-cycle regression, **Mixscape**, FindMarkers (DEG), Signac LR (DA), chromVAR /
  JASPAR motifs + footprinting, **Pando GRN**, validated against K562 H3K27ac/P300
  ChIP + Hi-C. Implication: emit fragments + **raw guide×cell counts** as first-class
  outputs; leave assignment + GRN/DA/motif to the recipe layer.

## Open questions

- nf-core single-binary pipeline shape (one module vs a thin pipeline) and how the
  determinism/provenance metadata surfaces in nf-core's schema.
- Canonical recipe schema details (the five properties) — to be drafted before the
  MCP layer is generalized.
- Whether to validate translation fidelity against a small curated workflow set in
  Phase I (recommended) and which workflows.
- Paper venue final call (Patterns vs Genome Biology vs NMeth Comment) and author list.

## Constraints

Per `AGENTS.md` / repo norms: independent buildability of Chromap-suite; no
STAR-specific assumptions in Chromap defaults; benchmarks serial by default; commit
only when asked and branch first (repo is under concurrent multi-agent development);
no AI co-authorship attribution in commits/PRs. These docs are untracked working
files — not committed.
