# Composable, deterministic bioinformatics for the agent era

*Single-binary pipelines, a vendor-neutral recipe layer, and what reproducibility
must mean when the analyst is a human–AI team.*

Status: **design / perspective writeup (working draft).** Some components described
here are implemented (Chromap-suite, STAR-suite, the recipe registry, the MCP
server, Bwb); others are planned or in progress (nf-core single-binary pipelines,
the multi-target translation/projection layer, the CAT-ATAC end-to-end
demonstration). Claims tied to results are marked **[pending]** until the run that
backs them exists. This document is the conceptual + architecture basis for a
planned implementation-with-perspective paper; the actionable roadmap lives in
[`plans/2026-06-09-agent-era-strategy-and-future-directions-handoff.md`](../plans/2026-06-09-agent-era-strategy-and-future-directions-handoff.md).

---

## 1. The argument in one paragraph

Reproducibility was never the goal; it was a *mechanism* for the real goals —
verification, reuse, error-correction, accountability. Pre-AI, "re-run the exact
pipeline and get the same bits" was the cheapest available mechanism because the
pipeline *was* the analysis. In a world where the analysis engine is increasingly
a human–AI team working adaptively, that mechanism breaks — but the goals do not.
The replacement is **deterministic, composable primitives whose actual invocation
is faithfully captured**: the binary is reproducible, the recipe is the provenance
marker, and the human–AI team composes. The same properties that make this
reproducible (determinism, captured provenance, standard interfaces, low footprint)
are exactly the properties that make tools *adoptable by agents* — so adoption and
reproducibility converge on one design.

## 2. What "reproducible" is for — goal vs. mechanism, and the layer it applies to

The word fuses several things (computational reproducibility, replicability,
robustness — the terms are contested; anchor to NAS 2019 / Goodman 2016). The
clarifying cut is **goal vs. mechanism**: computational reproducibility ("re-run
the exact pipeline") is a proxy that got mistaken for the goal. The second cut is
**which layer** you mean, because each needs a different mechanism:

| Layer | Right reproducibility mechanism |
|---|---|
| Data (raw measurements) | archival + standard format |
| **Processing** (raw → peaks/matrices) | **determinism** — bit-identical, achievable and valuable |
| **Analysis** (matrices → DE/GRN) | **captured provenance** of an exploratory human–AI process |
| Interpretation (→ claims) | replication |

Processing *should* be frozen-deterministic; analysis *cannot* be (it is
exploratory) and should not pretend to be. Conflating these layers is most of the
confusion in the field.

## 3. Two human-era answers: IGVF and MorPhiC

The field's two leading functional-genomics consortia solve *different aspects* of
reproducibility, and — tellingly — **neither does it by standardizing
analysis-execution**:

| | IGVF | MorPhiC |
|---|---|---|
| Structure | Federated, many labs/assays/cell types | Coordinated, uniform production |
| Aspect solved | **Interoperability / reuse** across heterogeneous producers | **Internal consistency / comparability** of one corpus |
| Reproducibility lever | Standardize **interfaces** (seqspec, metadata, formats) | Standardize **production** (defined lines, perturbation) + **deterministic processing** |
| Why that lever | *Cannot* control execution across labs → describe it uniformly | *Can* control execution → make it deterministic |
| What it does **not** buy | cross-lab conclusion consistency | breadth / generalizability |

IGVF is the instructive case: its centrally-pushed (Kallisto-based) workflows did
not achieve broad adoption — data-production centers had principled reasons to keep
their own processing — so IGVF's durable reproducibility lever became **seqspec, an
interface**, not an enforced pipeline. The lesson generalizes: *you can standardize
what people will adopt (interfaces), not what you wish they would run (execution).*
Both consortia locate reproducibility **below the analysis layer** — IGVF in the
interface substrate, MorPhiC in the production substrate — and leave analysis to
flexible composition. That convergence is the empirical heart of the argument here.

## 4. Implementation

### 4.1 Deterministic single-binary pipelines

Chromap-suite and STAR-suite collapse the common multi-step path into single
deterministic binaries that emit standard interchange formats:

- **Chromap-suite** — alignment → fragments → in-process narrow peaks
  (libMACS3, **byte-identical to MACS3 v3.0.3**) → matrices; BED/BAM/CRAM/TagAlign,
  and Hi-C `.pairs` via the same binary; an embeddable `libchromap` API; native CBQ
  (BINSEQ) input as a FASTQ drop-in; an AEV1 ATAC fragment sidecar.
- **STAR-suite** — STARsolo RNA (intron-inclusive `GeneFull` for single-nucleus),
  feature-barcoding / guide counts, and MuData assembly.

Determinism is the load-bearing property: it is *because* the primitives are
bit-identical and versioned that loose, adaptive composition above them remains
reproducible. The byte-identical MACS3 work is therefore not a correctness
footnote — it is the foundation that licenses everything in §6.

### 4.2 nf-core participation: discoverability for humans *and* agents

Contributing the transcriptomics/ATAC capabilities as **single-binary nf-core
pipelines** [pending] reaches the existing human/workflow community on its own
terms and makes the tools discoverable in the channel that community already
trusts. A single-binary pipeline is also the right shape for an agent: the whole
pipeline is one input→output unit, the multi-step complexity executing *below* the
interface. nf-core participation is treated as **one generated surface** of the
canonical recipe (§4.4), not a separate fork.

### 4.3 A vendor-neutral recipe layer + MCP server

For the agent side: an **MCP server exposing a small set of major composable
recipes**, deliberately **agnostic to workflow engine and scripting language** —
the recipe declares intent, typed I/O contract (composability), determinism/version,
and footprint; the LLM handles glue and exceptions. MCP is the lead surface because
it is the *converging cross-vendor protocol* (adopted across Anthropic, OpenAI,
Google, Microsoft, Cursor, and others). Discovery at scale uses progressive /
deferred tool loading so a large catalog stays within an agent's finite budget;
execution is the MCP tool that calls the scheduler.

### 4.4 The translation engine as a multi-target projector

The canonical artifact is a **simple, vendor-neutral recipe** (the source of truth).
A translation/projection engine — building on the Bwb workflow lineage — maps:

- *workflows → recipes* (bootstrap content from already-curated corpora such as
  nf-core, riding their curation instead of re-curating), and
- *recipes → surfaces*: MCP tool descriptors (lead), nf-core modules/pipelines
  (human reach), and per-vendor "skills" formats **only when one consolidates**.

Two design rules make this work rather than become "nf-core + a wrapper":
1. **Granularity hierarchy.** Translate each workflow into a *single* input→output
   unit by default (fits an agent's budget; DAG executes below the interface), and
   *also* expose the key intermediate I/O seams for cross-workflow composition.
   Footprint metadata tells the agent which granularity to use.
2. **Inherit determinism/provenance** from the source (e.g., nf-core versioning +
   containers) so translated units carry the reproducibility guarantee, and
   validate **fidelity** (translated unit reproduces the original's output).

This makes the stack **format-agnostic by construction**: bet on the converging
protocol (MCP) and a vendor-neutral base; treat divergent per-vendor packaging
(Claude/Cursor/Gemini "skills", GPTs) as cheap, deferred projections.

## 5. Demonstration: CAT-ATAC multiome Perturb-seq, end-to-end **[pending run]**

The demonstration substrate is **GSE288996 (CAT-ATAC)** — an open, complete,
cell-line **multiome Perturb-seq** (CRISPR + RNA + ATAC) in K562 (CRISPRi,
dasatinib resistance) and iPSC (CRISPRa). The reviewer-runnable K562 subset
(DMSO rep1 + Dasatinib rep1; ATAC verified complete `R1/R2-barcode/R3` via
`fasterq-dump --include-technical`) is downloaded to `/mnt/pikachu/catatac_gse288996`.

The planned demonstration:
- Chromap-suite maps the ATAC arm → fragments + libMACS3 narrow peaks (the dataset
  itself used MACS2 via Signac `CallPeaks` — peaks equivalent to MACS3).
- STAR-suite maps RNA (`GeneFull`) + guide counts → MuData with the peak matrix.
- Hi-C and ChIP **validation arms** are added *through recipes* using legacy
  Chromap mapping (`.pairs` → cooler) and public K562 H3K27ac/P300 tracks (or full
  MACS3), exercising the cross-modal "predict (GRN) and validate (3D/marks)" loop.
- An agent composes the whole run by **discovering recipes via MCP within its
  budget**, where the equivalent raw multi-module DAG overruns it — the
  self-instantiating evidence that findability + composability drive completion.

Results to be populated once the run is executed; until then this section states a
plan, not a finding.

## 6. The opinionated perspective

### 6.1 Binary-not-workflow; interface-not-execution
Collapsing the common path into a deterministic binary deletes most of the
orchestration surface, leaving a short, legible recipe an agent can author and adapt.
Standardize the interface (formats/metadata) and the canonical recipe; generate the
surfaces; do not enforce one execution.

### 6.2 Adoption at the blast radius — feasibility, not preference
An agent has a finite blast radius: bounded search, context, and iterations. It
returns the solution that *fits*, not the one that is marginally cheaper. A
40-module workflow that the agent must orchestrate is often *outside* the radius;
a findable binary + composable recipes behind a single input→output interface fits.
This is a feasibility threshold, and it **binds tighter as tasks get harder** — so
the composable substrate is not merely competitive at scale, it becomes the only
feasible one for complex agent-driven analysis. Acceptance is thus driven by what
the agent can complete, which routes around — rather than fights — the human trust
barrier. Discovery is the *first* gate of the blast radius: the MCP registry is what
brings a unit *into* the radius at all.

### 6.3 Determinism + captured provenance = post-AI reproducibility
Loose, adaptive composition is reproducible *only* because the primitives are
deterministic and versioned, and the actual invocation is captured. Provenance must
be **auto-captured from execution**, not hand-reconstructed — which is where the
human–AI model can be *stronger* than a curated methods section.

### 6.4 Adoption ≡ reproducibility
The properties that win the agent (composable, deterministic, standard-I/O, low
footprint, discoverable) are the same properties that keep the science reproducible
and make frictionless agent adoption *safe*. One design serves both.

### 6.5 A discovery layer for agents — "Pedro's page for the agent era"
Early curated indexes (Pedro's BioMolecular Research Tools) won by being first and
useful, and died on the manual-curation treadmill. The agent-era index escapes that
fate by being **machine-readable** (agents query contracts, not click links),
**maintained by an AI/usage-driven flywheel** (translated and validated as agents
use it), and **pointing at deterministic units** (which do not rot). nf-core
supplies the *quality/curation model* to borrow; the index is the *role* to own.

### 6.6 Protocol converges, packaging diverges
MCP is the converging cross-vendor *protocol*; vendor "skills" are diverging
*packaging*. Lead with MCP and a vendor-neutral canonical recipe; defer skills to a
generated projection. The translation engine being a *multi-target projector* turns
vendor divergence from a risk into a non-issue by design.

### 6.7 Credibility, risks, and honest limits
- **Credential.** These claims come from practitioners who built and sustained a
  funded workflow engine (Bwb): the critique of execution-standardization is an
  *evolution of our own work*, not an outsider's dismissal. Bwb's goals
  (reproducibility, accessibility) persist; the mechanisms evolve (GUI→AI;
  container-replay→deterministic-primitives+captured-provenance), and Bwb remains a
  backend behind the common scheduler interface.
- **The new model's own failure mode** is provenance drift / automation bias —
  answered by deterministic primitives + auto-captured provenance + human-in-the-loop
  verification. A perspective that critiques the old model while hand-waving the new
  model's risk is the one thing that sinks it.
- **Bootstrap / fidelity.** A registry is worthless empty; we seed from our own
  tools and from translated, already-curated corpora, start narrow, and let the
  blast-radius advantage pull adoption. Translation fidelity must be validated.
- **Scope honestly.** The blast-radius argument is strongest for *agent-driven*
  end-to-end analysis; if a human pre-wires a pipeline, the constraint does not bite.
  Workflow engines still earn their place at standardized high-throughput scale —
  hence the scheduler's common interface (recipes today, any engine when warranted).

## 7. Why this matters now

The composable, deterministic, captured-provenance design is simultaneously the
public-good case (reproducibility preserved as the mode of production changes) and
the market case (own the discovery layer agents route through). It is also a
durable posture: built on a converging open protocol and a vendor-neutral base, with
a translation engine that reaches every surface — infrastructure, not a bet on one
agent vendor.
