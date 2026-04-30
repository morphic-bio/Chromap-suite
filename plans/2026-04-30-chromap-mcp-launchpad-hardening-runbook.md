# Runbook: Chromap MCP, Agents, and Launchpad Hardening

Date: 2026-04-30
Status: Stage 0-4 implemented; Stage 5+ proposed
Prerequisites:

- MCP/Launchpad Stage 1 is present.
- `make test-libchromap-core-smoke` passes.
- Optional ENCODE S1 smoke harness is available for real-data checks.

## Goal

Make Chromap's MCP server, Launchpad, and agent documentation a robust
operator surface for reproducible local Chromap work.

The intended end state is:

- MCP and Launchpad read from one recipe registry.
- Recipes have typed inputs, preflight rules, expected outputs, smoke coverage,
  and runtime class metadata.
- Every run emits a durable run manifest under `plans/artifacts/`.
- Agents can discover where docs, recipes, tests, and sister-repo handoffs live.
- Unsafe execution paths are minimized by allowlisting recipes and validating
  parameters before command construction.

## Non-Goals

- Do not change Chromap mapping behavior.
- Do not turn Launchpad into a scheduler.
- Do not run production benchmarks from Launchpad by default.
- Do not make ENCODE downloads part of the default smoke gate.
- Do not merge STAR Suite execution logic into Chromap-suite.
- Do not add downstream Hi-C TAD callers; Hi-C stops at `.pairs`.

## Design Principles

- One source of truth for recipe metadata.
- Serial benchmark policy by default.
- Generated outputs stay under `plans/artifacts/`.
- Default tests remain cheap and hermetic.
- Real-data tests remain explicit and environment-gated.
- Agents do not claim co-authorship in commits or pushes.

## Stage 0: Baseline Audit

Purpose: record the current MCP/Launchpad surface before changing structure.

Estimated effort: 0.5 day.

### Tasks

1. Inventory current workflows:
   - `mcp_server/workflows/*.yaml`
   - Launchpad-visible recipe list
   - MCP tools exposed by `mcp_server/app.py`
   - tests under `mcp_server/tests/`

2. Record current gates:
   - `python3 -m pytest mcp_server/tests -q`
   - `make test-libchromap-core-smoke`
   - optional `tests/run_encode_downsample_smoke.sh`

3. Identify drift:
   - workflow fields duplicated in docs,
   - Launchpad UI assumptions not represented in YAML,
   - preflight checks that live outside recipe metadata,
   - undocumented output artifacts.

### Deliverables

- Short audit note appended to this runbook or placed under
  `plans/artifacts/mcp_launchpad_hardening/<timestamp>/audit.md`.
- No behavior changes.

### Gates

```bash
python3 -m pytest mcp_server/tests -q
make test-libchromap-core-smoke
```

### Implementation Note

Stage 0 audit result:

- Current public workflow YAMLs:
  - `chromap_index`
  - `chromap_atac_bed`
  - `chromap_atac_bam_fragments`
  - `chromap_hic_pairs`
- Existing MCP/Launchpad tests live under `mcp_server/tests/`.
- Current Launchpad behavior is workflow-YAML driven; recipe intent, output
  handoff metadata, benchmark policy, and smoke coverage were not represented
  in one shared metadata layer.
- Stage 1 addresses that drift by adding `mcp_server/recipes/registry.yaml`.

## Stage 1: Recipe Registry

Purpose: define a durable recipe registry that MCP, Launchpad, tests, and docs
can consume.

Estimated effort: 1-2 days.

### Proposed Files

Add:

- `mcp_server/recipes/registry.yaml`
- `mcp_server/recipes/schema.yaml` or a Python schema model
- `mcp_server/tests/test_recipe_registry.py`

Update:

- `mcp_server/workflows/AUTHORING.md`
- `mcp_server/README.md`
- `tests/README.md`

### Registry Fields

Each recipe should declare:

- `id`: stable recipe id, e.g. `chromap_atac_bam_fragments`.
- `title`: human-readable label.
- `purpose`: concise workflow description.
- `runtime_class`: `smoke`, `interactive`, `long`, or `benchmark`.
- `benchmark_policy`: `serial_required`, `parallel_safe`, or `not_benchmark`.
- `command_template`: argv template, not shell string.
- `inputs`: typed input definitions.
- `outputs`: expected output artifacts.
- `preflight`: list of validation rule ids.
- `smoke_coverage`: S0/S1/S2 coverage reference.
- `docs`: README/runbook/doc anchors.
- `handoff_artifacts`: artifacts consumed by STAR Suite or downstream tools.

### Initial Recipes

Move or mirror current workflow YAML metadata into registry entries for:

- `chromap_index`
- `chromap_atac_bed`
- `chromap_atac_bam_fragments`
- `chromap_hic_pairs`

Add metadata-only entries for future recipes:

- `chromap_chip_tagalign`
- `chromap_sorted_bam`
- `chromap_y_noy_split`
- `chromap_macs3_frag_peaks`
- `chromap_lib_runner_parity`

### Acceptance Criteria

- Registry validates in tests.
- Every public Launchpad workflow has a registry entry.
- Every registry entry has expected outputs and smoke coverage metadata.
- No generated artifacts are tracked.

### Implementation Note

Stage 1 added:

- `mcp_server/recipes/registry.yaml`
- `mcp_server/schemas/recipe.py`
- `mcp_server/tools/recipes.py`
- `mcp_server/tests/test_recipe_registry.py`

The registry mirrors the current public workflows and adds metadata-only
planned recipes for ChIP TagAlign, sorted BAM, Y/noY split, MACS3 FRAG peaks,
and `chromap_lib_runner` parity.

### Gates

```bash
python3 -m pytest mcp_server/tests -q
python3 -m mcp_server.app --help
```

## Stage 2: Preflight Framework

Purpose: make validation explicit, reusable, and recipe-driven.

Estimated effort: 1-2 days.

### Proposed Files

Add or update:

- `mcp_server/tools/preflight.py`
- `mcp_server/tests/test_chromap_preflight.py`
- `mcp_server/recipes/preflight_rules.yaml` if rule definitions are separated

### Core Rules

Implement reusable rules for:

- path exists,
- parent directory exists or can be created,
- path is under trusted roots,
- FASTQ read counts are pair-compatible for a bounded sample,
- read1/read2 list lengths match,
- barcode list length matches read1 list length,
- reference FASTA exists,
- Chromap index exists,
- output path differs from secondary artifact paths,
- `--write-index` requires sorted BAM/CRAM output,
- `--atac-fragments` requires paired-end barcoded BAM/CRAM workflow,
- Hi-C recipe stops at `.pairs`,
- S1 ENCODE smoke requires `CHROMAP_GRCH38_REF` and
  `CHROMAP_GRCH38_INDEX`.

### Preflight Output

Return structured results:

- `status`: `pass`, `warn`, or `fail`
- `rule_id`
- `message`
- `path` when applicable
- `suggested_fix` when practical

### Acceptance Criteria

- MCP exposes a preflight call for a recipe and parameter set.
- Launchpad can display pass/warn/fail before command rendering.
- Tests cover invalid paths, missing index, mismatched read lists, output
  collisions, and trusted-root rejection.

### Implementation Note

Stage 2 added:

- `preflight_recipe(recipe_id, params)` in `mcp_server/tools/preflight.py`
- MCP tool `preflight_recipe`
- structured `RecipePreflightCheck` and `RecipePreflightResponse` models
- `mcp_server/tests/test_chromap_preflight.py`

Implemented rule ids include reference/index existence, output parent trust and
writability, read-pair lane count matching, barcode lane count matching,
primary/secondary output collision checks, ATAC fragments requirements, Hi-C
`.pairs` boundary validation, Y/noY format validation, MACS3 peak option
validation, and binary presence checks for planned runner parity recipes.

### Gates

```bash
python3 -m pytest mcp_server/tests/test_chromap_preflight.py -q
python3 -m pytest mcp_server/tests -q
```

## Stage 3: Run Manifests

Purpose: every MCP/Launchpad execution should leave a reproducible record.

Estimated effort: 1 day.

### Proposed Files

Add or update:

- `mcp_server/tools/run_manifest.py`
- `mcp_server/tests/test_run_manifest.py`
- `docs/chromap_launchpad.md`

### Manifest Location

Default:

```text
plans/artifacts/mcp_runs/<YYYYMMDDTHHMMSSZ>_<recipe_id>/run.json
```

### Manifest Fields

Record:

- recipe id and registry version,
- rendered argv,
- shell preview,
- git commit and dirty-worktree flag,
- binary paths and versions,
- input paths,
- reference/index paths,
- output paths,
- artifact directory,
- stdout/stderr log paths,
- start/end UTC timestamps,
- exit code,
- host/user,
- relevant environment variables,
- preflight results,
- smoke or benchmark classification.

### Acceptance Criteria

- Dry runs can emit a manifest with `execution_status: dry_run`.
- Real runs emit logs plus `run.json`.
- Manifest records serial benchmark policy when `runtime_class=benchmark`.
- Tests assert required fields.

### Gates

```bash
python3 -m pytest mcp_server/tests/test_run_manifest.py -q
python3 -m pytest mcp_server/tests -q
```

### Stage 3 Implementation Note

Implemented in `mcp_server/tools/run_manifest.py` with MCP exposure through
`write_recipe_run_manifest` and Launchpad launch integration through
`mcp_server/launchpad/api.py`.

Run manifests are written under `plans/artifacts/mcp_runs/.../run.json` and
record the recipe id, registry version, rendered argv, shell preview, git
state, binary metadata, input/reference/output paths, artifact/log paths,
host/user context, relevant environment variables, preflight results, and
smoke/benchmark classification. Launchpad now captures stdout/stderr into the
same artifact directory instead of discarding them.

Validation:

```bash
python3 -m pytest mcp_server/tests/test_run_manifest.py -q
python3 -m pytest mcp_server/tests -q
```

## Stage 4: Launchpad From Registry

Purpose: make Launchpad forms generated from recipe metadata so docs, MCP, and
UI do not drift.

Estimated effort: 2-3 days.

### Proposed Files

Update:

- `mcp_server/launchpad/api.py`
- `mcp_server/launchpad/static/app.js`
- `mcp_server/launchpad/static/index.html`
- `mcp_server/launchpad/static/style.css`
- `mcp_server/tests/test_launchpad_chromap.py`

### UI Behavior

Launchpad should:

- list recipes from registry,
- show runtime class and benchmark policy,
- render typed fields from recipe inputs,
- show expected outputs,
- run preflight before execution,
- show command argv and shell preview,
- write or link to run manifests,
- make S1 ENCODE and benchmark recipes visibly opt-in.

### Acceptance Criteria

- No recipe-specific hardcoding is needed for the existing public recipes.
- Adding a metadata-only recipe does not expose it for execution until marked
  `enabled: true`.
- Launchpad API tests cover recipe listing, form schema, render, preflight, and
  dry-run manifest creation.

### Gates

```bash
python3 -m pytest mcp_server/tests/test_launchpad_chromap.py -q
bash scripts/launchpad_server.sh up
curl -fsS http://127.0.0.1:8765/launchpad/api/recipes >/tmp/chromap_recipes.json
bash scripts/launchpad_server.sh down
```

### Stage 4 Implementation Note

Implemented recipe-native Launchpad endpoints in `mcp_server/launchpad/api.py`:

- `GET /launchpad/api/recipes`
- `GET /launchpad/api/recipes/{recipe_id}/schema`
- `GET /launchpad/api/recipes/{recipe_id}/describe`
- `POST /launchpad/api/recipes/{recipe_id}/preflight`
- `POST /launchpad/api/recipes/{recipe_id}/render`
- `POST /launchpad/api/recipes/{recipe_id}/manifest`

The browser now populates its recipe selector from the registry, renders typed
fields from `RecipeEntry.inputs`, shows runtime/benchmark metadata and expected
outputs, runs registry preflight before command rendering/execution, and writes
dry-run manifests from the recipe route. Metadata-only recipes can be listed
through the explicit planned/long opt-in but their schema/render/preflight
routes reject execution until `enabled: true`.

Validation:

```bash
python3 -m pytest mcp_server/tests/test_launchpad_chromap.py -q
python3 -m pytest mcp_server/tests -q
python3 -m mcp_server.app --help
git diff --check
```

## Stage 5: MCP Execution Hardening

Purpose: keep MCP useful while avoiding arbitrary shell execution and unsafe
paths.

Estimated effort: 1-2 days.

### Tasks

1. Ensure executable actions use allowlisted recipe ids.
2. Construct commands as argv arrays.
3. Reject shell metacharacter-only command construction paths.
4. Enforce trusted roots for inputs and outputs.
5. Require preflight pass or explicit override for warnings.
6. Block `runtime_class=benchmark` unless the caller explicitly opts in.
7. Add timeouts for local execution.
8. Store logs under the run manifest artifact directory.

### Acceptance Criteria

- MCP cannot execute an arbitrary shell string through recipe APIs.
- Invalid paths fail before command execution.
- Long/benchmark recipes require explicit opt-in.
- Tests cover command injection attempts and untrusted output locations.

### Gates

```bash
python3 -m pytest mcp_server/tests -q
```

## Stage 6: Agent and User Documentation

Purpose: make the repository self-explanatory for future agents and users.

Estimated effort: 0.5-1 day.

### Update

- `AGENTS.md`
- `mcp_server/README.md`
- `docs/chromap_launchpad.md`
- `docs/index.md`
- `tests/README.md`

### Required Content

Document:

- recipe registry location,
- how to add a recipe,
- required tests for new recipes,
- artifact sandbox rules,
- benchmark policy,
- S0/S1/S2 test tiers,
- STAR Suite sister-repo handoff points,
- no co-authorship policy for agents,
- Launchpad/MCP execution safety model.

### Acceptance Criteria

- A new agent can find the registry, docs, tests, and artifact policy from
  `AGENTS.md`.
- Recipe authoring docs list every required registry field.
- Docs clearly state that S1 ENCODE downloads are opt-in.

## Stage 7: STAR Suite Handoff Metadata

Purpose: document integration boundaries without coupling Chromap execution to
STAR Suite.

Estimated effort: 0.5-1 day.

### Add Metadata For

- ATAC BAM
- fragments TSV or binary sidecar
- MACS3 peaks
- ATAC evidence TSV
- Hi-C `.pairs`

### Acceptance Criteria

- Registry records which artifacts are Chromap-owned handoff points.
- Docs say STAR Suite lives at `/mnt/pikachu/STAR-suite`.
- No STAR Suite workflow becomes a Chromap Launchpad default.

## Stage 8: CI and Local Gates

Purpose: keep the system stable as recipes grow.

Estimated effort: 1 day.

### Default Local Gate

```bash
python3 -m pytest mcp_server/tests -q
make test-libchromap-core-smoke
git diff --check
```

### Optional Real-Data Gate

```bash
CHROMAP_GRCH38_REF=/path/to/genome.fa \
CHROMAP_GRCH38_INDEX=/path/to/genome.index \
ENCODE_SKIP_PREPARE=1 \
tests/run_encode_downsample_smoke.sh
```

### Optional Launchpad Gate

```bash
bash scripts/launchpad_server.sh up
curl -fsS http://127.0.0.1:8765/launchpad/ >/tmp/chromap_launchpad.html
curl -fsS http://127.0.0.1:8765/launchpad/api/recipes >/tmp/chromap_recipes.json
bash scripts/launchpad_server.sh down
```

### Acceptance Criteria

- Default gate remains cheap and hermetic.
- S1 ENCODE gate is documented but opt-in.
- Benchmark recipes are excluded from default CI/local gates.

## Suggested Implementation Order

1. Stage 0: audit current state.
2. Stage 1: recipe registry and validation tests.
3. Stage 2: preflight framework.
4. Stage 3: run manifests.
5. Stage 4: Launchpad generation from registry.
6. Stage 5: MCP execution hardening.
7. Stage 6: docs and `AGENTS.md`.
8. Stage 7: STAR Suite handoff metadata.
9. Stage 8: CI/local gate cleanup.

The first practical pull request should include Stages 0-1 only. That keeps the
change reviewable and gives the later preflight/UI work a stable metadata
contract.
