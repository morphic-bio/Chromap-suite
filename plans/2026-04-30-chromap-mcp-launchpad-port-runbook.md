# Runbook: Chromap MCP Server + Launchpad Port

Date: 2026-04-30
Status: proposed two-stage implementation plan
Source reference: STAR-suite `mcp_server/` and Launchpad

## Goal

Bring the useful STAR-suite MCP server and Launchpad functionality into
Chromap-suite so agents and users can discover Chromap workflows, validate
parameters, render commands, and use a browser recipe builder.

The port should start with a small, reliable Chromap-specific surface rather
than copying every STAR workflow and then pruning later. STAR-suite's MCP
framework is reusable; the workflow catalog, build tools, docs, and validation
rules need Chromap-specific adaptation.

## Non-Goals

- Do not change Chromap mapping behavior.
- Do not make Launchpad execute long production benchmarks in the MVP.
- Do not copy STAR-only workflows, modules, or labels into Chromap docs.
- Do not require libMACS3 or STAR-suite to be installed for basic Chromap
  Launchpad use.
- Do not add remote execution beyond the existing local MCP server model.

## Existing STAR-suite Surface To Reuse

STAR-suite `mcp_server/` is approximately 22.6k lines across 68 source/config/UI
files. Reusable pieces:

- MCP app and transport: `app.py`
- Config and schema models: `config.py`, `schemas/`
- Workflow validation/rendering: `tools/workflows.py`
- Generic discovery/preflight/execution scaffolding: `tools/`
- Browser Launchpad static app: `launchpad/static/`
- Launchpad API layer: `launchpad/api.py`
- Docker and local server scripts, after rebranding
- Workflow authoring documentation and tests, after Chromap-specific fixtures

Chromap-specific pieces to write:

- Workflow YAMLs for Chromap commands.
- Build tool target mapping for `make`, `make libchromap.a`, and optionally
  `chromap_lib_runner`.
- Preflight checks for FASTA, index, FASTQs, paired inputs, barcode inputs,
  output collisions, temp dirs, and optional MACS3 peak outputs.
- README/docs language and Launchpad UI labels.

## Stage 1: MVP Port

Target effort: 2-3 focused days.

Purpose: land a working Chromap MCP + Launchpad recipe renderer with a small
workflow catalog and tests. This should be useful immediately, but conservative:
render and validate commands first; execution can remain limited to safe local
smokes.

### Stage 1 File Scope

Add:

- `mcp_server/`
- `mcp_server/workflows/chromap_index.yaml`
- `mcp_server/workflows/chromap_atac_bed.yaml`
- `mcp_server/workflows/chromap_atac_bam_fragments.yaml`
- `mcp_server/workflows/chromap_hic_pairs.yaml`
- `scripts/launchpad_server.sh`
- `docs/chromap_launchpad.md`

Update:

- `README.md`
- `docs/index.md`

Optional if needed:

- `mcp_server/workflows/chromap_lib_runner_atac.yaml`

### Stage 1 Tasks

1. Copy STAR-suite `mcp_server/` into Chromap-suite.

2. Remove generated/cache artifacts:
   - `__pycache__/`
   - `.pytest_cache/`
   - stale STAR-specific temporary files

3. Rebrand the server:
   - "STAR-suite" -> "Chromap-suite"
   - "STAR Launchpad" -> "Chromap Launchpad"
   - default public recipe filter from `star_*` to `chromap_*`
   - Docker labels and README text

4. Reduce `mcp_server/config.yaml` to Chromap defaults:
   - `repo_root: /mnt/pikachu/Chromap-suite`
   - `artifact_log_root: /mnt/pikachu/Chromap-suite/plans/artifacts`
   - trusted roots: repo root, `/tmp`, `/storage`, `/mnt/pikachu`
   - no STAR datasets by default
   - only Chromap workflows listed under `workflows:`

5. Adapt build tools:
   - `build_chromap(target="chromap|libchromap|runner", clean=false)`
   - `check_build_status(target=...)`
   - `ensure_fresh_build(target=...)`
   - map targets to:
     - `make`
     - `make libchromap.a`
     - `make chromap_lib_runner` if the target exists locally

6. Add core public workflow YAMLs:
   - `chromap_index`: `chromap -i -r REF -o INDEX`
   - `chromap_atac_bed`: ATAC paired-end BED/fragments output
   - `chromap_atac_bam_fragments`: BAM/CRAM primary output plus
     `--atac-fragments`/binary sidecar path where applicable
   - `chromap_hic_pairs`: Hi-C pairs output

7. Add minimal Chromap preflight checks:
   - FASTA path exists for mapping/index workflows.
   - index path exists for mapping workflows.
   - read1/read2 list counts match when both are present.
   - barcode list count matches read1 count when barcode is present.
   - output path is not identical to secondary fragments path.
   - output parent directory exists or is creatable under trusted roots.
   - `--write-index` requires `--sort-bam`.
   - `--atac-fragments` requires paired-end barcoded BAM/CRAM workflow.

8. Trim tests to the Chromap MVP:
   - config load tests
   - workflow schema tests
   - workflow render tests for each public recipe
   - Launchpad API smoke
   - build tool dry-run tests

9. Update docs:
   - README: add "Chromap Launchpad" section.
   - `docs/chromap_launchpad.md`: quick start, workflow list, limitations.
   - `mcp_server/README.md`: Chromap-specific MCP tools and examples.

### Stage 1 Gates

Build and static checks:

```bash
cd /mnt/pikachu/Chromap-suite
python -m pytest mcp_server/tests -q
python -m mcp_server.app --transport stdio --help
```

Launchpad smoke:

```bash
bash scripts/launchpad_server.sh up
curl -fsS http://127.0.0.1:8765/launchpad/ >/tmp/chromap_launchpad.html
curl -fsS http://127.0.0.1:8765/launchpad/api/workflows >/tmp/chromap_workflows.json
bash scripts/launchpad_server.sh down
```

Workflow render checks:

- `chromap_index` renders one argv containing `chromap -i -r ... -o ...`.
- `chromap_atac_bed` renders paired reads and `--preset atac`.
- `chromap_atac_bam_fragments` renders BAM/CRAM plus secondary fragments.
- `chromap_hic_pairs` renders `--preset hic --pairs`.

Acceptance:

- Launchpad loads in browser.
- Public recipes are Chromap-only by default.
- Rendered commands are valid shell previews.
- No STAR-suite workflow remains exposed by default.
- No generated cache files are tracked.

## Stage 2: Production Hardening

Target effort: 5-7 focused days after Stage 1.

Purpose: turn the MVP into a reliable Chromap workflow assistant with stronger
validation, fixture-backed tests, optional local execution, benchmark-oriented
recipes, and complete documentation.

### Stage 2 File Scope

Add or expand:

- `mcp_server/workflows/chromap_atac_macs3_frag_peaks.yaml`
- `mcp_server/workflows/chromap_atac_binary_sidecar.yaml`
- `mcp_server/workflows/chromap_lib_runner_atac.yaml`
- `mcp_server/workflows/chromap_bam_sort_index.yaml`
- `mcp_server/tests/test_chromap_preflight.py`
- `mcp_server/tests/test_chromap_workflow_e2e.py`
- `docs/chromap_mcp_workflows.md`
- `docs/chromap_launchpad_screenshots/` if screenshots are added

Update:

- `README.md`
- `docs/index.md`
- `docs/chromap.html` only if generated docs are intentionally refreshed
- `mcp_server/README.md`
- `mcp_server/workflows/AUTHORING.md`

### Stage 2 Tasks

1. Add benchmark/advanced workflow recipes:
   - ATAC BAM + binary sidecar
   - ATAC + libMACS3 FRAG peaks
   - low-memory ATAC workflow
   - `chromap_lib_runner` parity/debug workflow
   - BAM sort/index workflow

2. Strengthen preflight:
   - inspect file suffixes and output format compatibility
   - estimate disk usage for BAM, fragments/sidecar, and temp spill files
   - warn on `.gz` secondary fragments when binary sidecar is expected
   - validate MACS3 peak options:
     - peaks output path required
     - summits output path required
     - p-value > 0
     - min length > 0
     - max gap >= 0
   - validate low-memory temp directory and available space

3. Add optional local execution:
   - keep execution disabled or conservative by default
   - allow only allowlisted scripts/workflows
   - enforce trusted roots and timeout
   - support collecting logs and output manifests

4. Add fixture-backed tests:
   - synthetic tiny index build
   - tiny ATAC paired-end render/preflight
   - tiny Hi-C render/preflight
   - invalid path and invalid option tests
   - output collision tests

5. Improve Launchpad UX:
   - Chromap-specific workflow descriptions
   - grouped parameter panels: inputs, output, barcode, BAM/CRAM, MACS3, low-memory
   - examples for ATAC, Hi-C, index
   - save/load parameter JSON retained from STAR Launchpad

6. Documentation:
   - quick start
   - supported workflows table
   - authoring guide for adding new Chromap recipes
   - preflight rules
   - local execution policy
   - known limitations

7. Optional container support:
   - update Dockerfile labels and config
   - document `docker compose up`
   - smoke Launchpad in container

### Stage 2 Gates

Test suite:

```bash
cd /mnt/pikachu/Chromap-suite
python -m pytest mcp_server/tests -q
```

Build tools:

```bash
python -m pytest mcp_server/tests/test_build.py -q
make -j8
make -j8 libchromap.a
```

Workflow validation:

- Every public `chromap_*` workflow validates with a minimal valid parameter set.
- Every public `chromap_*` workflow rejects at least one known-invalid parameter set.
- Rendered argv preserves deterministic flag order.

Launchpad:

- Browser UI loads.
- Workflow dropdown defaults to Chromap recipes.
- Each public recipe renders a command.
- Save/load parameters round-trips for at least one ATAC workflow.

Execution, if enabled:

- Synthetic index build completes.
- One tiny mapping smoke completes.
- Output collection returns logs and generated files.

Acceptance:

- Docs clearly describe MCP and Launchpad usage from a fresh clone.
- Chromap recipes are usable without STAR-suite.
- Production benchmark workflows are discoverable but clearly marked as long-running.
- Tests protect command rendering and preflight behavior.

## Suggested Public Workflow Catalog

| Workflow ID | Purpose | Stage |
|---|---|---|
| `chromap_index` | Build a Chromap index from FASTA | 1 |
| `chromap_atac_bed` | Paired-end ATAC BED/fragments output | 1 |
| `chromap_atac_bam_fragments` | BAM/CRAM primary output plus fragments/sidecar | 1 |
| `chromap_hic_pairs` | Hi-C pairs output | 1 |
| `chromap_atac_macs3_frag_peaks` | ATAC mapping plus libMACS3 FRAG peaks | 2 |
| `chromap_atac_binary_sidecar` | ATAC BAM + fixed-record sidecar path | 2 |
| `chromap_lib_runner_atac` | libchromap runner debug/parity path | 2 |
| `chromap_bam_sort_index` | BAM sort/index focused workflow | 2 |

## Risks

### STAR Coupling Leakage

Risk: copied code keeps STAR names, paths, or workflow filters.

Mitigation: Stage 1 has an explicit no-STAR-default gate. Run `rg -n
"STAR|star_" mcp_server README.md docs` and allow only historical references
inside migration notes.

### Overfitting To Local Pikachu Paths

Risk: config works only on `/mnt/pikachu`.

Mitigation: keep local paths in `config.yaml`, but document how to override
`repo_root`, trusted roots, and datasets. Prefer relative workflow script paths.

### Execution Safety

Risk: Launchpad or MCP runs arbitrary commands.

Mitigation: preserve allowlisted workflow/script execution only. Keep Stage 1
render-only unless a small synthetic smoke is explicitly allowlisted.

### Drift From Chromap CLI

Risk: hand-written workflow YAMLs go stale as flags change.

Mitigation: keep workflow count small, test rendered commands, and add a later
follow-up to generate schema hints from `chromap -h` if useful.

## Rollback

Stage 1 is isolated to `mcp_server/`, `scripts/launchpad_server.sh`, and docs.
If the MVP is not useful, remove those paths and leave Chromap code untouched.

Stage 2 should be committed separately from Stage 1 so production-hardening
changes can be reverted without removing the basic Launchpad.

## Commit Plan

Stage 1:

1. `Add Chromap MCP server and Launchpad MVP`
2. `Add Chromap workflow recipes`
3. `Document Chromap Launchpad quick start`

Stage 2:

1. `Expand Chromap workflow preflight checks`
2. `Add advanced ATAC and libMACS3 Launchpad recipes`
3. `Add fixture-backed Chromap MCP tests`
4. `Document production Chromap MCP workflows`

