# Chromap MCP Server and Launchpad

This directory contains the Stage 1 Chromap port of the existing MCP server
and Launchpad. The server exposes a small, schema-driven Chromap command
surface for agents and for the browser UI at `/launchpad/`.

## Start Locally

```sh
bash scripts/launchpad_server.sh up
```

Open:

- Launchpad: `http://127.0.0.1:8765/launchpad/`
- Recipe API: `http://127.0.0.1:8765/launchpad/api/recipes`
- Workflow API: `http://127.0.0.1:8765/launchpad/api/workflows`

Stop it with:

```sh
bash scripts/launchpad_server.sh down
```

## Launchpad Recipes

The default public recipe list is registry-driven and intentionally Chromap-only:

- `chromap_index`
- `chromap_atac_bed`
- `chromap_atac_bam_fragments`
- `chromap_hic_pairs`

Each enabled recipe links to a workflow in `mcp_server/workflows/*.yaml` for
argv rendering. Launchpad forms are generated from
`mcp_server/recipes/registry.yaml`, so the browser shows recipe runtime class,
benchmark policy, typed inputs, expected outputs, registry preflight checks,
rendered argv, and dry-run or launch manifests from one metadata source.

Metadata-only recipes are hidden by default. The browser can opt into showing
planned/long recipes, but schema/render/preflight routes reject them until the
registry marks them `enabled: true`.

## Recipe Registry

Recipe metadata lives in:

```text
mcp_server/recipes/registry.yaml
```

The registry is the shared contract for MCP, Launchpad, agents, tests, and
future run manifests. Workflow YAML controls command rendering; recipe metadata
adds operator intent, expected outputs, preflight rule ids, smoke coverage,
runtime class, benchmark policy, docs, and STAR Suite/downstream handoff
artifacts.

Current handoff metadata covers Chromap-owned ATAC BAM, fragments TSV, optional
fragments binary sidecar, integrated MACS3 peak outputs, ATAC evidence TSV, and
Hi-C `.pairs`. STAR Suite lives at `/mnt/pikachu/STAR-suite` and may consume
these artifacts, but STAR orchestration is not a Chromap Launchpad default.

Launchpad recipe endpoints:

- `GET /launchpad/api/recipes`
- `GET /launchpad/api/recipes/{recipe_id}/schema`
- `GET /launchpad/api/recipes/{recipe_id}/describe`
- `POST /launchpad/api/recipes/{recipe_id}/preflight`
- `POST /launchpad/api/recipes/{recipe_id}/render`
- `POST /launchpad/api/recipes/{recipe_id}/manifest`

Recipe execution safety model:

- executable browser launches resolve through an enabled recipe id and linked
  workflow id,
- commands are constructed as argv arrays, never shell strings,
- recipe preflight runs again on the server immediately before execution,
- input/output paths must stay under configured trusted roots,
- path-like recipe parameters containing shell metacharacters are rejected,
- long and benchmark recipes require explicit opt-in,
- local launches are monitored with a timeout,
- stdout/stderr and `run.json` live in the same manifest artifact directory.

Validate it with:

```sh
python3 -m pytest mcp_server/tests/test_recipe_registry.py -q
```

## MCP Tools

Important tools exposed by `mcp_server.app`:

- `list_recipes`
- `describe_recipe`
- `preflight_recipe`
- `write_recipe_run_manifest`
- `list_workflows`
- `describe_workflow`
- `get_workflow_parameter_schema`
- `validate_workflow_parameters`
- `render_workflow_command`
- `build_chromap`
- `check_build_status`
- `ensure_fresh_build`

Build targets are `chromap`, `libchromap`, and `runner`. Build state is tracked
in `.chromap_build_state.json`.

`preflight_recipe` validates recipe parameters against the registry's preflight
rule ids before command rendering or execution. It returns per-rule
`pass`/`warn`/`fail` results with messages, paths, and suggested fixes.

`write_recipe_run_manifest` writes a dry-run reproducibility record under
`plans/artifacts/mcp_runs/.../run.json`. The manifest captures recipe metadata,
rendered argv, shell preview, git state, binary metadata, inputs, outputs,
preflight results, log paths, host/user context, and benchmark policy before any
real execution.

## Test Tiers

- S0: hermetic synthetic tests, including `python3 -m pytest mcp_server/tests -q`
  and `make test-libchromap-core-smoke`.
- S1: ENCODE downsample smoke tests. Downloads are opt-in only with
  `ENCODE_ALLOW_DOWNLOAD=1`; generated caches stay under `plans/artifacts/`.
- S2: longer integration or benchmark gates, including full-depth multiome and
  MACS3 parity runs. Benchmarks are serial by policy.

## Verification

```sh
python -m pytest mcp_server/tests -q
python -m mcp_server.app --help
```

The copied Launchpad still includes generic Script Lane/editor bridge plumbing
for later stages. Chromap recipe execution remains routed through allowlisted
registry/workflow entries.
