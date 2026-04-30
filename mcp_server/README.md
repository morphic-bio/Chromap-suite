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
- Workflow API: `http://127.0.0.1:8765/launchpad/api/workflows`

Stop it with:

```sh
bash scripts/launchpad_server.sh down
```

## Stage 1 Workflows

The default public workflow list is intentionally Chromap-only:

- `chromap_index`
- `chromap_atac_bed`
- `chromap_atac_bam_fragments`
- `chromap_hic_pairs`

Each workflow lives in `mcp_server/workflows/*.yaml` and renders a single
Chromap command. Browser validation checks schema types by default; local
Launchpad users can enable server-side path checks before launch.

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

Validate it with:

```sh
python3 -m pytest mcp_server/tests/test_recipe_registry.py -q
```

## MCP Tools

Important tools exposed by `mcp_server.app`:

- `list_recipes`
- `describe_recipe`
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

## Verification

```sh
python -m pytest mcp_server/tests -q
python -m mcp_server.app --help
```

The copied Launchpad still includes generic Script Lane/editor bridge plumbing
for later stages. Stage 1 does not port project-specific production recipes.
