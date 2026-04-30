# Chromap Launchpad

Chromap Launchpad is a local browser UI backed by the MCP recipe registry and
workflow renderer.

Start it from the repository root:

```sh
bash scripts/launchpad_server.sh up
```

Then open `http://127.0.0.1:8765/launchpad/`.

The default public recipe list is generated from
`mcp_server/recipes/registry.yaml` and includes:

- `chromap_index`
- `chromap_atac_bed`
- `chromap_atac_bam_fragments`
- `chromap_hic_pairs`

The UI renders typed fields from recipe inputs, shows runtime class, benchmark
policy, expected outputs, and preflight checks, then renders the exact Chromap
argv through the linked workflow. It can launch the command only from a loopback
browser on the same host as the server. Remote browsers stay in planning mode
unless a trusted local path check is available.

The recipe API starts at:

```text
http://127.0.0.1:8765/launchpad/api/recipes
```

Metadata-only planned/long recipes are hidden by default. Launchpad can opt into
displaying them, but they are not executable until the registry marks them
`enabled: true`.

Execution is intentionally narrow. Launchpad resolves a selected recipe to its
allowlisted workflow, builds an argv array rather than a shell string, runs
recipe preflight again on the server, rejects untrusted or shell-like path
parameters, requires explicit confirmation for long or benchmark recipes, and
monitors launched commands with a timeout. Stdout, stderr, and the run manifest
are written together under `plans/artifacts/mcp_runs/`.

Every local launch writes a run manifest before the process starts:

```text
plans/artifacts/mcp_runs/<YYYYMMDDTHHMMSSZ>_<recipe_id>/run.json
```

The same manifest format is available without execution through the MCP
`write_recipe_run_manifest` tool, which records `execution_status: dry_run`.
Manifests capture the recipe id and registry version, rendered argv and shell
preview, git commit and dirty-worktree state, binary metadata, inputs, outputs,
artifact/log paths, host/user context, preflight results, and smoke or benchmark
classification. Benchmark-class records keep the repository policy explicit:
benchmarks are serial-only unless a recipe says otherwise.

Stop the server with:

```sh
bash scripts/launchpad_server.sh down
```
