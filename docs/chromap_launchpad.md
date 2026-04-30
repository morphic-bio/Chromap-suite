# Chromap Launchpad

Chromap Launchpad is a local browser UI backed by the same workflow schemas
served through the MCP server.

Start it from the repository root:

```sh
bash scripts/launchpad_server.sh up
```

Then open `http://127.0.0.1:8765/launchpad/`.

Stage 1 includes four public Chromap CLI recipes:

- `chromap_index`
- `chromap_atac_bed`
- `chromap_atac_bam_fragments`
- `chromap_hic_pairs`

The UI validates parameters, renders the exact Chromap argv, and can launch the
command only from a loopback browser on the same host as the server. Remote
browsers stay in planning mode unless a trusted local path check is available.

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
