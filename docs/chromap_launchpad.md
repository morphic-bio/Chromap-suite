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

Stop the server with:

```sh
bash scripts/launchpad_server.sh down
```
