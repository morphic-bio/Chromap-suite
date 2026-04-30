# Agent Notes

## Documentation Map

Start with these local docs before making broad changes:

- `README.md`: user-facing Chromap usage, fork notes, and common commands.
- `docs/index.md`: documentation entry point.
- `docs/chromap.html`: generated/manpage-style CLI option reference.
- `docs/chromap_launchpad.md`: MCP-backed Launchpad usage.
- `mcp_server/README.md`: Chromap MCP server and Launchpad implementation notes.
- `plans/2026-04-30-chromap-mcp-launchpad-port-runbook.md`: MCP/Launchpad port plan.
- `plans/2026-04-30-chromap-mcp-launchpad-hardening-runbook.md`: staged MCP,
  Launchpad, recipe registry, and agent-doc hardening plan.
- `plans/2026-04-30-chromap-core-libchromap-smoke-runbook.md`: core smoke and libchromap regression plan.
- `tests/README.md`: existing Y/noY test documentation.

The MCP recipe registry lives at `mcp_server/recipes/registry.yaml`. Recipe
preflight rules are implemented in `mcp_server/tools/preflight.py` and exposed
through `preflight_recipe`. Recipe run manifests are implemented in
`mcp_server/tools/run_manifest.py`, exposed through `write_recipe_run_manifest`,
and stored under `plans/artifacts/mcp_runs/`. Launchpad recipe forms are served
from registry metadata through `mcp_server/launchpad/api.py`; workflow YAML is
still responsible for argv rendering. When adding or changing a Launchpad
workflow, update the workflow YAML, config entry, recipe registry metadata,
preflight coverage, manifest expectations, and registry/preflight/Launchpad
tests together.

When documenting a new workflow, add the durable reference to `docs/` or
`mcp_server/README.md`, and keep long experimental notes under `plans/`.

## Sister Repositories

STAR Suite is a sister repository used for multiomics integration work:

```text
/mnt/pikachu/STAR-suite
```

Chromap-suite should remain independently buildable and testable, but agents
should know that STAR Suite owns the STAR orchestration side of multiome runs
and may consume `libchromap` or Chromap artifacts. Do not copy STAR-specific
workflow assumptions into Chromap defaults unless a runbook explicitly calls for
that integration.

## MCP Server

Chromap-suite includes an MCP server and browser Launchpad under:

```text
mcp_server/
scripts/launchpad_server.sh
```

Useful commands:

```bash
python3 -m mcp_server.app --help
python3 -m pytest mcp_server/tests -q
bash scripts/launchpad_server.sh up
bash scripts/launchpad_server.sh down
```

The Launchpad is a command builder and local runner for Chromap recipes. Keep
generated Launchpad logs and pid files under `plans/artifacts/`.

Execution safety model:

- Launchpad execution must resolve through an enabled recipe id and linked
  workflow id.
- Commands must be argv arrays, never shell strings.
- Server-side recipe preflight must pass immediately before execution.
- Input and output paths must stay under trusted roots and must not smuggle
  shell syntax through path parameters.
- Long and benchmark recipes require explicit opt-in.
- Local launches write stdout, stderr, and `run.json` under
  `plans/artifacts/mcp_runs/` and are monitored with a timeout.

## Benchmarking Policy

Benchmarks are serial by default. Do not run multiple benchmark or large smoke
jobs concurrently unless a runbook or the user explicitly says that is safe.

Rationale:

- Chromap, STAR, libMACS3, htslib compression, and sort/spill paths compete for
  CPU, RAM, disk bandwidth, and temporary storage.
- Parallel benchmarks make wall time and RSS comparisons hard to interpret.
- Several regression gates compare byte-for-byte outputs and should not share
  mutable output directories.

Smoke tests may be parallelized only after the specific smoke tests are
established as independent, resource-light, and isolated under distinct
artifact directories. Until then, run smokes serially.

For every benchmark or timing run, record:

- command line,
- git commit or dirty-worktree note,
- input fixture and reference/index path,
- thread counts,
- output directory under `plans/artifacts/`,
- wall/user/sys/max RSS when available.

## Git Authorship

Agents must not claim co-authorship in commits, pull requests, or pushes. Do not
add `Co-authored-by` trailers, generated-by trailers, agent signatures, or bot
authorship claims unless the user explicitly requests a specific attribution.

## Artifact Hygiene

Keep generated smoke-test and agent artifacts out of tracked source paths.

Default sandbox:

```text
plans/artifacts/
```

This directory is git-ignored. Use it for:

- generated synthetic FASTA/FASTQ fixtures,
- ENCODE download/downsample caches,
- Chromap indexes created only for smoke tests,
- BAM/CRAM/BED/pairs/fragments outputs,
- Launchpad logs and pid files,
- timing files and test summary TSVs.

Recommended layout:

```text
plans/artifacts/chromap_core_smoke/<timestamp>/
plans/artifacts/encode_fixture_cache/
plans/artifacts/encode_downsample_smoke/<timestamp>/
plans/artifacts/launchpad_server.log
```

Tests and runbooks should accept `CHROMAP_ARTIFACT_ROOT` and default it to
`plans/artifacts`. Do not commit generated FASTQs, indexes, alignment outputs,
or downloaded ENCODE data. Commit only scripts, manifests, checksums, and
documentation needed to reproduce those artifacts.
