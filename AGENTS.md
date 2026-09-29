# Agent Notes

## Clean-room rule for 10x Genomics software (read first; settled 2026-09-27)

**We never look at 10x code.** This is settled; do not re-check or re-argue it unless the licences change.

- Never read, fetch, grep, quote, paraphrase or summarise 10x Genomics source code: Cell Ranger, Space Ranger or any other 10x repository or package. This includes code files inside installed tarballs (`lib/python/...`, `.mro`, Rust sources) and commit history.
- Never port, follow or "match" vendor code. The Flex gDNA diagnostic was removed from STAR Suite on 2026-09-27 because it had been written from Cell Ranger source.
- Never cite vendor source files or line numbers in code, docs, notes or papers. Cite published documentation.
- Allowed: published 10x documentation, licence and EULA texts, repository READMEs, and outputs of running the binaries. The Visium HD barcode contract uses reference data shipped with Space Ranger (barcode whitelists, the Visium HD slide-design file, the chemistry definition), disclosed in the Visium paper's Methods. Ask before opening any other file inside a vendor package.
- If a question can only be answered from vendor source, stop and ask.

Licence facts, verified 2026-09-27 against Wayback Machine copies (the live 10x pages return 403/429 to scripted requests):
- The `LICENSE` file in the Cell Ranger and Space Ranger binary packages says use is governed by support.10xgenomics.com/license, the 10x End User Software License. The Space Ranger and Cell Ranger EULA texts are identical.
- EULA §3.4 "No Source Code": "Licensee is granted no rights to, and Licensee agrees that it will not, and will not authorize or permit any third party to decompile, disassemble, or reverse engineer, or attempt to derive the source code for the Licensed Software, in whole or in part." EULA §4.1(c) bars modifying the software and derivative software based on it.
- The source on GitHub (`10XGenomics/spaceranger`, `10XGenomics/cellranger`) has a separate licence: use, copying and modification solely for internal use with 10x products; no redistribution; derivative works bound by the same terms. It has no reverse-engineering clause.
- In writing: describe vendor output neutrally (what it does; never "bug" or imputed intent), and never call the software "closed".

## Documentation Map

Start with these local docs before making broad changes:

- `README.md`: user-facing Chromap Suite usage, project history, and common commands.
- `HISTORY.md`: project lineage (spinoff from `haowenz/chromap` in 2026) and pre-spinoff fork notes.
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

- Chromap, STAR, librapidmacs, htslib compression, and sort/spill paths compete for
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

## Releases

Versioning mirrors STAR Suite: an annotated `vX.Y.Z` git tag is the release, with
a matching `docs/RELEASE_NOTES_vX.Y.Z.md` and a rolled `CHANGELOG.md` entry. The
suite version lives in `src/version.h` (`CHROMAP_SUITE_VERSION`); `chromap
--version` reports it, and `chromap --upstream-version` reports the chromap engine
version (also the BAM/SAM `@PG VN`).

To cut a release:

1. Bump `CHROMAP_SUITE_VERSION` in `src/version.h`.
2. Add `docs/RELEASE_NOTES_vX.Y.Z.md` (see prior releases for the format).
3. Roll `CHANGELOG.md`: move `[Unreleased]` entries under a new `[X.Y.Z] - <date>`.
4. Commit, then `git tag -a vX.Y.Z -m "Chromap Suite vX.Y.Z"` and
   `git push origin master --follow-tags`.

Pushing the tag triggers `.github/workflows/release.yml`. The shared
`scripts/release/build_release.sh` entry point builds, runs the 16 hermetic
release test targets, checks staged binaries and the SDK, and only then creates
tarballs and Debian packages for Ubuntu 22.04 / 24.04. Publication also requires
clean-container tarball and Debian install/runtime/purge checks, a Debian source
rebuild with the same tests, and a successful Docker image smoke test. Manual
workflow runs validate everything without publishing. Docker publication is
enabled only by `RELEASE_PUSH_IMAGE=true`; configured image publication failures
also block the GitHub release.

See [the release procedure](docs/releasing.md) for commands, dependencies,
artifacts, and the validation matrix. `make test-release` runs the shared test
gate. `make test-smoke` additionally needs local 100K biological fixtures and is
not the portable CI gate; retain those fixture-based checks for local release
qualification. The build requires system `libhts-dev` for RapidMACS and linking.

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
