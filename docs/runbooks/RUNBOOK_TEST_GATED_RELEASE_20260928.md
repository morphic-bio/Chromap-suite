# Test-gated release pipeline: 2026-09-28 validation

## Implementation

This implements the requested STAR Suite-style release automation for Chromap
Suite. The reference was STAR Suite `482454873cc0322737d2b2d1745fa9736f5f8cc5`,
particularly its release workflow and Debian build/install checks. STAR source
was not modified. The durable operating procedure is [releasing.md](../releasing.md).

The shared build command runs the portable test gate and staged runtime/SDK
checks before creating tarballs and binary Debian packages. CI additionally
requires the clean-runtime matrix, a source-bundle rebuild with the same tests,
and a tested Docker image. Configured Docker publication must succeed before
the GitHub release can publish. Manual workflow runs do not publish.

## Validated source and environment

- Chromap implementation: `bccfd55d472754cb514b31c3d92cc4d58190fa37` on
  `feat/test-gated-release-20260928`, based on integration `0a4c71c`.
- RapidMACS: `34df44818853a1beb59e98b0848813e38a4a8c51`.
- Worktree: `/mnt/pikachu/Chromap-suite-release-pipeline-20260928`.
- Artifact root: `plans/artifacts/release-pipeline-20260928/` in that worktree.
- Package builds used committed source snapshots, Ubuntu 22.04/24.04 Docker
  builders, eight compilation jobs and the shared serial test runner. The
  Debian source rebuild used four compilation jobs. Runtime containers were
  capped at four CPUs. Untimed builds/tests ran at nice level 10; no biological
  throughput benchmark was run for this packaging change.
- The final documentation commit updates installation instructions and records
  these results; it does not alter the validated implementation. Tagged CI will
  rebuild artifacts with that tag's exact source revision and documentation.

## Results

| Gate | Result | Evidence under artifact root |
| --- | --- | --- |
| Failure and publication regression tests | 13 passed | `pipeline-tests.log` |
| GitHub Actions syntax/dependencies | actionlint v1.7.7 passed | `actionlint.log` |
| Ubuntu 22.04 build, tests, stage, SDK and packages | 14 targets passed; packaging passed | `ubuntu22-attempt3.log`, `ubuntu22-attempt3/tests/`, `ubuntu22-attempt3/packages/` |
| Ubuntu 24.04 build, tests, stage, SDK and packages | 14 targets passed; packaging passed | `ubuntu24-attempt1.log`, `ubuntu24-attempt1/tests/`, `ubuntu24-attempt1/packages/` |
| Ubuntu 22.04 artifacts on Ubuntu 22.04 | passed | `runtime-22-on22/checks.tsv` |
| Ubuntu 22.04 artifacts on Ubuntu 24.04 | passed | `runtime-22-on24/checks.tsv` |
| Ubuntu 24.04 artifacts on Ubuntu 24.04 | passed | `runtime-24-on24/checks.tsv` |
| Published-format source bundle extraction and Debian rebuild | 14 targets and SDK passed | `source-attempt3/`, `source-rebuild.log`, `source-rebuild/source-check/plans/artifacts/debian-release/` |
| Source-rebuilt Debian package on Ubuntu 24.04 | passed | `source-runtime/checks.tsv` |
| Docker build, installed SDK and image smoke | passed; 16 identical CLI/library alignment rows | `docker-build.log`, `docker-smoke/summary.tsv` |
| Artifact assembly and SHA256 verification | passed | `validation.json`, `validation-summary.log`, `validated-assets/SHA256SUMS` |

Each runtime row checks tarball execution, Debian installation, Debian
execution and purge. It verifies five executable dependencies, version, source
revision and binary hashes; notices and SDK files; plain FASTQ/gzip/BGZF mapping
with 16 expected fragment rows; and CLI/library sidecar parity. These containers
have the declared runtime dependencies and Python, without a compiler or source
checkout. The source package includes the pinned RapidMACS source and rebuilds
without a Git worktree.

The test runner rejects a skipped required test. The existing optional external
BQTools path in input-format tests can remain skipped; native CBQ tests ran and
passed. Larger fixture-based biological qualification remains the previously
recorded v1.1.0 integration validation, documented in the
[reader handoff](../handoffs/HANDOFF_STAR_MIRRORED_FASTQ_READER_20260928.md).

The first clean-build attempt exposed a missing GNU `time` dependency; a later
staged check exposed a CLI-only option in the new library smoke invocation.
Both were fixed and the affected pipeline rerun successfully. Neither failed
attempt produced release packages. Regression tests also inject compilation,
test, skipped-test and SDK failures and verify that gates stop the pipeline.

Local candidate image:

```text
local/chromap-suite:release-pipeline-check
sha256:a722bb96078401f83f862569117c376033a18378c075469a773267d383d18953
```

## Reproduction and handoff

The exact successful local commands are preserved in
`plans/artifacts/release-pipeline-20260928/run_remaining.sh`; the Ubuntu 22.04
build was:

```bash
nice -n 10 bash scripts/release/run_build_container.sh v1.1.0 ubuntu:22.04 \
  plans/artifacts/release-pipeline-20260928/ubuntu22-attempt3 8
```

Use fresh output directories when repeating a build. The source rebuild
extracts `Chromap-suite-v1.1.0-debian-source.tar.gz`, uses `dpkg-source -x`, and
runs `dpkg-buildpackage -b -us -uc -j4` inside the Ubuntu 24.04 builder. Source
packages are unsigned. PPA publication and additional architectures are outside
the current workflow.

The local candidate assets are in `validated-assets/`: both baseline tarballs
and `.deb` files, their build manifests, the source bundle and SHA256SUMS.
`validation.json` records their hashes and all gate results. CI regenerates
these from the final release tag, retaining logs when a gate fails.

The release branch is integrated locally into `master` and
`integration/chromap-v1.1.0-20260928`, and the unpublished annotated `v1.1.0` tag
is advanced to include the pipeline and documentation. Existing reader commits
and the primary checkout's unrelated files are preserved. No Git ref, release
asset or image has been pushed. GitHub execution itself remains pending the
first manual run or tag push; its workflow has been checked locally with
actionlint and the dependency regression tests.
