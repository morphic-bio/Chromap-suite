# Building and releasing Chromap Suite

The release pipeline follows STAR Suite's build, package, clean-install and
runtime validation sequence. All publishing jobs depend on successful tests.
The implementation is independent of STAR Suite and requires no STAR checkout.

## Required gates

1. Check that the requested tag matches `CHROMAP_SUITE_VERSION` and has release
   notes. Run regression tests for failure propagation and workflow dependencies.
2. Build on Ubuntu 22.04 and 24.04 (amd64). On each baseline, run the 15 targets in
   `scripts/release/run_release_tests.sh` serially, then check staged binaries and
   compile/link/run a consumer against the packaged SDK. Only after these pass
   does `scripts/release/build_release.sh` create the tarball and binary `.deb`.
3. Extract the tarball and install the `.deb` in clean Ubuntu containers. Test
   mapping, CLI/library parity, plain FASTQ, gzip, BGZF and sidecar output; check
   version, source revision, binary hashes, licenses and dynamic dependencies.
   Purge the Debian package and verify that commands and package files disappear.
4. Export the committed source and pinned RapidMACS submodule into a Debian
   source package. Rebuild it on Ubuntu 24.04 with the same test gate and SDK
   check, then test its installed `.deb` in a clean Ubuntu 24.04 container.
5. Build and smoke-test the Docker image. If image publication is enabled, push
   that exact tested image. Publish the GitHub release only after all these jobs
   succeed, including any configured image push.

| Built on | Tested on |
| --- | --- |
| Ubuntu 22.04 / glibc 2.35 | Ubuntu 22.04 and 24.04 |
| Ubuntu 24.04 / glibc 2.39 | Ubuntu 24.04 |

The test gate covers unit tests, barcode sampling, reference sidecars, fragment
storage and peak calling, spill records/materialization, FASTQ intake, input
formats, ATAC sidecars, libchromap and CBQ. Python 3, samtools, bgzip and GNU time are
required; test inputs are generated locally. Optional external BQTools paths can
be skipped by their existing tests. The local `make test-smoke` additionally
uses out-of-tree 100K fixtures for peak and low-memory qualification; those
checks remain part of local release qualification, rather than a CI dependency
on a workstation's datasets. S1/S2 biological validation remains opt-in.

## Reproduce the builds locally

Commit the release changes first: the container helper and source exporter use
`git archive HEAD`, including the recorded RapidMACS revision. Initialize that
submodule with `git submodule update --init third_party/rapidmacs`. Docker and
Python 3 are the host requirements. Run heavyweight builds/tests serially.

```bash
bash scripts/release/run_build_container.sh v1.1.0 ubuntu:22.04 \
  plans/artifacts/release-check/ubuntu22 4
bash scripts/release/run_build_container.sh v1.1.0 ubuntu:24.04 \
  plans/artifacts/release-check/ubuntu24 4
```

Each output directory contains `packages/` and `tests/`; use a fresh output
directory for each attempt. The builder image installs its dependencies from
`scripts/release/docker/Dockerfile.build`. To build directly on a supported
Ubuntu host with those dependencies installed:

```bash
bash scripts/release/build_release.sh --version v1.1.0 --jobs 4 \
  --out-dir dist/release --test-dir plans/artifacts/release
```

This command tests the working tree. Use a clean committed checkout for release
provenance. A compilation, test or SDK failure exits nonzero before packaging;
CI never uploads packages from a failed build. Output directories must be empty
to prevent old artifacts from being mistaken for new results. The legacy
`build_release_tarball.sh` entry point now forwards to this same tested pipeline
and produces both formats; it cannot bypass the test gate.

For each row of the runtime matrix, run:

```bash
bash scripts/release/check_artifacts_in_container.sh \
  path/to/Chromap-suite-v1.1.0-linux-amd64-glibc235.tar.gz \
  path/to/chromap-suite_1.1.0-1.ubuntu22.04.1_amd64.deb \
  ubuntu:24.04 1.1.0 "$(git rev-parse HEAD)" \
  plans/artifacts/release-check/ubuntu22-on24
```

This uses runtime packages declared by `dpkg-shlibdeps`, without a compiler,
development headers or the source checkout inside the validation container.
Tarballs contain shared-library requirements in
`share/chromap-suite/release.json`; they are not fully static executables.

## Artifacts and installation

Each baseline produces a tarball, binary `.deb`, build manifest and SHA256SUMS.
Both formats include all five executables, static libchromap/librapidmacs
archives, public headers, documentation and dependency license notices.

Extract a tarball and add its `bin/` to `PATH`, installing its declared runtime
dependencies first. Install a Debian package with `sudo apt install ./<file>.deb`;
executables appear under `/usr/bin`, with the SDK under `/usr/lib/chromap-suite`.
The package's suggested development dependencies are needed for SDK consumers.

The release also includes `Chromap-suite-vX.Y.Z-debian-source.tar.gz`, containing
the `.dsc`, `.orig.tar.gz` and `.debian.tar.xz` with their original names. This
bundle preserves Debian prerelease filenames when downloaded from GitHub.
Extract the bundle, use `dpkg-source -x <file>.dsc`, install the build dependencies
from `debian/control`, and run `dpkg-buildpackage -b -us -uc -j4` in the extracted
source. Source packages are unsigned; PPA upload and signing are outside this
workflow. The current package matrix supports amd64 only.

## CI and publication

Use the **Release** workflow's manual action with a version (for example,
`v1.1.0`) to build and validate a committed branch without publishing. Test logs
are retained as workflow artifacts on failure. The source-export validation
requires a complete checkout with the pinned submodule.

For a release, update `src/version.h`, `CHANGELOG.md` and
`docs/RELEASE_NOTES_vX.Y.Z.md`, commit and create an annotated `vX.Y.Z` tag. Pushing
the tag starts the same validation and then publishes the release assets and
notes. A failed dependency skips publication; there is no `always()` bypass.
Never move a tag that has already been published.

Set repository variable `RELEASE_PUSH_IMAGE=true` to enable Docker Hub
publication, with secrets `DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN`. Optional
`DOCKER_IMAGE_REPO` defaults to `biodepot/chromap-suite`. Image builds and smoke
checks run even when publishing is disabled. Stable releases starting at v1
also update `latest`; prereleases and v0 tags do not. Manual runs never push
images or publish GitHub releases.
