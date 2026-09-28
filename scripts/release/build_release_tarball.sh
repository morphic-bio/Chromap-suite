#!/usr/bin/env bash
# Low-level packager for an existing build. For the mandatory test gate use
# scripts/release/build_release.sh, which also produces the Debian package.
set -euo pipefail
exec python3 "$(dirname "$0")/release.py" tarball "$@"
