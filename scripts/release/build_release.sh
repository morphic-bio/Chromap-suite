#!/usr/bin/env bash
# Canonical build -> tests -> staging/SDK checks -> tarball + .deb entrypoint.
set -euo pipefail
exec python3 "$(dirname "$0")/release.py" build "$@"
