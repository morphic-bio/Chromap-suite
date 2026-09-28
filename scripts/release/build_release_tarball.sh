#!/usr/bin/env bash
# Compatibility entry point: release packaging always runs the shared test gate.
set -euo pipefail
exec bash "$(dirname "$0")/build_release.sh" "$@"
