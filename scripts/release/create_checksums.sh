#!/usr/bin/env bash
# Write a SHA256SUMS manifest over every file in the release directory.
set -euo pipefail
DIR="${1:-dist/release}"
cd "$DIR"
find . -type f ! -name 'SHA256SUMS' -print0 | sort -z \
  | xargs -0 sha256sum | sed 's#\./##' > SHA256SUMS
echo "wrote ${DIR}/SHA256SUMS"
cat SHA256SUMS
