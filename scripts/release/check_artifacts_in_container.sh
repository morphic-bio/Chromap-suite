#!/usr/bin/env bash
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
if [[ $# -ne 6 ]]; then
  echo "usage: $0 TARBALL DEB ubuntu:22.04|ubuntu:24.04 VERSION REVISION OUTPUT_DIR" >&2
  exit 2
fi
tarball="$(realpath "$1")"
deb="$(realpath "$2")"
base="$3"
version="$4"
revision="$5"
out="$(realpath -m "$6")"
[[ "$base" == ubuntu:22.04 || "$base" == ubuntu:24.04 ]]
mkdir -p "$out"
docker run --rm --cpus "${RELEASE_CHECK_CPUS:-4}" \
  --mount "type=bind,src=$tarball,dst=/artifacts/release.tar.gz,readonly" \
  --mount "type=bind,src=$deb,dst=/artifacts/release.deb,readonly" \
  --mount "type=bind,src=$repo_root/tests,dst=/checks/tests,readonly" \
  --mount "type=bind,src=$repo_root/scripts/release/docker,dst=/checks/docker,readonly" \
  --mount "type=bind,src=$out,dst=/results" \
  "$base" nice -n 10 bash /checks/docker/check_artifacts.sh \
  /artifacts/release.tar.gz /artifacts/release.deb "$version" "$revision"
