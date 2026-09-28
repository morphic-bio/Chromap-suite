#!/usr/bin/env bash
# Local reproduction of one CI build/test/package baseline from committed source.
set -euo pipefail
repo_root="$(cd "$(dirname "$0")/../.." && pwd)"
if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "usage: $0 VERSION ubuntu:22.04|ubuntu:24.04 OUTPUT_DIR [JOBS]" >&2
  exit 2
fi
version="$1"
base="$2"
out="$(realpath -m "$3")"
jobs="${4:-4}"
[[ "$base" == ubuntu:22.04 || "$base" == ubuntu:24.04 ]]
[[ "$jobs" =~ ^[1-9][0-9]*$ ]]
mkdir -p "$out"
snapshot="$(mktemp -d "$out/snapshot.XXXXXX")"
trap 'rm -rf "$snapshot"' EXIT
python3 "$repo_root/scripts/release/release.py" snapshot --version "$version" --out-dir "$snapshot/source"
image="local/chromap-release-build:${base#ubuntu:}"
docker buildx build --builder default --load --build-arg "BASE_IMAGE=$base" \
  -f "$repo_root/scripts/release/docker/Dockerfile.build" -t "$image" "$repo_root"
docker run --rm --cpus "$jobs" \
  --mount "type=bind,src=$snapshot/source,dst=/source,readonly" \
  --mount "type=bind,src=$out,dst=/output" \
  "$image" nice -n 10 bash -c '
    set -euo pipefail
    cp -a /source/. /build/
    bash scripts/release/build_release.sh --version "$1" --out-dir /output/packages --test-dir /output/tests --jobs "$2"
  ' bash "$version" "$jobs"
