#!/usr/bin/env bash
# Package a Chromap Suite release tarball from already-built binaries.
#
# Assumes `make chromap chromap_callpeaks chromap_lib_runner libchromap.a` has
# already run (the release workflow builds first). Produces
#   <out-dir>/Chromap-suite-<version>-linux-amd64-<glibc-label>.tar.gz
# with the binaries, libchromap.a + header, docs, and a VERSION manifest.
set -euo pipefail

VERSION=""
OUT_DIR="dist/release"
LABEL=""
while [ $# -gt 0 ]; do
  case "$1" in
    --version) VERSION="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --label)   LABEL="$2";   shift 2 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
done
[ -n "$VERSION" ] || { echo "--version is required" >&2; exit 2; }

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$REPO_ROOT"

# glibc baseline label (e.g. glibc235) from the build host, unless overridden.
if [ -z "$LABEL" ]; then
  glibc="$(ldd --version 2>/dev/null | head -1 | grep -oE '[0-9]+\.[0-9]+' | head -1 || true)"
  LABEL="glibc${glibc//./}"
  [ "$LABEL" = "glibc" ] && LABEL="glibc-unknown"
fi

NAME="Chromap-suite-${VERSION}-linux-amd64-${LABEL}"
WORK="$(mktemp -d)"
STAGE="${WORK}/${NAME}"
trap 'rm -rf "${WORK}"' EXIT
mkdir -p "${STAGE}/bin" "${STAGE}/lib" "${STAGE}/include" "${STAGE}/docs"

for b in chromap chromap_callpeaks chromap_lib_runner; do
  [ -x "$b" ] || { echo "missing built binary: $b (run make first)" >&2; exit 1; }
  cp "$b" "${STAGE}/bin/"
done
[ -f libchromap.a ]   && cp libchromap.a   "${STAGE}/lib/"
[ -f src/libchromap.h ] && cp src/libchromap.h "${STAGE}/include/"
cp README.md LICENSE CHANGELOG.md "${STAGE}/"
[ -f "docs/RELEASE_NOTES_${VERSION}.md" ] && cp "docs/RELEASE_NOTES_${VERSION}.md" "${STAGE}/docs/"

engine="$(./chromap --upstream-version 2>&1 || echo unknown)"   # printed to stderr
cat > "${STAGE}/VERSION" <<EOF
chromap-suite ${VERSION}
chromap-engine ${engine}
glibc-baseline ${LABEL}
platform linux-amd64
EOF

mkdir -p "${OUT_DIR}"
tar -C "${WORK}" -czf "${OUT_DIR}/${NAME}.tar.gz" "${NAME}"
echo "wrote ${OUT_DIR}/${NAME}.tar.gz"
