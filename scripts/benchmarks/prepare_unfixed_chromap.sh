#!/usr/bin/env bash
# Materialize and build the unfixed upstream Chromap baseline used by the
# publication benchmark.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

SOURCE_REPO="${SOURCE_REPO:-/mnt/pikachu/chromap}"
UPSTREAM_REF="${UPSTREAM_REF:-upstream/master}"
ARTIFACT_ROOT="${ARTIFACT_ROOT:-${REPO_ROOT}/plans/artifacts/upstream_chromap_unfixed}"

if [[ ! -d "${SOURCE_REPO}/.git" ]]; then
  echo "ERROR: SOURCE_REPO is not a git repo: ${SOURCE_REPO}" >&2
  exit 2
fi

commit="$(git -C "${SOURCE_REPO}" rev-parse "${UPSTREAM_REF}")"
worktree="${ARTIFACT_ROOT}/${commit}"

mkdir -p "${ARTIFACT_ROOT}"
if [[ ! -d "${worktree}/.git" ]]; then
  git -C "${SOURCE_REPO}" worktree add --detach "${worktree}" "${commit}" >&2
fi

make -C "${worktree}" >&2

{
  printf 'source_repo\t%s\n' "${SOURCE_REPO}"
  printf 'upstream_ref\t%s\n' "${UPSTREAM_REF}"
  printf 'commit\t%s\n' "${commit}"
  printf 'worktree\t%s\n' "${worktree}"
  printf 'chromap\t%s/chromap\n' "${worktree}"
  printf 'version\t'
  "${worktree}/chromap" --version 2>&1 || true
  printf 'sha256\t%s\n' "$(sha256sum "${worktree}/chromap" | awk '{print $1}')"
} >"${worktree}/build_manifest.tsv"

printf '%s/chromap\n' "${worktree}"
