#!/usr/bin/env bash
# Runs as root only inside a disposable validation container.
set -euo pipefail
tarball="$1"
deb="$2"
version="$3"
revision="$4"
export DEBIAN_FRONTEND=noninteractive
mkdir -p /work/tarball /results
apt-get update
apt-get install -y --no-install-recommends ca-certificates python3
python3 - "$tarball" <<'PY'
from pathlib import PurePosixPath
import sys
import tarfile
with tarfile.open(sys.argv[1]) as archive:
    for entry in archive.getmembers():
        path = PurePosixPath(entry.name)
        if path.is_absolute() or '..' in path.parts or not (entry.isfile() or entry.isdir()):
            raise SystemExit(f'unsafe tarball member: {entry.name}')
    archive.extractall('/work/tarball')
PY
prefixes=(/work/tarball/Chromap-suite-*)
test "${#prefixes[@]}" -eq 1
prefix="${prefixes[0]}"
# Resolve the tarball's exact dependency relationships in this clean OS, using
# an empty temporary package. This tests declared runtime deps without adding
# compilers, HTS development packages or unrelated tools that could mask omissions.
python3 - "$prefix" <<'PY'
import json
from pathlib import Path
import sys
meta = json.loads((Path(sys.argv[1]) / 'share/chromap-suite/release.json').read_text())
control = Path('/work/runtime-deps/DEBIAN/control')
control.parent.mkdir(parents=True)
control.write_text('Package: chromap-release-runtime-deps\nVersion: 1\nArchitecture: all\n'
                  'Maintainer: Chromap Suite\nDescription: temporary release validation dependencies\n'
                  f"Depends: {meta['runtime_dependencies']}\n")
PY
dpkg-deb --build /work/runtime-deps /work/runtime-deps.deb
apt-get install -y --no-install-recommends /work/runtime-deps.deb
bash /checks/tests/run_release_artifact_smoke.sh "$prefix" "$version" "$revision" /results/tarball
apt-get purge -y chromap-release-runtime-deps
apt-get autoremove -y --purge

apt-get install -y --no-install-recommends "$deb"
for tool in chromap rapidmacs chromap_callpeaks chromap_lib_runner chromap_atac_spill_materializer; do
  test "$(command -v "$tool")" = "/usr/bin/$tool"
done
bash /checks/tests/run_release_artifact_smoke.sh /usr/lib/chromap-suite "$version" "$revision" /results/deb
apt-get purge -y chromap-suite
apt-get autoremove -y --purge
hash -r
test ! -e /usr/lib/chromap-suite
for tool in chromap rapidmacs chromap_callpeaks chromap_lib_runner chromap_atac_spill_materializer; do
  test ! -e "/usr/bin/$tool"
  test ! -L "/usr/bin/$tool"
done
if dpkg-query -W -f='${Status}' chromap-suite 2>/dev/null | grep -q 'install ok installed'; then exit 1; fi
printf 'tarball\tPASS\ndeb-install\tPASS\ndeb-runtime\tPASS\ndeb-purge\tPASS\n' > /results/checks.tsv
echo 'PASS: clean-container tarball runtime and Debian install/runtime/purge'
