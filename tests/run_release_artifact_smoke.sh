#!/usr/bin/env bash
# Test only the supplied installed/extracted prefix; no source-tree binaries.
set -euo pipefail
prefix="$(realpath "$1")"
version="$2"
revision="$3"
out="$(realpath -m "$4")"
mkdir -p "$out"
for tool in chromap rapidmacs chromap_callpeaks chromap_lib_runner chromap_atac_spill_materializer; do
  test -x "$prefix/bin/$tool"
  ldd "$prefix/bin/$tool" > "$out/$tool.ldd"
  if grep -q 'not found' "$out/$tool.ldd"; then cat "$out/$tool.ldd" >&2; exit 1; fi
done
test "$("$prefix/bin/chromap" --version 2>&1)" = "$version"
python3 - "$prefix" "$version" "$revision" "$out" <<'PY'
import gzip
import hashlib
import json
from pathlib import Path
import random
import struct
import sys
import zlib

prefix, version, revision, out = Path(sys.argv[1]), sys.argv[2], sys.argv[3], Path(sys.argv[4])
meta = json.loads((prefix / 'share/chromap-suite/release.json').read_text())
assert meta['suite_version'] == version
assert meta['source_revision'] == revision
assert len(meta['binaries']) == 5
for name, checksum in meta['binaries'].items():
    assert hashlib.sha256((prefix / 'bin' / name).read_bytes()).hexdigest() == checksum, name
for name in ('LICENSE', 'STAR-input-LICENSE', 'RapidMACS-LICENSE', 'HTSlib-LICENSE'):
    assert (prefix / 'share/licenses/chromap-suite' / name).stat().st_size > 0, name
for path in ('include/chromap-suite/libchromap.h', 'include/chromap-suite/mapping_parameters.h',
             'lib/libchromap.a', 'lib/librapidmacs.a'):
    assert (prefix / path).stat().st_size > 0, path
rng = random.Random(1701)
genome = ''.join(rng.choice('ACGT') for _ in range(12000))
(out / 'ref.fa').write_text('>chrSynthetic\n' + genome + '\n')
rc = str.maketrans('ACGT', 'TGCA')
barcode = 'ACGTACGTACGTACGT'
(out / 'whitelist.txt').write_text(barcode + '\n')
records = [[], [], []]
for i in range(16):
    start = 500 + i * 311
    seqs = [genome[start:start + 90], genome[start + 130:start + 220].translate(rc)[::-1], barcode]
    for mate, seq in enumerate(seqs):
        records[mate].append(f'@pair{i}/{mate + 1}\n{seq}\n+\n' + 'I' * len(seq) + '\n')
for name, rows in zip(('R1', 'R2', 'BC'), records):
    data = ''.join(rows).encode()
    (out / f'{name}.fq').write_bytes(data)
    (out / f'{name}.fq.gz').write_bytes(gzip.compress(data, mtime=0))
    # A small synthetic BGZF member plus the standard empty terminator.
    def member(raw):
        compressor = zlib.compressobj(wbits=-15)
        payload = compressor.compress(raw) + compressor.flush()
        header = bytes.fromhex('1f8b08040000000000ff060042430200')
        return header + struct.pack('<H', 26 + len(payload) - 1) + payload + struct.pack('<II', zlib.crc32(raw), len(raw))
    (out / f'{name}.fq.bgz').write_bytes(member(data) + member(b''))
PY
"$prefix/bin/chromap" --build-index -r "$out/ref.fa" -o "$out/ref.index" -k 11 -w 5 \
  > "$out/index.stdout" 2> "$out/index.stderr"
common=(-x "$out/ref.index" -r "$out/ref.fa" -t 1 --BED)
for format in fq fq.gz fq.bgz; do
  for tool in chromap chromap_lib_runner; do
    "$prefix/bin/$tool" "${common[@]}" -1 "$out/R1.$format" -2 "$out/R2.$format" \
      -o "$out/$tool.$format.bed" > "$out/$tool.$format.stdout" 2> "$out/$tool.$format.stderr"
    LC_ALL=C sort "$out/$tool.$format.bed" > "$out/$tool.$format.sorted"
    test "$(wc -l < "$out/$tool.$format.sorted")" -eq 16
    cmp "$out/chromap.fq.sorted" "$out/$tool.$format.sorted"
  done
done
for tool in chromap chromap_lib_runner; do
  "$prefix/bin/$tool" --preset atac -x "$out/ref.index" -r "$out/ref.fa" -t 1 \
    -1 "$out/R1.fq.bgz" -2 "$out/R2.fq.bgz" -b "$out/BC.fq.bgz" \
    --barcode-whitelist "$out/whitelist.txt" --deterministic-mapping \
    --input-bgzf-mode on --atac-sidecar-only \
    --atac-fragment-binary-output "$out/$tool.fragments.bin" \
    > "$out/$tool.sidecar.stdout" 2> "$out/$tool.sidecar.stderr"
  test -s "$out/$tool.fragments.bin"
done
cmp "$out/chromap.fragments.bin" "$out/chromap_lib_runner.fragments.bin"
cmp "$out/chromap.fragments.bin.chroms.tsv" "$out/chromap_lib_runner.fragments.bin.chroms.tsv"
printf 'result\tPASS\nversion\t%s\nsource_revision\t%s\nalignment_rows\t16\n' "$version" "$revision" > "$out/summary.tsv"
echo "PASS: packaged binaries, dependencies, notices, FASTQ/gzip/BGZF CLI/library parity and sidecars ($out)"
