# Mirrored STAR Suite FASTQ input reader

This directory is a copy of STAR Suite's in-tree BGZF FASTQ reader, so that
Chromap Suite reads FASTQ the way STAR Suite does. It is a stopgap: a shared,
open FASTQ library will replace it. Until then, change these files only by
re-mirroring them from STAR Suite.

## Source

- Repository: STAR Suite, <https://github.com/morphic-bio/STAR-suite>
- Commit: `6e8385375942d59b25b90b3aa9db0a87b9fefc7a`
  (branch `feat/atac-sidecar-without-bam`), directory
  `core/legacy/source/input/`.
- The four files are byte-identical at STAR Suite `master` (`4824548`) and at
  the public release tag `v1.9.5.a` (`f0d9f27de4105bce5e1513cca53f36d5d797ea95`).
- Mirrored on 2026-09-28.

| Mirrored file | STAR Suite source | sha256 of the source file | sha256 of the mirrored file |
|---|---|---|---|
| `BgzfBlockReader.h` | `core/legacy/source/input/BgzfBlockReader.h` | `fcc8cdf08da9bb9b7ad1dc6368de501110aa6ab1d221e97a294bbe63398f5fd1` | `e943ac972838b0b4db6a62fe7329bbbea183faa94153dc29eab0ce1fa5186edd` |
| `BgzfBlockReader.cpp` | `core/legacy/source/input/BgzfBlockReader.cpp` | `1755f003ebe9de491c0a5e9fd267452df73f3840c8d00b78905e5ece16620c6e` | `ede519a7b939ecb380fd01839c3ddef2e29e8125248b30279b629109fb65b6cd` |
| `BgzfRangeReader.h` | `core/legacy/source/input/BgzfRangeReader.h` | `3b8d306001c223d64c5bc2eddcaf73fba51379e49c2fec3dabb33b6afe6e6e19` | `4731a151508fe4d4e7bba6dd63a466962acca8cd6ba24ff1a547353abbf1809d` |
| `BgzfRangeReader.cpp` | `core/legacy/source/input/BgzfRangeReader.cpp` | `d140bc6b50ad7940f220447e950ce96a768eed4b601b775382aa1aaa302f0718` | `429fa8984fc2aa89a08ef8c9bf14b5870bf0592d2f89ab903e0e97a338da52bd` |

`LICENSE` reproduces STAR Suite's `LICENSE` and `core/legacy/LICENSE`
verbatim; both are MIT. The source files carry no per-file notices.

## Local edits

Only namespaces, include paths and include guards changed. Each edit, in every
file where it applies:

1. Namespaces: `namespace star {` -> `namespace chromap {`,
   `namespace input {` -> `namespace star_input {`, and the matching closing
   comments. STAR Suite links `libchromap.a` into its own binary, which also
   contains `star::input`; the new namespace keeps the two copies apart.
2. Include paths: `#include "input/BgzfBlockReader.h"` and
   `#include "input/BgzfRangeReader.h"` -> `"star_input/..."`
   (`BgzfBlockReader.cpp`, `BgzfRangeReader.h`, `BgzfRangeReader.cpp`).
3. Include guards: `CODE_input_BgzfBlockReader` ->
   `CHROMAP_STAR_INPUT_BgzfBlockReader` and `CODE_input_BgzfRangeReader` ->
   `CHROMAP_STAR_INPUT_BgzfRangeReader`, so a translation unit that sees both
   copies cannot silently skip one.

`diff` against the source shows nothing else. To check:

```bash
git -C STAR-suite show 6e83853:core/legacy/source/input/BgzfRangeReader.cpp |
  diff - src/star_input/BgzfRangeReader.cpp
```

## Consulted, not copied

- `BgzfPipeGroup.h` feeds ordered inflated bytes through FIFOs into STAR's own
  chunk parser; Chromap parses records directly and has no use for it.
- `FastxInputModule.*` and `InputContract.h` are STAR's generic line reader
  and record contract. They depend on STAR's `IncludeDefine.h`, and Chromap's
  zlib path keeps its own reader (kseq).
- `BgzfStarAdapter.*` (two ordered streams) and
  `core/features/process_features/src/pf_bgzf_input.cpp` (two or three
  ordered streams) define how STAR Suite keeps mates synchronized: each file has
  its own reader and inflate workers; records are taken from the files in
  lockstep; the records of one pair must have the same ordinal and the same
  read-name stem, a trailing `/1`, `/2` or `/3` removed; a stream that ends
  early is an error; and inflate workers are split evenly across the files.
  `src/fastq_bgzf_input.cc` implements the same rules for Chromap rather than
  copying these files, which depend on STAR types. It consumes each file on
  its own thread, as STAR's core path consumes each mate through its own
  `BgzfPipeGroup` producer, and checks the pairing rules for every record of
  each batch.

## Where Chromap uses it

Only `src/fastq_bgzf_input.{h,cc}` includes these headers. Swapping in the
shared library therefore means replacing this directory, the
`star_input_cpp_source` line in the Makefile, and the implementation in
`src/fastq_bgzf_input.cc`; the rest of Chromap sees only
`BgzfPairedEndReadProvider`, `BgzfFastqStream`, `SelectBgzfFastqInput` and
`BgzfInflateWorkersPerStream`.
