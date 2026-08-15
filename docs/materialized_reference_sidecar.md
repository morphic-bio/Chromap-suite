# Materialized reference sidecar

Chromap normally parses the complete reference FASTA at the start of every
mapping process. The minimizer index contains candidate locations but not the
reference bases, names, or lengths needed for candidate verification and
output.

Chromap Suite can optionally materialize that reference information once while
building the minimizer index and load it as binary data during mapping.

## Build

The feature is opt-in. Supply `--reference-sidecar` to index construction:

```sh
chromap --build-index \
  --ref genome.fa \
  --output genome.index \
  --reference-sidecar genome.chromapref \
  --num-threads 32
```

This produces:

- `genome.index`: the ordinary Chromap minimizer index plus a 40-byte binding
  footer;
- `genome.chromapref`: the materialized reference bases and contig dictionary.

An index built without `--reference-sidecar` retains the historical byte
layout exactly.

## Map

Use the same sidecar with the index that created it:

```sh
chromap --preset atac \
  --index genome.index \
  --reference-sidecar genome.chromapref \
  --read1 R1.fastq.gz --read2 R2.fastq.gz \
  --output fragments.bed \
  --num-threads 32
```

The FASTA is not opened by the mapping engine in this mode. The sidecar loader
uses aligned 256 MiB blocks, one `O_DIRECT` file descriptor, and static
parallel positioned reads into per-thread staging buffers. Each completed
block is copied into the existing contig allocations; the sidecar layout is
unchanged. `--num-threads` caps the number of concurrent block reads. If a
filesystem does not support direct I/O, Chromap falls back to balanced 256 MiB
buffered `pread` tasks.

CRAM output remains the exception: HTSlib needs the original FASTA for CRAM
reference encoding. Supply both options for CRAM; Chromap still uses the binary
sidecar for alignment:

```sh
chromap --CRAM \
  --ref genome.fa \
  --index genome.index \
  --reference-sidecar genome.chromapref \
  --read1 R1.fastq.gz --read2 R2.fastq.gz \
  --output alignments.cram
```

## Contract and validation

Sidecar format version 1 contains:

- a fixed header with magic, version, endian marker, file dimensions, and
  reference fingerprint;
- one directory entry per contig with name/base offsets, length, and CRC32;
- concatenated contig names;
- concatenated uncompressed reference bases.

The index footer records the same reference fingerprint, sequence count, and
total base count. Chromap rejects:

- a sidecar used with an historical/unbound index;
- a sidecar generated for a different reference;
- unsupported versions or byte order;
- malformed offsets, truncation, or trailing data.

Production loading deliberately does not checksum reference blocks. Artifact
integrity is validated once when the reference bundle is downloaded or staged;
workers perform only constant-time header, dimension, offset, and index-binding
checks plus short-read detection. The sequence CRC values used to construct the
binding fingerprint are read from the directory; workers do not recompute them
over the base payload. This keeps the hot path as pure parallel I/O.

Legacy FASTA mapping remains the default and works with both historical and
sidecar-bound indexes. The standalone ATAC spill materializer does not use this
sidecar because its spill headers already contain the contig dictionary and it
does not perform alignment verification.

## Test

The hermetic smoke exercises FASTA/sidecar mapping parity, legacy fallback,
index binding, mismatch rejection, truncated-file rejection, and the direct
parallel loader diagnostic:

```sh
make test-materialized-reference
```

To isolate the exact production sidecar-load path from mapping, build and run
the focused probe:

```sh
make tests/materialized_reference_load_probe
tests/materialized_reference_load_probe genome.chromapref 32
```

Generated test data and outputs remain under
`plans/artifacts/materialized_reference_smoke/`.
