# Direct minimizer-index loading

Chromap's historical index is a native binary image containing a small prefix,
the khash flags/keys/values arrays, and the repeated-minimizer occurrence
array. Chromap Suite preserves that file format and its lookup representation.

On platforms that support `O_DIRECT`, mapping loads the file through one
page-aligned 256 MiB staging buffer. Each block is copied into the existing
khash or occurrence destination according to its historical file offset. This
bypasses the filesystem client page cache without rebuilding or rehashing any
entry. The final unaligned file tail is read through the ordinary descriptor.

If direct open, aligned allocation, or direct reading is unsupported, Chromap
automatically rereads the four destination spans with positioned buffered
reads. Startup logs identify the selected path:

```text
Index input path: direct I/O block reader, block size 268435456 bytes.
```

or:

```text
Index input path: buffered positioned-read fallback (...).
```

Neither path checksums index blocks during production loading. Reference and
index artifacts should be digest-validated once when downloaded or staged.
The loader still rejects inconsistent dimensions, offsets, trailing data, and
invalid materialized-reference binding footers using constant-time metadata
checks.

For an index-only timing gate, build `tests/index_load_probe` and pass it one
index path. The probe invokes the production `Index::Load()` implementation but
does not load a FASTA or map reads, keeping unrelated validation work outside
the measured interval.
