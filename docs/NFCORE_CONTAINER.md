# nf-core Container Package

Chromap Suite `v1.0.1` is the source basis for the nf-core package. The build is
pinned to:

- Chromap Suite commit `98a4da086f81b7cb159d8fe44efff2fb168e0785`.
- RapidMACS `v1.0.1` commit `c22ff439cc66337e1b7fc9b87b547478eb5c0d3b`.
- Ubuntu 22.04 linux/amd64 manifest
  `sha256:79676deb51ebb02885b0b9d33788e78a37cf1045ad79d1bb04c6a222c3556b3d`.

The runtime image contains these entrypoints:

- `chromap`
- `rapidmacs`
- `chromap_callpeaks` (compatibility name)
- `chromap_lib_runner`
- `chromap_atac_spill_materializer`

`libchromap.a` and `librapidmacs.a` are retained under
`/opt/chromap-suite/lib`, with their Chromap, RapidMACS, and htslib headers under
`/opt/chromap-suite/include`. A downstream multi-stage container build can copy
that prefix as an exact SDK/runtime composition boundary. STAR Suite executables
that use the direct libchromap contract remain statically linked; they do not
require these archives at runtime.

## Build and test

```bash
git submodule update --init --recursive
docker build --pull=false -t local/chromap-suite:v1.0.1 .
tests/run_container_smoke.sh local/chromap-suite:v1.0.1
```

The container smoke verifies OCI provenance, checks every executable for
missing shared libraries, builds a synthetic index, and requires exact BED-row
parity between the container's Chromap CLI and `libchromap` runner.

Do not publish an nf-core module with only a mutable image tag. After the image
is pushed to its approved registry, record the registry manifest digest and use
that digest in the submission evidence. The candidate human-readable tag is
`v1.0.1`; the digest is the immutable identity.

## STAR Suite 1.7.1 composition gate

The exact STAR Suite `v1.7.1` source (`b523c1f58c7f99eb7bc3d3f1b418ac4ab59112a4`)
was clean-built against Chromap Suite `v1.0.1`. Validation on 2026-08-27
established:

- Chromap's hermetic CLI/libchromap parity suite passed for bulk, ChIP, ATAC,
  low-memory ATAC, Hi-C, single-cell barcode, sorted BAM, dual BAM/fragments,
  and Y/noY cases.
- STAR's libchromap contract tools and full `WITH_CHROMAP=1` core linked.
- STAR's synthetic FASTQ-versus-native-CBQ contract smoke passed.
- The 100K concurrent GEX+ATAC multiome run completed with 320,017 retained
  ATAC mappings; GEX and ATAC BAMs passed `samtools quickcheck`.
- The in-process and sidecar callers produced the same 2,047 peak and summit
  records. Coordinates and all numeric fields were identical, with normalized
  peak SHA-256
  `56af66d75d4d53934ea00f7e7428f6631238101c7ac43bd41c35992c837c3e08`
  and summit SHA-256
  `42843cc33d3d5d8bfa36e221fd127651882849c52a729864d4a794dec038f43f`.

The existing STAR 1.7.1 byte-comparison script reports a naming-only mismatch:
the legacy in-process path uses MACS3's default `NA_peak_N`, while the production
sidecar path uses `peak_N`. STAR 1.7.1 documentation designates the sidecar
peak/MEX path as the production boundary, so this does not alter mapping, peak
coordinates, peak scores, or MEX construction. A future STAR maintenance
change should give both paths the same explicit name prefix and update the
byte-identity gate.
