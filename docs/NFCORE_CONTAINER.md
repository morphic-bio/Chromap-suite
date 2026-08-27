# nf-core Container Package

Chromap Suite `v1.0.1` is the source basis for the nf-core package. The build is
pinned to:

- Chromap Suite commit `98a4da086f81b7cb159d8fe44efff2fb168e0785`.
- RapidMACS `v1.0.1` commit `c22ff439cc66337e1b7fc9b87b547478eb5c0d3b`.
- Ubuntu 22.04 multi-architecture index
  `sha256:2edbbc5dc405e9612ba3584ce95480277e3eb374407b5505fe26f17df77c7dbc`.

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
docker buildx build --platform linux/amd64 --load \
  -t local/chromap-suite:v1.0.1 .
tests/run_container_smoke.sh local/chromap-suite:v1.0.1
docker buildx build --platform linux/arm64 --load \
  -t local/chromap-suite:v1.0.1-arm64 .
tests/run_container_smoke.sh local/chromap-suite:v1.0.1-arm64
```

The container smoke verifies OCI provenance, checks every executable for
missing shared libraries, builds a synthetic index, and requires exact BED-row
parity between the container's Chromap CLI and `libchromap` runner.

The image supports linux/amd64 and linux/arm64. Chromap 0.3.3-r519 includes x86
SSE2/SSE4.1 intrinsics in `ksw.cc` and `alignment.cc`. The amd64 build keeps
native `-msse4.1`; the arm64 build places two packaging-only compatibility
headers earlier on the include path and resolves the same intrinsic API through
Ubuntu's MIT-licensed SIMDe 0.7.2 headers. SIMDe lowers supported operations to
NEON without changing the Chromap source or scalar algorithm. RapidMACS already
selects architecture flags independently and contains no SSE intrinsics.

Both descriptors must pass the container smoke before publication. The smoke
requires CLI/libchromap output parity, executable linkage, index construction
and a synthetic alignment. The installed SDK also compiles and runs a downstream
`libchromap` consumer on each target architecture. STAR Suite 1.7.1's full 100K
composition gate was run on linux/amd64; the arm64 SDK and mapping gates cover
the binary composition boundary pending a future full arm64 100K run.

Pre-publication evidence from 2026-08-27:

| Architecture | Container smoke | Alignment rows | Sorted alignment SHA-256 |
|---|---:|---:|---|
| linux/amd64 | PASS | 16 | `e9f1721d56f6b0a166a1ea4f358a959cbfcbb0e0844fe1f1f2edfbfa8f43dbed` |
| linux/arm64 | PASS | 16 | `e9f1721d56f6b0a166a1ea4f358a959cbfcbb0e0844fe1f1f2edfbfa8f43dbed` |

Within each architecture, the CLI and direct `libchromap` runner outputs are
byte-identical after sorting. They are also byte-identical across architectures.
The local evidence summaries are under `plans/artifacts/container_smoke/`; that
directory is intentionally ignored. Re-run the smoke from the committed source
before publishing the registry manifest.

Do not publish an nf-core module with only a mutable image tag. After the image
is pushed to its approved registry, record the registry manifest digest and use
that digest in the submission evidence. The candidate human-readable tag is
`v1.0.1`; the digest is the immutable identity.

The owner-only publication step is:

```bash
docker buildx build --platform linux/amd64,linux/arm64 \
  --tag docker.io/biodepot/chromap-suite:v1.0.1 \
  --push .
docker buildx imagetools inspect docker.io/biodepot/chromap-suite:v1.0.1
```

Run it from the reviewed, committed package branch. The nf-core module files
must then replace the human-readable tag with the registry manifest digest and
pass Docker plus Singularity/Apptainer tests against that public reference.

## nf-core component boundary

The first nf-core component packet uses the following composable path:

1. `chromapsuite/index` builds the genome index.
2. `chromapsuite/align` writes the standards-compliant native read-level,
   coordinate-sorted BAM, BAI and mapping summary.
3. `rapidmacs/callpeak` consumes that BAM as `BAMPE` and writes narrowPeak and
   summit files.

The first packet does not expose `chromap --atac-fragments` together with native
BAM output. That dual-output mode deliberately emits a compact fragment carrier
whose records omit read sequence, qualities and CIGAR. `samtools quickcheck`
accepts the container framing, but strict SAM/BAM readers reject the mapped
records as `INVALID_CIGAR`, so it is not suitable for an nf-core BAM channel.
Fragment export will use a separate component contract after that representation
is either made standards-compliant or explicitly modelled as a non-BAM
intermediate.

The standalone materialiser is also deferred from the first packet. It requires
a versioned set of mergeable shards plus exact sample, input and ordinal
identity; a dedicated component and fixture will be clearer than folding that
distributed gather contract into the ordinary aligner.

RapidMACS v1.0.1 does not provide a standalone `--version` option. The nf-core
module therefore reports the enclosing Chromap Suite release with `chromap
--version`; this is `1.0.1`, and the image pins RapidMACS to the matching v1.0.1
commit listed above.

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
