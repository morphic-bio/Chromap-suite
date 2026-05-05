---
name: Bulk-ATAC inline FRAG peaks (pointer)
description: Pointer to libMACS3 runbook describing the Chromap-suite change to allow bulk ATAC through --call-macs3-frag-peaks
type: pointer
---

# Pointer: Bulk-ATAC inline MACS3 FRAG peaks runbook

Date: 2026-05-05
Status: design — under review

The runbook for this change lives in the libMACS3 repo because the
companion script work (the new `runChromap_libmacs3.sh` for
morphic-atac-seq) is described there too. The Chromap-suite-side code
changes are summarized in the "Code-touch summary" section of that
runbook.

## Where to read

- libMACS3 main: `plans/2026-05-05-bulk-atac-inline-peaks-runbook.md`
  (commit `b43767d`).
- GitHub: https://github.com/morphic-bio/libMACS3/blob/main/plans/2026-05-05-bulk-atac-inline-peaks-runbook.md

## What the Chromap-suite reviewer needs to look at

Three small spots, ~30 LOC total:

1. `src/chromap_driver.cc:1109-1121` — drop the `barcode_file_paths`
   requirement from the `be_bed` gate; update the error string to
   stop saying "barcoded".
2. `src/libchromap.cc:117-121` — the same gate, same edit.
3. `src/mapping_writer.cc:229-244` — add `OutputHeader` and bucket-push
   for the unbarcoded BED writer
   (`MappingWriter<PairedEndMappingWithoutBarcode>`), mirroring the
   barcoded versions at `mapping_writer.cc:281-334`.

The structural justification: `macs3::FragmentRecord` has no barcode
field; the FRAG peak pipeline does not consume barcodes. The current
gates are upstream-plumbing relics from when bucket capture was first
wired into the barcoded writer.

## Out of scope here

- libMACS3 changes (none).
- BAM dual-ATAC support for bulk (deferred — Phase 1.4 in the runbook).
- BigWig subsumption into chromap (deferred; bigwig stays in shell).
