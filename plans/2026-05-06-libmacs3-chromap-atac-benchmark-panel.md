# libMACS3 + Chromap ATAC Benchmark Panel

Date: 2026-05-06

Purpose: define a small public benchmark panel for comparing the integrated
Chromap + libMACS3 FRAG peak-calling path against established Chromap plus
standalone MACS2/MACS3-style ATAC workflows.

## Important Scope Fix

The existing 10x PBMC multiome fixture is useful as a continuity or compatibility
check, but it should not be counted as either core scATAC or bulk ATAC in the
paper benchmark. The core panel below uses ATAC-only ENCODE bulk ATAC datasets
and ATAC-only 10x Single Cell ATAC PBMC datasets.

## Recommended Core Panel

### Bulk ATAC: moderate ENCODE human sample

- Dataset: ENCSR803FKU
- Assay: bulk ATAC-seq
- Biosample: Homo sapiens T-helper 17 cell, male adult, 38 years
- FASTQs:
  - R1: ENCFF013DXQ, 41,980,719 reads, 100 bp, 2.36 GB
  - R2: ENCFF473SZG, 41,980,719 reads, 100 bp, 2.48 GB
- Role: first full public bulk run; already appears in local ENCODE smoke
  planning, small enough to run repeatedly.
- Source:
  - https://www.encodeproject.org/experiments/ENCSR803FKU/
  - https://www.encodeproject.org/files/ENCFF013DXQ/@@download/ENCFF013DXQ.fastq.gz
  - https://www.encodeproject.org/files/ENCFF473SZG/@@download/ENCFF473SZG.fastq.gz

### Bulk ATAC: canonical larger ENCODE human sample

- Dataset: ENCSR868FGK
- Assay: bulk ATAC-seq
- Biosample: Homo sapiens K562
- Recommended slice: biological replicate 2 only, technical replicate 1
- FASTQ pairs:
  - ENCFF098UCE / ENCFF703BGR, 29,189,692 pairs
  - ENCFF819HSY / ENCFF168MQE, 29,441,473 pairs
  - ENCFF774CPY / ENCFF160WRX, 28,956,167 pairs
  - ENCFF903FIK / ENCFF924YGE, 29,086,283 pairs
- Approximate total: 116,673,615 read pairs, 100 bp.
- Role: canonical ENCODE ATAC stress case without using the entire experiment.
- Source:
  - https://www.encodeproject.org/experiments/ENCSR868FGK/

### scATAC: Chromap-paper 10x PBMC

- Dataset: atac_v1_pbmc_10k
- Assay: 10x Genomics Single Cell ATAC, not multiome.
- Organism/sample: human PBMC from a healthy donor; 9,668 detected nuclei in
  the Cell Ranger ATAC 1.2.0 public page.
- FASTQ tar:
  - https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fastqs.tar
  - Size from HEAD check: 33.6 GB.
- Fragment file:
  - https://cf.10xgenomics.com/samples/cell-atac/1.2.0/atac_v1_pbmc_10k/atac_v1_pbmc_10k_fragments.tsv.gz
  - Size from HEAD check: 1.96 GB.
- Role: main scATAC benchmark. This is directly aligned with the Chromap paper,
  which benchmarked a 10k PBMC 10x Genomics scATAC-seq dataset, and is also used
  in Signac documentation.
- Sources:
  - https://www.nature.com/articles/s41467-021-26865-w
  - https://www.10xgenomics.com/datasets/10-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-1-2-0
  - https://stuartlab.org/signac/1.2.0/articles/pbmc_vignette

### scATAC: small true scATAC development run

- Dataset: atac_pbmc_1k_nextgem
- Assay: 10x Genomics Single Cell ATAC, not multiome.
- Organism/sample: human PBMC from a healthy donor; 1,004 detected nuclei in
  the Cell Ranger ATAC 2.0.0 public page.
- FASTQ tar:
  - https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_fastqs.tar
  - Size from HEAD check: 4.69 GB.
- Fragment file:
  - https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_1k_nextgem/atac_pbmc_1k_nextgem_fragments.tsv.gz
  - Size from HEAD check: 217 MB.
- Role: quick end-to-end true scATAC FASTQ run before spending cycles on the
  full 10k PBMC panel member.
- Source:
  - https://www.10xgenomics.com/datasets/1-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-2-0-0

## Deferred / Optional

- 10x PBMC 3K Multiome ATAC+GEX fixture: keep for regression continuity only.
  Label it as multiome ATAC, not core scATAC and not bulk ATAC.
- 10x mouse cortex ATAC-only data: biologically useful, but defer until the
  integrated wrapper exposes a non-human effective genome size. The current
  production wrapper accepts human-only genome sizes for the libMACS3 path, so
  mouse would not be an apples-to-apples MACS comparison today.
- GM12878 ENCODE ATAC: a reasonable backup if K562 download/runtime is too high,
  but it is less diverse relative to the PBMC/Th17 immune samples.

## Completed Compatibility Benchmark

The 10x PBMC 3K Multiome ATAC fixture was used as a regression and publication
workflow compatibility benchmark, not as a core scATAC or bulk ATAC panel member.

- Run ID: `20260506T_full_pbmc3k_16t`
- Dataset root:
  `/mnt/pikachu/atac-seq/benchmarks/pbmc_unsorted_3k_100k/extracted/pbmc_unsorted_3k/atac`
- Output root:
  `plans/artifacts/libmacs3_chromap_atac_panel/20260506T_full_pbmc3k_16t/full_pbmc3k/`
- Baseline: unfixed upstream Chromap `0.3.3-r519` from commit
  `c9d8ae058bfdf04d45bc5f99a164a460f759a6a7`, emitting SAM through
  `samtools view | samtools sort` for BAM and emitting FRAG rows for standalone
  MACS3 `callpeak -f FRAG`.
- Integrated: Chromap-suite `0.3.3-r519` emitting sorted BAM, BAM index,
  fragments, and libMACS3 FRAG peaks in one command.
- Threading: 16 Chromap threads and 16 samtools threads.
- Result summary:
  - normal mode: PASS; 53,969,811 baseline fragments and 53,969,811 integrated
    fragments; sorted five-column fragment MD5 matched exactly; 50,274 MACS3
    peaks and 50,274 integrated libMACS3 peaks; BED3 peak Jaccard 1.0; summit
    Jaccard 1.0; all 50,274 summits at distance 0 bp.
  - low-memory mode: PASS with the same fragment, peak, and summit parity
    metrics as normal mode.
- BAM note: standalone SAM/samtools BAM and the integrated dual-fragment BAM had
  different record counts after sorting, but the validated biological outputs
  were identical. Treat this as variation in non-biological read-level fields,
  especially read names and tie-selected duplicate representatives, rather than
  a fragment or peak failure.

## Comparison Lanes

Run benchmarks serially unless a runbook explicitly relaxes that policy.

Established workflow anchors:

- Chromap paper: https://www.nature.com/articles/s41467-021-26865-w
- ENCODE ATAC-seq pipeline: https://www.encodeproject.org/atac-seq/
- Signac PBMC scATAC vignette: https://stuartlab.org/signac/1.2.0/articles/pbmc_vignette

For each dataset:

1. Chromap mapping/preprocessing only
   - Emit BAM or fragments as appropriate.
   - Use this as the fragment identity source.

2. Chromap + standalone MACS2/MACS3-style peak calling
   - Bulk field baseline: call peaks from Chromap BAM/BAMPE using MACS2 or MACS3
     with ATAC-style parameters.
   - FRAG parity baseline: call MACS3 FRAG peaks from the same Chromap-derived
     fragments when the format is compatible. This isolates peak-caller
     integration from mapping and shifting differences.
   - scATAC baseline: call aggregate all-cell peaks from fragments. Optional
     cell-type pseudobulk peaks can be reported separately if trustworthy cell
     labels are available, but they are not the direct comparison for the current
     integrated aggregate peak caller.

3. Integrated Chromap + libMACS3
   - Use the production flags:
     - `--Tn5-shift-mode symmetric`
     - `--call-macs3-frag-peaks`
     - `--macs3-frag-peaks-source memory`
     - `--macs3-frag-pvalue 0.01`
     - `--macs3-frag-min-length 200`
     - `--macs3-frag-max-gap 30`
     - `--low-mem`
   - Add `--macs3-frag-low-mem` on large samples if needed.

## Metrics

- Wall time, user time, system time, and max RSS from `/usr/bin/time -v`.
- Output and temporary disk usage.
- Peak count and total bp covered.
- BED3 Jaccard and peak overlap.
- Reciprocal 50% peak overlap where useful.
- Summit nearest-neighbor distance.
- FRiP on the same fragment set, if feasible.
- Fragment identity between no-peaks Chromap and integrated Chromap outputs.
- BAM alignment counts or `idxstats` where BAMs are emitted.
- BAM comparisons should be interpreted after sorting and normalization. Variation
  in non-biological read-level fields, especially read names and tie-selected
  duplicate representatives, should be reported separately rather than used as a
  failure when fragment coordinates, barcode, duplicate count, and peak geometry
  are identical.

## Operational Notes

- Use `--low-mem` for production-scale runs. Avoid fully in-memory runs on large
  samples; previous production attempts OOMed around 126-129 GB RSS.
- Develop on downsampled FASTQs first, then run the full panel once commands and
  metrics are stable.
- Put benchmark outputs under:
  - `plans/artifacts/libmacs3_chromap_atac_panel/<timestamp>/<dataset>/`
- Keep large FASTQ/BAM/bigWig outputs in per-sample large-file directories and
  copy or transfer them separately.
- Record Chromap-suite, libMACS3, and wrapper git commits or dirty-worktree
  notes in every benchmark manifest.
