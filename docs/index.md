## Documentation

* [README](https://github.com/morphic-bio/Chromap-suite/blob/master/README.md): platform overview, additions over Chromap, sample commands
* [HISTORY](https://github.com/morphic-bio/Chromap-suite/blob/master/HISTORY.md): project lineage (spun off from `haowenz/chromap` in 2026) and pre-spinoff fork notes
* [Manpage](chromap.html): legacy CLI option reference (inherited from Chromap; the README covers Chromap Suite's new flags)
* [BAM sort specification](sort_spec.md): coordinate-sort key, `samtools sort` differences, indexing rules
* [ATAC runtime spill schema runbook](atac_runtime_spill_schema_runbook.md): design, implementation status, and harness plan for one low-memory ATAC spill buffer across fragments, sidecar, and BAM output
* [Mergeable ATAC spill and materializer](mergeable_atac_spill.md): opt-in post-alignment shard contract, raw barcode evidence, merged correction histograms, ordinal validation, and standalone gather commands
* [Materialized reference sidecar](materialized_reference_sidecar.md): optional index-time binary reference generation, parallel mapping-time load, binding validation, and FASTA fallback
* [Direct minimizer-index loading](direct_index_loading.md): backward-compatible aligned direct I/O with automatic buffered fallback
* [Chromap Launchpad](chromap_launchpad.md): browser-based recipe builder served from the MCP server
* [MCP server](https://github.com/morphic-bio/Chromap-suite/blob/master/mcp_server/README.md): recipe registry, Launchpad API, preflight, run manifests, test tiers
* [GitHub Issues](https://github.com/morphic-bio/Chromap-suite/issues): report bugs, request features, ask questions
* [Chromap Suite preprint](https://github.com/morphic-bio/chromap_suite_paper): methodology and headline benchmarks


## Acquiring Chromap Suite

* `git clone https://github.com/morphic-bio/Chromap-suite.git`
* See [Building & Installing](https://github.com/morphic-bio/Chromap-suite/blob/master/README.md#building--installing) in the README

## Companion projects

* [RapidMACS](https://github.com/morphic-bio/rapidmacs): portable C++ implementation (`librapidmacs`) of MACS3's narrow peak caller (byte-identical to MACS3 v3.0.3)
* [STAR Suite](https://github.com/morphic-bio/STAR-suite): RNA-side counterpart; embeds `libchromap` + `librapidmacs` for the end-to-end multiomic single-binary pipeline
