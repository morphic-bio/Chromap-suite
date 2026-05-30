## Documentation

* [README](https://github.com/morphic-bio/Chromap-suite/blob/master/README.md): platform overview, additions over Chromap, sample commands
* [HISTORY](https://github.com/morphic-bio/Chromap-suite/blob/master/HISTORY.md): project lineage (spun off from `haowenz/chromap` in 2026) and pre-spinoff fork notes
* [Manpage](chromap.html): legacy CLI option reference (inherited from Chromap; the README covers Chromap Suite's new flags)
* [BAM sort specification](sort_spec.md): coordinate-sort key, `samtools sort` differences, indexing rules
* [ATAC runtime spill schema runbook](atac_runtime_spill_schema_runbook.md): design, implementation status, and harness plan for one low-memory ATAC spill buffer across fragments, sidecar, and BAM output
* [Chromap Launchpad](chromap_launchpad.md): browser-based recipe builder served from the MCP server
* [MCP server](https://github.com/morphic-bio/Chromap-suite/blob/master/mcp_server/README.md): recipe registry, Launchpad API, preflight, run manifests, test tiers
* [GitHub Issues](https://github.com/morphic-bio/Chromap-suite/issues): report bugs, request features, ask questions
* [Chromap Suite preprint](https://github.com/morphic-bio/chromap_suite_paper): methodology and headline benchmarks

## Acquiring Chromap Suite

* `git clone https://github.com/morphic-bio/Chromap-suite.git`
* See [Building & Installing](https://github.com/morphic-bio/Chromap-suite/blob/master/README.md#building--installing) in the README

## Companion projects

* [libMACS3](https://github.com/morphic-bio/libMACS3): portable C++ implementation of MACS3's narrow peak caller (byte-identical to MACS3 v3.0.3)
* [STAR Suite](https://github.com/morphic-bio/STAR-suite): RNA-side counterpart; embeds `libchromap` + `libMACS3` for the end-to-end multiomic single-binary pipeline
