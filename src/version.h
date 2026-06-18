#ifndef CHROMAP_VERSION_H_
#define CHROMAP_VERSION_H_

// Chromap Suite release version (semantic versioning).
//
// This is the version reported by `chromap --version`. It versions the suite as
// a whole (ATAC alignment + native BAM, libMACS3 peaks, libchromap API, MCP
// server / Launchpad), independent of the underlying chromap engine.
//
// The chromap engine / upstream lineage version is CHROMAP_VERSION (defined in
// chromap.h); it is reported by `chromap --upstream-version` and written to the
// BAM/SAM @PG VN tag for provenance.
#define CHROMAP_SUITE_VERSION "1.0.0"

#endif  // CHROMAP_VERSION_H_
