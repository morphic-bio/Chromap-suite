#ifndef LIBCHROMAP_H_
#define LIBCHROMAP_H_

#include <string>

#include "mapping_parameters.h"

namespace chromap {

struct ChromapRunResult {
  bool ok = false;
  int exit_code = 1;
  std::string message;
  std::string output_path;
  std::string summary_path;
};

// Runs Chromap mapping through the callable core. This first library boundary
// preserves Chromap's file-output semantics so it can be validated directly
// against the existing CLI before any STAR-suite integration.
ChromapRunResult RunMapping(const MappingParameters &mapping_parameters);

// Narrow convenience entrypoint for the multiome ATAC path: paired-end reads
// with a barcode read and BED/TagAlign-style fragment output.
ChromapRunResult RunAtacMapping(const MappingParameters &mapping_parameters);

}  // namespace chromap

#endif  // LIBCHROMAP_H_
