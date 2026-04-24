#include "libchromap.h"

#include <exception>
#include <string>

#include "chromap.h"
#include "utils.h"

namespace chromap {
namespace {

ChromapRunResult MakeSuccess(const MappingParameters &mapping_parameters) {
  ChromapRunResult result;
  result.ok = true;
  result.exit_code = 0;
  result.message = "ok";
  result.output_path = mapping_parameters.mapping_output_file_path;
  result.summary_path = mapping_parameters.summary_metadata_file_path;
  return result;
}

ChromapRunResult MakeFailure(const MappingParameters &mapping_parameters,
                             const std::string &message) {
  ChromapRunResult result;
  result.ok = false;
  result.exit_code = 1;
  result.message = message;
  result.output_path = mapping_parameters.mapping_output_file_path;
  result.summary_path = mapping_parameters.summary_metadata_file_path;
  return result;
}

}  // namespace

ChromapRunResult RunMapping(const MappingParameters &mapping_parameters) {
  try {
    Chromap chromap_for_mapping(mapping_parameters);

    if (mapping_parameters.read_file2_paths.empty()) {
      switch (mapping_parameters.mapping_output_format) {
        case MAPPINGFORMAT_PAF:
          chromap_for_mapping.MapSingleEndReads<PAFMapping>();
          break;
        case MAPPINGFORMAT_SAM:
        case MAPPINGFORMAT_BAM:
        case MAPPINGFORMAT_CRAM:
          chromap_for_mapping.MapSingleEndReads<SAMMapping>();
          break;
        case MAPPINGFORMAT_PAIRS:
          return MakeFailure(mapping_parameters,
                             "single-end PAIRS output is not supported");
        case MAPPINGFORMAT_BED:
        case MAPPINGFORMAT_TAGALIGN:
          if (!mapping_parameters.barcode_file_paths.empty()) {
            chromap_for_mapping.MapSingleEndReads<MappingWithBarcode>();
          } else {
            chromap_for_mapping.MapSingleEndReads<MappingWithoutBarcode>();
          }
          break;
        default:
          return MakeFailure(mapping_parameters, "unknown mapping output format");
      }
    } else {
      switch (mapping_parameters.mapping_output_format) {
        case MAPPINGFORMAT_PAF:
          chromap_for_mapping.MapPairedEndReads<PairedPAFMapping>();
          break;
        case MAPPINGFORMAT_SAM:
        case MAPPINGFORMAT_BAM:
        case MAPPINGFORMAT_CRAM:
          chromap_for_mapping.MapPairedEndReads<SAMMapping>();
          break;
        case MAPPINGFORMAT_PAIRS:
          chromap_for_mapping.MapPairedEndReads<PairsMapping>();
          break;
        case MAPPINGFORMAT_BED:
        case MAPPINGFORMAT_TAGALIGN:
          if (!mapping_parameters.barcode_file_paths.empty()) {
            chromap_for_mapping.MapPairedEndReads<PairedEndMappingWithBarcode>();
          } else {
            chromap_for_mapping.MapPairedEndReads<
                PairedEndMappingWithoutBarcode>();
          }
          break;
        default:
          return MakeFailure(mapping_parameters, "unknown mapping output format");
      }
    }
  } catch (const std::exception &error) {
    return MakeFailure(mapping_parameters, error.what());
  } catch (...) {
    return MakeFailure(mapping_parameters, "unknown Chromap mapping failure");
  }

  return MakeSuccess(mapping_parameters);
}

ChromapRunResult RunAtacMapping(const MappingParameters &mapping_parameters) {
  if (mapping_parameters.read_file1_paths.empty()) {
    return MakeFailure(mapping_parameters, "ATAC mapping requires read 1 input");
  }
  if (mapping_parameters.read_file2_paths.empty()) {
    return MakeFailure(mapping_parameters, "ATAC mapping requires read 2 input");
  }
  if (mapping_parameters.barcode_file_paths.empty()) {
    return MakeFailure(mapping_parameters,
                       "ATAC mapping requires barcode input");
  }
  if (mapping_parameters.mapping_output_format != MAPPINGFORMAT_BED &&
      mapping_parameters.mapping_output_format != MAPPINGFORMAT_TAGALIGN) {
    return MakeFailure(mapping_parameters,
                       "ATAC mapping requires BED or TagAlign output");
  }

  return RunMapping(mapping_parameters);
}

}  // namespace chromap
