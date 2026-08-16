#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "atac_spill_materializer.h"

namespace {

void Usage() {
  std::cerr
      << "Usage: chromap_atac_spill_materializer --spill SHARD [--spill "
         "SHARD ...] ((--output-bam BAM | --output-cram CRAM) --fragments "
         "TSV[.gz] | --output-bed BED) [options]\n"
      << "Options:\n"
      << "  --evidence FILE             Final AEV1 fragment sidecar (BAM mode)\n"
      << "  --materialized-binary FILE  Preserve canonical binary fragment blocks\n"
      << "  --summary FILE              Gathered Chromap alignment summary\n"
      << "  --reference FASTA           Reference FASTA (required for CRAM)\n"
      << "  --output-noY FILE           Secondary BAM/CRAM without Y-hit reads\n"
      << "  --output-Y FILE             Secondary BAM/CRAM with Y-hit reads\n"
      << "  --barcode-translate FILE    Barcode translation table\n"
      << "  --barcode-translate-from-first\n"
      << "  --threads N                 Output/terminal-export threads [1]\n"
      << "  --temp-dir DIR              BAM sort spill directory\n"
      << "  --sort-bam | --no-sort-bam  Coordinate-sort BAM [sort]\n"
      << "  --write-index               Write BAM index (requires sort)\n"
      << "  --read-group ID             Optional BAM read group\n";
}

std::string RequireValue(int argc, char **argv, int *index) {
  if (*index + 1 >= argc) {
    throw std::runtime_error(std::string("missing value for ") + argv[*index]);
  }
  ++(*index);
  return argv[*index];
}

}  // namespace

int main(int argc, char **argv) {
  try {
    std::vector<std::string> spills;
    chromap::MappingParameters parameters;
    parameters.mapping_output_format = chromap::MAPPINGFORMAT_UNKNOWN;
    parameters.sort_bam = true;
    parameters.num_threads = 1;

    for (int i = 1; i < argc; ++i) {
      const std::string option = argv[i];
      if (option == "--spill") {
        spills.push_back(RequireValue(argc, argv, &i));
      } else if (option == "--output-bam") {
        parameters.mapping_output_format = chromap::MAPPINGFORMAT_BAM;
        parameters.mapping_output_file_path = RequireValue(argc, argv, &i);
      } else if (option == "--output-cram") {
        parameters.mapping_output_format = chromap::MAPPINGFORMAT_CRAM;
        parameters.mapping_output_file_path = RequireValue(argc, argv, &i);
      } else if (option == "--output-bed") {
        parameters.mapping_output_format = chromap::MAPPINGFORMAT_BED;
        parameters.mapping_output_file_path = RequireValue(argc, argv, &i);
      } else if (option == "--fragments") {
        parameters.atac_fragment_output_file_path =
            RequireValue(argc, argv, &i);
      } else if (option == "--evidence") {
        parameters.atac_fragment_binary_output_file_path =
            RequireValue(argc, argv, &i);
      } else if (option == "--materialized-binary") {
        parameters.atac_materialized_binary_output_file_path =
            RequireValue(argc, argv, &i);
      } else if (option == "--summary") {
        parameters.summary_metadata_file_path = RequireValue(argc, argv, &i);
      } else if (option == "--reference") {
        parameters.reference_file_path = RequireValue(argc, argv, &i);
      } else if (option == "--output-noY") {
        parameters.emit_noY_stream = true;
        parameters.noY_output_path = RequireValue(argc, argv, &i);
      } else if (option == "--output-Y") {
        parameters.emit_Y_stream = true;
        parameters.Y_output_path = RequireValue(argc, argv, &i);
      } else if (option == "--barcode-translate") {
        parameters.barcode_translate_table_file_path =
            RequireValue(argc, argv, &i);
      } else if (option == "--barcode-translate-from-first") {
        parameters.barcode_translate_from_first_column = true;
      } else if (option == "--threads") {
        parameters.num_threads =
            std::stoi(RequireValue(argc, argv, &i));
        if (parameters.num_threads < 1) {
          throw std::runtime_error("--threads must be positive");
        }
      } else if (option == "--temp-dir") {
        parameters.temp_directory_path = RequireValue(argc, argv, &i);
      } else if (option == "--sort-bam") {
        parameters.sort_bam = true;
      } else if (option == "--no-sort-bam") {
        parameters.sort_bam = false;
      } else if (option == "--write-index") {
        parameters.write_index = true;
      } else if (option == "--read-group") {
        parameters.read_group_id = RequireValue(argc, argv, &i);
      } else if (option == "--help" || option == "-h") {
        Usage();
        return 0;
      } else {
        throw std::runtime_error("unknown option: " + option);
      }
    }

    if (spills.empty() ||
        parameters.mapping_output_format == chromap::MAPPINGFORMAT_UNKNOWN ||
        parameters.mapping_output_file_path.empty()) {
      Usage();
      return 2;
    }
    if ((parameters.mapping_output_format == chromap::MAPPINGFORMAT_BAM ||
         parameters.mapping_output_format == chromap::MAPPINGFORMAT_CRAM) &&
        parameters.atac_fragment_output_file_path.empty()) {
      throw std::runtime_error("BAM/CRAM output requires --fragments");
    }
    if (parameters.mapping_output_format == chromap::MAPPINGFORMAT_CRAM &&
        parameters.reference_file_path.empty()) {
      throw std::runtime_error("--output-cram requires --reference");
    }
    if ((parameters.emit_noY_stream && parameters.noY_output_path.empty()) ||
        (parameters.emit_Y_stream && parameters.Y_output_path.empty())) {
      throw std::runtime_error("Y-routing output path is empty");
    }
    if (parameters.write_index && !parameters.sort_bam) {
      throw std::runtime_error("--write-index requires coordinate sorting");
    }

    const chromap::AtacSpillMaterializationResult result =
        chromap::MaterializeAtacSpillRecords(spills, parameters);
    if (!result.ok) {
      std::cerr << "ATAC spill materialization failed: " << result.message
                << "\n";
      return 1;
    }
    std::cout << "ATAC spill materialization complete\n"
              << "sample_id\t" << result.sample_id << "\n"
              << "input_id\t" << result.input_id << "\n"
              << "shards\t" << result.shard_count << "\n"
              << "input_records\t" << result.input_record_count << "\n"
              << "spill_records\t" << result.spill_record_count << "\n"
              << "corrected_barcode_records\t"
              << result.corrected_barcode_record_count << "\n"
              << "rejected_barcode_records\t"
              << result.rejected_barcode_record_count << "\n"
              << "output_fragments\t" << result.output_fragment_count << "\n"
              << "parallel_hot_spill\t"
              << (result.used_parallel_hot_spill ? "yes" : "no") << "\n"
              << "merge_output_seconds\t" << result.merge_output_seconds
              << "\n"
              << "terminal_bed_export_seconds\t"
              << result.terminal_bed_export_seconds << "\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "ATAC spill materialization failed: " << error.what() << "\n";
    return 2;
  }
}
