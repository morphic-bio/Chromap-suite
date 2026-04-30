#include <cstdint>
#include <algorithm>
#include <cctype>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "chromap.h"
#include "cxxopts.hpp"
#include "libchromap.h"
#include "utils.h"

namespace {

void ApplyPreset(const std::string &preset,
                 chromap::MappingParameters *mapping_parameters) {
  if (preset == "atac") {
    mapping_parameters->max_insert_size = 2000;
    mapping_parameters->trim_adapters = true;
    mapping_parameters->remove_pcr_duplicates = true;
    mapping_parameters->remove_pcr_duplicates_at_bulk_level = false;
    mapping_parameters->Tn5_shift = true;
    mapping_parameters->mapping_output_format = chromap::MAPPINGFORMAT_BED;
    mapping_parameters->low_memory_mode = true;
  } else if (preset == "chip") {
    mapping_parameters->max_insert_size = 2000;
    mapping_parameters->remove_pcr_duplicates = true;
    mapping_parameters->low_memory_mode = true;
    mapping_parameters->mapping_output_format = chromap::MAPPINGFORMAT_BED;
  } else if (preset == "hic") {
    mapping_parameters->error_threshold = 4;
    mapping_parameters->mapq_threshold = 1;
    mapping_parameters->split_alignment = true;
    mapping_parameters->low_memory_mode = true;
    mapping_parameters->mapping_output_format = chromap::MAPPINGFORMAT_PAIRS;
  } else {
    chromap::ExitWithMessage("Unrecognized preset parameters " + preset);
  }
}

uint64_t ParseSizeString(const std::string &value,
                         const std::string &option_name) {
  if (value.empty()) {
    chromap::ExitWithMessage(option_name + " requires a non-empty value");
  }

  size_t suffix_pos = value.size();
  while (suffix_pos > 0 &&
         ((value[suffix_pos - 1] >= 'A' && value[suffix_pos - 1] <= 'Z') ||
          (value[suffix_pos - 1] >= 'a' && value[suffix_pos - 1] <= 'z'))) {
    --suffix_pos;
  }

  const std::string number_part = value.substr(0, suffix_pos);
  const std::string suffix = value.substr(suffix_pos);
  if (number_part.empty()) {
    chromap::ExitWithMessage(option_name + " requires a numeric value");
  }

  std::stringstream parser(number_part);
  double number = 0.0;
  parser >> number;
  if (!parser || !parser.eof() || number < 0.0) {
    chromap::ExitWithMessage("Invalid size for " + option_name + ": " + value);
  }

  uint64_t multiplier = 1;
  if (suffix.empty() || suffix == "B" || suffix == "b") {
    multiplier = 1;
  } else if (suffix == "K" || suffix == "KB" || suffix == "k" ||
             suffix == "kb") {
    multiplier = 1024ULL;
  } else if (suffix == "M" || suffix == "MB" || suffix == "m" ||
             suffix == "mb") {
    multiplier = 1024ULL * 1024ULL;
  } else if (suffix == "G" || suffix == "GB" || suffix == "g" ||
             suffix == "gb") {
    multiplier = 1024ULL * 1024ULL * 1024ULL;
  } else {
    chromap::ExitWithMessage("Invalid suffix for " + option_name + ": " +
                             suffix);
  }

  return static_cast<uint64_t>(number * multiplier);
}

std::string DeriveSecondaryOutputPath(const std::string &primary_path,
                                      const std::string &suffix) {
  if (primary_path == "/dev/stdout" || primary_path == "/dev/stderr" ||
      primary_path == "-") {
    chromap::ExitWithMessage(
        "Y/noY stream outputs require explicit secondary output paths when "
        "primary output is stdout");
    return "";
  }

  const size_t slash_pos = primary_path.rfind('/');
  const size_t search_start =
      slash_pos == std::string::npos ? 0 : slash_pos + 1;

  std::string lower_path = primary_path;
  std::transform(lower_path.begin(), lower_path.end(), lower_path.begin(),
                 [](unsigned char c) { return std::tolower(c); });

  if (lower_path.size() > 7 &&
      lower_path.substr(lower_path.size() - 7) == ".sam.gz") {
    return primary_path.substr(0, primary_path.size() - 7) + suffix +
           ".sam.gz";
  }
  if (lower_path.size() > 4 &&
      lower_path.substr(lower_path.size() - 4) == ".bam") {
    return primary_path.substr(0, primary_path.size() - 4) + suffix + ".bam";
  }
  if (lower_path.size() > 5 &&
      lower_path.substr(lower_path.size() - 5) == ".cram") {
    return primary_path.substr(0, primary_path.size() - 5) + suffix + ".cram";
  }
  if (lower_path.size() > 4 &&
      lower_path.substr(lower_path.size() - 4) == ".sam") {
    return primary_path.substr(0, primary_path.size() - 4) + suffix + ".sam";
  }

  const size_t dot_pos = primary_path.rfind('.');
  if (dot_pos != std::string::npos && dot_pos > search_start) {
    return primary_path.substr(0, dot_pos) + suffix +
           primary_path.substr(dot_pos);
  }
  return primary_path + suffix + ".sam";
}

std::string DeriveYReadNamesOutputPath(const std::string &primary_path) {
  if (primary_path == "/dev/stdout" || primary_path == "/dev/stderr" ||
      primary_path == "-") {
    chromap::ExitWithMessage(
        "--emit-Y-read-names requires --Y-read-names-output when primary "
        "output is stdout");
    return "";
  }

  const size_t dot_pos = primary_path.rfind('.');
  if (dot_pos == std::string::npos) {
    return primary_path + ".Y.names.txt";
  }
  return primary_path.substr(0, dot_pos) + ".Y.names.txt";
}

void ValidateInputs(const chromap::MappingParameters &mapping_parameters) {
  if (mapping_parameters.reference_file_path.empty()) {
    chromap::ExitWithMessage("No reference specified");
  }
  if (mapping_parameters.index_file_path.empty()) {
    chromap::ExitWithMessage("No index specified");
  }
  if (mapping_parameters.mapping_output_file_path.empty()) {
    chromap::ExitWithMessage("No output specified");
  }
  if (mapping_parameters.read_file1_paths.empty()) {
    chromap::ExitWithMessage("No read 1 input specified");
  }
  if (!mapping_parameters.read_file2_paths.empty() &&
      mapping_parameters.read_file1_paths.size() !=
          mapping_parameters.read_file2_paths.size()) {
    chromap::ExitWithMessage("Read 1 and read 2 input counts differ");
  }
  if (!mapping_parameters.barcode_file_paths.empty() &&
      mapping_parameters.read_file1_paths.size() !=
          mapping_parameters.barcode_file_paths.size()) {
    chromap::ExitWithMessage("Read 1 and barcode input counts differ");
  }
  if (mapping_parameters.barcode_file_paths.empty() &&
      !mapping_parameters.barcode_whitelist_file_path.empty()) {
    chromap::ExitWithMessage(
        "Barcode whitelist was supplied without barcode reads");
  }
  if (mapping_parameters.write_index && !mapping_parameters.sort_bam) {
    chromap::ExitWithMessage("--write-index requires --sort-bam");
  }
  const bool primary_is_stdout =
      mapping_parameters.mapping_output_file_path == "-" ||
      mapping_parameters.mapping_output_file_path == "/dev/stdout" ||
      mapping_parameters.mapping_output_file_path == "/dev/stderr";
  if (primary_is_stdout && mapping_parameters.emit_noY_stream &&
      mapping_parameters.noY_output_path.empty()) {
    chromap::ExitWithMessage(
        "--emit-noY-bam requires --noY-output when primary output is stdout");
  }
  if (primary_is_stdout && mapping_parameters.emit_Y_stream &&
      mapping_parameters.Y_output_path.empty()) {
    chromap::ExitWithMessage(
        "--emit-Y-bam requires --Y-output when primary output is stdout");
  }
}

void PrintRunSummary(const chromap::MappingParameters &mapping_parameters) {
  std::cerr << "chromap_lib_runner launching libchromap\n";
  std::cerr << "Reference file: " << mapping_parameters.reference_file_path
            << "\n";
  std::cerr << "Index file: " << mapping_parameters.index_file_path << "\n";
  std::cerr << "Output file: " << mapping_parameters.mapping_output_file_path
            << "\n";
  if (!mapping_parameters.summary_metadata_file_path.empty()) {
    std::cerr << "Summary file: "
              << mapping_parameters.summary_metadata_file_path << "\n";
  }
  std::cerr << "Number of threads: " << mapping_parameters.num_threads << "\n";
}

}  // namespace

int main(int argc, char **argv) {
  cxxopts::Options options(
      "chromap_lib_runner",
      "Parity harness for running Chromap through libchromap");

  options.set_width(120).add_options()
      ("h,help", "Print help")
      ("v,version", "Print Chromap version")
      ("preset", "Preset parameters: atac, chip, or hic",
       cxxopts::value<std::string>(), "STR")
      ("t,num-threads", "Number of threads",
       cxxopts::value<int>(), "INT")
      ("x,index", "Chromap index file",
       cxxopts::value<std::string>(), "FILE")
      ("r,ref", "Reference FASTA",
       cxxopts::value<std::string>(), "FILE")
      ("1,read1", "Read 1 FASTQ file(s), comma separated",
       cxxopts::value<std::vector<std::string>>(), "FILE[,FILE]")
      ("2,read2", "Read 2 FASTQ file(s), comma separated",
       cxxopts::value<std::vector<std::string>>(), "FILE[,FILE]")
      ("b,barcode", "Barcode FASTQ file(s), comma separated",
       cxxopts::value<std::vector<std::string>>(), "FILE[,FILE]")
      ("barcode-whitelist", "Barcode whitelist file",
       cxxopts::value<std::string>(), "FILE")
      ("barcode-translate", "Barcode translation table",
       cxxopts::value<std::string>(), "FILE")
      ("barcode-translate-from-first",
       "Treat the translation table as <from_bc>\\t<to_bc> (col1 is the "
       "hash key / source). Default is the historical Chromap convention "
       "<to_bc>\\t<from_bc> (col2 is the hash key).",
       cxxopts::value<bool>()->default_value("false"))
      ("o,output", "Mapping output path",
       cxxopts::value<std::string>(), "FILE")
      ("atac-fragments", "Secondary ATAC fragments output path",
       cxxopts::value<std::string>(), "FILE")
      ("summary", "Summary metadata output path",
       cxxopts::value<std::string>(), "FILE")
      ("temp-dir", "Temporary directory",
       cxxopts::value<std::string>(), "DIR")
      ("read-format", "Read format string",
       cxxopts::value<std::string>(), "STR")
      ("chr-order", "Custom chromosome order file",
       cxxopts::value<std::string>(), "FILE")
      ("pairs-natural-chr-order", "Pairs natural chromosome order file",
       cxxopts::value<std::string>(), "FILE")
      ("e,error-threshold", "Maximum read mapping errors",
       cxxopts::value<int>(), "INT")
      ("s,min-num-seeds", "Minimum seeds required for mapping",
       cxxopts::value<int>(), "INT")
      ("f,max-seed-frequencies", "Maximum seed frequencies",
       cxxopts::value<std::vector<int>>(), "INT[,INT]")
      ("l,max-insert-size", "Maximum insert size",
       cxxopts::value<int>(), "INT")
      ("q,MAPQ-threshold", "Minimum MAPQ",
       cxxopts::value<int>(), "INT")
      ("min-read-length", "Minimum read length",
       cxxopts::value<int>(), "INT")
      ("bc-error-threshold", "Barcode correction error threshold",
       cxxopts::value<int>(), "INT")
      ("bc-probability-threshold", "Barcode correction probability threshold",
       cxxopts::value<double>(), "FLOAT")
      ("drop-repetitive-reads", "Drop reads with more mappings than this",
       cxxopts::value<int>(), "INT")
      ("trim-adapters", "Trim adapters")
      ("remove-pcr-duplicates", "Remove PCR duplicates")
      ("remove-pcr-duplicates-at-bulk-level",
       "Remove PCR duplicates at bulk level")
      ("remove-pcr-duplicates-at-cell-level",
       "Remove PCR duplicates at cell level")
      ("allocate-multi-mappings", "Allocate multi-mappings")
      ("output-mappings-not-in-whitelist",
       "Output mappings even when barcode is not in whitelist")
      ("Tn5-shift", "Apply Tn5 shift")
      ("Tn5-shift-mode", "Tn5 shift mode: classical or symmetric",
       cxxopts::value<std::string>(), "STR")
      ("split-alignment", "Allow split alignment")
      ("BED", "Output BED/BEDPE")
      ("TagAlign", "Output TagAlign/PairedTagAlign")
      ("PAF", "Output PAF")
      ("SAM", "Output SAM")
      ("BAM", "Output BAM")
      ("CRAM", "Output CRAM")
      ("pairs", "Output pairs")
      ("low-mem", "Use low memory mode")
      ("low-mem-ram", "Low-memory spill threshold",
       cxxopts::value<std::string>(), "SIZE")
      ("hts-threads", "Htslib worker threads",
       cxxopts::value<int>(), "INT")
      ("read-group", "Read group ID",
       cxxopts::value<std::string>(), "STR")
      ("sort-bam", "Coordinate-sort BAM/CRAM output")
      ("sort-bam-ram", "Sort RAM limit",
       cxxopts::value<std::string>(), "SIZE")
      ("write-index", "Write BAM/CRAM index")
      ("emit-noY-bam", "Emit additional stream excluding Y-chromosome reads")
      ("noY-output", "Explicit noY output path",
       cxxopts::value<std::string>(), "FILE")
      ("emit-Y-bam", "Emit additional stream with only Y-chromosome reads")
      ("Y-output", "Explicit Y-only output path",
       cxxopts::value<std::string>(), "FILE")
      ("emit-Y-read-names", "Emit Y read names")
      ("Y-read-names-output", "Explicit Y read names output path",
       cxxopts::value<std::string>(), "FILE")
      ("emit-Y-noY-fastq", "Emit Y/noY FASTQ files")
      ("emit-Y-noY-fastq-compression", "Y/noY FASTQ compression: gz or none",
       cxxopts::value<std::string>(), "STR")
      ("Y-fastq-output-prefix", "Prefix for Y FASTQ outputs",
       cxxopts::value<std::string>(), "PREFIX")
      ("noY-fastq-output-prefix", "Prefix for noY FASTQ outputs",
       cxxopts::value<std::string>(), "PREFIX")
      ("skip-barcode-check", "Skip barcode compatibility check");

  chromap::MappingParameters mapping_parameters;

  try {
    const auto result = options.parse(argc, argv);

    if (result.count("help")) {
      std::cerr << options.help() << "\n";
      return 0;
    }
    if (result.count("version")) {
      std::cerr << CHROMAP_VERSION << "\n";
      return 0;
    }

    if (result.count("preset")) {
      ApplyPreset(result["preset"].as<std::string>(), &mapping_parameters);
    }
    if (result.count("num-threads")) {
      mapping_parameters.num_threads = result["num-threads"].as<int>();
    }
    if (result.count("index")) {
      mapping_parameters.index_file_path = result["index"].as<std::string>();
    }
    if (result.count("ref")) {
      mapping_parameters.reference_file_path = result["ref"].as<std::string>();
    }
    if (result.count("read1")) {
      mapping_parameters.read_file1_paths =
          result["read1"].as<std::vector<std::string>>();
    }
    if (result.count("read2")) {
      mapping_parameters.read_file2_paths =
          result["read2"].as<std::vector<std::string>>();
    }
    if (result.count("barcode")) {
      mapping_parameters.is_bulk_data = false;
      mapping_parameters.barcode_file_paths =
          result["barcode"].as<std::vector<std::string>>();
    }
    if (result.count("barcode-whitelist")) {
      mapping_parameters.barcode_whitelist_file_path =
          result["barcode-whitelist"].as<std::string>();
    }
    if (result.count("barcode-translate")) {
      mapping_parameters.barcode_translate_table_file_path =
          result["barcode-translate"].as<std::string>();
    }
    mapping_parameters.barcode_translate_from_first_column =
        result["barcode-translate-from-first"].as<bool>();
    if (result.count("output")) {
      mapping_parameters.mapping_output_file_path =
          result["output"].as<std::string>();
    }
    if (result.count("atac-fragments")) {
      mapping_parameters.atac_fragment_output_file_path =
          result["atac-fragments"].as<std::string>();
    }
    if (result.count("summary")) {
      mapping_parameters.summary_metadata_file_path =
          result["summary"].as<std::string>();
    }
    if (result.count("temp-dir")) {
      mapping_parameters.temp_directory_path =
          result["temp-dir"].as<std::string>();
    }
    if (result.count("read-format")) {
      mapping_parameters.read_format = result["read-format"].as<std::string>();
    }
    if (result.count("chr-order")) {
      mapping_parameters.custom_rid_order_file_path =
          result["chr-order"].as<std::string>();
    }
    if (result.count("pairs-natural-chr-order")) {
      mapping_parameters.pairs_flipping_custom_rid_order_file_path =
          result["pairs-natural-chr-order"].as<std::string>();
    }
    if (result.count("error-threshold")) {
      mapping_parameters.error_threshold =
          result["error-threshold"].as<int>();
    }
    if (result.count("min-num-seeds")) {
      mapping_parameters.min_num_seeds_required_for_mapping =
          result["min-num-seeds"].as<int>();
    }
    if (result.count("max-seed-frequencies")) {
      mapping_parameters.max_seed_frequencies =
          result["max-seed-frequencies"].as<std::vector<int>>();
    }
    if (result.count("max-insert-size")) {
      mapping_parameters.max_insert_size =
          result["max-insert-size"].as<int>();
    }
    if (result.count("MAPQ-threshold")) {
      const int mapq = result["MAPQ-threshold"].as<int>();
      if (mapq < 0 || mapq > 255) {
        chromap::ExitWithMessage("--MAPQ-threshold must be in [0, 255]");
      }
      mapping_parameters.mapq_threshold = static_cast<uint8_t>(mapq);
    }
    if (result.count("min-read-length")) {
      mapping_parameters.min_read_length =
          result["min-read-length"].as<int>();
    }
    if (result.count("bc-error-threshold")) {
      mapping_parameters.barcode_correction_error_threshold =
          result["bc-error-threshold"].as<int>();
    }
    if (result.count("bc-probability-threshold")) {
      mapping_parameters.barcode_correction_probability_threshold =
          result["bc-probability-threshold"].as<double>();
    }
    if (result.count("drop-repetitive-reads")) {
      mapping_parameters.drop_repetitive_reads =
          result["drop-repetitive-reads"].as<int>();
    }
    if (result.count("trim-adapters")) {
      mapping_parameters.trim_adapters = true;
    }
    if (result.count("remove-pcr-duplicates")) {
      mapping_parameters.remove_pcr_duplicates = true;
    }
    if (result.count("remove-pcr-duplicates-at-bulk-level")) {
      mapping_parameters.remove_pcr_duplicates_at_bulk_level = true;
    }
    if (result.count("remove-pcr-duplicates-at-cell-level")) {
      mapping_parameters.remove_pcr_duplicates_at_bulk_level = false;
    }
    if (result.count("allocate-multi-mappings")) {
      mapping_parameters.allocate_multi_mappings = true;
      mapping_parameters.only_output_unique_mappings = false;
    }
    if (result.count("output-mappings-not-in-whitelist")) {
      mapping_parameters.output_mappings_not_in_whitelist = true;
    }
    if (result.count("Tn5-shift")) {
      mapping_parameters.Tn5_shift = true;
    }
    if (result.count("Tn5-shift-mode")) {
      const std::string mode = result["Tn5-shift-mode"].as<std::string>();
      if (mode == "classical") {
        mapping_parameters.Tn5_forward_shift = 4;
        mapping_parameters.Tn5_reverse_shift = -5;
      } else if (mode == "symmetric") {
        mapping_parameters.Tn5_forward_shift = 4;
        mapping_parameters.Tn5_reverse_shift = -4;
      } else {
        chromap::ExitWithMessage(
            "Unrecognized --Tn5-shift-mode '" + mode + "'");
      }
      mapping_parameters.Tn5_shift = true;
    }
    if (result.count("split-alignment")) {
      mapping_parameters.split_alignment = true;
    }
    if (result.count("BED")) {
      mapping_parameters.mapping_output_format = chromap::MAPPINGFORMAT_BED;
    }
    if (result.count("TagAlign")) {
      mapping_parameters.mapping_output_format =
          chromap::MAPPINGFORMAT_TAGALIGN;
    }
    if (result.count("PAF")) {
      mapping_parameters.mapping_output_format = chromap::MAPPINGFORMAT_PAF;
    }
    if (result.count("SAM")) {
      mapping_parameters.mapping_output_format = chromap::MAPPINGFORMAT_SAM;
    }
    if (result.count("BAM")) {
      mapping_parameters.mapping_output_format = chromap::MAPPINGFORMAT_BAM;
    }
    if (result.count("CRAM")) {
      mapping_parameters.mapping_output_format = chromap::MAPPINGFORMAT_CRAM;
    }
    if (result.count("pairs")) {
      mapping_parameters.mapping_output_format = chromap::MAPPINGFORMAT_PAIRS;
    }
    if (result.count("low-mem")) {
      mapping_parameters.low_memory_mode = true;
    }
    if (result.count("low-mem-ram")) {
      mapping_parameters.low_mem_ram_limit =
          ParseSizeString(result["low-mem-ram"].as<std::string>(),
                          "--low-mem-ram");
    }
    if (result.count("hts-threads")) {
      mapping_parameters.hts_threads = result["hts-threads"].as<int>();
    }
    if (result.count("read-group")) {
      mapping_parameters.read_group_id =
          result["read-group"].as<std::string>();
    }
    if (result.count("sort-bam")) {
      mapping_parameters.sort_bam = true;
    }
    if (result.count("sort-bam-ram")) {
      mapping_parameters.sort_bam_ram_limit =
          ParseSizeString(result["sort-bam-ram"].as<std::string>(),
                          "--sort-bam-ram");
    }
    if (result.count("write-index")) {
      mapping_parameters.write_index = true;
    }
    if (result.count("emit-noY-bam")) {
      mapping_parameters.emit_noY_stream = true;
    }
    if (result.count("noY-output")) {
      mapping_parameters.noY_output_path =
          result["noY-output"].as<std::string>();
    } else if (mapping_parameters.emit_noY_stream) {
      mapping_parameters.noY_output_path = DeriveSecondaryOutputPath(
          mapping_parameters.mapping_output_file_path, ".noY");
    }
    if (result.count("emit-Y-bam")) {
      mapping_parameters.emit_Y_stream = true;
    }
    if (result.count("Y-output")) {
      mapping_parameters.Y_output_path =
          result["Y-output"].as<std::string>();
    } else if (mapping_parameters.emit_Y_stream) {
      mapping_parameters.Y_output_path = DeriveSecondaryOutputPath(
          mapping_parameters.mapping_output_file_path, ".Y");
    }
    if (result.count("emit-Y-read-names")) {
      mapping_parameters.emit_y_read_names = true;
    }
    if (result.count("Y-read-names-output")) {
      mapping_parameters.y_read_names_output_path =
          result["Y-read-names-output"].as<std::string>();
    } else if (mapping_parameters.emit_y_read_names) {
      mapping_parameters.y_read_names_output_path =
          DeriveYReadNamesOutputPath(
              mapping_parameters.mapping_output_file_path);
    }
    if (result.count("emit-Y-noY-fastq")) {
      mapping_parameters.emit_y_noy_fastq = true;
    }
    if (result.count("emit-Y-noY-fastq-compression")) {
      mapping_parameters.y_noy_fastq_compression =
          result["emit-Y-noY-fastq-compression"].as<std::string>();
    }
    if (result.count("Y-fastq-output-prefix")) {
      mapping_parameters.y_fastq_output_prefix =
          result["Y-fastq-output-prefix"].as<std::string>();
    }
    if (result.count("noY-fastq-output-prefix")) {
      mapping_parameters.noy_fastq_output_prefix =
          result["noY-fastq-output-prefix"].as<std::string>();
    }
    if (result.count("skip-barcode-check")) {
      mapping_parameters.skip_barcode_check = true;
    }

    ValidateInputs(mapping_parameters);
    PrintRunSummary(mapping_parameters);

    const chromap::ChromapRunResult run_result =
        chromap::RunMapping(mapping_parameters);
    if (!run_result.ok) {
      std::cerr << "libchromap failed: " << run_result.message << "\n";
      return run_result.exit_code == 0 ? 1 : run_result.exit_code;
    }
  } catch (const cxxopts::OptionException &error) {
    std::cerr << "chromap_lib_runner option error: " << error.what() << "\n";
    return 1;
  }

  return 0;
}
