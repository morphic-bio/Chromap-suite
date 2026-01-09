#include "chromap_driver.h"

#include <glob.h>
#include <cctype>

#include <algorithm>
#include <array>
#include <cassert>
#include <iomanip>
#include <regex>
#include <string>
#include <vector>

#include "chromap.h"
#include "cxxopts.hpp"

namespace chromap {
namespace {

void AddIndexingOptions(cxxopts::Options &options) {
  options.add_options("Indexing")("i,build-index", "Build index")(
      "min-frag-length",
      "Min fragment length for choosing k and w automatically [30]",
      cxxopts::value<int>(),
      "INT")("k,kmer", "Kmer length [17]", cxxopts::value<int>(), "INT")(
      "w,window", "Window size [7]", cxxopts::value<int>(), "INT");
}

void AddMappingOptions(cxxopts::Options &options) {
  options.set_width(120).add_options("Mapping")(
      "preset",
      "Preset parameters for mapping reads (always applied before other "
      "options) []\natac: mapping ATAC-seq/scATAC-seq reads\nchip: mapping "
      "ChIP-seq reads\nhic: mapping Hi-C reads",
      cxxopts::value<std::string>(),
      "STR")("split-alignment", "Allow split alignments")(
      "e,error-threshold", "Max # errors allowed to map a read [8]",
      cxxopts::value<int>(), "INT")
      //("A,match-score", "Match score [1]", cxxopts::value<int>(), "INT")
      //("B,mismatch-penalty", "Mismatch penalty [4]", cxxopts::value<int>(),
      //"INT")
      //("O,gap-open-penalties", "Gap open penalty [6,6]",
      // cxxopts::value<std::vector<int>>(), "INT[,INT]")
      //("E,gap-extension-penalties", "Gap extension penalty [1,1]",
      // cxxopts::value<std::vector<int>>(), "INT[,INT]")
      ("s,min-num-seeds", "Min # seeds to try to map a read [2]",
       cxxopts::value<int>(),
       "INT")("f,max-seed-frequencies",
              "Max seed frequencies for a seed to be selected [500,1000]",
              cxxopts::value<std::vector<int>>(), "INT[,INT]")
      //("n,max-num-best-mappings", "Only report n best mappings [1]",
      // cxxopts::value<int>(), "INT")
      ("l,max-insert-size",
       "Max insert size, only for paired-end read mapping [1000]",
       cxxopts::value<int>(),
       "INT")("q,MAPQ-threshold",
              "Min MAPQ in range [0, 60] for mappings to be output [30]",
              cxxopts::value<uint8_t>(),
              "INT")("min-read-length", "Min read length [30]",
                     cxxopts::value<int>(), "INT")
      //("multi-mapping-allocation-distance", "Uni-mappings within this distance
      // from any end of multi-mappings are used for allocation [0]",
      // cxxopts::value<int>(), "INT")
      //("multi-mapping-allocation-seed", "Seed for random number generator in
      // multi-mapping allocation [11]", cxxopts::value<int>(), "INT")
      //("drop-repetitive-reads", "Drop reads with too many best mappings
      //[500000]", cxxopts::value<int>(), "INT")
      ("trim-adapters", "Try to trim adapters on 3'")("remove-pcr-duplicates",
                                                      "Remove PCR duplicates")(
          "remove-pcr-duplicates-at-bulk-level",
          "Remove PCR duplicates at bulk level for single cell data")(
          "remove-pcr-duplicates-at-cell-level",
          "Remove PCR duplicates at cell level for single cell data")
      //("allocate-multi-mappings", "Allocate multi-mappings")
      ("Tn5-shift", "Perform Tn5 shift")("low-mem", "Use low memory mode")(
          "low-mem-ram",
          "Max RAM for low-mem spill buffer before writing temp files [default: 1G; 512M for SAM/PAF/PAIRS]",
          cxxopts::value<std::string>(), "SIZE")(
          "bc-error-threshold",
          "Max Hamming distance allowed to correct a barcode [1]",
          cxxopts::value<int>(),
          "INT")("bc-probability-threshold",
                 "Min probability to correct a barcode [0.9]",
                 cxxopts::value<double>(),
                 "FLT")("t,num-threads", "# threads for mapping [1]",
                        cxxopts::value<int>(), "INT")
      ("frip-est-params", "coefficients used for frip est calculation, separated by semi-colons",
      cxxopts::value<std::string>(), "STR")
      ("turn-off-num-uniq-cache-slots", "turn off the output of number of cache slots in summary file")
      ("temp-dir", "Directory for temporary files [system temp directory]",
      cxxopts::value<std::string>(), "DIR");
}

void AddInputOptions(cxxopts::Options &options) {
  options.add_options("Input")("r,ref", "Reference file",
                               cxxopts::value<std::string>(), "FILE")(
      "x,index", "Index file", cxxopts::value<std::string>(), "FILE")(
      "1,read1", "Single-end read files or paired-end read files 1",
      cxxopts::value<std::vector<std::string>>(),
      "FILE")("2,read2", "Paired-end read files 2",
              cxxopts::value<std::vector<std::string>>(),
              "FILE")("b,barcode", "Cell barcode files",
                      cxxopts::value<std::vector<std::string>>(), "FILE")(
      "barcode-whitelist", "Cell barcode whitelist file",
      cxxopts::value<std::string>(),
      "FILE")("read-format",
              "Format for read files and barcode files  [\"r1:0:-1,bc:0:-1\" "
              "as 10x Genomics single-end format]",
              cxxopts::value<std::string>(), "STR");
}

void AddOutputOptions(cxxopts::Options &options) {
  options.add_options("Output")("o,output", "Output file",
                                cxxopts::value<std::string>(), "FILE")
      //("p,matrix-output-prefix", "Prefix of matrix output files",
      // cxxopts::value<std::string>(), "FILE")
      ("output-mappings-not-in-whitelist",
       "Output mappings with barcode not in the whitelist")(
          "chr-order",
          "Custom chromosome order file. If not specified, the order of "
          "reference sequences will be used",
          cxxopts::value<std::string>(),
          "FILE")("BED", "Output mappings in BED/BEDPE format")(
          "TagAlign", "Output mappings in TagAlign/PairedTagAlign format")(
          "SAM", "Output mappings in SAM format")(
          "BAM", "Output mappings in BAM format")(
          "CRAM", "Output mappings in CRAM format (requires -r reference)")(
          "pairs",
          "Output mappings in pairs format (defined by 4DN for HiC data)")(
          "pairs-natural-chr-order",
          "Custom chromosome order file for pairs flipping. If not specified, "
          "the custom chromosome order will be used",
          cxxopts::value<std::string>(),
          "FILE")(          "barcode-translate",
                  "Convert barcode to the specified sequences during output",
                  cxxopts::value<std::string>(), "FILE")(
          "summary",
          "Summarize the mapping statistics at bulk or barcode level",
          cxxopts::value<std::string>(), "FILE")(
          "emit-noY-bam",
          "Emit additional stream excluding Y-chromosome reads (requires --SAM/--BAM/--CRAM)")(
          "noY-output",
          "Explicit path for noY output [default: <output>.noY.sam/.bam/.cram]",
          cxxopts::value<std::string>(), "FILE")(
          "emit-Y-bam",
          "Emit additional stream with only Y-chromosome reads (requires --SAM/--BAM/--CRAM)")(
          "Y-output",
          "Explicit path for Y-only output [default: <output>.Y.sam/.bam/.cram]",
          cxxopts::value<std::string>(), "FILE")(
          "emit-Y-read-names",
          "Emit list of read names with Y-chromosome alignments")(
          "Y-read-names-output",
          "Explicit path for Y read names output [default: <output>.Y.names.txt]",
          cxxopts::value<std::string>(), "FILE")(
          "emit-Y-noY-fastq",
          "Emit Y/noY FASTQ files split by Y-chromosome alignment status")(
          "emit-Y-noY-fastq-compression",
          "Compression for Y/noY FASTQ files [gz|none, default: gz]",
          cxxopts::value<std::string>(), "STR")(
          "Y-fastq-output-prefix",
          "Prefix for Y FASTQ output files (overrides auto-naming)",
          cxxopts::value<std::string>(), "PREFIX")(
          "noY-fastq-output-prefix",
          "Prefix for noY FASTQ output files (overrides auto-naming)",
          cxxopts::value<std::string>(), "PREFIX")(
          "hts-threads",
          "Number of threads for BAM/CRAM compression [default: min(num_threads, 4)]",
          cxxopts::value<int>(), "INT")(
          "read-group",
          "Read group ID, or 'auto' to generate from input filenames",
          cxxopts::value<std::string>(), "ID")(
          "write-index",
          "Write .bai/.crai index (requires --sort-bam for coordinate-sorted output)")(
          "sort-bam",
          "Enable coordinate sorting for BAM/CRAM output. Sets @HD SO:coordinate. Required for --write-index. Sort key: (tid,pos,flag,mtid,mpos,isize,read_id). Note: differs from samtools sort ordering.")(
          "sort-bam-ram",
          "Max RAM for sorting before spilling to disk [8G]",
          cxxopts::value<std::string>(), "SIZE");
  //("PAF", "Output mappings in PAF format (only for test)");
}

void AddDevelopmentOptions(cxxopts::Options &options) {
  options.add_options("Development options")("A,match-score", "Match score [1]",
                                             cxxopts::value<int>(), "INT")(
      "B,mismatch-penalty", "Mismatch penalty [4]", cxxopts::value<int>(),
      "INT")("O,gap-open-penalties", "Gap open penalty [6,6]",
             cxxopts::value<std::vector<int>>(), "INT[,INT]")(
      "E,gap-extension-penalties", "Gap extension penalty [1,1]",
      cxxopts::value<std::vector<int>>(),
      "INT[,INT]")("n,max-num-best-mappings", "Only report n best mappings [1]",
                   cxxopts::value<int>(),
                   "INT")("multi-mapping-allocation-distance",
                          "Uni-mappings within this distance from any end of "
                          "multi-mappings are used for allocation [0]",
                          cxxopts::value<int>(), "INT")(
      "multi-mapping-allocation-seed",
      "Seed for random number generator in multi-mapping allocation [11]",
      cxxopts::value<int>(), "INT")(
      "drop-repetitive-reads",
      "Drop reads with too many best mappings [500000]", cxxopts::value<int>(),
      "INT")("allocate-multi-mappings", "Allocate multi-mappings")(
      "PAF", "Output mappings in PAF format (only for test)")(
      "skip-barcode-check",
      "Do not check whether too few barcodes are in the whitelist")
      ("cache-size", "number of cache entries [4000003]", cxxopts::value<int>(), "INT")
      ("cache-update-param", "value used to control number of reads sampled [0.01]", cxxopts::value<double>(), "FLT")
      ("debug-cache", "verbose output for debugging cache used in chromap")
      ("k-for-minhash", "number of values stored in each MinHash sketch [250]", cxxopts::value<int>(), "INT");
}

void AddPeakOptions(cxxopts::Options &options) {
  options.add_options("Peak")("cell-by-bin", "Generate cell-by-bin matrix")(
      "bin-size", "Bin size to generate cell-by-bin matrix [5000]",
      cxxopts::value<int>(),
      "INT")("depth-cutoff", "Depth cutoff for peak calling [3]",
             cxxopts::value<int>(),
             "INT")("peak-min-length", "Min length of peaks to report [30]",
                    cxxopts::value<int>(), "INT")(
      "peak-merge-max-length", "Peaks within this length will be merged [30]",
      cxxopts::value<int>(), "INT");
}

// Parse size string like "8G", "512M", "1024K" to bytes
uint64_t ParseSizeString(const std::string& sizeStr, const std::string& option_name) {
  if (sizeStr.empty()) {
    chromap::ExitWithMessage("Empty size string for " + option_name);
  }
  
  size_t endPos = sizeStr.length() - 1;
  char unit = std::toupper(sizeStr[endPos]);
  std::string numStr = sizeStr.substr(0, endPos);
  
  uint64_t multiplier = 1;
  if (unit == 'K') {
    multiplier = 1024ULL;
  } else if (unit == 'M') {
    multiplier = 1024ULL * 1024;
  } else if (unit == 'G') {
    multiplier = 1024ULL * 1024 * 1024;
  } else if (unit == 'T') {
    multiplier = 1024ULL * 1024 * 1024 * 1024;
  } else {
    // No unit, assume bytes
    numStr = sizeStr;
    multiplier = 1;
  }
  
  try {
    uint64_t num = std::stoull(numStr);
    return num * multiplier;
  } catch (const std::exception& e) {
    chromap::ExitWithMessage("Invalid size string for " + option_name + ": " + sizeStr);
  }
  return 0;  // Never reached
}

// Return all file paths that match the input pattern.
std::vector<std::string> GetMatchedFilePaths(const std::string &pattern) {
  glob_t glob_result;
  memset(&glob_result, 0, sizeof(glob_result));

  const int return_value =
      glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);

  if (return_value != 0) {
    globfree(&glob_result);
    chromap::ExitWithMessage("glob() failed with return value " +
                             std::to_string(return_value) + "\n");
  }

  std::vector<std::string> matched_file_paths;
  matched_file_paths.reserve(glob_result.gl_pathc);
  for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
    matched_file_paths.push_back(std::string(glob_result.gl_pathv[i]));
    std::cerr << matched_file_paths.back() << "\n";
  }
  globfree(&glob_result);

  return matched_file_paths;
}

// Return all file paths that match the input patterns.
std::vector<std::string> GetMatchedFilePaths(
    const std::vector<std::string> &patterns) {
  std::vector<std::string> all_matched_file_paths;
  for (const auto &pattern : patterns) {
    std::vector<std::string> matched_file_paths = GetMatchedFilePaths(pattern);
    all_matched_file_paths.reserve(all_matched_file_paths.size() +
                                   matched_file_paths.size());
    all_matched_file_paths.insert(
        std::end(all_matched_file_paths),
        std::make_move_iterator(std::begin(matched_file_paths)),
        std::make_move_iterator(std::end(matched_file_paths)));
  }
  return all_matched_file_paths;
}

// Derive secondary output path from primary, handling /dev/stdout and compound extensions
std::string DeriveSecondaryOutputPath(const std::string &primary_path,
                                       const std::string &suffix) {
  // Handle special device paths and stdout indicator
  if (primary_path == "/dev/stdout" || primary_path == "/dev/stderr" || primary_path == "-") {
    // For stdout, cannot derive secondary paths - user must specify explicitly
    // Return a default name but this should be caught by validation
    return "chromap_output" + suffix + ".sam";
  }
  
  size_t slash_pos = primary_path.rfind('/');
  size_t search_start = (slash_pos == std::string::npos) ? 0 : slash_pos + 1;
  
  // Check for compound extensions: .sam.gz, .bam
  std::string lower_path = primary_path;
  std::transform(lower_path.begin(), lower_path.end(), lower_path.begin(), ::tolower);
  
  // Handle .sam.gz -> output.noY.sam.gz
  if (lower_path.length() > 7 && lower_path.substr(lower_path.length() - 7) == ".sam.gz") {
    return primary_path.substr(0, primary_path.length() - 7) + suffix + ".sam.gz";
  }
  // Handle .bam -> output.noY.bam (preserve extension)
  if (lower_path.length() > 4 && lower_path.substr(lower_path.length() - 4) == ".bam") {
    return primary_path.substr(0, primary_path.length() - 4) + suffix + ".bam";
  }
  // Handle .cram -> output.noY.cram (preserve extension)
  if (lower_path.length() > 5 && lower_path.substr(lower_path.length() - 5) == ".cram") {
    return primary_path.substr(0, primary_path.length() - 5) + suffix + ".cram";
  }
  // Handle .sam -> output.noY.sam
  if (lower_path.length() > 4 && lower_path.substr(lower_path.length() - 4) == ".sam") {
    return primary_path.substr(0, primary_path.length() - 4) + suffix + ".sam";
  }
  
  // Find last dot for other extensions
  size_t dot_pos = primary_path.rfind('.');
  if (dot_pos != std::string::npos && dot_pos > search_start) {
    return primary_path.substr(0, dot_pos) + suffix + primary_path.substr(dot_pos);
  }
  
  // No extension: output -> output.noY.sam
  return primary_path + suffix + ".sam";
}

// Extract directory from a file path
std::string GetDirectoryFromPath(const std::string &path) {
  size_t slash_pos = path.rfind('/');
  if (slash_pos == std::string::npos) {
    return ".";
  }
  return path.substr(0, slash_pos);
}

// Determine base read extension (.fastq/.fq/.fasta/.fa/.fna) from input filename
std::string GetFastqBaseExtension(const std::string &input_path) {
  size_t slash_pos = input_path.rfind('/');
  std::string filename =
      (slash_pos == std::string::npos) ? input_path : input_path.substr(slash_pos + 1);
  std::string lower = filename;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

  if (lower.size() >= 9 && lower.substr(lower.size() - 9) == ".fastq.gz") {
    return ".fastq";
  }
  if (lower.size() >= 9 && lower.substr(lower.size() - 9) == ".fasta.gz") {
    return ".fasta";
  }
  if (lower.size() >= 7 && lower.substr(lower.size() - 7) == ".fna.gz") {
    return ".fna";
  }
  if (lower.size() >= 6 && lower.substr(lower.size() - 6) == ".fq.gz") {
    return ".fq";
  }
  if (lower.size() >= 6 && lower.substr(lower.size() - 6) == ".fa.gz") {
    return ".fa";
  }
  if (lower.size() >= 6 && lower.substr(lower.size() - 6) == ".fastq") {
    return ".fastq";
  }
  if (lower.size() >= 6 && lower.substr(lower.size() - 6) == ".fasta") {
    return ".fasta";
  }
  if (lower.size() >= 4 && lower.substr(lower.size() - 4) == ".fna") {
    return ".fna";
  }
  if (lower.size() >= 3 && lower.substr(lower.size() - 3) == ".fq") {
    return ".fq";
  }
  if (lower.size() >= 3 && lower.substr(lower.size() - 3) == ".fa") {
    return ".fa";
  }
  if (lower.size() >= 3 && lower.substr(lower.size() - 3) == ".gz") {
    return ".fastq";
  }
  return ".fastq";
}

std::string BuildFastqExtension(const std::string &base_ext,
                                const std::string &compression) {
  if (compression == "gz") {
    return base_ext + ".gz";
  }
  return base_ext;
}

// Derive FASTQ output path from input filename
// Inserts _Y or _noY before the last _R[0-9]+ token
// Preserves base read extension; applies compression setting
// file_index: 1-based index for multiple input files (0 = no index)
// mate_index: 1 or 2 (used for fallback naming)
std::string DeriveFastqOutputPath(const std::string &input_path,
                                   const std::string &output_dir,
                                   const std::string &suffix,  // "_Y" or "_noY"
                                   const std::string &compression,
                                   int file_index,
                                   int mate_index) {
  size_t slash_pos = input_path.rfind('/');
  std::string filename =
      (slash_pos == std::string::npos) ? input_path : input_path.substr(slash_pos + 1);

  std::string base_ext = GetFastqBaseExtension(filename);
  std::string output_ext = BuildFastqExtension(base_ext, compression);

  // Strip extension (including optional .gz) to get stem
  std::string lower = filename;
  std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
  std::string stem = filename;
  if (lower.size() >= 3 && lower.substr(lower.size() - 3) == ".gz") {
    stem = filename.substr(0, filename.size() - 3);
    lower = lower.substr(0, lower.size() - 3);
  }
  std::string base_ext_lower = base_ext;
  std::transform(base_ext_lower.begin(), base_ext_lower.end(),
                 base_ext_lower.begin(), ::tolower);
  if (lower.size() >= base_ext_lower.size() &&
      lower.substr(lower.size() - base_ext_lower.size()) == base_ext_lower) {
    stem = stem.substr(0, stem.size() - base_ext.size());
  }

  std::string base_name;
  std::regex r_token(R"(_R(\d+))");
  std::sregex_iterator it(stem.begin(), stem.end(), r_token);
  std::sregex_iterator end;
  if (it != end) {
    std::sregex_iterator last = it;
    for (; it != end; ++it) {
      last = it;
    }
    size_t pos = last->position(0);
    base_name = stem;
    base_name.insert(pos, suffix);
  } else {
    base_name = (suffix == "_Y" ? "Y_reads" : "noY_reads") + std::string(".mate") +
                std::to_string(mate_index);
  }

  std::string index_suffix = (file_index > 0) ? (".f" + std::to_string(file_index)) : "";

  std::string result = output_dir;
  if (!result.empty() && result.back() != '/') {
    result += "/";
  }
  result += base_name + index_suffix + output_ext;

  return result;
}

}  // namespace

void ChromapDriver::ParseArgsAndRun(int argc, char *argv[]) {
  cxxopts::Options options(
      "chromap", "Fast alignment and preprocessing of chromatin profiles");

  options.add_options()("v,version", "Print version")("h,help", "Print help");

  AddIndexingOptions(options);
  AddMappingOptions(options);

  // We don't support peak options for now.
  // AddPeakOptions(options);

  AddInputOptions(options);
  AddOutputOptions(options);

  AddDevelopmentOptions(options);

  auto result = options.parse(argc, argv);
  if (result.count("h")) {
    std::cerr << options.help(
        {"", "Indexing", "Mapping", "Peak", "Input", "Output"});
    return;
  }
  if (result.count("v")) {
    std::cerr << CHROMAP_VERSION << "\n";
    return;
  }
  // Parameters and their default
  IndexParameters index_parameters;
  MappingParameters mapping_parameters;

  if (result.count("preset")) {
    std::string read_type = result["preset"].as<std::string>();
    if (read_type == "atac") {
      std::cerr << "Preset parameters for ATAC-seq/scATAC-seq are used.\n";
      mapping_parameters.max_insert_size = 2000;
      mapping_parameters.trim_adapters = true;
      mapping_parameters.remove_pcr_duplicates = true;
      mapping_parameters.remove_pcr_duplicates_at_bulk_level = false;
      mapping_parameters.Tn5_shift = true;
      mapping_parameters.mapping_output_format = MAPPINGFORMAT_BED;
      mapping_parameters.low_memory_mode = true;
    } else if (read_type == "chip") {
      std::cerr << "Preset parameters for ChIP-seq are used.\n";
      mapping_parameters.max_insert_size = 2000;
      mapping_parameters.remove_pcr_duplicates = true;
      mapping_parameters.low_memory_mode = true;
      mapping_parameters.mapping_output_format = MAPPINGFORMAT_BED;
    } else if (read_type == "hic") {
      std::cerr << "Preset parameters for Hi-C are used.\n";
      mapping_parameters.error_threshold = 4;
      mapping_parameters.mapq_threshold = 1;
      mapping_parameters.split_alignment = true;
      mapping_parameters.low_memory_mode = true;
      mapping_parameters.mapping_output_format = MAPPINGFORMAT_PAIRS;
    } else {
      chromap::ExitWithMessage("Unrecognized preset parameters " + read_type +
                               "\n");
    }
  }
  // Optional parameters
  if (result.count("min-frag-length")) {
    int min_fragment_length = result["min-frag-length"].as<int>();
    if (min_fragment_length <= 60) {
      index_parameters.kmer_size = 17;
      index_parameters.window_size = 7;
    } else if (min_fragment_length <= 80) {
      index_parameters.kmer_size = 19;
      index_parameters.window_size = 10;
    } else {
      index_parameters.kmer_size = 23;
      index_parameters.window_size = 11;
    }
  }
  if (result.count("k")) {
    index_parameters.kmer_size = result["kmer"].as<int>();
  }
  if (result.count("w")) {
    index_parameters.window_size = result["window"].as<int>();
  }
  if (result.count("e")) {
    mapping_parameters.error_threshold = result["error-threshold"].as<int>();
  }
  if (result.count("A")) {
    mapping_parameters.match_score = result["match-score"].as<int>();
  }
  if (result.count("B")) {
    mapping_parameters.mismatch_penalty = result["mismatch-penalty"].as<int>();
  }
  if (result.count("O")) {
    mapping_parameters.gap_open_penalties =
        result["gap-open-penalties"].as<std::vector<int>>();
  }
  if (result.count("E")) {
    mapping_parameters.gap_extension_penalties =
        result["gap-extension-penalties"].as<std::vector<int>>();
  }
  if (result.count("s")) {
    mapping_parameters.min_num_seeds_required_for_mapping =
        result["min-num-seeds"].as<int>();
  }
  if (result.count("f")) {
    mapping_parameters.max_seed_frequencies =
        result["max-seed-frequencies"].as<std::vector<int>>();
  }
  if (result.count("n")) {
    mapping_parameters.max_num_best_mappings =
        result["max-num-best-mappings"].as<int>();
  }
  if (result.count("l")) {
    mapping_parameters.max_insert_size = result["max-insert-size"].as<int>();
  }
  if (result.count("q")) {
    mapping_parameters.mapq_threshold = result["MAPQ-threshold"].as<uint8_t>();
  }
  if (result.count("t")) {
    mapping_parameters.num_threads = result["num-threads"].as<int>();
  }


  // check cache-related parameters
  if (result.count("cache-update-param")) {
    mapping_parameters.cache_update_param = result["cache-update-param"].as<double>();
    if (mapping_parameters.cache_update_param < 0.0 || mapping_parameters.cache_update_param > 1.0){
      chromap::ExitWithMessage("cache update param is not approriate, must be in this range (0, 1]");
    }
  } 
  if (result.count("cache-size")) {
    mapping_parameters.cache_size = result["cache-size"].as<int>();
    if (mapping_parameters.cache_size < 2000000 || mapping_parameters.cache_size > 15000000) {
        chromap::ExitWithMessage("cache size is not in appropriate range\n");
    }
  }
  if (result.count("debug-cache")) {
    mapping_parameters.debug_cache = true;
  }
  if (result.count("frip-est-params")) {
    mapping_parameters.frip_est_params = result["frip-est-params"].as<std::string>();
  }
  if (result.count("turn-off-num-uniq-cache-slots")) {
    mapping_parameters.output_num_uniq_cache_slots = false;
  } 
  if (result.count("temp-dir")) {
    mapping_parameters.temp_directory_path = result["temp-dir"].as<std::string>();
  }
  if (result.count("k-for-minhash")) {
    mapping_parameters.k_for_minhash = result["k-for-minhash"].as<int>();
    if (mapping_parameters.k_for_minhash < 1 || mapping_parameters.k_for_minhash >= 2000) {
      chromap::ExitWithMessage("Invalid paramter for size of MinHash sketch (--k-for-minhash)");
    }
  }


  if (result.count("min-read-length")) {
    mapping_parameters.min_read_length = result["min-read-length"].as<int>();
  }
  if (result.count("bc-error-threshold")) {
    mapping_parameters.barcode_correction_error_threshold =
        result["bc-error-threshold"].as<int>();
  }
  if (result.count("bc-probability-threshold")) {
    mapping_parameters.barcode_correction_probability_threshold =
        result["bc-probability-threshold"].as<double>();
  }
  if (result.count("multi-mapping-allocation-distance")) {
    mapping_parameters.multi_mapping_allocation_distance =
        result["multi-mapping-allocation-distance"].as<int>();
  }
  if (result.count("multi-mapping-allocation-seed")) {
    mapping_parameters.multi_mapping_allocation_seed =
        result["multi-mapping-allocation-seed"].as<int>();
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
  if (result.count("Tn5-shift")) {
    mapping_parameters.Tn5_shift = true;
  }
  if (result.count("split-alignment")) {
    mapping_parameters.split_alignment = true;
  }
  if (result.count("output-mappings-not-in-whitelist")) {
    mapping_parameters.output_mappings_not_in_whitelist = true;
  }
  if (result.count("BED")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_BED;
  }
  if (result.count("TagAlign")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_TAGALIGN;
  }
  if (result.count("PAF")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_PAF;
  }
  if (result.count("pairs")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_PAIRS;
  }
  if (result.count("SAM")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_SAM;
  }
  if (result.count("BAM")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_BAM;
  }
  if (result.count("CRAM")) {
    mapping_parameters.mapping_output_format = MAPPINGFORMAT_CRAM;
  }
  if (result.count("hts-threads")) {
    mapping_parameters.hts_threads = result["hts-threads"].as<int>();
    if (mapping_parameters.hts_threads < 0) {
      chromap::ExitWithMessage("--hts-threads must be >= 0");
    }
  }
  if (result.count("read-group")) {
    mapping_parameters.read_group_id = result["read-group"].as<std::string>();
  }
  if (result.count("write-index")) {
    mapping_parameters.write_index = true;
  }
  if (result.count("sort-bam")) {
    mapping_parameters.sort_bam = true;
  }
  if (result.count("sort-bam-ram")) {
    std::string sizeStr = result["sort-bam-ram"].as<std::string>();
    mapping_parameters.sort_bam_ram_limit = ParseSizeString(sizeStr, "--sort-bam-ram");
  }
  if (result.count("low-mem")) {
    mapping_parameters.low_memory_mode = true;
  }
  if (result.count("low-mem-ram")) {
    std::string sizeStr = result["low-mem-ram"].as<std::string>();
    mapping_parameters.low_mem_ram_limit = ParseSizeString(sizeStr, "--low-mem-ram");
  }
  if (result.count("cell-by-bin")) {
    mapping_parameters.cell_by_bin = true;
  }
  if (result.count("bin-size")) {
    mapping_parameters.bin_size = result["bin-size"].as<int>();
  }
  if (result.count("depth-cutoff")) {
    mapping_parameters.depth_cutoff_to_call_peak =
        result["depth-cutoff"].as<uint16_t>();
  }
  if (result.count("peak-min-length")) {
    mapping_parameters.peak_min_length = result["peak-min-length"].as<int>();
  }
  if (result.count("peak-merge-max-length")) {
    mapping_parameters.peak_merge_max_length =
        result["peak-merge-max-length"].as<int>();
  }

  std::cerr << std::setprecision(2) << std::fixed;
  if (result.count("i")) {
    if (result.count("r")) {
      index_parameters.reference_file_path = result["ref"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No reference specified!");
    }
    if (result.count("o")) {
      index_parameters.index_output_file_path =
          result["output"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No output file specified!");
    }
    std::cerr << "Build index for the reference.\n";
    std::cerr << "Kmer length: " << index_parameters.kmer_size
              << ", window size: " << index_parameters.window_size << "\n";
    std::cerr << "Reference file: " << index_parameters.reference_file_path
              << "\n";
    std::cerr << "Output file: " << index_parameters.index_output_file_path
              << "\n";
    chromap::Chromap chromap_for_indexing(index_parameters);
    chromap_for_indexing.ConstructIndex();
  } else if (result.count("1")) {
    std::cerr << "Start to map reads.\n";
    if (result.count("r")) {
      mapping_parameters.reference_file_path = result["ref"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No reference specified!");
    }
    if (result.count("o")) {
      mapping_parameters.mapping_output_file_path =
          result["output"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No output file specified!");
    }
    if (result.count("x")) {
      mapping_parameters.index_file_path = result["index"].as<std::string>();
    } else {
      chromap::ExitWithMessage("No index file specified!");
    }
    if (result.count("1")) {
      mapping_parameters.read_file1_paths =
          GetMatchedFilePaths(result["read1"].as<std::vector<std::string>>());
    } else {
      chromap::ExitWithMessage("No read file specified!");
    }
    if (result.count("2")) {
      mapping_parameters.read_file2_paths =
          GetMatchedFilePaths(result["read2"].as<std::vector<std::string>>());
    }

    if (result.count("b")) {
      mapping_parameters.is_bulk_data = false;
      mapping_parameters.barcode_file_paths =
          GetMatchedFilePaths(result["barcode"].as<std::vector<std::string>>());
      if (result.count("barcode-whitelist") == 0) {
        std::cerr << "WARNING: there are input barcode files but a barcode "
                     "whitelist file is missing!\n";
      }
    }

    if (result.count("barcode-whitelist")) {
      if (mapping_parameters.is_bulk_data) {
        chromap::ExitWithMessage(
            "No barcode file specified but the barcode whitelist file is "
            "given!");
      }
      mapping_parameters.barcode_whitelist_file_path =
          result["barcode-whitelist"].as<std::string>();
    }

    if (result.count("p")) {
      mapping_parameters.matrix_output_prefix =
          result["matrix-output-prefix"].as<std::string>();
      if (mapping_parameters.is_bulk_data) {
        chromap::ExitWithMessage(
            "No barcode file specified but asked to output matrix files!");
      }
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

    if (result.count("barcode-translate")) {
      mapping_parameters.barcode_translate_table_file_path =
          result["barcode-translate"].as<std::string>();
    }

    if (result.count("summary")) {
      mapping_parameters.summary_metadata_file_path =
          result["summary"].as<std::string>();
    }

    if (result.count("skip-barcode-check")) {
      mapping_parameters.skip_barcode_check = true;
    }

    // Y-chromosome stream filtering
    if (result.count("emit-noY-bam")) {
      mapping_parameters.emit_noY_stream = true;
    }
    if (result.count("emit-Y-bam")) {
      mapping_parameters.emit_Y_stream = true;
    }
    
    // Validate: Y-filtering requires SAM/BAM/CRAM mode
    if ((mapping_parameters.emit_noY_stream || mapping_parameters.emit_Y_stream) &&
        mapping_parameters.mapping_output_format != MAPPINGFORMAT_SAM &&
        mapping_parameters.mapping_output_format != MAPPINGFORMAT_BAM &&
        mapping_parameters.mapping_output_format != MAPPINGFORMAT_CRAM) {
      chromap::ExitWithMessage(
          "--emit-noY-bam and --emit-Y-bam require --SAM, --BAM, or --CRAM output format");
    }
    
    // Validate: Y-filtering with stdout requires explicit output paths
    if ((mapping_parameters.emit_noY_stream || mapping_parameters.emit_Y_stream) &&
        (mapping_parameters.mapping_output_file_path == "-" ||
         mapping_parameters.mapping_output_file_path == "/dev/stdout" ||
         mapping_parameters.mapping_output_file_path == "/dev/stderr")) {
      bool has_noY_path = result.count("noY-output");
      bool has_Y_path = result.count("Y-output");
      if (mapping_parameters.emit_noY_stream && !has_noY_path) {
        chromap::ExitWithMessage("--emit-noY-bam requires --noY-output when primary output is stdout");
      }
      if (mapping_parameters.emit_Y_stream && !has_Y_path) {
        chromap::ExitWithMessage("--emit-Y-bam requires --Y-output when primary output is stdout");
      }
    }
    
    // Validate: CRAM requires reference
    if (mapping_parameters.mapping_output_format == MAPPINGFORMAT_CRAM &&
        mapping_parameters.reference_file_path.empty()) {
      chromap::ExitWithMessage("--CRAM requires --ref/-r reference file");
    }
    
    // Validate: --sort-bam constraints
    if (mapping_parameters.sort_bam) {
      if (mapping_parameters.mapping_output_format != MAPPINGFORMAT_BAM &&
          mapping_parameters.mapping_output_format != MAPPINGFORMAT_CRAM) {
        chromap::ExitWithMessage("--sort-bam requires --BAM or --CRAM output format");
      }
    }
    
    // Validate: --write-index constraints
    if (mapping_parameters.write_index) {
      if (!mapping_parameters.sort_bam) {
        chromap::ExitWithMessage("--write-index requires --sort-bam for coordinate-sorted output");
      }
      if (mapping_parameters.mapping_output_file_path == "-" ||
          mapping_parameters.mapping_output_file_path == "/dev/stdout" ||
          mapping_parameters.mapping_output_file_path == "/dev/stderr") {
        chromap::ExitWithMessage("--write-index is incompatible with stdout output (-o - or -o /dev/stdout).");
      }
      if (mapping_parameters.mapping_output_format != MAPPINGFORMAT_BAM &&
          mapping_parameters.mapping_output_format != MAPPINGFORMAT_CRAM) {
        chromap::ExitWithMessage("--write-index only works with --BAM or --CRAM output format");
      }
    }
    
    // Derive or set explicit paths
    if (mapping_parameters.emit_noY_stream) {
      if (result.count("noY-output")) {
        mapping_parameters.noY_output_path = result["noY-output"].as<std::string>();
      } else {
        mapping_parameters.noY_output_path = DeriveSecondaryOutputPath(
            mapping_parameters.mapping_output_file_path, ".noY");
      }
      std::cerr << "noY output file: " << mapping_parameters.noY_output_path << "\n";
    }
    
    if (mapping_parameters.emit_Y_stream) {
      if (result.count("Y-output")) {
        mapping_parameters.Y_output_path = result["Y-output"].as<std::string>();
      } else {
        mapping_parameters.Y_output_path = DeriveSecondaryOutputPath(
            mapping_parameters.mapping_output_file_path, ".Y");
      }
      std::cerr << "Y-only output file: " << mapping_parameters.Y_output_path << "\n";
    }

    // Y read names emission
    if (result.count("emit-Y-read-names")) {
      mapping_parameters.emit_y_read_names = true;
    }
    if (mapping_parameters.emit_y_read_names) {
      if (result.count("Y-read-names-output")) {
        mapping_parameters.y_read_names_output_path = result["Y-read-names-output"].as<std::string>();
      } else {
        // Derive from output path: <output>.Y.names.txt
        std::string output_path = mapping_parameters.mapping_output_file_path;
        if (output_path == "-" || output_path == "/dev/stdout" || output_path == "/dev/stderr") {
          chromap::ExitWithMessage("--emit-Y-read-names requires --Y-read-names-output when primary output is stdout");
        }
        // Find last dot or use full path
        size_t dot_pos = output_path.rfind('.');
        if (dot_pos != std::string::npos) {
          mapping_parameters.y_read_names_output_path = output_path.substr(0, dot_pos) + ".Y.names.txt";
        } else {
          mapping_parameters.y_read_names_output_path = output_path + ".Y.names.txt";
        }
      }
      std::cerr << "Y read names output file: " << mapping_parameters.y_read_names_output_path << "\n";
    }

    // Y/noY FASTQ emission
    if (result.count("emit-Y-noY-fastq")) {
      mapping_parameters.emit_y_noy_fastq = true;
    }
    if (result.count("emit-Y-noY-fastq-compression")) {
      std::string compression = result["emit-Y-noY-fastq-compression"].as<std::string>();
      if (compression != "gz" && compression != "none") {
        chromap::ExitWithMessage("--emit-Y-noY-fastq-compression must be 'gz' or 'none'");
      }
      mapping_parameters.y_noy_fastq_compression = compression;
    }
    if (result.count("Y-fastq-output-prefix")) {
      mapping_parameters.y_fastq_output_prefix = result["Y-fastq-output-prefix"].as<std::string>();
    }
    if (result.count("noY-fastq-output-prefix")) {
      mapping_parameters.noy_fastq_output_prefix = result["noY-fastq-output-prefix"].as<std::string>();
    }
    
    // Validate: FASTQ emission with stdout requires explicit prefixes
    if (mapping_parameters.emit_y_noy_fastq &&
        (mapping_parameters.mapping_output_file_path == "-" ||
         mapping_parameters.mapping_output_file_path == "/dev/stdout" ||
         mapping_parameters.mapping_output_file_path == "/dev/stderr")) {
      if (mapping_parameters.y_fastq_output_prefix.empty() ||
          mapping_parameters.noy_fastq_output_prefix.empty()) {
        chromap::ExitWithMessage("--emit-Y-noY-fastq requires --Y-fastq-output-prefix and --noY-fastq-output-prefix when primary output is stdout");
      }
    }
    
    // Derive FASTQ output paths if FASTQ emission is enabled
    if (mapping_parameters.emit_y_noy_fastq) {
      std::string output_dir =
          GetDirectoryFromPath(mapping_parameters.mapping_output_file_path);
      const std::string compression = mapping_parameters.y_noy_fastq_compression;

      const bool is_paired = !mapping_parameters.read_file2_paths.empty();
      const size_t num_files = mapping_parameters.read_file1_paths.size();

      if (num_files > 1) {
        std::cerr << "WARNING: Multiple input files detected. FASTQ outputs will include file index suffixes.\n";
      }

      mapping_parameters.y_fastq_output_paths_per_file.resize(num_files);
      mapping_parameters.noy_fastq_output_paths_per_file.resize(num_files);

      for (size_t file_index = 0; file_index < num_files; ++file_index) {
        const int output_file_index = (num_files > 1) ? static_cast<int>(file_index + 1) : 0;

        if (is_paired) {
          mapping_parameters.y_fastq_output_paths_per_file[file_index].resize(2);
          mapping_parameters.noy_fastq_output_paths_per_file[file_index].resize(2);

          const std::string read1_ext =
              BuildFastqExtension(GetFastqBaseExtension(mapping_parameters.read_file1_paths[file_index]),
                                  compression);
          const std::string read2_ext =
              BuildFastqExtension(GetFastqBaseExtension(mapping_parameters.read_file2_paths[file_index]),
                                  compression);

          if (!mapping_parameters.y_fastq_output_prefix.empty()) {
            const std::string index_suffix =
                (output_file_index > 0) ? (".f" + std::to_string(output_file_index)) : "";
            mapping_parameters.y_fastq_output_paths_per_file[file_index][0] =
                mapping_parameters.y_fastq_output_prefix + "mate1" + index_suffix + read1_ext;
            mapping_parameters.y_fastq_output_paths_per_file[file_index][1] =
                mapping_parameters.y_fastq_output_prefix + "mate2" + index_suffix + read2_ext;
          } else {
            mapping_parameters.y_fastq_output_paths_per_file[file_index][0] =
                DeriveFastqOutputPath(mapping_parameters.read_file1_paths[file_index],
                                      output_dir, "_Y", compression, output_file_index, 1);
            mapping_parameters.y_fastq_output_paths_per_file[file_index][1] =
                DeriveFastqOutputPath(mapping_parameters.read_file2_paths[file_index],
                                      output_dir, "_Y", compression, output_file_index, 2);
          }

          if (!mapping_parameters.noy_fastq_output_prefix.empty()) {
            const std::string index_suffix =
                (output_file_index > 0) ? (".f" + std::to_string(output_file_index)) : "";
            mapping_parameters.noy_fastq_output_paths_per_file[file_index][0] =
                mapping_parameters.noy_fastq_output_prefix + "mate1" + index_suffix + read1_ext;
            mapping_parameters.noy_fastq_output_paths_per_file[file_index][1] =
                mapping_parameters.noy_fastq_output_prefix + "mate2" + index_suffix + read2_ext;
          } else {
            mapping_parameters.noy_fastq_output_paths_per_file[file_index][0] =
                DeriveFastqOutputPath(mapping_parameters.read_file1_paths[file_index],
                                      output_dir, "_noY", compression, output_file_index, 1);
            mapping_parameters.noy_fastq_output_paths_per_file[file_index][1] =
                DeriveFastqOutputPath(mapping_parameters.read_file2_paths[file_index],
                                      output_dir, "_noY", compression, output_file_index, 2);
          }
        } else {
          mapping_parameters.y_fastq_output_paths_per_file[file_index].resize(1);
          mapping_parameters.noy_fastq_output_paths_per_file[file_index].resize(1);

          const std::string read_ext =
              BuildFastqExtension(GetFastqBaseExtension(mapping_parameters.read_file1_paths[file_index]),
                                  compression);

          if (!mapping_parameters.y_fastq_output_prefix.empty()) {
            const std::string index_suffix =
                (output_file_index > 0) ? (".f" + std::to_string(output_file_index)) : "";
            mapping_parameters.y_fastq_output_paths_per_file[file_index][0] =
                mapping_parameters.y_fastq_output_prefix + "reads" + index_suffix + read_ext;
          } else {
            mapping_parameters.y_fastq_output_paths_per_file[file_index][0] =
                DeriveFastqOutputPath(mapping_parameters.read_file1_paths[file_index],
                                      output_dir, "_Y", compression, output_file_index, 1);
          }

          if (!mapping_parameters.noy_fastq_output_prefix.empty()) {
            const std::string index_suffix =
                (output_file_index > 0) ? (".f" + std::to_string(output_file_index)) : "";
            mapping_parameters.noy_fastq_output_paths_per_file[file_index][0] =
                mapping_parameters.noy_fastq_output_prefix + "reads" + index_suffix + read_ext;
          } else {
            mapping_parameters.noy_fastq_output_paths_per_file[file_index][0] =
                DeriveFastqOutputPath(mapping_parameters.read_file1_paths[file_index],
                                      output_dir, "_noY", compression, output_file_index, 1);
          }
        }
      }

      std::cerr << "Y FASTQ output files:\n";
      for (size_t file_index = 0; file_index < num_files; ++file_index) {
        const auto &paths = mapping_parameters.y_fastq_output_paths_per_file[file_index];
        const std::string file_label =
            (num_files > 1) ? ("  file" + std::to_string(file_index + 1) + ": ") : "  ";
        if (is_paired) {
          std::cerr << file_label << "mate1: " << paths[0] << "\n";
          std::cerr << file_label << "mate2: " << paths[1] << "\n";
        } else {
          std::cerr << file_label << paths[0] << "\n";
        }
      }
      std::cerr << "noY FASTQ output files:\n";
      for (size_t file_index = 0; file_index < num_files; ++file_index) {
        const auto &paths = mapping_parameters.noy_fastq_output_paths_per_file[file_index];
        const std::string file_label =
            (num_files > 1) ? ("  file" + std::to_string(file_index + 1) + ": ") : "  ";
        if (is_paired) {
          std::cerr << file_label << "mate1: " << paths[0] << "\n";
          std::cerr << file_label << "mate2: " << paths[1] << "\n";
        } else {
          std::cerr << file_label << paths[0] << "\n";
        }
      }
    }

    // std::cerr << "Parameters: error threshold: " << error_threshold << ",
    // match score: " << match_score << ", mismatch_penalty: " <<
    // mismatch_penalty << ", gap open penalties for deletions and insertions: "
    // << gap_open_penalties[0] << "," << gap_open_penalties[1] << ", gap
    // extension penalties for deletions and insertions: " <<
    // gap_extension_penalties[0] << "," << gap_extension_penalties[1] << ",
    // min-num-seeds: " << min_num_seeds_required_for_mapping << ",
    // max-seed-frequency: " << max_seed_frequencies[0] << "," <<
    // max_seed_frequencies[1] << ", max-num-best-mappings: " <<
    // max_num_best_mappings << ", max-insert-size: " << max_insert_size << ",
    // MAPQ-threshold: " << (int)mapq_threshold << ", min-read-length: " <<
    // min_read_length << ", multi-mapping-allocation-distance: " <<
    // multi_mapping_allocation_distance << ", multi-mapping-allocation-seed: "
    // << multi_mapping_allocation_seed << ", drop-repetitive-reads: " <<
    // drop_repetitive_reads << "\n";
    std::cerr << "Parameters: error threshold: "
              << mapping_parameters.error_threshold << ", min-num-seeds: "
              << mapping_parameters.min_num_seeds_required_for_mapping
              << ", max-seed-frequency: "
              << mapping_parameters.max_seed_frequencies[0] << ","
              << mapping_parameters.max_seed_frequencies[1]
              << ", max-num-best-mappings: "
              << mapping_parameters.max_num_best_mappings
              << ", max-insert-size: " << mapping_parameters.max_insert_size
              << ", MAPQ-threshold: " << (int)mapping_parameters.mapq_threshold
              << ", min-read-length: " << mapping_parameters.min_read_length
              << ", bc-error-threshold: "
              << mapping_parameters.barcode_correction_error_threshold
              << ", bc-probability-threshold: "
              << mapping_parameters.barcode_correction_probability_threshold
              << "\n";
    std::cerr << "Number of threads: " << mapping_parameters.num_threads
              << "\n";
    if (mapping_parameters.is_bulk_data) {
      std::cerr << "Analyze bulk data.\n";
    } else {
      std::cerr << "Analyze single-cell data.\n";
    }
    if (mapping_parameters.trim_adapters) {
      std::cerr << "Will try to remove adapters on 3'.\n";
    } else {
      std::cerr << "Won't try to remove adapters on 3'.\n";
    }
    if (mapping_parameters.remove_pcr_duplicates) {
      std::cerr << "Will remove PCR duplicates after mapping.\n";
    } else {
      std::cerr << "Won't remove PCR duplicates after mapping.\n";
    }
    if (mapping_parameters.remove_pcr_duplicates_at_bulk_level) {
      std::cerr << "Will remove PCR duplicates at bulk level.\n";
    } else {
      std::cerr << "Will remove PCR duplicates at cell level.\n";
    }
    if (mapping_parameters.allocate_multi_mappings) {
      std::cerr << "Will allocate multi-mappings after mapping.\n";
    } else {
      std::cerr << "Won't allocate multi-mappings after mapping.\n";
    }
    if (mapping_parameters.only_output_unique_mappings) {
      std::cerr << "Only output unique mappings after mapping.\n";
    }
    if (!mapping_parameters.output_mappings_not_in_whitelist) {
      std::cerr << "Only output mappings of which barcodes are in whitelist.\n";
    } else {
      std::cerr << "No filtering of mappings based on whether their barcodes "
                   "are in whitelist.\n";
    }
    // if (allocate_multi_mappings && only_output_unique_mappings) {
    //  std::cerr << "WARNING: you want to output unique mappings only but you
    //  ask to allocate multi-mappings! In this case, it won't allocate
    //  multi-mappings and will only output unique mappings.\n";
    //  allocate_multi_mappings = false;
    //}
    if (mapping_parameters.max_num_best_mappings >
        mapping_parameters.drop_repetitive_reads) {
      std::cerr << "WARNING: you want to drop mapped reads with more than "
                << mapping_parameters.drop_repetitive_reads
                << " mappings. But you want to output top "
                << mapping_parameters.max_num_best_mappings
                << " best mappings. In this case, only reads with <="
                << mapping_parameters.drop_repetitive_reads
                << " best mappings will be output.\n";
      mapping_parameters.max_num_best_mappings =
          mapping_parameters.drop_repetitive_reads;
    }
    if (mapping_parameters.Tn5_shift) {
      std::cerr << "Perform Tn5 shift.\n";
    }
    if (mapping_parameters.split_alignment) {
      std::cerr << "Allow split alignment.\n";
    }

    switch (mapping_parameters.mapping_output_format) {
      case MAPPINGFORMAT_BED:
        std::cerr << "Output mappings in BED/BEDPE format.\n";
        break;
      case MAPPINGFORMAT_TAGALIGN:
        std::cerr << "Output mappings in TagAlign/PairedTagAlign format.\n";
        break;
      case MAPPINGFORMAT_PAF:
        std::cerr << "Output mappings in PAF format.\n";
        break;
      case MAPPINGFORMAT_SAM:
        std::cerr << "Output mappings in SAM format.\n";
        break;
      case MAPPINGFORMAT_BAM:
        std::cerr << "Output mappings in BAM format.\n";
        break;
      case MAPPINGFORMAT_CRAM:
        std::cerr << "Output mappings in CRAM format.\n";
        break;
      case MAPPINGFORMAT_PAIRS:
        std::cerr << "Output mappings in pairs format.\n";
        break;
      default:
        chromap::ExitWithMessage("Unknown mapping output format!");
        break;
    }

    std::cerr << "Reference file: " << mapping_parameters.reference_file_path
              << "\n";
    std::cerr << "Index file: " << mapping_parameters.index_file_path << "\n";
    for (size_t i = 0; i < mapping_parameters.read_file1_paths.size(); ++i) {
      std::cerr << i + 1
                << "th read 1 file: " << mapping_parameters.read_file1_paths[i]
                << "\n";
    }
    if (result.count("2") != 0) {
      for (size_t i = 0; i < mapping_parameters.read_file2_paths.size(); ++i) {
        std::cerr << i + 1 << "th read 2 file: "
                  << mapping_parameters.read_file2_paths[i] << "\n";
      }
    }
    if (result.count("b") != 0) {
      for (size_t i = 0; i < mapping_parameters.barcode_file_paths.size();
           ++i) {
        std::cerr << i + 1 << "th cell barcode file: "
                  << mapping_parameters.barcode_file_paths[i] << "\n";
      }
    }
    if (result.count("barcode-whitelist") != 0) {
      std::cerr << "Cell barcode whitelist file: "
                << mapping_parameters.barcode_whitelist_file_path << "\n";
    }
    std::cerr << "Output file: " << mapping_parameters.mapping_output_file_path
              << "\n";
    if (result.count("matrix-output-prefix") != 0) {
      std::cerr << "Matrix output prefix: "
                << mapping_parameters.matrix_output_prefix << "\n";
    }

    chromap::Chromap chromap_for_mapping(mapping_parameters);

    if (result.count("2") == 0) {
      // Single-end reads.
      switch (mapping_parameters.mapping_output_format) {
        case MAPPINGFORMAT_PAF: {
          chromap_for_mapping.MapSingleEndReads<chromap::PAFMapping>();
          break;
        }
        case MAPPINGFORMAT_SAM:
        case MAPPINGFORMAT_BAM:
        case MAPPINGFORMAT_CRAM: {
          chromap_for_mapping.MapSingleEndReads<chromap::SAMMapping>();
          break;
        }
        case MAPPINGFORMAT_PAIRS:
          chromap::ExitWithMessage("No support for single-end HiC yet!");
          break;
        case MAPPINGFORMAT_BED:
        case MAPPINGFORMAT_TAGALIGN:
          if (result.count("b") != 0) {
            chromap_for_mapping
                .MapSingleEndReads<chromap::MappingWithBarcode>();
          } else {
            chromap_for_mapping
                .MapSingleEndReads<chromap::MappingWithoutBarcode>();
          }
          break;
        default:
          chromap::ExitWithMessage("Unknown mapping output format!");
          break;
      }
    } else {
      // Paired-end reads.
      switch (mapping_parameters.mapping_output_format) {
        case MAPPINGFORMAT_PAF: {
          chromap_for_mapping.MapPairedEndReads<chromap::PairedPAFMapping>();
          break;
        }
        case MAPPINGFORMAT_SAM:
        case MAPPINGFORMAT_BAM:
        case MAPPINGFORMAT_CRAM: {
          chromap_for_mapping.MapPairedEndReads<chromap::SAMMapping>();
          break;
        }
        case MAPPINGFORMAT_PAIRS: {
          chromap_for_mapping.MapPairedEndReads<chromap::PairsMapping>();
          break;
        }
        case MAPPINGFORMAT_BED:
        case MAPPINGFORMAT_TAGALIGN:
          if (result.count("b") != 0) {
            chromap_for_mapping
                .MapPairedEndReads<chromap::PairedEndMappingWithBarcode>();
          } else {
            chromap_for_mapping
                .MapPairedEndReads<chromap::PairedEndMappingWithoutBarcode>();
          }
          break;
        default:
          chromap::ExitWithMessage("Unknown mapping output format!");
          break;
      }
    }
  } else {
    std::cerr << options.help(
        {"", "Indexing", "Mapping", "Peak", "Input", "Output"});
  }
}

}  // namespace chromap

int main(int argc, char *argv[]) {
  chromap::ChromapDriver chromap_driver;
  chromap_driver.ParseArgsAndRun(argc, argv);
  return 0;
}
