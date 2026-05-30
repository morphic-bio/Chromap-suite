#include "y_noy_path_utils.h"

#include <sys/stat.h>
#include <sys/types.h>

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <cstring>
#include <string>

#include "utils.h"

namespace chromap {
namespace {

bool IsStdStreamPath(const std::string &path) {
  return path == "-" || path == "/dev/stdout" || path == "/dev/stderr";
}

std::string NormalizeDirectoryPath(std::string path) {
  while (path.size() > 1 && path[path.size() - 1] == '/') {
    path.resize(path.size() - 1);
  }
  return path;
}

bool IsDirectory(const std::string &path) {
  struct stat st;
  return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}

void EnsureDirectoryExists(const std::string &path) {
  if (path.empty()) {
    return;
  }
  const std::string normalized = NormalizeDirectoryPath(path);
  if (normalized.empty() || normalized == ".") {
    return;
  }
  if (IsDirectory(normalized)) {
    return;
  }

  const size_t slash_pos = normalized.rfind('/');
  if (slash_pos != std::string::npos && slash_pos > 0) {
    EnsureDirectoryExists(normalized.substr(0, slash_pos));
  }

  if (mkdir(normalized.c_str(), 0775) != 0 && errno != EEXIST) {
    ExitWithMessage("Failed to create Y/noY FASTQ output directory " +
                    normalized + ": " + std::strerror(errno));
  }
  if (!IsDirectory(normalized)) {
    ExitWithMessage("Y/noY FASTQ output path is not a directory: " +
                    normalized);
  }
}

std::string PathBasename(const std::string &path) {
  const size_t slash_pos = path.find_last_of("/\\");
  if (slash_pos == std::string::npos) {
    return path;
  }
  return path.substr(slash_pos + 1);
}

std::string OutputDirectoryFromPrimary(const std::string &primary_path) {
  if (IsStdStreamPath(primary_path)) {
    return "y_separated";
  }
  const size_t slash_pos = primary_path.find_last_of("/\\");
  if (slash_pos == std::string::npos) {
    return "y_separated";
  }
  return primary_path.substr(0, slash_pos + 1) + "y_separated";
}

std::string LowerCopy(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  return value;
}

bool EndsWith(const std::string &value, const std::string &suffix) {
  return value.size() >= suffix.size() &&
         value.compare(value.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::string GetFastxBaseExtension(const std::string &input_path) {
  const std::string lower = LowerCopy(PathBasename(input_path));
  if (EndsWith(lower, ".fastq.gz")) return ".fastq";
  if (EndsWith(lower, ".fasta.gz")) return ".fasta";
  if (EndsWith(lower, ".fna.gz")) return ".fna";
  if (EndsWith(lower, ".fq.gz")) return ".fq";
  if (EndsWith(lower, ".fa.gz")) return ".fa";
  if (EndsWith(lower, ".fastq")) return ".fastq";
  if (EndsWith(lower, ".fasta")) return ".fasta";
  if (EndsWith(lower, ".fna")) return ".fna";
  if (EndsWith(lower, ".fq")) return ".fq";
  if (EndsWith(lower, ".fa")) return ".fa";
  return ".fastq";
}

std::string BuildFastxOutputExtension(const std::string &base_ext,
                                      const std::string &compression) {
  return compression == "gz" ? base_ext + ".gz" : base_ext;
}

std::string StripFastxExtension(const std::string &filename) {
  std::string stem = filename;
  std::string lower = LowerCopy(filename);
  if (EndsWith(lower, ".gz")) {
    stem = stem.substr(0, stem.size() - 3);
    lower = lower.substr(0, lower.size() - 3);
  }
  const std::string base_ext = GetFastxBaseExtension(filename);
  const std::string base_ext_lower = LowerCopy(base_ext);
  if (EndsWith(lower, base_ext_lower)) {
    stem = stem.substr(0, stem.size() - base_ext.size());
  }
  return stem;
}

bool InsertBeforeLastReadToken(const std::string &stem,
                               const std::string &tag,
                               std::string *tagged) {
  size_t last_pos = std::string::npos;
  for (size_t i = 0; i + 2 < stem.size(); ++i) {
    if (stem[i] != '_' || (stem[i + 1] != 'R' && stem[i + 1] != 'r') ||
        !std::isdigit(static_cast<unsigned char>(stem[i + 2]))) {
      continue;
    }
    size_t j = i + 3;
    while (j < stem.size() &&
           std::isdigit(static_cast<unsigned char>(stem[j]))) {
      ++j;
    }
    last_pos = i;
    i = j - 1;
  }
  if (last_pos == std::string::npos) {
    return false;
  }
  *tagged = stem;
  tagged->insert(last_pos, tag);
  return true;
}

std::string JoinPath(const std::string &dir, const std::string &filename) {
  if (dir.empty()) {
    return filename;
  }
  if (dir[dir.size() - 1] == '/' || dir[dir.size() - 1] == '\\') {
    return dir + filename;
  }
  return dir + "/" + filename;
}

std::string DeriveFastqOutputPath(const std::string &input_path,
                                  const std::string &output_dir,
                                  const std::string &tag,
                                  const std::string &compression,
                                  int file_index,
                                  int mate_index) {
  const std::string filename = PathBasename(input_path);
  const std::string base_ext = GetFastxBaseExtension(filename);
  const std::string output_ext =
      BuildFastxOutputExtension(base_ext, compression);
  const std::string stem = StripFastxExtension(filename);

  std::string base_name;
  if (!InsertBeforeLastReadToken(stem, tag, &base_name)) {
    base_name = (tag == "_Y" ? "Y_reads" : "noY_reads") +
                std::string(".mate") + std::to_string(mate_index);
  }

  const std::string index_suffix =
      file_index > 0 ? ".f" + std::to_string(file_index) : "";
  return JoinPath(output_dir, base_name + index_suffix + output_ext);
}

const std::string &InputSourcePathForMate(
    const MappingParameters &mapping_parameters,
    size_t file_index,
    int mate_index) {
  if (mapping_parameters.UsesCbqInput()) {
    return mapping_parameters.read_pair_cbq_paths[file_index];
  }
  if (mate_index == 2) {
    return mapping_parameters.read_file2_paths[file_index];
  }
  return mapping_parameters.read_file1_paths[file_index];
}

std::string PrefixOutputPath(const std::string &prefix,
                             const std::string &label,
                             int output_file_index,
                             const std::string &extension) {
  const std::string index_suffix =
      output_file_index > 0 ? ".f" + std::to_string(output_file_index) : "";
  return prefix + label + index_suffix + extension;
}

}  // namespace

void InitializeYNoYFastqOutputPaths(
    MappingParameters *mapping_parameters) {
  if (mapping_parameters == nullptr ||
      !mapping_parameters->emit_y_noy_fastq) {
    return;
  }
  if (mapping_parameters->y_noy_fastq_compression != "gz" &&
      mapping_parameters->y_noy_fastq_compression != "none") {
    ExitWithMessage("--emit-Y-noY-fastq-compression must be 'gz' or 'none'");
  }

  const size_t num_files = mapping_parameters->NumInputLanes();
  if (num_files == 0) {
    ExitWithMessage("--emit-Y-noY-fastq requires read input");
  }

  std::string output_dir = mapping_parameters->y_noy_fastq_output_dir;
  if (output_dir.empty()) {
    output_dir =
        OutputDirectoryFromPrimary(mapping_parameters->mapping_output_file_path);
  }
  const bool create_default_dir =
      mapping_parameters->y_fastq_output_prefix.empty() ||
      mapping_parameters->noy_fastq_output_prefix.empty();
  if (create_default_dir) {
    EnsureDirectoryExists(output_dir);
  }

  const bool is_paired = mapping_parameters->HasPairedEndInput();
  mapping_parameters->y_fastq_output_paths_per_file.assign(
      num_files, std::vector<std::string>());
  mapping_parameters->noy_fastq_output_paths_per_file.assign(
      num_files, std::vector<std::string>());

  for (size_t file_index = 0; file_index < num_files; ++file_index) {
    const int output_file_index =
        num_files > 1 ? static_cast<int>(file_index + 1) : 0;
    const size_t mate_count = is_paired ? 2 : 1;
    mapping_parameters->y_fastq_output_paths_per_file[file_index].resize(
        mate_count);
    mapping_parameters->noy_fastq_output_paths_per_file[file_index].resize(
        mate_count);

    for (size_t mate = 0; mate < mate_count; ++mate) {
      const int mate_index = static_cast<int>(mate + 1);
      const std::string &source_path =
          InputSourcePathForMate(*mapping_parameters, file_index, mate_index);
      const std::string output_ext = BuildFastxOutputExtension(
          GetFastxBaseExtension(source_path),
          mapping_parameters->y_noy_fastq_compression);

      if (!mapping_parameters->y_fastq_output_prefix.empty()) {
        mapping_parameters->y_fastq_output_paths_per_file[file_index][mate] =
            PrefixOutputPath(mapping_parameters->y_fastq_output_prefix,
                             is_paired ? "mate" + std::to_string(mate_index)
                                       : "reads",
                             output_file_index,
                             output_ext);
      } else {
        mapping_parameters->y_fastq_output_paths_per_file[file_index][mate] =
            DeriveFastqOutputPath(source_path, output_dir, "_Y",
                                  mapping_parameters->y_noy_fastq_compression,
                                  output_file_index, mate_index);
      }

      if (!mapping_parameters->noy_fastq_output_prefix.empty()) {
        mapping_parameters->noy_fastq_output_paths_per_file[file_index][mate] =
            PrefixOutputPath(mapping_parameters->noy_fastq_output_prefix,
                             is_paired ? "mate" + std::to_string(mate_index)
                                       : "reads",
                             output_file_index,
                             output_ext);
      } else {
        mapping_parameters->noy_fastq_output_paths_per_file[file_index][mate] =
            DeriveFastqOutputPath(source_path, output_dir, "_noY",
                                  mapping_parameters->y_noy_fastq_compression,
                                  output_file_index, mate_index);
      }
    }
  }
}

}  // namespace chromap
