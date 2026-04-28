#include "libchromap.h"

#include <exception>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

#include "atac_dual_mapping.h"
#include "chromap.h"
#include "peak_caller/fragment_input.h"
#include "peak_caller/macs3_frag_peak_pipeline.h"
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

bool HasEmptyPath(const std::vector<std::string> &paths) {
  for (size_t i = 0; i < paths.size(); ++i) {
    if (paths[i].empty()) {
      return true;
    }
  }
  return false;
}

bool HasEmptyPath(const std::vector<std::vector<std::string>> &paths) {
  for (size_t i = 0; i < paths.size(); ++i) {
    if (HasEmptyPath(paths[i])) {
      return true;
    }
  }
  return false;
}

void WriteMacs3FragPeakSidecar(
    const MappingParameters &mapping_parameters, const std::string &work_used,
    const std::string &fragments_source,
    const std::string &memory_storage_mode_label) {
  if (mapping_parameters.summary_metadata_file_path.empty()) {
    return;
  }
  const std::string sidecar =
      mapping_parameters.summary_metadata_file_path + ".macs3_frag_peaks.tsv";
  FILE *fp = std::fopen(sidecar.c_str(), "w");
  if (fp == nullptr) {
    return;
  }
  std::fprintf(fp, "# Opt-in MACS3-compatible FRAG peaks (same C++ path as chromap_callpeaks).\n");
  std::fprintf(fp, "# Note: qValue and signalValue columns are not byte-for-byte MACS3 parity yet.\n");
  std::fprintf(fp, "key\tvalue\n");
  std::fprintf(fp, "fragments_file\t%s\n",
               mapping_parameters.atac_fragment_output_file_path.c_str());
  std::fprintf(fp, "fragments_source\t%s\n", fragments_source.c_str());
  if (fragments_source == "memory" && !memory_storage_mode_label.empty()) {
    std::fprintf(fp, "memory_storage_mode\t%s\n",
                 memory_storage_mode_label.c_str());
  }
  std::fprintf(fp, "narrowpeak_out\t%s\n",
               mapping_parameters.macs3_frag_peaks_narrowpeak_path.c_str());
  std::fprintf(fp, "summits_out\t%s\n",
               mapping_parameters.macs3_frag_peaks_summits_path.c_str());
  std::fprintf(fp, "macs3_frag_pvalue\t%.10g\n",
               mapping_parameters.macs3_frag_pvalue);
  std::fprintf(fp, "bdgpeakcall_cutoff_neg_log10_p\t%.10g\n",
               static_cast<double>(peaks::BdgPeakCallCutoffFromPValue(
                   mapping_parameters.macs3_frag_pvalue)));
  std::fprintf(fp, "macs3_frag_min_length\t%d\n",
               mapping_parameters.macs3_frag_min_length);
  std::fprintf(fp, "macs3_frag_max_gap\t%d\n",
               mapping_parameters.macs3_frag_max_gap);
  std::fprintf(fp, "macs3_uint8_counts\t%d\n",
               mapping_parameters.macs3_frag_uint8_counts ? 1 : 0);
  std::fprintf(fp, "effective_genome_size\t%lld\n",
               static_cast<long long>(2913022398LL));
  std::fprintf(fp, "llocal_bp\t%d\n", 10000);
  if (!mapping_parameters.macs3_frag_keep_intermediates_dir.empty()) {
    std::fprintf(fp, "keep_intermediates_dir\t%s\n",
                 mapping_parameters.macs3_frag_keep_intermediates_dir.c_str());
  } else if (!work_used.empty()) {
    std::fprintf(fp, "temp_work_dir_removed\t%s\n", work_used.c_str());
  }
  std::fclose(fp);
}

}  // namespace

ChromapRunResult ValidateMappingParameters(
    const MappingParameters &mapping_parameters) {
  if (mapping_parameters.reference_file_path.empty()) {
    return MakeFailure(mapping_parameters, "reference file is required");
  }
  if (mapping_parameters.index_file_path.empty()) {
    return MakeFailure(mapping_parameters, "index file is required");
  }
  if (mapping_parameters.mapping_output_file_path.empty()) {
    return MakeFailure(mapping_parameters, "mapping output path is required");
  }
  if (mapping_parameters.read_file1_paths.empty()) {
    return MakeFailure(mapping_parameters, "read 1 input is required");
  }
  if (HasEmptyPath(mapping_parameters.read_file1_paths)) {
    return MakeFailure(mapping_parameters, "read 1 input paths must be non-empty");
  }
  if (HasEmptyPath(mapping_parameters.read_file2_paths)) {
    return MakeFailure(mapping_parameters, "read 2 input paths must be non-empty");
  }
  if (HasEmptyPath(mapping_parameters.barcode_file_paths)) {
    return MakeFailure(mapping_parameters, "barcode input paths must be non-empty");
  }
  if (!mapping_parameters.read_file2_paths.empty() &&
      mapping_parameters.read_file1_paths.size() !=
          mapping_parameters.read_file2_paths.size()) {
    return MakeFailure(mapping_parameters,
                       "read 1 and read 2 input counts differ");
  }
  if (!mapping_parameters.barcode_file_paths.empty() &&
      mapping_parameters.read_file1_paths.size() !=
          mapping_parameters.barcode_file_paths.size()) {
    return MakeFailure(mapping_parameters,
                       "read 1 and barcode input counts differ");
  }
  if (mapping_parameters.mapping_output_format == MAPPINGFORMAT_UNKNOWN) {
    return MakeFailure(mapping_parameters, "mapping output format is required");
  }
  if (mapping_parameters.num_threads < 1) {
    return MakeFailure(mapping_parameters, "number of threads must be >= 1");
  }
  if (mapping_parameters.hts_threads < 0) {
    return MakeFailure(mapping_parameters, "hts threads must be >= 0");
  }
  if (mapping_parameters.max_seed_frequencies.size() != 2) {
    return MakeFailure(mapping_parameters,
                       "max seed frequencies must contain two values");
  }
  if (mapping_parameters.gap_open_penalties.size() != 2 ||
      mapping_parameters.gap_extension_penalties.size() != 2) {
    return MakeFailure(mapping_parameters,
                       "gap open and extension penalties must contain two values");
  }
  if (mapping_parameters.cache_update_param < 0.0 ||
      mapping_parameters.cache_update_param > 1.0) {
    return MakeFailure(mapping_parameters,
                       "cache update parameter must be in [0, 1]");
  }
  if (mapping_parameters.cache_size < 2000000 ||
      mapping_parameters.cache_size > 15000000) {
    return MakeFailure(mapping_parameters,
                       "cache size must be in [2000000, 15000000]");
  }
  if (mapping_parameters.k_for_minhash < 1 ||
      mapping_parameters.k_for_minhash >= 2000) {
    return MakeFailure(mapping_parameters,
                       "MinHash sketch size must be in [1, 2000)");
  }
  if (!mapping_parameters.barcode_whitelist_file_path.empty() &&
      mapping_parameters.barcode_file_paths.empty()) {
    return MakeFailure(mapping_parameters,
                       "barcode whitelist requires barcode reads");
  }
  if (!mapping_parameters.matrix_output_prefix.empty() &&
      mapping_parameters.barcode_file_paths.empty()) {
    return MakeFailure(mapping_parameters,
                       "matrix output requires barcode reads");
  }
  if (mapping_parameters.mapping_output_format == MAPPINGFORMAT_CRAM &&
      mapping_parameters.reference_file_path.empty()) {
    return MakeFailure(mapping_parameters,
                       "CRAM output requires a reference file");
  }
  if (mapping_parameters.emit_noY_stream || mapping_parameters.emit_Y_stream) {
    if (mapping_parameters.mapping_output_format != MAPPINGFORMAT_SAM &&
        mapping_parameters.mapping_output_format != MAPPINGFORMAT_BAM &&
        mapping_parameters.mapping_output_format != MAPPINGFORMAT_CRAM) {
      return MakeFailure(mapping_parameters,
                         "Y/noY stream output requires SAM, BAM, or CRAM output");
    }
  }
  if (mapping_parameters.emit_noY_stream &&
      mapping_parameters.noY_output_path.empty()) {
    return MakeFailure(mapping_parameters, "noY output path is required");
  }
  if (mapping_parameters.emit_Y_stream &&
      mapping_parameters.Y_output_path.empty()) {
    return MakeFailure(mapping_parameters, "Y output path is required");
  }
  if (mapping_parameters.emit_noY_stream && mapping_parameters.emit_Y_stream &&
      mapping_parameters.noY_output_path == mapping_parameters.Y_output_path) {
    return MakeFailure(mapping_parameters,
                       "Y and noY output paths must differ");
  }
  if (mapping_parameters.emit_y_read_names &&
      mapping_parameters.y_read_names_output_path.empty()) {
    return MakeFailure(mapping_parameters,
                       "Y read names output path is required");
  }
  if (mapping_parameters.emit_y_noy_fastq) {
    if (mapping_parameters.y_noy_fastq_compression != "gz" &&
        mapping_parameters.y_noy_fastq_compression != "none") {
      return MakeFailure(mapping_parameters,
                         "Y/noY FASTQ compression must be 'gz' or 'none'");
    }
    const size_t n_files = mapping_parameters.read_file1_paths.size();
    const size_t n_mates = mapping_parameters.read_file2_paths.empty() ? 1 : 2;
    if (mapping_parameters.y_fastq_output_paths_per_file.size() != n_files ||
        mapping_parameters.noy_fastq_output_paths_per_file.size() != n_files) {
      return MakeFailure(mapping_parameters,
                         "Y/noY FASTQ paths must be provided for every input file");
    }
    for (size_t i = 0; i < n_files; ++i) {
      if (mapping_parameters.y_fastq_output_paths_per_file[i].size() != n_mates ||
          mapping_parameters.noy_fastq_output_paths_per_file[i].size() != n_mates) {
        return MakeFailure(mapping_parameters,
                           "Y/noY FASTQ paths must match input mate count");
      }
    }
    if (HasEmptyPath(mapping_parameters.y_fastq_output_paths_per_file) ||
        HasEmptyPath(mapping_parameters.noy_fastq_output_paths_per_file)) {
      return MakeFailure(mapping_parameters,
                         "Y/noY FASTQ output paths must be non-empty");
    }
  }
  if (mapping_parameters.sort_bam &&
      mapping_parameters.mapping_output_format != MAPPINGFORMAT_BAM &&
      mapping_parameters.mapping_output_format != MAPPINGFORMAT_CRAM) {
    return MakeFailure(mapping_parameters,
                       "coordinate sorting requires BAM or CRAM output");
  }
  if (mapping_parameters.write_index) {
    if (!mapping_parameters.sort_bam) {
      return MakeFailure(mapping_parameters,
                         "writing an index requires coordinate sorting");
    }
    if (mapping_parameters.mapping_output_file_path == "-" ||
        mapping_parameters.mapping_output_file_path == "/dev/stdout" ||
        mapping_parameters.mapping_output_file_path == "/dev/stderr") {
      return MakeFailure(mapping_parameters,
                         "writing an index is incompatible with stdout output");
    }
    if (mapping_parameters.mapping_output_format != MAPPINGFORMAT_BAM &&
        mapping_parameters.mapping_output_format != MAPPINGFORMAT_CRAM) {
      return MakeFailure(mapping_parameters,
                         "writing an index requires BAM or CRAM output");
    }
  }
  if (!mapping_parameters.atac_fragment_output_file_path.empty()) {
    if (mapping_parameters.read_file2_paths.empty()) {
      return MakeFailure(mapping_parameters,
                         "ATAC fragments require paired-end reads");
    }
    if (mapping_parameters.barcode_file_paths.empty()) {
      return MakeFailure(mapping_parameters,
                         "ATAC fragments require barcode reads");
    }
    if (mapping_parameters.mapping_output_format != MAPPINGFORMAT_BAM &&
        mapping_parameters.mapping_output_format != MAPPINGFORMAT_CRAM) {
      return MakeFailure(mapping_parameters,
                         "ATAC fragments require BAM or CRAM primary output");
    }
    if (mapping_parameters.low_memory_mode) {
      return MakeFailure(mapping_parameters,
                         "ATAC fragments are not supported with low-memory mode");
    }
    if (mapping_parameters.atac_fragment_output_file_path ==
        mapping_parameters.mapping_output_file_path) {
      return MakeFailure(mapping_parameters,
                         "ATAC fragments path must differ from mapping output");
    }
  }
  if (mapping_parameters.macs3_frag_peaks_source ==
          Macs3FragPeaksSource::kMemory &&
      !mapping_parameters.call_macs3_frag_peaks) {
    return MakeFailure(
        mapping_parameters,
        "MACS3 FRAG memory source requires peak calling to be enabled");
  }
  if (mapping_parameters.macs3_frag_compact_min_count_bits < 1 ||
      mapping_parameters.macs3_frag_compact_min_count_bits > 31) {
    return MakeFailure(mapping_parameters,
                       "MACS3 FRAG compact count bits must be in [1, 31]");
  }
  if (!mapping_parameters.call_macs3_frag_peaks) {
    return MakeSuccess(mapping_parameters);
  }
  if (!mapping_parameters.AtacDualFragmentAndBam()) {
    return MakeFailure(
        mapping_parameters,
        "--call-macs3-frag-peaks requires paired-end barcoded reads, BAM/CRAM, and --atac-fragments");
  }
  if (mapping_parameters.macs3_frag_peaks_narrowpeak_path.empty() ||
      mapping_parameters.macs3_frag_peaks_summits_path.empty()) {
    return MakeFailure(
        mapping_parameters,
        "--call-macs3-frag-peaks requires narrowPeak and summits output paths");
  }
  if (mapping_parameters.macs3_frag_pvalue <= 0.0 ||
      mapping_parameters.macs3_frag_pvalue > 1.0 ||
      peaks::BdgPeakCallCutoffFromPValue(
          mapping_parameters.macs3_frag_pvalue) <= 0.f) {
    return MakeFailure(mapping_parameters,
                       "invalid MACS3 FRAG p-value for peak calling");
  }
  if (mapping_parameters.macs3_frag_min_length < 1 ||
      mapping_parameters.macs3_frag_max_gap < 0) {
    return MakeFailure(mapping_parameters,
                       "invalid MACS3 FRAG min length or max gap");
  }
  return MakeSuccess(mapping_parameters);
}

namespace {

ChromapRunResult RunMacs3FragPeaksAfterMapping(
    MappingParameters *mapping_parameters) {
  if (mapping_parameters == nullptr ||
      !mapping_parameters->call_macs3_frag_peaks) {
    return MakeSuccess(mapping_parameters == nullptr ? MappingParameters()
                                                    : *mapping_parameters);
  }
  std::string fragments_source = "file";
  std::string mem_mode;
  std::vector<peaks::ChromFragments> chs;
  if (mapping_parameters->macs3_frag_peaks_source ==
      Macs3FragPeaksSource::kMemory) {
    fragments_source = "memory";
    mem_mode = "workspace_events";
    if (!mapping_parameters->macs3_frag_workspace) {
      return MakeFailure(*mapping_parameters,
                         "MACS3 FRAG peaks: missing memory workspace");
    }
  } else if (!peaks::LoadFragmentsFromTsv(
                 mapping_parameters->atac_fragment_output_file_path, &chs)) {
    return MakeFailure(
        *mapping_parameters,
        "MACS3 FRAG peaks: failed to read fragments file " +
            mapping_parameters->atac_fragment_output_file_path);
  }

  peaks::Macs3FragPeakPipelineParams pr;
  pr.bdgpeakcall_cutoff =
      peaks::BdgPeakCallCutoffFromPValue(mapping_parameters->macs3_frag_pvalue);
  pr.min_length = mapping_parameters->macs3_frag_min_length;
  pr.max_gap = mapping_parameters->macs3_frag_max_gap;
  pr.macs3_uint8_counts = mapping_parameters->macs3_frag_uint8_counts;
  pr.peak_caller_threads = mapping_parameters->num_threads;

  std::string err;
  std::string work_used;
  const std::string &keep =
      mapping_parameters->macs3_frag_keep_intermediates_dir;
  const std::string parent = mapping_parameters->temp_directory_path;
  if (mapping_parameters->macs3_frag_peaks_source ==
      Macs3FragPeaksSource::kMemory) {
    if (!peaks::RunMacs3FragPeakPipelineFromWorkspace(
            mapping_parameters->macs3_frag_workspace.get(), pr,
            peaks::Macs3FragPeakPipelinePaths(),
            mapping_parameters->macs3_frag_peaks_narrowpeak_path,
            mapping_parameters->macs3_frag_peaks_summits_path, keep, parent,
            &work_used, &err)) {
      return MakeFailure(*mapping_parameters, "MACS3 FRAG peaks: " + err);
    }
  } else {
    if (!peaks::RunMacs3FragPeakPipelineFromFragments(
            &chs, pr, peaks::Macs3FragPeakPipelinePaths(),
            mapping_parameters->macs3_frag_peaks_narrowpeak_path,
            mapping_parameters->macs3_frag_peaks_summits_path, keep, parent,
            &work_used, &err, nullptr)) {
      return MakeFailure(*mapping_parameters, "MACS3 FRAG peaks: " + err);
    }
  }
  WriteMacs3FragPeakSidecar(*mapping_parameters, work_used, fragments_source,
                            mem_mode);
  return MakeSuccess(*mapping_parameters);
}

}  // namespace

ChromapRunResult RunMapping(const MappingParameters &mapping_parameters) {
  MappingParameters params = mapping_parameters;
  const ChromapRunResult validation = ValidateMappingParameters(params);
  if (!validation.ok) {
    return validation;
  }
  if (params.call_macs3_frag_peaks &&
      params.macs3_frag_peaks_source == Macs3FragPeaksSource::kMemory &&
      !params.macs3_frag_workspace) {
    params.macs3_frag_workspace =
        std::make_shared<peaks::Macs3FragPeakWorkspace>();
  }
  try {
    Chromap chromap_for_mapping(params);

    if (params.read_file2_paths.empty()) {
      switch (params.mapping_output_format) {
        case MAPPINGFORMAT_PAF:
          chromap_for_mapping.MapSingleEndReads<PAFMapping>();
          break;
        case MAPPINGFORMAT_SAM:
        case MAPPINGFORMAT_BAM:
        case MAPPINGFORMAT_CRAM:
          chromap_for_mapping.MapSingleEndReads<SAMMapping>();
          break;
        case MAPPINGFORMAT_PAIRS:
          return MakeFailure(params,
                             "single-end PAIRS output is not supported");
        case MAPPINGFORMAT_BED:
        case MAPPINGFORMAT_TAGALIGN:
          if (!params.barcode_file_paths.empty()) {
            chromap_for_mapping.MapSingleEndReads<MappingWithBarcode>();
          } else {
            chromap_for_mapping.MapSingleEndReads<MappingWithoutBarcode>();
          }
          break;
        default:
          return MakeFailure(params, "unknown mapping output format");
      }
    } else {
      if (params.AtacDualFragmentAndBam()) {
        chromap_for_mapping.MapPairedEndReads<PairedEndAtacDualMapping>();
      } else {
        switch (params.mapping_output_format) {
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
            if (!params.barcode_file_paths.empty()) {
              chromap_for_mapping
                  .MapPairedEndReads<PairedEndMappingWithBarcode>();
            } else {
              chromap_for_mapping.MapPairedEndReads<
                  PairedEndMappingWithoutBarcode>();
            }
            break;
          default:
            return MakeFailure(params,
                               "unknown mapping output format");
        }
      }
    }
  } catch (const std::exception &error) {
    return MakeFailure(params, error.what());
  } catch (...) {
    return MakeFailure(params, "unknown Chromap mapping failure");
  }

  return RunMacs3FragPeaksAfterMapping(&params);
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
      mapping_parameters.mapping_output_format != MAPPINGFORMAT_TAGALIGN &&
      mapping_parameters.mapping_output_format != MAPPINGFORMAT_SAM &&
      mapping_parameters.mapping_output_format != MAPPINGFORMAT_BAM &&
      mapping_parameters.mapping_output_format != MAPPINGFORMAT_CRAM) {
    return MakeFailure(mapping_parameters,
                       "ATAC mapping requires BED, TagAlign, SAM, BAM, or CRAM output");
  }
  if (mapping_parameters.AtacDualFragmentAndBam() &&
      mapping_parameters.low_memory_mode) {
    return MakeFailure(mapping_parameters,
                       "dual ATAC BAM/fragments output is not supported with low-memory mode");
  }

  return RunMapping(mapping_parameters);
}

}  // namespace chromap
