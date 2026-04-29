#include "libchromap.h"

#include <exception>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atac_dual_mapping.h"
#include "chromap.h"
#include "libmacs3/fragment_input.h"
#include "libmacs3/fragments.h"
#include "libmacs3/macs3_frag_peak_pipeline.h"
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
      if (mapping_parameters.AtacDualFragmentAndBam()) {
        chromap_for_mapping.MapPairedEndReads<PairedEndAtacDualMapping>();
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
              chromap_for_mapping
                  .MapPairedEndReads<PairedEndMappingWithBarcode>();
            } else {
              chromap_for_mapping.MapPairedEndReads<
                  PairedEndMappingWithoutBarcode>();
            }
            break;
          default:
            return MakeFailure(mapping_parameters,
                               "unknown mapping output format");
        }
      }
    }
  } catch (const std::exception &error) {
    return MakeFailure(mapping_parameters, error.what());
  } catch (...) {
    return MakeFailure(mapping_parameters, "unknown Chromap mapping failure");
  }

  return MakeSuccess(mapping_parameters);
}

namespace {

ChromapRunResult RunMacs3FragPeaksFromMappingParameters(
    MappingParameters &mapping_parameters) {
  const bool be_dual = mapping_parameters.AtacDualFragmentAndBam();
  const bool be_bed =
      (mapping_parameters.mapping_output_format == MAPPINGFORMAT_BED &&
       !mapping_parameters.read_file2_paths.empty() &&
       !mapping_parameters.barcode_file_paths.empty());
  if (!be_dual && !be_bed) {
    return MakeFailure(mapping_parameters,
                       "MACS3 FRAG peaks require BAM/CRAM dual + --atac-fragments "
                       "or paired-end barcoded BED output");
  }
  if (mapping_parameters.macs3_frag_peaks_narrowpeak_path.empty() ||
      mapping_parameters.macs3_frag_peaks_summits_path.empty()) {
    return MakeFailure(mapping_parameters,
                       "MACS3 FRAG peaks require narrowPeak and summits output paths");
  }

  peaks::Macs3FragPeakPipelineParams pr;
  pr.bdgpeakcall_cutoff =
      peaks::BdgPeakCallCutoffFromPValue(mapping_parameters.macs3_frag_pvalue);
  if (pr.bdgpeakcall_cutoff <= 0.f) {
    return MakeFailure(mapping_parameters,
                       "invalid --macs3-frag-pvalue for bdgpeakcall cutoff");
  }
  pr.min_length = mapping_parameters.macs3_frag_min_length;
  pr.max_gap = mapping_parameters.macs3_frag_max_gap;
  pr.macs3_uint8_counts = mapping_parameters.macs3_frag_uint8_counts;
  pr.peak_caller_threads = mapping_parameters.num_threads;

  const std::string &keep =
      mapping_parameters.macs3_frag_keep_intermediates_dir;
  const std::string &parent = mapping_parameters.temp_directory_path;
  std::string err;
  std::string work_used;

  if (mapping_parameters.macs3_frag_peaks_source ==
      Macs3FragPeaksSource::kMemory) {
    if (!mapping_parameters.macs3_frag_buffer ||
        !mapping_parameters.macs3_frag_chrom_names) {
      return MakeFailure(
          mapping_parameters,
          "MACS3 FRAG peaks (memory source): missing in-memory buffer or chrom_names");
    }
    auto &buckets = *mapping_parameters.macs3_frag_buffer;
    auto &chrom_names = *mapping_parameters.macs3_frag_chrom_names;
    if (mapping_parameters.macs3_frag_low_mem) {
      std::vector<macs3::FragmentRecord> flat;
      size_t total = 0;
      for (const auto &b : buckets) total += b.size();
      flat.reserve(total);
      for (auto &b : buckets) {
        for (auto &rec : b) flat.push_back(rec);
        std::vector<macs3::FragmentRecord>().swap(b);
      }
      std::vector<std::vector<macs3::FragmentRecord>>().swap(buckets);
      auto iter = macs3::WrapVectorFragmentIterator(
          std::move(flat), std::move(chrom_names));
      mapping_parameters.macs3_frag_buffer.reset();
      mapping_parameters.macs3_frag_chrom_names.reset();
      if (!peaks::RunMacs3FragPeakPipelineFromSortedIterator(
              *iter, pr, peaks::Macs3FragPeakPipelinePaths(),
              mapping_parameters.macs3_frag_peaks_narrowpeak_path,
              mapping_parameters.macs3_frag_peaks_summits_path, keep, parent,
              &work_used, &err)) {
        return MakeFailure(mapping_parameters, "MACS3 FRAG peaks: " + err);
      }
    } else {
      std::vector<peaks::ChromFragments> per_chrom;
      per_chrom.reserve(buckets.size());
      for (size_t i = 0; i < buckets.size(); ++i) {
        peaks::ChromFragments cf;
        cf.name = (i < chrom_names.size()) ? chrom_names[i] : std::string();
        cf.frags.reserve(buckets[i].size());
        for (const auto &rec : buckets[i]) {
          peaks::Fragment f;
          f.start = rec.start;
          f.end = rec.end;
          f.count = static_cast<int32_t>(rec.count);
          cf.frags.push_back(f);
        }
        std::vector<macs3::FragmentRecord>().swap(buckets[i]);
        per_chrom.push_back(std::move(cf));
      }
      std::vector<std::vector<macs3::FragmentRecord>>().swap(buckets);
      mapping_parameters.macs3_frag_buffer.reset();
      mapping_parameters.macs3_frag_chrom_names.reset();
      if (!peaks::RunMacs3FragPeakPipelineFromFragments(
              &per_chrom, pr, peaks::Macs3FragPeakPipelinePaths(),
              mapping_parameters.macs3_frag_peaks_narrowpeak_path,
              mapping_parameters.macs3_frag_peaks_summits_path, keep, parent,
              &work_used, &err, nullptr)) {
        return MakeFailure(mapping_parameters, "MACS3 FRAG peaks: " + err);
      }
    }
  } else {
    const std::string &fragments_path =
        be_dual ? mapping_parameters.atac_fragment_output_file_path
                : mapping_parameters.mapping_output_file_path;
    std::vector<peaks::ChromFragments> chs;
    if (!peaks::LoadFragmentsFromTsv(fragments_path, &chs)) {
      return MakeFailure(mapping_parameters,
                         "MACS3 FRAG peaks: failed to read fragments file " +
                             fragments_path);
    }
    if (!peaks::RunMacs3FragPeakPipelineFromFragments(
            &chs, pr, peaks::Macs3FragPeakPipelinePaths(),
            mapping_parameters.macs3_frag_peaks_narrowpeak_path,
            mapping_parameters.macs3_frag_peaks_summits_path, keep, parent,
            &work_used, &err, nullptr)) {
      return MakeFailure(mapping_parameters, "MACS3 FRAG peaks: " + err);
    }
  }
  return MakeSuccess(mapping_parameters);
}

}  // namespace

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

  // Local mutable copy so we can allocate the in-memory fragment buckets
  // before constructing Chromap (the constructor copies the parameters,
  // so the buckets must already be present on the input).
  MappingParameters params = mapping_parameters;
  if (params.call_macs3_frag_peaks &&
      params.macs3_frag_peaks_source == Macs3FragPeaksSource::kMemory) {
    params.macs3_frag_buffer =
        std::make_shared<std::vector<std::vector<macs3::FragmentRecord>>>();
    params.macs3_frag_chrom_names =
        std::make_shared<std::vector<std::string>>();
  }

  const ChromapRunResult mapping_result = RunMapping(params);
  if (!mapping_result.ok || !params.call_macs3_frag_peaks) {
    return mapping_result;
  }

  return RunMacs3FragPeaksFromMappingParameters(params);
}

}  // namespace chromap
