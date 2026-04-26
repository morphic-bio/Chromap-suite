#ifndef MAPPING_WRITER_H_
#define MAPPING_WRITER_H_

#include <assert.h>

#include <cinttypes>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <unordered_set>

#include "barcode_translator.h"
#include "bed_mapping.h"
#include "atac_dual_mapping.h"
#include "mapping.h"
#include "mapping_parameters.h"
#include "paf_mapping.h"
#include "pairs_mapping.h"
#include "sam_mapping.h"
#include "sequence_batch.h"
#include "temp_mapping.h"
#include "utils.h"
#include "summary_metadata.h"

#ifndef LEGACY_OVERFLOW
#include "overflow_writer.h"
#include "overflow_reader.h"
#include <memory>
#include <mutex>
#include <vector>
#endif

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <zlib.h>
#include <memory>

// Include bam_sorter.h before namespace chromap to avoid namespace collision
// We need complete type for unique_ptr destructor
#include "bam_sorter.h"

namespace chromap {

template <typename MappingRecord>
class MappingWriter {
 public:
  MappingWriter() = delete;

  MappingWriter(const MappingParameters mapping_parameters,
                const uint32_t cell_barcode_length,
                const std::vector<int> &pairs_custom_rid_rank)
      : mapping_parameters_(mapping_parameters),
        cell_barcode_length_(cell_barcode_length),
        pairs_custom_rid_rank_(pairs_custom_rid_rank) {
    if (!mapping_parameters_.barcode_translate_table_file_path.empty()) {
      barcode_translator_.SetTranslateTable(
          mapping_parameters_.barcode_translate_table_file_path,
          mapping_parameters_.barcode_translate_from_first_column);
    }
    summary_metadata_.SetBarcodeLength(cell_barcode_length);
    if (mapping_parameters_.AtacDualFragmentAndBam()) {
      const std::string &af =
          mapping_parameters_.atac_fragment_output_file_path;
      if (af.size() >= 3 && af.compare(af.size() - 3, 3, ".gz") == 0) {
        atac_fragment_gz_output_ = gzopen(af.c_str(), "wb");
        if (!atac_fragment_gz_output_) {
          ExitWithMessage("Failed to open --atac-fragments output: " + af);
        }
      } else {
        atac_fragment_output_file_ = fopen(af.c_str(), "w");
        if (!atac_fragment_output_file_) {
          ExitWithMessage("Failed to open --atac-fragments output: " + af);
        }
      }
    } else if (mapping_parameters_.mapping_output_format != MAPPINGFORMAT_BAM &&
               mapping_parameters_.mapping_output_format != MAPPINGFORMAT_CRAM) {
      mapping_output_file_ =
          fopen(mapping_parameters_.mapping_output_file_path.c_str(), "w");
      assert(mapping_output_file_ != nullptr);
    }
  }

  ~MappingWriter() {
    if (atac_fragment_gz_output_) {
      gzclose(atac_fragment_gz_output_);
      atac_fragment_gz_output_ = nullptr;
    }
    if (atac_fragment_output_file_) {
      fclose(atac_fragment_output_file_);
      atac_fragment_output_file_ = nullptr;
    }
    if (mapping_output_file_) {
      fclose(mapping_output_file_);
    }
    // Finalize sorted output before closing streams (if sorting was enabled)
    FinalizeSortedOutput();  // No-op for non-SAMMapping, specialized for SAMMapping
    CloseYFilterStreams();
    CloseHtsOutput();  // No-op for non-SAMMapping, specialized for SAMMapping
  }

  void OutputTempMappings(
      uint32_t num_reference_sequences,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs,
      std::vector<TempMappingFileHandle<MappingRecord>>
          &temp_mapping_file_handles);

#ifndef LEGACY_OVERFLOW
  void OutputTempMappingsToOverflow(
      uint32_t num_reference_sequences,
      std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs);
      
  void ProcessAndOutputMappingsInLowMemoryFromOverflow(
      uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
      const SequenceBatch &reference,
      const khash_t(k64_seq) * barcode_whitelist_lookup_table);
  
  // Helper to close thread-local writer and collect paths
  static void CloseThreadOverflowWriter();
  
  // Helper to rotate thread-local writer after each flush
  // Closes current writer (collecting file paths) and resets it so next spill creates fresh files
  static void RotateThreadOverflowWriter();
#endif

  void OutputMappings(uint32_t num_reference_sequences,
                      const SequenceBatch &reference,
                      const std::vector<std::vector<MappingRecord>> &mappings);

  void OutputHeader(uint32_t num_reference_sequences,
                    const SequenceBatch &reference);

  void ProcessAndOutputMappingsInLowMemory(
      uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
      const SequenceBatch &reference,
      const khash_t(k64_seq) * barcode_whitelist_lookup_table,
      std::vector<TempMappingFileHandle<MappingRecord>>
          &temp_mapping_file_handles);

  void OutputSummaryMetadata(std::vector<double> frip_est_coeffs = {0.0, 0.0, 0.0, 0.0, 0.0}, bool output_num_cache_slots_info = true);
  void UpdateSummaryMetadata(uint64_t barcode, int type, int change);
  void UpdateSpeicalCategorySummaryMetadata(int category, int type, int change);
  void AdjustSummaryPairedEndOverCount();

  // Set the Y-hit read IDs for filtering (call before any output)
  void SetYHitFilter(const std::unordered_set<uint32_t> *reads_with_y_hit) {
    reads_with_y_hit_ = reads_with_y_hit;
  }
  
  // Open secondary output streams for Y-filtering
  // MUST be called BEFORE OutputHeader() so headers are mirrored
  void OpenYFilterStreams() {
    if (mapping_parameters_.emit_noY_stream && !noY_output_file_) {
      noY_output_file_ = fopen(mapping_parameters_.noY_output_path.c_str(), "w");
      if (!noY_output_file_) {
        ExitWithMessage("Failed to open noY output file: " + 
                        mapping_parameters_.noY_output_path);
      }
    }
    if (mapping_parameters_.emit_Y_stream && !Y_output_file_) {
      Y_output_file_ = fopen(mapping_parameters_.Y_output_path.c_str(), "w");
      if (!Y_output_file_) {
        ExitWithMessage("Failed to open Y output file: " + 
                        mapping_parameters_.Y_output_path);
      }
    }
  }
  
  // Close secondary streams (call at end of mapping)
  void CloseYFilterStreams() {
    if (noY_output_file_) {
      fclose(noY_output_file_);
      noY_output_file_ = nullptr;
    }
    if (Y_output_file_) {
      fclose(Y_output_file_);
      Y_output_file_ = nullptr;
    }
  }
  
  // Finalize sorted output (call before closing streams if sorting was enabled)
  // Default no-op, specialized for SAMMapping
  void FinalizeSortedOutput() {}

 protected:
  void AppendMapping(uint32_t rid, const SequenceBatch &reference,
                     const MappingRecord &mapping);

  inline void AppendMappingOutput(const std::string &line) {
    if (mapping_output_file_) {
      (void)fwrite(line.data(), 1, line.size(), mapping_output_file_);
    }
  }

  // Secondary ATAC fragment stream when AtacDualFragmentAndBam() is active.
  inline void AppendAtacFragmentOutput(const std::string &line) {
    if (!mapping_parameters_.AtacDualFragmentAndBam()) {
      return;
    }
    if (atac_fragment_gz_output_) {
      if (gzwrite(atac_fragment_gz_output_, line.data(),
                  static_cast<unsigned>(line.size())) <= 0) {
        ExitWithMessage("gzwrite failed for --atac-fragments output");
      }
    } else if (atac_fragment_output_file_) {
      if (fwrite(line.data(), 1, line.size(), atac_fragment_output_file_) !=
          line.size()) {
        ExitWithMessage("fwrite failed for --atac-fragments output");
      }
    }
  }
  
  // Helper methods for htslib BAM/CRAM output (only implemented for SAMMapping)
  void OpenHtsOutput() {}  // Default no-op, specialized for SAMMapping
  void CloseHtsOutput() {}  // Default no-op, specialized for SAMMapping
  void BuildHtsHeader(uint32_t num_reference_sequences,
                      const SequenceBatch &reference) {}  // Default no-op
  bam1_t* ConvertToHtsBam(uint32_t rid,
                           const SequenceBatch &reference,
                           const MappingRecord &mapping) { return nullptr; }  // Default no-op
  std::string DeriveReadGroupFromFilename(const std::string &filename) { return ""; }  // Default no-op

  size_t FindBestMappingIndexFromDuplicates(
      const khash_t(k64_seq) * barcode_whitelist_lookup_table,
      const std::vector<MappingRecord> &duplicates);

  void OutputMappingsInVector(
      uint8_t mapq_threshold, uint32_t num_reference_sequences,
      const SequenceBatch &reference,
      const std::vector<std::vector<MappingRecord>> &mappings);

  // Output the mappings in a temp file.
  inline void OutputTempMapping(
      const std::string &temp_mapping_output_file_path,
      uint32_t num_reference_sequences,
      const std::vector<std::vector<MappingRecord>> &mappings) {
    FILE *temp_mapping_output_file =
        fopen(temp_mapping_output_file_path.c_str(), "wb");
    assert(temp_mapping_output_file != NULL);
    for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
      // make sure mappings[ri] exists even if its size is 0
      size_t num_mappings = mappings[ri].size();
      fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
      if (mappings[ri].size() > 0) {
        fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
               temp_mapping_output_file);
      }
    }
    fclose(temp_mapping_output_file);
  }

  // TODO(Haowen): use mapping_output_format_ variable to decide output in BED
  // or TagAlign. It should be removed later.
  const MappingParameters mapping_parameters_;
  const uint32_t cell_barcode_length_;
  FILE *mapping_output_file_ = nullptr;
  BarcodeTranslator barcode_translator_;
  SummaryMetadata summary_metadata_;

  // for pairs
  const std::vector<int> pairs_custom_rid_rank_;

  FILE *atac_fragment_output_file_ = nullptr;
  gzFile atac_fragment_gz_output_ = nullptr;

  // Y-chromosome filtering (SAM/BAM/CRAM mode)
  FILE *noY_output_file_ = nullptr;
  FILE *Y_output_file_ = nullptr;
  const std::unordered_set<uint32_t> *reads_with_y_hit_ = nullptr;
  
  // htslib handles for BAM/CRAM output
  samFile *hts_out_ = nullptr;
  bam_hdr_t *hts_hdr_ = nullptr;
  samFile *noY_hts_out_ = nullptr;
  samFile *Y_hts_out_ = nullptr;
  // Persistent index paths (must outlive htslib file handles)
  std::string hts_index_path_;
  std::string noY_index_path_;
  std::string Y_index_path_;
  
  // BAM sorter for coordinate sorting (only for SAMMapping)
  // bam_sorter.h is included above, so BamSorter is a complete type
  std::unique_ptr<chromap::BamSorter> bam_sorter_;

#ifndef LEGACY_OVERFLOW
  // Thread-local overflow writer (one per OpenMP worker thread)
  static thread_local std::unique_ptr<OverflowWriter> tls_overflow_writer_;
  
  // Shared collection of overflow file paths from all threads
  static std::vector<std::string> shared_overflow_file_paths_;
  static std::mutex overflow_paths_mutex_;
#endif
};

template <typename MappingRecord>
size_t MappingWriter<MappingRecord>::FindBestMappingIndexFromDuplicates(
    const khash_t(k64_seq) * barcode_whitelist_lookup_table,
    const std::vector<MappingRecord> &duplicates) {
  // Find the best barcode, break ties first by the number of the
  // barcodes in the dups, then by the barcode abundance.
  size_t best_mapping_index = 0;

  khiter_t barcode_whitelist_lookup_table_iterator =
      kh_get(k64_seq, barcode_whitelist_lookup_table,
             duplicates[best_mapping_index].GetBarcode());

  double best_mapping_barcode_abundance = kh_value(
      barcode_whitelist_lookup_table,
      barcode_whitelist_lookup_table_iterator);  /// (double)num_sample_barcodes_;

  for (size_t bulk_dup_i = 1; bulk_dup_i < duplicates.size(); ++bulk_dup_i) {
    barcode_whitelist_lookup_table_iterator =
        kh_get(k64_seq, barcode_whitelist_lookup_table,
               duplicates[bulk_dup_i].GetBarcode());

    const double current_mapping_barcode_abundance = kh_value(
        barcode_whitelist_lookup_table,
        barcode_whitelist_lookup_table_iterator);  /// (double)num_sample_barcodes_;

    const bool same_num_dups_with_higer_barcode_abundance =
        duplicates[bulk_dup_i].num_dups_ ==
            duplicates[best_mapping_index].num_dups_ &&
        current_mapping_barcode_abundance > best_mapping_barcode_abundance;

    if (duplicates[bulk_dup_i].num_dups_ >
            duplicates[best_mapping_index].num_dups_ ||
        same_num_dups_with_higer_barcode_abundance) {
      best_mapping_index = bulk_dup_i;
      best_mapping_barcode_abundance = current_mapping_barcode_abundance;
    }
  }
  return best_mapping_index;
}

template <typename MappingRecord>
void MappingWriter<MappingRecord>::ProcessAndOutputMappingsInLowMemory(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table,
    std::vector<TempMappingFileHandle<MappingRecord>>
        &temp_mapping_file_handles) {
  if (temp_mapping_file_handles.size() == 0) {
    return;
  }

  double sort_and_dedupe_start_time = GetRealTime();

  // Calculate block size and initialize
  uint64_t max_mem_size = 10 * ((uint64_t)1 << 30);
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM ||
      mapping_parameters_.mapping_output_format == MAPPINGFORMAT_PAIRS ||
      mapping_parameters_.mapping_output_format == MAPPINGFORMAT_PAF ||
      mapping_parameters_.AtacDualFragmentAndBam()) {
    max_mem_size = (uint64_t)1 << 30;
  }
  for (size_t hi = 0; hi < temp_mapping_file_handles.size(); ++hi) {
    const uint32_t temp_mapping_block_size =
        max_mem_size / temp_mapping_file_handles.size() / sizeof(MappingRecord);

    temp_mapping_file_handles[hi].InitializeTempMappingLoading(
        temp_mapping_block_size);
    temp_mapping_file_handles[hi].LoadTempMappingBlock(num_reference_sequences);
  }

  // Merge and dedupe.
  uint32_t last_rid = std::numeric_limits<uint32_t>::max();
  MappingRecord last_mapping = MappingRecord();
  uint32_t num_last_mapping_dups = 0;
  uint64_t num_uni_mappings = 0;
  uint64_t num_multi_mappings = 0;
  uint64_t num_mappings_passing_filters = 0;
  uint64_t num_total_mappings = 0;
  std::vector<MappingRecord> temp_dups_for_bulk_level_dedup;
  temp_dups_for_bulk_level_dedup.reserve(255);

  const bool deduplicate_at_bulk_level_for_single_cell_data =
      mapping_parameters_.remove_pcr_duplicates &&
      !mapping_parameters_.is_bulk_data &&
      mapping_parameters_.remove_pcr_duplicates_at_bulk_level;

  while (true) {
    // Merge, dedupe and output.
    // Find min first (sorted by rid and then barcode and then positions).
    size_t min_handle_index = temp_mapping_file_handles.size();
    uint32_t min_rid = std::numeric_limits<uint32_t>::max();

    for (size_t hi = 0; hi < temp_mapping_file_handles.size(); ++hi) {
      const TempMappingFileHandle<MappingRecord> &current_handle =
          temp_mapping_file_handles[hi];
      if (current_handle.HasMappings()) {
        const bool rid_is_smaller = current_handle.current_rid < min_rid;
        const bool same_rid_smaller_mapping =
            current_handle.current_rid == min_rid &&
            current_handle.GetCurrentMapping() <
                temp_mapping_file_handles[min_handle_index].GetCurrentMapping();

        if (rid_is_smaller || same_rid_smaller_mapping) {
          min_handle_index = hi;
          min_rid = current_handle.current_rid;
        }
      }
    }

    // All mappings are merged. We only have to handle the case when the last
    // mapping is a duplicate.
    if (min_handle_index == temp_mapping_file_handles.size()) {
      break;
    }

    ++num_total_mappings;

    // Output the current min mapping if it is not a duplicate.
    const MappingRecord &current_min_mapping =
        temp_mapping_file_handles[min_handle_index].GetCurrentMapping();

    const bool is_first_iteration = num_total_mappings == 1;
    const bool current_mapping_is_duplicated_at_cell_level =
        !is_first_iteration && current_min_mapping == last_mapping;
    const bool current_mapping_is_duplicated_at_bulk_level =
        !is_first_iteration && deduplicate_at_bulk_level_for_single_cell_data &&
        current_min_mapping.IsSamePosition(last_mapping);
    const bool current_mapping_is_duplicated =
        last_rid == min_rid && (current_mapping_is_duplicated_at_cell_level ||
                                current_mapping_is_duplicated_at_bulk_level);
    if (mapping_parameters_.remove_pcr_duplicates &&
        current_mapping_is_duplicated) {
      ++num_last_mapping_dups;
      if (deduplicate_at_bulk_level_for_single_cell_data) {
        if (!temp_dups_for_bulk_level_dedup.empty() &&
            current_min_mapping == temp_dups_for_bulk_level_dedup.back()) {
          // Merge if their barcodes are the same. Be careful of "==" here!
          temp_dups_for_bulk_level_dedup.back() = current_min_mapping;
          temp_dups_for_bulk_level_dedup.back().num_dups_ += 1;
        } else {
          temp_dups_for_bulk_level_dedup.push_back(current_min_mapping);
          temp_dups_for_bulk_level_dedup.back().num_dups_ = 1;
        }
      }
      // FIX: Always update to current mapping to match normal mode's behavior
      // Normal mode unconditionally does "last_it = it", keeping the LAST duplicate
      // in sort order (which has highest MAPQ due to ascending sort, and for equal
      // MAPQ, highest read_id). The previous ">" check kept the FIRST duplicate
      // for equal MAPQ, causing parity differences.
      last_mapping = current_min_mapping;
    } else {
      if (!is_first_iteration) {
        if (deduplicate_at_bulk_level_for_single_cell_data) {
          size_t best_mapping_index = FindBestMappingIndexFromDuplicates(
              barcode_whitelist_lookup_table, temp_dups_for_bulk_level_dedup);
          last_mapping = temp_dups_for_bulk_level_dedup[best_mapping_index];

          temp_dups_for_bulk_level_dedup.clear();
        }

        if (last_mapping.mapq_ >= mapping_parameters_.mapq_threshold) {
          last_mapping.num_dups_ =
              std::min((uint32_t)std::numeric_limits<uint8_t>::max(),
                       num_last_mapping_dups);
          if (mapping_parameters_.Tn5_shift) {
            last_mapping.Tn5Shift(mapping_parameters_.Tn5_forward_shift,
                                  mapping_parameters_.Tn5_reverse_shift);
          }

          AppendMapping(last_rid, reference, last_mapping);
          ++num_mappings_passing_filters;
          if (!mapping_parameters_.summary_metadata_file_path.empty())
            summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_DUP,
              num_last_mapping_dups - 1);
        } else {
          if (!mapping_parameters_.summary_metadata_file_path.empty())
            summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_LOWMAPQ, 
                      num_last_mapping_dups);
        }
        if (!mapping_parameters_.summary_metadata_file_path.empty())
          summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_MAPPED, 
                num_last_mapping_dups);

        if (last_mapping.is_unique_ == 1) {
          ++num_uni_mappings;
        } else {
          ++num_multi_mappings;
        }
      }

      last_mapping = current_min_mapping;
      last_rid = min_rid;
      num_last_mapping_dups = 1;

      if (deduplicate_at_bulk_level_for_single_cell_data) {
        temp_dups_for_bulk_level_dedup.push_back(current_min_mapping);
        temp_dups_for_bulk_level_dedup.back().num_dups_ = 1;
      }
    }

    temp_mapping_file_handles[min_handle_index].Next(num_reference_sequences);
  }

  if (last_mapping.mapq_ >= mapping_parameters_.mapq_threshold) {
    if (deduplicate_at_bulk_level_for_single_cell_data) {
      size_t best_mapping_index = FindBestMappingIndexFromDuplicates(
          barcode_whitelist_lookup_table, temp_dups_for_bulk_level_dedup);
      last_mapping = temp_dups_for_bulk_level_dedup[best_mapping_index];

      temp_dups_for_bulk_level_dedup.clear();
    }

    last_mapping.num_dups_ = std::min(
        (uint32_t)std::numeric_limits<uint8_t>::max(), num_last_mapping_dups);
    if (mapping_parameters_.Tn5_shift) {
      last_mapping.Tn5Shift(mapping_parameters_.Tn5_forward_shift,
                            mapping_parameters_.Tn5_reverse_shift);
    }
    AppendMapping(last_rid, reference, last_mapping);
    ++num_mappings_passing_filters;
    
    if (!mapping_parameters_.summary_metadata_file_path.empty())
      summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_DUP,
          num_last_mapping_dups - 1);
  } else {
    if (!mapping_parameters_.summary_metadata_file_path.empty())
      summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_LOWMAPQ, 
          num_last_mapping_dups);
  }
  if (!mapping_parameters_.summary_metadata_file_path.empty())
    summary_metadata_.UpdateCount(last_mapping.GetBarcode(), SUMMARY_METADATA_MAPPED, 
          num_last_mapping_dups);

  if (last_mapping.is_unique_ == 1) {
    ++num_uni_mappings;
  } else {
    ++num_multi_mappings;
  }

  // Delete temp files.
  for (size_t hi = 0; hi < temp_mapping_file_handles.size(); ++hi) {
    temp_mapping_file_handles[hi].FinalizeTempMappingLoading();
    remove(temp_mapping_file_handles[hi].file_path.c_str());
  }

  if (mapping_parameters_.remove_pcr_duplicates) {
    std::cerr << "Sorted, deduped and outputed mappings in "
              << GetRealTime() - sort_and_dedupe_start_time << "s.\n";
  } else {
    std::cerr << "Sorted and outputed mappings in "
              << GetRealTime() - sort_and_dedupe_start_time << "s.\n";
  }
  std::cerr << "# uni-mappings: " << num_uni_mappings
            << ", # multi-mappings: " << num_multi_mappings
            << ", total: " << num_uni_mappings + num_multi_mappings << ".\n";
  std::cerr << "Number of output mappings (passed filters): "
            << num_mappings_passing_filters << "\n";
}

template <typename MappingRecord>
void MappingWriter<MappingRecord>::OutputTempMappings(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs,
    std::vector<TempMappingFileHandle<MappingRecord>>
        &temp_mapping_file_handles) {
  TempMappingFileHandle<MappingRecord> temp_mapping_file_handle;
  temp_mapping_file_handle.file_path =
      mapping_parameters_.mapping_output_file_path + ".temp" +
      std::to_string(temp_mapping_file_handles.size());
  if (mapping_parameters_.mapping_output_file_path == "/dev/stdout"
      || mapping_parameters_.mapping_output_file_path == "/dev/stderr")
  {
    temp_mapping_file_handle.file_path = "chromap_output.temp" +
      std::to_string(temp_mapping_file_handles.size());
  }
  temp_mapping_file_handles.emplace_back(temp_mapping_file_handle);

  OutputTempMapping(temp_mapping_file_handle.file_path, num_reference_sequences,
                    mappings_on_diff_ref_seqs);

  for (uint32_t i = 0; i < num_reference_sequences; ++i) {
    mappings_on_diff_ref_seqs[i].clear();
  }
}

template <typename MappingRecord>
void MappingWriter<MappingRecord>::OutputMappingsInVector(
    uint8_t mapq_threshold, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const std::vector<std::vector<MappingRecord>> &mappings) {
  uint64_t num_mappings_passing_filters = 0;
  for (uint32_t ri = 0; ri < num_reference_sequences; ++ri) {
    for (auto it = mappings[ri].begin(); it != mappings[ri].end(); ++it) {
      uint8_t mapq = (it->mapq_);
      // uint8_t is_unique = (it->is_unique);
      if (mapq >= mapq_threshold) {
        // if (allocate_multi_mappings_ || (only_output_unique_mappings_ &&
        // is_unique == 1)) {
        AppendMapping(ri, reference, *it);
        ++num_mappings_passing_filters;
        //}
        //it->num_dups_ is capped by 255 here, so the count might be different in the
        //  low-mem mode.
        if (!mapping_parameters_.summary_metadata_file_path.empty())
          summary_metadata_.UpdateCount(it->GetBarcode(), SUMMARY_METADATA_DUP,
              it->num_dups_ - 1);
      } else {
        if (!mapping_parameters_.summary_metadata_file_path.empty())
          summary_metadata_.UpdateCount(it->GetBarcode(), SUMMARY_METADATA_LOWMAPQ,
              it->num_dups_);
      }
      if (!mapping_parameters_.summary_metadata_file_path.empty())
        summary_metadata_.UpdateCount(it->GetBarcode(), SUMMARY_METADATA_MAPPED,
            it->num_dups_);
    }
  }
  std::cerr << "Number of output mappings (passed filters): "
            << num_mappings_passing_filters << "\n";
}

template <typename MappingRecord>
void MappingWriter<MappingRecord>::OutputMappings(
    uint32_t num_reference_sequences, const SequenceBatch &reference,
    const std::vector<std::vector<MappingRecord>> &mappings) {
  // if (only_output_unique_mappings_ && mapq_threshold_ < 4)
  //  mapq_threshold_ = 4;
  OutputMappingsInVector(mapping_parameters_.mapq_threshold,
                         num_reference_sequences, reference, mappings);
}

template <typename MappingRecord>
void MappingWriter<MappingRecord>::OutputSummaryMetadata(std::vector<double> frip_est_coeffs, bool output_num_cache_slots_info) {
  if (!mapping_parameters_.summary_metadata_file_path.empty())
  {
    summary_metadata_.Output(
        mapping_parameters_.summary_metadata_file_path.c_str(),
        !mapping_parameters_.barcode_whitelist_file_path.empty() && !mapping_parameters_.output_mappings_not_in_whitelist,
        frip_est_coeffs,
        output_num_cache_slots_info
        );
  }
}

template <typename MappingRecord>
  void MappingWriter<MappingRecord>::UpdateSummaryMetadata(uint64_t barcode, int type, int change) {
  if (!mapping_parameters_.summary_metadata_file_path.empty())
    summary_metadata_.UpdateCount(barcode, type, change);
}

// category: 0: non-whitelist barcode
template <typename MappingRecord>
  void MappingWriter<MappingRecord>::UpdateSpeicalCategorySummaryMetadata(int category, int type, int change) {
  if (!mapping_parameters_.summary_metadata_file_path.empty()) {
    if (category == 0)
      summary_metadata_.UpdateNonWhitelistCount(type, change);
  }
}

template <typename MappingRecord>
  void MappingWriter<MappingRecord>::AdjustSummaryPairedEndOverCount() {
  if (!mapping_parameters_.summary_metadata_file_path.empty()
      && mapping_parameters_.mapping_output_format == MAPPINGFORMAT_SAM)
      summary_metadata_.AdjustPairedEndOverCount() ; 
}


// Specialization for BED format.
template <>
void MappingWriter<MappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<MappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithBarcode &mapping);

template <>
void MappingWriter<MappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<MappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithoutBarcode &mapping);

// Specialization for BEDPE format.
template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithoutBarcode &mapping);

template <>
void MappingWriter<PairedEndMappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<PairedEndMappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithBarcode &mapping);

// Specialization for PAF format.
template <>
void MappingWriter<PAFMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference);

template <>
void MappingWriter<PAFMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const PAFMapping &mapping);

template <>
void MappingWriter<PAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PAFMapping>> &mappings);

// Specialization for PairedPAF format.
template <>
void MappingWriter<PairedPAFMapping>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<PairedPAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairedPAFMapping>> &mappings);

template <>
void MappingWriter<PairedPAFMapping>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedPAFMapping &mapping);

// Specialization for SAM format.
template <>
void MappingWriter<SAMMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference);

template <>
void MappingWriter<SAMMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const SAMMapping &mapping);

template <>
void MappingWriter<SAMMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<SAMMapping>> &mappings);

// Specialization for htslib helper methods
template <>
void MappingWriter<SAMMapping>::OpenHtsOutput();

template <>
void MappingWriter<SAMMapping>::CloseHtsOutput();

template <>
void MappingWriter<SAMMapping>::BuildHtsHeader(uint32_t num_reference_sequences,
                                               const SequenceBatch &reference);

template <>
bam1_t* MappingWriter<SAMMapping>::ConvertToHtsBam(uint32_t rid,
                                                    const SequenceBatch &reference,
                                                    const SAMMapping &mapping);

template <>
std::string MappingWriter<SAMMapping>::DeriveReadGroupFromFilename(
    const std::string &filename);

// Specialization for Y-filter streams with htslib support
template <>
void MappingWriter<SAMMapping>::OpenYFilterStreams();

template <>
void MappingWriter<SAMMapping>::CloseYFilterStreams();

// Specialization for sorted output finalization
template <>
void MappingWriter<SAMMapping>::FinalizeSortedOutput();

// ATAC dual output: BED fragments + BAM/CRAM (paired barcoded reads only).
template <>
void MappingWriter<PairedEndAtacDualMapping>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<PairedEndAtacDualMapping>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndAtacDualMapping &mapping);

template <>
void MappingWriter<PairedEndAtacDualMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairedEndAtacDualMapping>> &mappings);

template <>
void MappingWriter<PairedEndAtacDualMapping>::OpenHtsOutput();

template <>
void MappingWriter<PairedEndAtacDualMapping>::CloseHtsOutput();

template <>
void MappingWriter<PairedEndAtacDualMapping>::BuildHtsHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference);

template <>
void MappingWriter<PairedEndAtacDualMapping>::OpenYFilterStreams();

template <>
void MappingWriter<PairedEndAtacDualMapping>::CloseYFilterStreams();

template <>
void MappingWriter<PairedEndAtacDualMapping>::FinalizeSortedOutput();

// Specialization for pairs format.
template <>
void MappingWriter<PairsMapping>::OutputHeader(uint32_t num_reference_sequences,
                                               const SequenceBatch &reference);

template <>
void MappingWriter<PairsMapping>::AppendMapping(uint32_t rid,
                                                const SequenceBatch &reference,
                                                const PairsMapping &mapping);

template <>
void MappingWriter<PairsMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairsMapping>> &mappings);

}  // namespace chromap

#endif  // MAPPING_WRITER_H_
