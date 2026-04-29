#ifndef MAPPING_PARAMETERS_H_
#define MAPPING_PARAMETERS_H_

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "libmacs3/frag_compact_store.h"
#include "libmacs3/fragments.h"

namespace chromap {

enum MappingOutputFormat {
  MAPPINGFORMAT_UNKNOWN,
  MAPPINGFORMAT_BED,
  MAPPINGFORMAT_TAGALIGN,
  MAPPINGFORMAT_PAF,
  MAPPINGFORMAT_SAM,
  MAPPINGFORMAT_PAIRS,
  MAPPINGFORMAT_BAM,
  MAPPINGFORMAT_CRAM
};

// Source for --call-macs3-frag-peaks fragment rows (file reread vs in-memory).
enum class Macs3FragPeaksSource { kFile, kMemory };

struct MappingParameters {
  int error_threshold = 8;
  int match_score = 1;
  int mismatch_penalty = 4;
  std::vector<int> gap_open_penalties = {6, 6};
  std::vector<int> gap_extension_penalties = {1, 1};
  int min_num_seeds_required_for_mapping = 2;
  std::vector<int> max_seed_frequencies = {500, 1000};

  double cache_update_param = 0.01;
  int cache_size = 4000003;
  bool debug_cache = false;
  std::string frip_est_params = "-1.0996;4.2391;3.0164e-05;-2.1087e-04;-5.5825e-05";
  bool output_num_uniq_cache_slots = true;
  int k_for_minhash = 250;

  // Read with # best mappings greater than it will have this number of best
  // mappings reported.
  int max_num_best_mappings = 1;
  int max_insert_size = 1000;
  uint8_t mapq_threshold = 30;
  int num_threads = 1;
  int min_read_length = 30;
  int barcode_correction_error_threshold = 1;
  double barcode_correction_probability_threshold = 0.9;
  int multi_mapping_allocation_distance = 0;
  int multi_mapping_allocation_seed = 11;
  // Read with more than this number of mappings will be dropped.
  int drop_repetitive_reads = 500000;
  bool trim_adapters = false;
  bool remove_pcr_duplicates = false;
  bool remove_pcr_duplicates_at_bulk_level = true;
  bool is_bulk_data = true;
  bool allocate_multi_mappings = false;
  bool only_output_unique_mappings = true;
  bool output_mappings_not_in_whitelist = false;
  bool Tn5_shift = false;
  // Tn5 cut-site offsets applied to fragment 5' ends when Tn5_shift=true.
  // Defaults match ARC / Cell Ranger ATAC / Buenrostro 2013 ("classical" +4/-5).
  // Set to (+4, -4) for the "symmetric" convention used by ChromBPNet.
  // rs is signed; -5 means the reverse-strand 5' end moves by 5 bases.
  int Tn5_forward_shift = 4;
  int Tn5_reverse_shift = -5;
  bool split_alignment = false;
  MappingOutputFormat mapping_output_format = MAPPINGFORMAT_BED;
  bool low_memory_mode = false;
  uint64_t low_mem_ram_limit = 0;  // 0 = use default spill threshold
  bool cell_by_bin = false;
  int bin_size = 5000;
  uint16_t depth_cutoff_to_call_peak = 3;
  int peak_min_length = 30;
  int peak_merge_max_length = 30;
  std::string reference_file_path;
  std::string index_file_path;
  std::vector<std::string> read_file1_paths;
  std::vector<std::string> read_file2_paths;
  std::vector<std::string> barcode_file_paths;
  std::string barcode_whitelist_file_path;
  std::string read_format;
  std::string mapping_output_file_path;
  // When set alongside paired-end barcoded reads and BAM/CRAM primary output,
  // emit fragment lines (BED/TagAlign) to this path in the same alignment pass.
  std::string atac_fragment_output_file_path;
  // Optional fixed-record binary replacement for the secondary ATAC fragments
  // TSV. Records store chrom_id/start/end/count/packed barcode key; chrom
  // names are written to <path>.chroms.tsv.
  std::string atac_fragment_binary_output_file_path;
  std::string matrix_output_prefix;
  // The order for general sorting.
  std::string custom_rid_order_file_path;
  // The order for pairs format flipping.
  std::string pairs_flipping_custom_rid_order_file_path;
  std::string barcode_translate_table_file_path;
  // If true, the barcode translate table is read in natural source<TAB>dest
  // order (col1 = from_bc / hash key, col2 = to_bc / value). The default
  // (false) preserves the historical Chromap convention where col1 = to_bc
  // and col2 = from_bc.
  bool barcode_translate_from_first_column = false;
  std::string summary_metadata_file_path;
  std::string temp_directory_path;  // Directory for temporary files (NEW_OVERFLOW)
  bool skip_barcode_check = false;

  // Y-chromosome stream filtering (SAM/BAM/CRAM mode)
  bool emit_noY_stream = false;
  bool emit_Y_stream = false;
  std::string noY_output_path;
  std::string Y_output_path;

  // Y read names and FASTQ emission
  bool emit_y_read_names = false;
  std::string y_read_names_output_path;
  bool emit_y_noy_fastq = false;
  std::string y_noy_fastq_compression = "gz";  // "gz" or "none"
  std::string y_fastq_output_prefix;
  std::string noy_fastq_output_prefix;
  std::vector<std::vector<std::string>> y_fastq_output_paths_per_file;  // [file][mate]
  std::vector<std::vector<std::string>> noy_fastq_output_paths_per_file;  // [file][mate]

  // htslib BAM/CRAM output parameters
  int hts_threads = 0;           // 0 = auto: min(num_threads, 4)
  std::string read_group_id;     // "" = none, "auto" = generate from filenames, else literal RG ID
  bool write_index = false;
  bool sort_bam = false;          // Enable coordinate sorting for BAM/CRAM output
  uint64_t sort_bam_ram_limit = 8ULL * 1024 * 1024 * 1024;  // 8GB default

  // Opt-in MACS3-compatible FRAG narrowPeak (post-process fragments file; not default).
  bool call_macs3_frag_peaks = false;
  // file: reread --atac-fragments (default). memory: in-memory rows from mapping.
  Macs3FragPeaksSource macs3_frag_peaks_source = Macs3FragPeaksSource::kFile;
  int macs3_frag_compact_min_count_bits = 16;
  std::shared_ptr<peaks::FragPeakMemoryAccumulator> macs3_frag_memory_accumulator;
  // kMemory peak source: writer pushes FragmentRecord into per-chrom
  // buckets (indexed by chrom_id) instead of one flat vector. Fragments
  // arrive chrom-grouped + start-sorted from OutputMappingsInVector
  // (single-threaded loop), so each bucket ends up sorted with no extra
  // work. The driver feeds the buckets into either:
  //   - the events workspace (RunMacs3FragPeakPipelineFromFragments,
  //     parallel per-chrom finalize) when macs3_frag_low_mem == false
  //     (default — fast, ~+3.8 GB RAM for the events buffer), or
  //   - the sweep workspace (WrapVectorFragmentIterator +
  //     RunMacs3FragPeakPipelineFromSortedIterator) when
  //     macs3_frag_low_mem == true (~tens of KB peak heap; pick this
  //     when sharing a host with STAR or another large-RAM consumer).
  // Memory profile of the buckets themselves: ~16 B/fragment regardless
  // of mode.
  std::shared_ptr<std::vector<std::vector<macs3::FragmentRecord>>>
      macs3_frag_buffer;
  std::shared_ptr<std::vector<std::string>> macs3_frag_chrom_names;

  std::string macs3_frag_peaks_narrowpeak_path;
  std::string macs3_frag_peaks_summits_path;
  double macs3_frag_pvalue = 1e-5;
  int macs3_frag_min_length = 200;
  int macs3_frag_max_gap = 30;
  bool macs3_frag_uint8_counts = true;
  // When true, force the sweep-line workspace for the in-memory peak
  // path (lower peak RSS, slightly slower wall). When false (default),
  // use the events workspace (parallel per-chrom finalize, faster wall,
  // ~+3.8 GB RAM at PBMC 3k scale). The kFile peak source always uses
  // the events workspace regardless of this flag.
  bool macs3_frag_low_mem = false;
  std::string macs3_frag_keep_intermediates_dir;

  // Optional STAR-managed dynamic-permit hooks. Signatures match
  // process_features' pf_permit_{acquire,release}_fn so the same STAR-side
  // shim shape can carry both domains. When both pointers are non-null
  // the chromap PE batch loop will call acquire() before each per-thread
  // mini-batch (~64 read pairs) and release() at the end of the batch
  // with telemetry (waitNs from acquire, work_units, work_bytes, work_ns).
  uint64_t (*permit_acquire_hook)(void *hook_ctx) = nullptr;
  void (*permit_release_hook)(void *hook_ctx,
                              uint64_t wait_ns,
                              uint64_t work_units,
                              uint64_t work_bytes,
                              uint64_t work_ns) = nullptr;
  void *permit_hook_ctx = nullptr;
  bool PermitHooksEnabled() const {
    return permit_acquire_hook != nullptr && permit_release_hook != nullptr;
  }

  // Dual ATAC: BAM/CRAM to mapping_output_file_path and fragments to
  // atac_fragment_output_file_path (one pass; not supported with --low-mem).
  bool AtacDualFragmentAndBam() const {
    return !read_file2_paths.empty() && !barcode_file_paths.empty() &&
           !atac_fragment_output_file_path.empty() &&
           (mapping_output_format == MAPPINGFORMAT_BAM ||
            mapping_output_format == MAPPINGFORMAT_CRAM);
  }

  int GetNumVPULanes() const {
    int NUM_VPU_LANES = 0;
    if (error_threshold < 8) {
      NUM_VPU_LANES = 8;
    } else if (error_threshold < 16) {
      NUM_VPU_LANES = 4;
    }
    return NUM_VPU_LANES;
  }
};

}  // namespace chromap

#endif  // MAPPING_PARAMETERS_H_
