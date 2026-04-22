#include "mapping_writer.h"

#include <queue>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <unistd.h>
#include <cstdlib>
#include "chromap.h"
#include "bam_sorter.h"

namespace chromap {

#ifndef LEGACY_OVERFLOW
// Static member definitions for thread-local and shared overflow handling
template <typename MappingRecord>
thread_local std::unique_ptr<OverflowWriter> MappingWriter<MappingRecord>::tls_overflow_writer_;

template <typename MappingRecord>
std::vector<std::string> MappingWriter<MappingRecord>::shared_overflow_file_paths_;

template <typename MappingRecord>
std::mutex MappingWriter<MappingRecord>::overflow_paths_mutex_;
#endif

// Specialization for BED format.
template <>
void MappingWriter<MappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<MappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithBarcode &mapping) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\t" +
                              barcode_translator_.Translate(
                                  mapping.cell_barcode_, cell_barcode_length_) +
                              "\t" + std::to_string(mapping.num_dups_) + "\n");
  } else {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\n");
  }
}

template <>
void MappingWriter<MappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<MappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const MappingWithoutBarcode &mapping) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\t" + std::to_string(mapping.num_dups_) + "\n");
  } else {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\t" + std::to_string(mapping.num_dups_) + "\n");
  }
}

// Specialization for BEDPE format.
template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithoutBarcode &mapping) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\tN\t" + std::to_string(mapping.mapq_) + "\t" +
                              strand + "\t" + std::to_string(mapping.num_dups_) + "\n");
  } else {
    bool positive_strand = mapping.IsPositiveStrand();
    uint32_t positive_read_end =
        mapping.fragment_start_position_ + mapping.positive_alignment_length_;
    uint32_t negative_read_end =
        mapping.fragment_start_position_ + mapping.fragment_length_;
    uint32_t negative_read_start =
        negative_read_end - mapping.negative_alignment_length_;
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    if (positive_strand) {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\t" + 
          std::to_string(mapping.num_dups_) + "\n");
    } else {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\t" +
          std::to_string(mapping.num_dups_) + "\n");
    }
  }
}

template <>
void MappingWriter<PairedEndMappingWithBarcode>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<PairedEndMappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithBarcode &mapping) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\t" +
                              barcode_translator_.Translate(
                                  mapping.cell_barcode_, cell_barcode_length_) +
                              "\t" + std::to_string(mapping.num_dups_) + "\n");
  } else {
    bool positive_strand = mapping.IsPositiveStrand();
    uint32_t positive_read_end =
        mapping.fragment_start_position_ + mapping.positive_alignment_length_;
    uint32_t negative_read_end =
        mapping.fragment_start_position_ + mapping.fragment_length_;
    uint32_t negative_read_start =
        negative_read_end - mapping.negative_alignment_length_;
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    if (positive_strand) {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\n");
    } else {
      this->AppendMappingOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(negative_read_start) + "\t" +
          std::to_string(negative_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t-\n" +
          std::string(reference_sequence_name) + "\t" +
          std::to_string(mapping.fragment_start_position_) + "\t" +
          std::to_string(positive_read_end) + "\tN\t" +
          std::to_string(mapping.mapq_) + "\t+\n");
    }
  }
}

// Specialization for PAF format.
template <>
void MappingWriter<PAFMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference) {}

template <>
void MappingWriter<PAFMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const PAFMapping &mapping) {
  const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
  uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
  std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
  uint32_t mapping_end_position =
      mapping.fragment_start_position_ + mapping.fragment_length_;
  this->AppendMappingOutput(
      mapping.read_name_ + "\t" + std::to_string(mapping.read_length_) + "\t" +
      std::to_string(0) + "\t" + std::to_string(mapping.read_length_) + "\t" +
      strand + "\t" + std::string(reference_sequence_name) + "\t" +
      std::to_string(reference_sequence_length) + "\t" +
      std::to_string(mapping.fragment_start_position_) + "\t" +
      std::to_string(mapping_end_position) + "\t" +
      std::to_string(mapping.read_length_) + "\t" +
      std::to_string(mapping.fragment_length_) + "\t" +
      std::to_string(mapping.mapq_) + "\n");
}

template <>
void MappingWriter<PAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PAFMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

// Specialization for PairedPAF format.
template <>
void MappingWriter<PairedPAFMapping>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {}

template <>
void MappingWriter<PairedPAFMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairedPAFMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

template <>
void MappingWriter<PairedPAFMapping>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedPAFMapping &mapping) {
  bool positive_strand = mapping.IsPositiveStrand();
  uint32_t positive_read_end =
      mapping.fragment_start_position_ + mapping.positive_alignment_length_;
  uint32_t negative_read_end =
      mapping.fragment_start_position_ + mapping.fragment_length_;
  uint32_t negative_read_start =
      negative_read_end - mapping.negative_alignment_length_;
  const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
  uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
  if (positive_strand) {
    this->AppendMappingOutput(
        mapping.read1_name_ + "\t" + std::to_string(mapping.read1_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" + "+" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(mapping.fragment_start_position_) + "\t" +
        std::to_string(positive_read_end) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" +
        std::to_string(mapping.positive_alignment_length_) + "\t" +
        std::to_string(mapping.mapq1_) + "\n");
    this->AppendMappingOutput(
        mapping.read2_name_ + "\t" + std::to_string(mapping.read2_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" + "-" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(negative_read_start) + "\t" +
        std::to_string(negative_read_end) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" +
        std::to_string(mapping.negative_alignment_length_) + "\t" +
        std::to_string(mapping.mapq2_) + "\n");
  } else {
    this->AppendMappingOutput(
        mapping.read1_name_ + "\t" + std::to_string(mapping.read1_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" + "-" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(negative_read_start) + "\t" +
        std::to_string(negative_read_end) + "\t" +
        std::to_string(mapping.read1_length_) + "\t" +
        std::to_string(mapping.negative_alignment_length_) + "\t" +
        std::to_string(mapping.mapq1_) + "\n");
    this->AppendMappingOutput(
        mapping.read2_name_ + "\t" + std::to_string(mapping.read2_length_) +
        "\t" + std::to_string(0) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" + "+" + "\t" +
        std::string(reference_sequence_name) + "\t" +
        std::to_string(reference_sequence_length) + "\t" +
        std::to_string(mapping.fragment_start_position_) + "\t" +
        std::to_string(positive_read_end) + "\t" +
        std::to_string(mapping.read2_length_) + "\t" +
        std::to_string(mapping.positive_alignment_length_) + "\t" +
        std::to_string(mapping.mapq2_) + "\n");
  }
}

// Specialization for SAM format.
template <>
void MappingWriter<SAMMapping>::OutputHeader(uint32_t num_reference_sequences,
                                             const SequenceBatch &reference) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM ||
      mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM) {
    // BAM/CRAM path: open htslib output first
    if (!hts_out_) {
      OpenHtsOutput();
    }
    // Build and write header
    BuildHtsHeader(num_reference_sequences, reference);
    
    // Write header to Y-filter streams if they exist
    if (noY_hts_out_ && hts_hdr_) {
      if (sam_hdr_write(noY_hts_out_, hts_hdr_) < 0) {
        ExitWithMessage("Failed to write header to noY BAM/CRAM output");
      }
      // Initialize indexing only if requested and output is not stdout
      // Store index path in persistent member to avoid dangling pointer
      if (mapping_parameters_.write_index &&
          mapping_parameters_.noY_output_path != "-" &&
          mapping_parameters_.noY_output_path != "/dev/stdout" &&
          mapping_parameters_.noY_output_path != "/dev/stderr") {
        this->noY_index_path_ = mapping_parameters_.noY_output_path;
        if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) {
          this->noY_index_path_ += ".bai";  // Use BAI for compatibility
        } else if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM) {
          this->noY_index_path_ += ".crai";
        }
        // Use BAI format (min_shift=0) for wide compatibility
        int min_shift = 0;
        if (sam_idx_init(noY_hts_out_, hts_hdr_, min_shift, this->noY_index_path_.c_str()) < 0) {
          ExitWithMessage("Failed to initialize index for noY BAM/CRAM output");
        }
      }
    }
    if (Y_hts_out_ && hts_hdr_) {
      if (sam_hdr_write(Y_hts_out_, hts_hdr_) < 0) {
        ExitWithMessage("Failed to write header to Y BAM/CRAM output");
      }
      // Initialize indexing only if requested and output is not stdout
      // Store index path in persistent member to avoid dangling pointer
      if (mapping_parameters_.write_index &&
          mapping_parameters_.Y_output_path != "-" &&
          mapping_parameters_.Y_output_path != "/dev/stdout" &&
          mapping_parameters_.Y_output_path != "/dev/stderr") {
        this->Y_index_path_ = mapping_parameters_.Y_output_path;
        if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) {
          this->Y_index_path_ += ".bai";  // Use BAI for compatibility
        } else if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM) {
          this->Y_index_path_ += ".crai";
        }
        // Use BAI format (min_shift=0) for wide compatibility
        int min_shift = 0;
        if (sam_idx_init(Y_hts_out_, hts_hdr_, min_shift, this->Y_index_path_.c_str()) < 0) {
          ExitWithMessage("Failed to initialize index for Y BAM/CRAM output");
        }
      }
    }
  } else {
    // SAM text path
    for (uint32_t rid = 0; rid < num_reference_sequences; ++rid) {
      const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
      uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
      std::string header_line = "@SQ\tSN:" + std::string(reference_sequence_name) +
                                "\tLN:" + std::to_string(reference_sequence_length) + "\n";
      
      // Write to primary output
      this->AppendMappingOutput(header_line);
      
      // Mirror to secondary streams (must be open before this call)
      if (noY_output_file_) {
        fwrite(header_line.data(), 1, header_line.size(), noY_output_file_);
      }
      if (Y_output_file_) {
        fwrite(header_line.data(), 1, header_line.size(), Y_output_file_);
      }
    }
  }
}

template <>
void MappingWriter<SAMMapping>::AppendMapping(uint32_t rid,
                                              const SequenceBatch &reference,
                                              const SAMMapping &mapping) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM ||
      mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM) {
    // === htslib BAM/CRAM path ===
    bam1_t *b = ConvertToHtsBam(rid, reference, mapping);
    
    if (mapping_parameters_.sort_bam && bam_sorter_) {
      // Buffer BAM record to sorter instead of writing directly
      bool hasY = reads_with_y_hit_ && reads_with_y_hit_->count(mapping.read_id_) > 0;
      
      // Pass bam1_t directly to sorter (sorter handles serialization internally)
      bam_sorter_->addRecord(b, mapping.read_id_, hasY);
      
      bam_destroy1(b);
    } else {
      // Direct write path (existing behavior)
      // Write to primary output
      if (sam_write1(hts_out_, hts_hdr_, b) < 0) {
        bam_destroy1(b);
        ExitWithMessage("Failed to write BAM/CRAM record");
      }
      
      // Route to Y-filter streams (existing --emit-noY-bam / --emit-Y-bam)
      // MUST write BEFORE bam_destroy1 to avoid use-after-free
      if (reads_with_y_hit_ && (noY_hts_out_ || Y_hts_out_)) {
        bool is_y_hit = reads_with_y_hit_->count(mapping.read_id_) > 0;
        if (Y_hts_out_ && is_y_hit) {
          if (sam_write1(Y_hts_out_, hts_hdr_, b) < 0) {
            bam_destroy1(b);
            ExitWithMessage("Failed to write BAM/CRAM record to Y stream");
          }
        }
        if (noY_hts_out_ && !is_y_hit) {
          if (sam_write1(noY_hts_out_, hts_hdr_, b) < 0) {
            bam_destroy1(b);
            ExitWithMessage("Failed to write BAM/CRAM record to noY stream");
          }
        }
      }
      
      // Destroy AFTER all writes complete
      bam_destroy1(b);
    }
  } else {
    // === SAM text path ===
    const char *reference_sequence_name =
        (mapping.flag_ & BAM_FUNMAP) > 0 ? "*" : reference.GetSequenceNameAt(rid);
    const char *mate_ref_sequence_name =
        mapping.mrid_ < 0 ? "*" : 
        ((uint32_t)mapping.mrid_ == rid ? "=" : reference.GetSequenceNameAt(mapping.mrid_));
    const uint32_t mapping_start_position = mapping.GetStartPosition();
    const uint32_t mate_mapping_start_position = mapping.mrid_ < 0 ? 0 : (mapping.mpos_ + 1);

    // Build the entire SAM record in one string to ensure atomic output
    std::string out;
    out.reserve(256 + mapping.sequence_.size() + mapping.sequence_qual_.size() + mapping.MD_.size());
    
    out.append(mapping.read_name_);
    out.push_back('\t');
    out.append(std::to_string(mapping.flag_));
    out.push_back('\t');
    out.append(reference_sequence_name);
    out.push_back('\t');
    out.append(std::to_string(mapping_start_position));
    out.push_back('\t');
    out.append(std::to_string(mapping.mapq_));
    out.push_back('\t');
    out.append(mapping.GenerateCigarString());
    out.push_back('\t');
    out.append(mate_ref_sequence_name);
    out.push_back('\t');
    out.append(std::to_string(mate_mapping_start_position));
    out.push_back('\t');
    out.append(std::to_string(mapping.tlen_));
    out.push_back('\t');
    out.append(mapping.sequence_);
    out.push_back('\t');
    out.append(mapping.sequence_qual_);
    out.push_back('\t');
    out.append(mapping.GenerateIntTagString("NM", mapping.NM_));
    out.append("\tMD:Z:");
    out.append(mapping.MD_);
    
    if (cell_barcode_length_ > 0) {
      out.append("\tCB:Z:");
      out.append(barcode_translator_.Translate(mapping.cell_barcode_, cell_barcode_length_));
    }
    
    out.push_back('\n');
    
    // Write to primary output
    this->AppendMappingOutput(out);
    
    // Route to Y-filter streams based on read ID
    if (reads_with_y_hit_ && (noY_output_file_ || Y_output_file_)) {
      bool is_y_hit = reads_with_y_hit_->count(mapping.read_id_) > 0;
      
      if (Y_output_file_ && is_y_hit) {
        fwrite(out.data(), 1, out.size(), Y_output_file_);
      }
      if (noY_output_file_ && !is_y_hit) {
        fwrite(out.data(), 1, out.size(), noY_output_file_);
      }
    }
  }
}

template <>
void MappingWriter<SAMMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<SAMMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

// Specialization for pairs format.
template <>
void MappingWriter<PairsMapping>::OutputHeader(uint32_t num_reference_sequences,
                                               const SequenceBatch &reference) {
  std::vector<uint32_t> rid_order;
  rid_order.resize(num_reference_sequences);
  uint32_t i;
  for (i = 0; i < num_reference_sequences; ++i) {
    rid_order[pairs_custom_rid_rank_[i]] = i;
  }
  this->AppendMappingOutput("## pairs format v1.0.0\n#shape: upper triangle\n");
  for (i = 0; i < num_reference_sequences; ++i) {
    uint32_t rid = rid_order[i];
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t reference_sequence_length = reference.GetSequenceLengthAt(rid);
    this->AppendMappingOutput(
        "#chromsize: " + std::string(reference_sequence_name) + " " +
        std::to_string(reference_sequence_length) + "\n");
  }
  this->AppendMappingOutput(
      "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n");
}

template <>
void MappingWriter<PairsMapping>::AppendMapping(uint32_t rid,
                                                const SequenceBatch &reference,
                                                const PairsMapping &mapping) {
  const char *reference_sequence_name1 =
      reference.GetSequenceNameAt(mapping.rid1_);
  const char *reference_sequence_name2 =
      reference.GetSequenceNameAt(mapping.rid2_);
  this->AppendMappingOutput(mapping.read_name_ + "\t" +
                            std::string(reference_sequence_name1) + "\t" +
                            std::to_string(mapping.GetPosition(1)) + "\t" +
                            std::string(reference_sequence_name2) + "\t" +
                            std::to_string(mapping.GetPosition(2)) + "\t" +
                            std::string(1, mapping.GetStrand(1)) + "\t" +
                            std::string(1, mapping.GetStrand(2)) + "\tUU\n");
}

template <>
void MappingWriter<PairsMapping>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<PairsMapping> > &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    // make sure mappings[ri] exists even if its size is 0
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
      // fwrite(mappings[ri].data(), sizeof(MappingRecord), mappings[ri].size(),
      // temp_mapping_output_file);
    }
  }
  fclose(temp_mapping_output_file);
}

#ifndef LEGACY_OVERFLOW
// Overflow writer methods
template <typename MappingRecord>
void MappingWriter<MappingRecord>::OutputTempMappingsToOverflow(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingRecord>> &mappings_on_diff_ref_seqs) {
  
  // Initialize thread-local overflow writer if not already done
  if (!tls_overflow_writer_) {
    // Use user-specified temp directory, or let OverflowWriter choose optimal location
    std::string base_dir = mapping_parameters_.temp_directory_path;
    tls_overflow_writer_ = std::unique_ptr<OverflowWriter>(new OverflowWriter(base_dir, "chromap"));
  }
  
  // Write all mappings to overflow files using thread-local writer
  for (uint32_t rid = 0; rid < num_reference_sequences; ++rid) {
    for (const auto& mapping : mappings_on_diff_ref_seqs[rid]) {
      tls_overflow_writer_->Write(rid, mapping);
    }
    mappings_on_diff_ref_seqs[rid].clear();
  }
}

template <typename MappingRecord>
void MappingWriter<MappingRecord>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  
  // All thread-local writers should already be closed by now
  // Use the shared collection of overflow file paths
  
  if (shared_overflow_file_paths_.empty()) {
    return;
  }
  
  std::cerr << "Processing " << shared_overflow_file_paths_.size() << " overflow files for k-way merge\n";

  double sort_and_dedupe_start_time = GetRealTime();

  // Structure to hold current record from each file for k-way merge
  struct FileRecord {
    uint32_t rid;
    MappingRecord mapping;
    size_t file_index;
    
    bool operator<(const FileRecord& other) const {
      if (rid != other.rid) {
        return rid > other.rid;  // Min-heap: smaller rid first
      }
      // FIX: Use strict weak ordering for equal mappings
      // For equal mappings (!(a<b) && !(b<a)), use file_index as tiebreaker
      // This ensures deterministic ordering and prevents undefined behavior in std::priority_queue
      bool a_less_b = mapping < other.mapping;
      bool b_less_a = other.mapping < mapping;
      if (!a_less_b && !b_less_a) {
        // Equal mappings: use file_index as tiebreaker for deterministic ordering
        return file_index > other.file_index;  // Smaller file_index first (stable)
      }
      return !a_less_b;  // Min-heap: smaller mapping first (inverted for max-heap semantics)
    }
  };

  // First, scan files to determine which rids they contain
  // This allows us to group files by rid and process rids in ascending order
  // Note: OverflowWriter creates one file per rid, but we scan to handle edge cases
  std::unordered_map<uint32_t, std::vector<size_t>> rid_to_files;
  std::unordered_set<uint32_t> all_rids;
  
  // Scan all files to build rid -> file_index mapping
  for (size_t fi = 0; fi < shared_overflow_file_paths_.size(); ++fi) {
    OverflowReader scanner(shared_overflow_file_paths_[fi]);
    if (!scanner.IsValid()) {
      continue;
    }
    
    // Scan file to find all rids it contains
    uint32_t rid;
    std::string payload;
    std::unordered_set<uint32_t> rids_in_file;
    while (scanner.ReadNext(rid, payload)) {
      rids_in_file.insert(rid);
      all_rids.insert(rid);
    }
    
    // Add this file to the list for each rid it contains
    for (uint32_t r : rids_in_file) {
      rid_to_files[r].push_back(fi);
    }
  }
  
  // Create readers for actual merge (will be opened per-rid to limit FDs)
  std::vector<std::unique_ptr<OverflowReader>> readers;
  readers.resize(shared_overflow_file_paths_.size());

  // Merge and dedupe (replicating logic from ProcessAndOutputMappingsInLowMemory)
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

  // Process rids in ascending order to preserve coordinate-sorted output
  std::vector<uint32_t> rids(all_rids.begin(), all_rids.end());
  std::sort(rids.begin(), rids.end());

  // For each rid, k-way merge all files containing that rid
  for (uint32_t current_rid : rids) {
    const auto& file_indices = rid_to_files[current_rid];
    
    // Priority queue for k-way merge within this rid
    std::priority_queue<FileRecord> merge_heap;
    
    // Initialize: read first record from each file for this rid
    // Note: OverflowWriter creates one file per rid, so each file contains only records for one rid
    for (size_t fi : file_indices) {
      readers[fi].reset(new OverflowReader(shared_overflow_file_paths_[fi]));
      if (!readers[fi] || !readers[fi]->IsValid()) {
        continue;
      }
      
      uint32_t rid;
      std::string payload;
      if (readers[fi]->ReadNext(rid, payload)) {
        // Verify rid matches (should always be true, but check for safety)
        if (rid == current_rid) {
          FILE* mem_file = fmemopen(const_cast<char*>(payload.data()), payload.size(), "rb");
          if (mem_file) {
            MappingRecord mapping;
            mapping.LoadFromFile(mem_file);
            fclose(mem_file);
            
            FileRecord rec;
            rec.rid = rid;
            rec.mapping = mapping;
            rec.file_index = fi;
            merge_heap.push(rec);
          }
        }
      }
    }

    // K-way merge for this rid
    while (!merge_heap.empty()) {
      FileRecord min_rec = merge_heap.top();
      merge_heap.pop();
      
      ++num_total_mappings;
      
      const MappingRecord& current_min_mapping = min_rec.mapping;
      const uint32_t min_rid = min_rec.rid;
      
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

      // Read next record from the file we just processed
      // Note: Since OverflowWriter creates one file per rid, all records in this file are for current_rid
      if (min_rec.file_index < readers.size() && readers[min_rec.file_index]) {
        uint32_t rid;
        std::string payload;
        if (readers[min_rec.file_index]->ReadNext(rid, payload)) {
          // Should always be current_rid, but verify for safety
          if (rid == current_rid) {
            FILE* mem_file = fmemopen(const_cast<char*>(payload.data()), payload.size(), "rb");
            if (mem_file) {
              MappingRecord mapping;
              mapping.LoadFromFile(mem_file);
              fclose(mem_file);
              
              FileRecord rec;
              rec.rid = rid;
              rec.mapping = mapping;
              rec.file_index = min_rec.file_index;
              merge_heap.push(rec);
            }
          }
        }
      }
    }
    
    // Close readers for this rid to free file descriptors before moving to next rid
    for (size_t fi : file_indices) {
      readers[fi].reset();
    }
  }

  // Output the last mapping if any
  if (num_total_mappings > 0) {
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
  }

  // Clean up overflow files
  for (const std::string& file_path : shared_overflow_file_paths_) {
    unlink(file_path.c_str());
  }
  shared_overflow_file_paths_.clear();

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

// Explicit template instantiations for overflow methods
// Only instantiate for types that have WriteToFile/LoadFromFile methods
template void MappingWriter<SAMMapping>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<SAMMapping>>&);
template void MappingWriter<SAMMapping>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t, uint32_t, const SequenceBatch&, const khash_t(k64_seq)*);

template void MappingWriter<PAFMapping>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<PAFMapping>>&);
template void MappingWriter<PAFMapping>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t, uint32_t, const SequenceBatch&, const khash_t(k64_seq)*);

template void MappingWriter<PairsMapping>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<PairsMapping>>&);
template void MappingWriter<PairsMapping>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t, uint32_t, const SequenceBatch&, const khash_t(k64_seq)*);

template void MappingWriter<PairedPAFMapping>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<PairedPAFMapping>>&);
template void MappingWriter<PairedPAFMapping>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t, uint32_t, const SequenceBatch&, const khash_t(k64_seq)*);

// Helper function to close thread-local writer and collect paths
template <typename MappingRecord>
void MappingWriter<MappingRecord>::CloseThreadOverflowWriter() {
  if (!tls_overflow_writer_) {
    return;
  }
  
  // Close the thread-local writer and get its file paths
  auto file_paths = tls_overflow_writer_->Close();
  tls_overflow_writer_.reset();
  
  // Add the paths to the shared collection under lock
  {
    std::lock_guard<std::mutex> lock(overflow_paths_mutex_);
    shared_overflow_file_paths_.insert(shared_overflow_file_paths_.end(),
                                       file_paths.begin(), file_paths.end());
  }
}

// Helper function to rotate thread-local writer after each flush
// Closes current writer (collecting file paths) and resets it so next spill creates fresh files
// This ensures each overflow file contains exactly one sorted run for correct k-way merge
template <typename MappingRecord>
void MappingWriter<MappingRecord>::RotateThreadOverflowWriter() {
  if (!tls_overflow_writer_) {
    return;
  }
  
  // Close the thread-local writer and get its file paths
  auto file_paths = tls_overflow_writer_->Close();
  tls_overflow_writer_.reset();  // Next spill will create new files
  
  // Add the paths to the shared collection under lock
  {
    std::lock_guard<std::mutex> lock(overflow_paths_mutex_);
    shared_overflow_file_paths_.insert(shared_overflow_file_paths_.end(),
                                       file_paths.begin(), file_paths.end());
  }
}

// MappingWithBarcode and MappingWithoutBarcode don't have WriteToFile/LoadFromFile methods
// so we provide fallback implementations that use the legacy temp file system
template <>
void MappingWriter<MappingWithBarcode>::OutputTempMappingsToOverflow(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingWithBarcode>> &mappings_on_diff_ref_seqs) {
  // Fallback to legacy temp file system for barcode types
  // This should not be called in practice since barcode types typically don't use low-memory mode
  // But we need the symbol to exist for linking
}

template <>
void MappingWriter<MappingWithBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  // Fallback - should not be called for barcode types
}

template <>
void MappingWriter<MappingWithoutBarcode>::OutputTempMappingsToOverflow(
    uint32_t num_reference_sequences,
    std::vector<std::vector<MappingWithoutBarcode>> &mappings_on_diff_ref_seqs) {
  // Fallback to legacy temp file system for barcode types
  // This should not be called in practice since barcode types typically don't use low-memory mode
  // But we need the symbol to exist for linking
}

template <>
void MappingWriter<MappingWithoutBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  // Fallback - should not be called for barcode types
}

template <>
void MappingWriter<PairedEndMappingWithBarcode>::OutputTempMappingsToOverflow(
    uint32_t num_reference_sequences,
    std::vector<std::vector<PairedEndMappingWithBarcode>> &mappings_on_diff_ref_seqs) {
  // Fallback - should not be called for paired-end barcode types
}

template <>
void MappingWriter<PairedEndMappingWithBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  // Fallback - should not be called for barcode types
}

template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::OutputTempMappingsToOverflow(
    uint32_t num_reference_sequences,
    std::vector<std::vector<PairedEndMappingWithoutBarcode>> &mappings_on_diff_ref_seqs) {
  // Fallback - should not be called for paired-end barcode types
}

template <>
void MappingWriter<PairedEndMappingWithoutBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  // Fallback - should not be called for barcode types
}

// Add explicit instantiation for CloseThreadOverflowWriter
template void MappingWriter<SAMMapping>::CloseThreadOverflowWriter();
template void MappingWriter<PAFMapping>::CloseThreadOverflowWriter();
template void MappingWriter<PairedPAFMapping>::CloseThreadOverflowWriter();
template void MappingWriter<PairsMapping>::CloseThreadOverflowWriter();
template void MappingWriter<MappingWithBarcode>::CloseThreadOverflowWriter();
template void MappingWriter<MappingWithoutBarcode>::CloseThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithBarcode>::CloseThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithoutBarcode>::CloseThreadOverflowWriter();

// Add explicit instantiation for RotateThreadOverflowWriter
template void MappingWriter<SAMMapping>::RotateThreadOverflowWriter();
template void MappingWriter<PAFMapping>::RotateThreadOverflowWriter();
template void MappingWriter<PairedPAFMapping>::RotateThreadOverflowWriter();
template void MappingWriter<PairsMapping>::RotateThreadOverflowWriter();
template void MappingWriter<MappingWithBarcode>::RotateThreadOverflowWriter();
template void MappingWriter<MappingWithoutBarcode>::RotateThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithBarcode>::RotateThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithoutBarcode>::RotateThreadOverflowWriter();

#endif  // LEGACY_OVERFLOW

// htslib helper methods for SAMMapping specialization
template <>
void MappingWriter<SAMMapping>::OpenHtsOutput() {
  if (mapping_parameters_.mapping_output_format != MAPPINGFORMAT_BAM &&
      mapping_parameters_.mapping_output_format != MAPPINGFORMAT_CRAM) {
    return;  // Not BAM/CRAM format
  }
  
  const char *output_path = mapping_parameters_.mapping_output_file_path.c_str();
  const char *hts_mode = (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) ? "wb" : "wc";
  
  hts_out_ = sam_open(output_path, hts_mode);
  if (!hts_out_) {
    ExitWithMessage("Failed to open BAM/CRAM output file: " + 
                    mapping_parameters_.mapping_output_file_path);
  }
  
  // Set compression threads
  int effective_hts_threads = mapping_parameters_.hts_threads;
  if (effective_hts_threads == 0) {
    effective_hts_threads = std::min(mapping_parameters_.num_threads, 4);
  }
  hts_set_threads(hts_out_, effective_hts_threads);
  
  // For CRAM, set reference and validate
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM &&
      !mapping_parameters_.reference_file_path.empty()) {
    int ret = hts_set_fai_filename(hts_out_, mapping_parameters_.reference_file_path.c_str());
    if (ret != 0) {
      ExitWithMessage("Failed to set CRAM reference file: " + 
                      mapping_parameters_.reference_file_path);
    }
  }
}

template <>
void MappingWriter<SAMMapping>::FinalizeSortedOutput() {
  if (!bam_sorter_) return;
  
  bam_sorter_->finalize();
  
  bam1_t* b;
  bool hasY;
  
  while (bam_sorter_->nextRecord(&b, &hasY)) {
    // b is a fully reconstructed bam1_t owned by the sorter
    // The sorter ensures all core fields (l_qname, n_cigar, l_qseq, etc.) are valid
    // and data buffer is properly allocated
    
    // Write to primary output (always)
    if (sam_write1(hts_out_, hts_hdr_, b) < 0) {
      ExitWithMessage("Failed to write sorted BAM/CRAM record");
    }
    
    // Route to Y/noY streams based on stored hasY flag
    // Note: hasY was computed at AppendMapping time using reads_with_y_hit_,
    // so we only need to check if the streams are open (not reads_with_y_hit_)
    if (hasY && Y_hts_out_) {
      if (sam_write1(Y_hts_out_, hts_hdr_, b) < 0) {
        ExitWithMessage("Failed to write sorted BAM/CRAM record to Y stream");
      }
    }
    if (!hasY && noY_hts_out_) {
      if (sam_write1(noY_hts_out_, hts_hdr_, b) < 0) {
        ExitWithMessage("Failed to write sorted BAM/CRAM record to noY stream");
      }
    }
    
    // Note: b is valid until next nextRecord() call; do not destroy here
  }
  
  // Cleanup sorter
  bam_sorter_.reset();
}

template <>
void MappingWriter<SAMMapping>::CloseHtsOutput() {
  // Finalize sorted output before closing (if sorting was enabled)
  if (bam_sorter_) {
    FinalizeSortedOutput();
  }
  
  if (hts_out_) {
    // Save index if requested and output is not stdout
    // Note: For CRAM, indexing is handled by cram_close() automatically
    if (mapping_parameters_.write_index &&
        mapping_parameters_.mapping_output_file_path != "-" &&
        mapping_parameters_.mapping_output_file_path != "/dev/stdout" &&
        mapping_parameters_.mapping_output_file_path != "/dev/stderr") {
      // For BAM, explicitly save index before close
      // For CRAM, let cram_close() handle indexing (sam_idx_save may cause issues)
      if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) {
        int ret = sam_idx_save(hts_out_);
        if (ret < 0) {
          std::cerr << "Warning: Failed to save BAM index (error code: " << ret << ")\n";
        }
      }
      // For CRAM, indexing is handled by cram_close() - don't call sam_idx_save
    }
    sam_close(hts_out_);
    hts_out_ = nullptr;
  }
  if (hts_hdr_) {
    bam_hdr_destroy(hts_hdr_);
    hts_hdr_ = nullptr;
  }
  // Note: noY_hts_out_ and Y_hts_out_ are closed by CloseYFilterStreams(),
  // which runs before this method. This is just a safety check in case
  // CloseYFilterStreams() wasn't called.
  if (noY_hts_out_) {
    // Save index only if requested and output is not stdout
    // For CRAM, indexing is handled by cram_close() - don't call sam_idx_save
    if (mapping_parameters_.write_index &&
        mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM &&
        mapping_parameters_.noY_output_path != "-" &&
        mapping_parameters_.noY_output_path != "/dev/stdout" &&
        mapping_parameters_.noY_output_path != "/dev/stderr") {
      int ret = sam_idx_save(noY_hts_out_);  // Ignore errors on cleanup
      (void)ret;
    }
    sam_close(noY_hts_out_);
    noY_hts_out_ = nullptr;
  }
  if (Y_hts_out_) {
    // Save index only if requested and output is not stdout
    // For CRAM, indexing is handled by cram_close() - don't call sam_idx_save
    if (mapping_parameters_.write_index &&
        mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM &&
        mapping_parameters_.Y_output_path != "-" &&
        mapping_parameters_.Y_output_path != "/dev/stdout" &&
        mapping_parameters_.Y_output_path != "/dev/stderr") {
      int ret = sam_idx_save(Y_hts_out_);  // Ignore errors on cleanup
      (void)ret;
    }
    sam_close(Y_hts_out_);
    Y_hts_out_ = nullptr;
  }
}

template <>
void MappingWriter<SAMMapping>::OpenYFilterStreams() {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM ||
      mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM) {
    // BAM/CRAM path: use htslib
    const char *hts_mode = (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) ? "wb" : "wc";
    int effective_hts_threads = mapping_parameters_.hts_threads;
    if (effective_hts_threads == 0) {
      effective_hts_threads = std::min(mapping_parameters_.num_threads, 4);
    }
    
    if (mapping_parameters_.emit_noY_stream && !noY_hts_out_) {
      noY_hts_out_ = sam_open(mapping_parameters_.noY_output_path.c_str(), hts_mode);
      if (!noY_hts_out_) {
        ExitWithMessage("Failed to open noY BAM/CRAM output file: " + 
                        mapping_parameters_.noY_output_path);
      }
      hts_set_threads(noY_hts_out_, effective_hts_threads);
      if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM &&
          !mapping_parameters_.reference_file_path.empty()) {
        int ret = hts_set_fai_filename(noY_hts_out_, mapping_parameters_.reference_file_path.c_str());
        if (ret != 0) {
          ExitWithMessage("Failed to set CRAM reference file for noY stream: " + 
                          mapping_parameters_.reference_file_path);
        }
      }
      // Indexing will be initialized in OutputHeader after header is written
    }
    
    if (mapping_parameters_.emit_Y_stream && !Y_hts_out_) {
      Y_hts_out_ = sam_open(mapping_parameters_.Y_output_path.c_str(), hts_mode);
      if (!Y_hts_out_) {
        ExitWithMessage("Failed to open Y BAM/CRAM output file: " + 
                        mapping_parameters_.Y_output_path);
      }
      hts_set_threads(Y_hts_out_, effective_hts_threads);
      if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM &&
          !mapping_parameters_.reference_file_path.empty()) {
        int ret = hts_set_fai_filename(Y_hts_out_, mapping_parameters_.reference_file_path.c_str());
        if (ret != 0) {
          ExitWithMessage("Failed to set CRAM reference file for Y stream: " + 
                          mapping_parameters_.reference_file_path);
        }
      }
      // Indexing will be initialized in OutputHeader after header is written
    }
  } else {
    // SAM text path: use FILE*
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
}

template <>
void MappingWriter<SAMMapping>::CloseYFilterStreams() {
  // Close htslib streams
  if (noY_hts_out_) {
    // Save index only if requested and output is not stdout
    // For CRAM, indexing is handled by cram_close() - don't call sam_idx_save
    if (mapping_parameters_.write_index &&
        mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM &&
        mapping_parameters_.noY_output_path != "-" &&
        mapping_parameters_.noY_output_path != "/dev/stdout" &&
        mapping_parameters_.noY_output_path != "/dev/stderr") {
      int ret = sam_idx_save(noY_hts_out_);  // Ignore errors on cleanup
      (void)ret;
    }
    sam_close(noY_hts_out_);
    noY_hts_out_ = nullptr;
  }
  if (Y_hts_out_) {
    // Save index only if requested and output is not stdout
    // For CRAM, indexing is handled by cram_close() - don't call sam_idx_save
    if (mapping_parameters_.write_index &&
        mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM &&
        mapping_parameters_.Y_output_path != "-" &&
        mapping_parameters_.Y_output_path != "/dev/stdout" &&
        mapping_parameters_.Y_output_path != "/dev/stderr") {
      int ret = sam_idx_save(Y_hts_out_);  // Ignore errors on cleanup
      (void)ret;
    }
    sam_close(Y_hts_out_);
    Y_hts_out_ = nullptr;
  }
  
  // Close FILE* streams
  if (noY_output_file_) {
    fclose(noY_output_file_);
    noY_output_file_ = nullptr;
  }
  if (Y_output_file_) {
    fclose(Y_output_file_);
    Y_output_file_ = nullptr;
  }
}

template <>
std::string MappingWriter<SAMMapping>::DeriveReadGroupFromFilename(
    const std::string &filename) {
  // Extract sample name from filename (e.g., "sample1_R1.fastq.gz" -> "sample1")
  size_t last_slash = filename.find_last_of("/\\");
  std::string basename = (last_slash == std::string::npos) ? filename : filename.substr(last_slash + 1);
  
  // Remove common suffixes
  std::vector<std::string> suffixes = {".fastq.gz", ".fq.gz", ".fastq", ".fq", ".gz"};
  for (const auto &suffix : suffixes) {
    if (basename.length() >= suffix.length() &&
        basename.substr(basename.length() - suffix.length()) == suffix) {
      basename = basename.substr(0, basename.length() - suffix.length());
      break;
    }
  }
  
  // Remove _R1/_R2 suffix
  size_t r_pos = basename.find("_R");
  if (r_pos != std::string::npos) {
    basename = basename.substr(0, r_pos);
  }
  
  return basename.empty() ? "default" : basename;
}

template <>
void MappingWriter<SAMMapping>::BuildHtsHeader(uint32_t num_ref_seqs,
                                                const SequenceBatch &reference) {
  std::string header_text = "@HD\tVN:1.6";
  
  // Set sort order based on --sort-bam flag
  if (mapping_parameters_.sort_bam) {
    header_text += "\tSO:coordinate";
  } else {
    header_text += "\tSO:unknown";
  }
  header_text += "\n";
  
  // @SQ lines
  for (uint32_t rid = 0; rid < num_ref_seqs; ++rid) {
    header_text += "@SQ\tSN:" + std::string(reference.GetSequenceNameAt(rid)) +
                   "\tLN:" + std::to_string(reference.GetSequenceLengthAt(rid)) + "\n";
  }
  
  // @PG line
  header_text += "@PG\tID:chromap\tPN:chromap\tVN:" + std::string(CHROMAP_VERSION) + "\n";
  
  // @RG line (if read group set)
  if (!mapping_parameters_.read_group_id.empty()) {
    std::string rg_id = mapping_parameters_.read_group_id;
    if (rg_id == "auto") {
      if (!mapping_parameters_.read_file1_paths.empty()) {
        rg_id = DeriveReadGroupFromFilename(mapping_parameters_.read_file1_paths[0]);
      } else {
        rg_id = "default";
      }
    }
    header_text += "@RG\tID:" + rg_id + "\tSM:" + rg_id + "\n";
  }
  
  hts_hdr_ = sam_hdr_parse(header_text.length(), header_text.c_str());
  if (!hts_hdr_) {
    ExitWithMessage("Failed to parse BAM/CRAM header");
  }
  
  // Write header to output
  if (sam_hdr_write(hts_out_, hts_hdr_) < 0) {
    ExitWithMessage("Failed to write BAM/CRAM header");
  }
  
  // For CRAM, validate that all header sequences exist in reference file
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM &&
      !mapping_parameters_.reference_file_path.empty()) {
    faidx_t *fai = fai_load3(mapping_parameters_.reference_file_path.c_str(), NULL, NULL, FAI_CREATE);
    if (!fai) {
      ExitWithMessage("Failed to load FASTA index for CRAM reference: " + 
                      mapping_parameters_.reference_file_path + 
                      ". Run: samtools faidx " + mapping_parameters_.reference_file_path);
    }
    
    std::vector<std::string> missing_contigs;
    for (int i = 0; i < hts_hdr_->n_targets; ++i) {
      const char *seq_name = hts_hdr_->target_name[i];
      if (!faidx_has_seq(fai, seq_name)) {
        missing_contigs.push_back(std::string(seq_name));
      }
    }
    
    fai_destroy(fai);
    
    if (!missing_contigs.empty()) {
      std::string error_msg = "CRAM reference file is missing " + 
                              std::to_string(missing_contigs.size()) + 
                              " contig(s) present in header:\n";
      for (size_t j = 0; j < missing_contigs.size() && j < 10; ++j) {
        error_msg += "  " + missing_contigs[j] + "\n";
      }
      if (missing_contigs.size() > 10) {
        error_msg += "  ... and " + std::to_string(missing_contigs.size() - 10) + " more\n";
      }
      error_msg += "Ensure the reference FASTA matches the index used for mapping.";
      ExitWithMessage(error_msg);
    }
  }
  
  // Initialize indexing if requested (after header is written)
  // Store index path in persistent member to avoid dangling pointer
  if (mapping_parameters_.write_index) {
    this->hts_index_path_ = mapping_parameters_.mapping_output_file_path;
    // Use BAI format for BAM (min_shift=0), CRAI for CRAM
    // BAI is more widely compatible with existing tools (samtools, IGV, etc.)
    // Note: min_shift=0 creates BAI, min_shift>0 creates CSI
    int min_shift = 0;  // 0 = BAI format (most compatible)
    if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) {
      this->hts_index_path_ += ".bai";  // BAI format for BAM (widely compatible)
    } else if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM) {
      this->hts_index_path_ += ".crai";  // CRAM uses .crai
    }
    int ret = sam_idx_init(hts_out_, hts_hdr_, min_shift, this->hts_index_path_.c_str());
    if (ret < 0) {
      ExitWithMessage("Failed to initialize BAM/CRAM index");
    }
  }
  
  // Initialize BamSorter if --sort-bam is enabled
  if (mapping_parameters_.sort_bam) {
    std::string tmpDir = mapping_parameters_.temp_directory_path;
    if (tmpDir.empty()) {
      // Use system temp directory
      const char* env_tmp = std::getenv("TMPDIR");
      if (env_tmp) {
        tmpDir = env_tmp;
      } else {
        tmpDir = "/tmp";
      }
    }
    bam_sorter_.reset(new BamSorter(mapping_parameters_.sort_bam_ram_limit, tmpDir));
  }
}

template <>
bam1_t* MappingWriter<SAMMapping>::ConvertToHtsBam(uint32_t rid,
                                                    const SequenceBatch &reference,
                                                    const SAMMapping &mapping) {
  bam1_t *b = bam_init1();
  if (!b) {
    ExitWithMessage("Failed to allocate bam1_t");
  }
  
  // Handle unmapped reads
  if (mapping.flag_ & BAM_FUNMAP) {
    b->core.tid = -1;
    b->core.pos = -1;
    b->core.qual = 0;
  } else {
    b->core.tid = rid;
    b->core.pos = mapping.pos_;  // already 0-based
    b->core.qual = mapping.mapq_;
  }
  
  // Handle unmapped mate
  if (mapping.mrid_ < 0 || (mapping.flag_ & BAM_FMUNMAP)) {
    b->core.mtid = -1;
    b->core.mpos = -1;
  } else {
    b->core.mtid = mapping.mrid_;
    b->core.mpos = mapping.mpos_;  // already 0-based
  }
  
  b->core.flag = mapping.flag_;
  b->core.isize = mapping.tlen_;
  b->core.n_cigar = (mapping.flag_ & BAM_FUNMAP) ? 0 : mapping.n_cigar_;
  
  // Set l_qname to actual length (including NUL), NOT padded length
  // BAM spec: l_qname is the length of the read name including trailing NUL
  b->core.l_qname = mapping.read_name_.length() + 1;  // Actual length including NUL
  
  // Set l_qseq: for unmapped reads with empty sequence, use 0
  // Otherwise use actual sequence length
  if ((mapping.flag_ & BAM_FUNMAP) && mapping.sequence_.empty()) {
    b->core.l_qseq = 0;
  } else {
    b->core.l_qseq = mapping.sequence_.length();
  }
  
  // Set bin for mapped reads (required for indexing)
  // BAM uses binning scheme with min_shift=14, n_lvls=5
  if (!(mapping.flag_ & BAM_FUNMAP)) {
    // Calculate end position: pos + alignment length on reference
    int32_t end = mapping.pos_ + mapping.GetAlignmentLength();
    // Ensure non-empty interval: hts_reg2bin requires end > pos
    if (end <= mapping.pos_) {
      end = mapping.pos_ + 1;
    }
    b->core.bin = hts_reg2bin(mapping.pos_, end, 14, 5);
  } else {
    b->core.bin = 0;  // Unmapped reads have bin = 0
  }
  
  // Calculate total data length (NO padding on qname - BAM expects CIGAR immediately after l_qname)
  int data_len = b->core.l_qname +                      // read name (no padding)
                 b->core.n_cigar * 4 +                   // CIGAR (4 bytes per op)
                 (b->core.l_qseq + 1) / 2 +              // sequence (4-bit packed, 2 bases per byte)
                 b->core.l_qseq;                         // qualities (1 byte per base)
  
  // Allocate data buffer
  if (static_cast<int>(b->m_data) < data_len) {
    b->data = (uint8_t*)realloc(b->data, data_len);
    if (!b->data) {
      bam_destroy1(b);
      ExitWithMessage("Failed to allocate bam1_t data buffer");
    }
    b->m_data = data_len;
  }
  b->l_data = data_len;
  
  uint8_t *p = b->data;
  
  // Pack read name (null-terminated, NO padding - CIGAR follows immediately)
  memcpy(p, mapping.read_name_.c_str(), mapping.read_name_.length() + 1);
  p += b->core.l_qname;  // Advance by actual length (no padding)
  
  // Pack CIGAR
  if (b->core.n_cigar > 0) {
    memcpy(p, mapping.cigar_.data(), b->core.n_cigar * 4);
    p += b->core.n_cigar * 4;
    
    // Sanity check: verify CIGAR consistency with sequence length
    hts_pos_t qlen = bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b));
    if (qlen != static_cast<hts_pos_t>(b->core.l_qseq)) {
      std::cerr << "WARNING: CIGAR/sequence length mismatch for read "
                << mapping.read_name_ << " at " << mapping.pos_
                << ": CIGAR qlen=" << qlen << " vs l_qseq=" << b->core.l_qseq
                << ". Writing as unmapped to avoid indexing issues.\n";
      std::vector<uint8_t> raw_qual(mapping.sequence_.size(), 0xFF);
      for (size_t qi = 0; qi < raw_qual.size(); ++qi) {
        if (qi < mapping.sequence_qual_.size() &&
            mapping.sequence_qual_[qi] >= 33) {
          int qual_val = mapping.sequence_qual_[qi] - 33;
          if (qual_val < 0) qual_val = 0;
          if (qual_val > 93) qual_val = 93;
          raw_qual[qi] = static_cast<uint8_t>(qual_val);
        }
      }
      const char *seq_ptr = mapping.sequence_.empty() ? nullptr : mapping.sequence_.c_str();
      const char *qual_ptr =
          raw_qual.empty() ? nullptr : reinterpret_cast<const char*>(raw_qual.data());
      uint16_t unmapped_flag =
          static_cast<uint16_t>(mapping.flag_ | BAM_FUNMAP | BAM_FMUNMAP);
      if (bam_set1(b,
                   mapping.read_name_.length(), mapping.read_name_.c_str(),
                   unmapped_flag, -1, -1, 0,
                   0, nullptr,
                   -1, -1, 0,
                   mapping.sequence_.size(), seq_ptr, qual_ptr,
                   0) < 0) {
        bam_destroy1(b);
        ExitWithMessage("Failed to rebuild BAM record after CIGAR/sequence mismatch");
      }
      // Re-add aux tags (RG, CB) for unmapped fallback record
      // These are important for downstream analysis even if read is unmapped
      if (cell_barcode_length_ > 0) {
        std::string cb = barcode_translator_.Translate(mapping.cell_barcode_, cell_barcode_length_);
        bam_aux_append(b, "CB", 'Z', cb.length() + 1, (const uint8_t*)cb.c_str());
      }
      if (!mapping_parameters_.read_group_id.empty()) {
        std::string rg_id = mapping_parameters_.read_group_id;
        if (rg_id == "auto") {
          if (!mapping_parameters_.read_file1_paths.empty()) {
            rg_id = DeriveReadGroupFromFilename(mapping_parameters_.read_file1_paths[0]);
          } else {
            rg_id = "default";
          }
        }
        bam_aux_append(b, "RG", 'Z', rg_id.length() + 1, (const uint8_t*)rg_id.c_str());
      }
      return b;
    }
  }
  
  // Pack sequence (4-bit encoding)
  // Helper function to convert base to 4-bit encoding (A=1, C=2, G=4, T=8, N=15)
  auto base_to_nt16 = [](char c) -> uint8_t {
    switch (c) {
      case 'A': case 'a': return 1;
      case 'C': case 'c': return 2;
      case 'G': case 'g': return 4;
      case 'T': case 't': return 8;
      case 'N': case 'n': return 15;
      default: return 15;  // Unknown bases as N
    }
  };
  for (int i = 0; i < b->core.l_qseq; i += 2) {
    uint8_t base1 = base_to_nt16(mapping.sequence_[i]);
    uint8_t base2 = (i + 1 < b->core.l_qseq) ? base_to_nt16(mapping.sequence_[i+1]) : 0;
    *p++ = (base1 << 4) | base2;
  }
  
  // Pack qualities (convert from ASCII Phred+33 to raw Phred)
  // Handle missing or short quality strings (e.g., FASTA input)
  // Clamp to valid range [0, 93] to prevent underflow/overflow
  for (int i = 0; i < b->core.l_qseq; ++i) {
    uint8_t raw_qual;
    if (i < static_cast<int>(mapping.sequence_qual_.size()) && 
        mapping.sequence_qual_[i] >= 33) {
      // Valid quality character: convert from ASCII Phred+33 to raw Phred
      int qual_val = mapping.sequence_qual_[i] - 33;
      // Clamp to valid Phred range [0, 93]
      if (qual_val < 0) qual_val = 0;
      if (qual_val > 93) qual_val = 93;
      raw_qual = static_cast<uint8_t>(qual_val);
    } else {
      // Missing quality: use 0xFF (BAM missing quality indicator)
      raw_qual = 0xFF;
    }
    *p++ = raw_qual;
  }
  
  // Append aux tags
  // NM_ is a bit-field, so copy to a variable before taking address
  int32_t nm_value = static_cast<int32_t>(mapping.NM_);
  bam_aux_append(b, "NM", 'i', sizeof(int32_t), (const uint8_t*)&nm_value);
  if (!mapping.MD_.empty()) {
    bam_aux_append(b, "MD", 'Z', mapping.MD_.length() + 1, (const uint8_t*)mapping.MD_.c_str());
  }
  if (cell_barcode_length_ > 0) {
    std::string cb = barcode_translator_.Translate(mapping.cell_barcode_, cell_barcode_length_);
    bam_aux_append(b, "CB", 'Z', cb.length() + 1, (const uint8_t*)cb.c_str());
  }
  if (!mapping_parameters_.read_group_id.empty()) {
    std::string rg_id = mapping_parameters_.read_group_id;
    if (rg_id == "auto") {
      if (!mapping_parameters_.read_file1_paths.empty()) {
        rg_id = DeriveReadGroupFromFilename(mapping_parameters_.read_file1_paths[0]);
      } else {
        rg_id = "default";
      }
    }
    bam_aux_append(b, "RG", 'Z', rg_id.length() + 1, (const uint8_t*)rg_id.c_str());
  }
  
  return b;
}

}  // namespace chromap
