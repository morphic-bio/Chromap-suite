#include "mapping_writer.h"

#include <map>
#include <queue>
#include <algorithm>
#include <atomic>
#include <unordered_set>
#include <unordered_map>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sstream>
#include "chromap.h"
#include "atac_hot_spill.h"
#include "atac_kway_spill.h"
#include "atac_lowmem_finalize.h"
#include "bam_sorter.h"
#include "rapidmacs/frag_compact_store.h"
#include "rapidmacs/fragments.h"

namespace chromap {

#ifndef LEGACY_OVERFLOW
static std::atomic<uint32_t> g_low_mem_mid_batch_overflow_flushes{0};

void RecordLowMemMidBatchOverflowFlush() {
  g_low_mem_mid_batch_overflow_flushes.fetch_add(1, std::memory_order_relaxed);
}

void ResetLowMemMidBatchOverflowFlushCount() {
  g_low_mem_mid_batch_overflow_flushes.store(0, std::memory_order_relaxed);
}

uint32_t LowMemMidBatchOverflowFlushCount() {
  return g_low_mem_mid_batch_overflow_flushes.load(std::memory_order_relaxed);
}
#endif

namespace {

// A low-memory spill file that cannot be opened stops the run: skipping it
// would drop its records without an error. The spill files are temporary, so
// they are removed before exiting. files_for_reference is the number of spill
// files the merge keeps open at once for the reference being merged (0 while
// scanning).
std::string UnopenableOverflowFileMessage(const std::string &path,
                                          int open_errno, uint32_t rid,
                                          size_t files_for_reference) {
  std::ostringstream message;
  message << "Cannot open low-memory overflow file " << path << ": "
          << (open_errno != 0 ? std::strerror(open_errno) : "unknown error");
  if (open_errno == EMFILE || open_errno == ENFILE) {
    std::string limit_text = "unknown";
    struct rlimit limit;
    if (getrlimit(RLIMIT_NOFILE, &limit) == 0) {
      limit_text = limit.rlim_cur == RLIM_INFINITY
                       ? std::string("unlimited")
                       : std::to_string(
                             static_cast<unsigned long long>(limit.rlim_cur));
    }
    message << ". The low-memory merge opens all spill files of a reference "
               "at once";
    if (files_for_reference > 0) {
      message << " (" << files_for_reference << " for reference " << rid
              << ")";
    }
    message << ", and the open-file limit is " << limit_text
            << ". Raise the limit (ulimit -n) or use a larger --low-mem-ram "
               "so that fewer spill files are written";
  }
  return message.str();
}

void ExitOnUnopenableOverflowFile(
    const std::vector<std::string> &overflow_file_paths,
    const std::string &path, int open_errno, uint32_t rid,
    size_t files_for_reference) {
  const std::string message = UnopenableOverflowFileMessage(
      path, open_errno, rid, files_for_reference);
  for (const std::string &spill_path : overflow_file_paths) {
    unlink(spill_path.c_str());
  }
  ExitWithMessage(message);
}

std::string DeriveReadGroupFromFilenameImpl(const std::string &filename) {
  size_t last_slash = filename.find_last_of("/\\");
  std::string basename =
      (last_slash == std::string::npos) ? filename
                                        : filename.substr(last_slash + 1);

  std::vector<std::string> suffixes = {".fastq.gz", ".fq.gz", ".fastq", ".fq",
                                       ".gz"};
  for (const auto &suffix : suffixes) {
    if (basename.length() >= suffix.length() &&
        basename.substr(basename.length() - suffix.length()) == suffix) {
      basename = basename.substr(0, basename.length() - suffix.length());
      break;
    }
  }

  size_t r_pos = basename.find("_R");
  if (r_pos != std::string::npos) {
    basename = basename.substr(0, r_pos);
  }

  return basename.empty() ? "default" : basename;
}

bam1_t *BuildBamRecordFromSamMappingFields(
    const MappingParameters &mapping_parameters,
    uint32_t cell_barcode_length, BarcodeTranslator &barcode_translator,
    uint32_t rid, const SequenceBatch &reference, const SAMMapping &mapping) {
  // `rid` is the per-chromosome bucket used when batching paired-end results
  // (e.g. read1's reference for BED/fragment rows). BAM @SQ index (tid) must
  // come from this read's own mapping.rid_ so read2 lands on the correct
  // contig for discordant / interchromosomal pairs.
  (void)rid;
  bam1_t *b = bam_init1();
  if (!b) {
    ExitWithMessage("Failed to allocate bam1_t");
  }

  if (mapping.n_cigar_ < 0 ||
      static_cast<size_t>(mapping.n_cigar_) > mapping.cigar_.size()) {
    bam_destroy1(b);
    ExitWithMessage("Invalid SAMMapping CIGAR state while writing BAM for read " +
                    mapping.read_name_ + ": n_cigar=" +
                    std::to_string(mapping.n_cigar_) + ", cigar_size=" +
                    std::to_string(mapping.cigar_.size()));
  }

  const size_t qname_len = mapping.read_name_.length();
  if (qname_len > 254) {
    bam_destroy1(b);
    ExitWithMessage("Read name is too long for BAM qname field while writing " +
                    mapping.read_name_);
  }

  if (mapping.sequence_.length() >
      static_cast<size_t>(std::numeric_limits<int32_t>::max())) {
    bam_destroy1(b);
    ExitWithMessage("Read sequence is too long for BAM while writing " +
                    mapping.read_name_);
  }

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

  uint16_t flag = static_cast<uint16_t>(mapping.flag_);
  int32_t tid = static_cast<int32_t>(mapping.rid_);
  hts_pos_t pos = mapping.pos_;
  uint8_t mapq = mapping.mapq_;
  size_t n_cigar = static_cast<size_t>(mapping.n_cigar_);
  const uint32_t *cigar = n_cigar == 0 ? nullptr : mapping.cigar_.data();
  int32_t mtid = mapping.mrid_;
  hts_pos_t mpos = mapping.mpos_;

  if (mapping.flag_ & BAM_FUNMAP) {
    tid = -1;
    pos = -1;
    mapq = 0;
    n_cigar = 0;
    cigar = nullptr;
  }
  if (mapping.mrid_ < 0 || (mapping.flag_ & BAM_FMUNMAP)) {
    mtid = -1;
    mpos = -1;
  }

  const char *seq_ptr =
      mapping.sequence_.empty() ? nullptr : mapping.sequence_.c_str();
  const char *qual_ptr =
      raw_qual.empty() ? nullptr : reinterpret_cast<const char *>(raw_qual.data());

  if (bam_set1(b, qname_len, mapping.read_name_.c_str(), flag, tid, pos, mapq,
               n_cigar, cigar, mtid, mpos, mapping.tlen_,
               mapping.sequence_.size(), seq_ptr, qual_ptr, 0) < 0) {
    bam_destroy1(b);
    ExitWithMessage("Failed to build BAM record for read " + mapping.read_name_);
  }

  int32_t nm_value = static_cast<int32_t>(mapping.NM_);
  bam_aux_append(b, "NM", 'i', sizeof(int32_t), (const uint8_t *)&nm_value);
  if (!mapping.MD_.empty()) {
    bam_aux_append(b, "MD", 'Z', mapping.MD_.length() + 1,
                   (const uint8_t *)mapping.MD_.c_str());
  }
  if (cell_barcode_length > 0) {
    std::string cb =
        barcode_translator.Translate(mapping.cell_barcode_, cell_barcode_length);
    bam_aux_append(b, "CB", 'Z', cb.length() + 1, (const uint8_t *)cb.c_str());
  }
  if (!mapping_parameters.read_group_id.empty()) {
    std::string rg_id = mapping_parameters.read_group_id;
    if (rg_id == "auto") {
      if (!mapping_parameters.ReadGroupSourcePath().empty()) {
        rg_id = DeriveReadGroupFromFilenameImpl(
            mapping_parameters.ReadGroupSourcePath());
      } else {
        rg_id = "default";
      }
    }
    bam_aux_append(b, "RG", 'Z', rg_id.length() + 1,
                   (const uint8_t *)rg_id.c_str());
  }

  (void)reference;
  return b;
}

}  // namespace

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
    uint32_t num_reference_sequences, const SequenceBatch &reference) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    if (mapping_parameters_.macs3_frag_chrom_names) {
      auto& names = *mapping_parameters_.macs3_frag_chrom_names;
      names.clear();
      names.reserve(num_reference_sequences);
      for (uint32_t i = 0; i < num_reference_sequences; ++i) {
        names.emplace_back(reference.GetSequenceNameAt(i));
      }
    }
    if (mapping_parameters_.macs3_frag_buffer) {
      mapping_parameters_.macs3_frag_buffer->assign(
          num_reference_sequences, std::vector<macs3::FragmentRecord>());
    }
  }
}

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
    if (mapping_parameters_.macs3_frag_buffer) {
      macs3::FragmentRecord rec;
      rec.chrom_id = static_cast<int32_t>(rid);
      rec.start = static_cast<int32_t>(mapping.GetStartPosition());
      rec.end = static_cast<int32_t>(mapping_end_position);
      // Match the emitted 7-column bulk FRAG row: column 7 is the collapsed
      // duplicate count consumed by standalone MACS3 and RapidMACS file input.
      rec.count = static_cast<uint32_t>(mapping.num_dups_);
      if (rec.end > rec.start && rec.count > 0) {
        auto& buckets = *mapping_parameters_.macs3_frag_buffer;
        if (rid >= buckets.size()) {
          buckets.resize(rid + 1);
        }
        buckets[rid].push_back(rec);
      }
    }
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
    uint32_t num_reference_sequences, const SequenceBatch &reference) {
  // Populate the macs3 chrom-name table + pre-size per-chrom buckets when
  // peak calling is requested with BED-only output (MAPPINGFORMAT_BED is
  // the BED-fragments-only path; TagAlign and other 6-col modes don't
  // emit fragments). Mirrors the dual writer's OutputHeader.
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    if (mapping_parameters_.macs3_frag_chrom_names) {
      auto& names = *mapping_parameters_.macs3_frag_chrom_names;
      names.clear();
      names.reserve(num_reference_sequences);
      for (uint32_t i = 0; i < num_reference_sequences; ++i) {
        names.emplace_back(reference.GetSequenceNameAt(i));
      }
    }
    if (mapping_parameters_.macs3_frag_buffer) {
      mapping_parameters_.macs3_frag_buffer->assign(
          num_reference_sequences, std::vector<macs3::FragmentRecord>());
    }
  }
}

template <>
void MappingWriter<PairedEndMappingWithBarcode>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const PairedEndMappingWithBarcode &mapping) {
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
    std::string strand = mapping.IsPositiveStrand() ? "+" : "-";
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    uint32_t mapping_end_position = mapping.GetEndPosition();
    const std::string translated_barcode = barcode_translator_.Translate(
        mapping.cell_barcode_, cell_barcode_length_);
    this->AppendMappingOutput(std::string(reference_sequence_name) + "\t" +
                              std::to_string(mapping.GetStartPosition()) +
                              "\t" + std::to_string(mapping_end_position) +
                              "\t" +
                              translated_barcode +
                              "\t" + std::to_string(mapping.num_dups_) + "\n");
    // BED-only peak calling: capture per-chrom (start, end, count) for the
    // events workspace (same shape as the dual writer's bucket capture).
    if (mapping_parameters_.macs3_frag_buffer) {
      macs3::FragmentRecord rec;
      rec.chrom_id = static_cast<int32_t>(rid);
      rec.start = static_cast<int32_t>(mapping.GetStartPosition());
      rec.end = static_cast<int32_t>(mapping_end_position);
      rec.count = static_cast<uint32_t>(mapping.num_dups_);
      if (rec.end > rec.start && rec.count > 0) {
        auto& buckets = *mapping_parameters_.macs3_frag_buffer;
        if (rid >= buckets.size()) {
          buckets.resize(rid + 1);
        }
        buckets[rid].push_back(rec);
      }
    }
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

// ---------------------------------------------------------------------------
// Per-reference low-memory finalisation of paired-end ATAC spills.
//
// Every spill file holds one reference, and duplicates never cross
// references, so each reference is merged and deduplicated by its own task.
// A task writes its output (AEV1 records, or BED/TagAlign rows) to a
// temporary partition file and its in-memory MACS3 fragments to its own
// bucket. The partitions are appended to the output in reference order.
//
// The serial merge carries its group state from one reference into the next
// and emits a reference's last group only when the next reference starts (or,
// for the last reference, after the loop, with the MAPQ test before the
// best-duplicate choice). Each task therefore returns its last group as a
// tail, and the assembly emits the tails in the serial merge's positions.
//
// Summary counts: while tasks run, no row is inserted into the summary table.
// Counts for barcodes that already have a row are added atomically; the other
// barcodes are logged per task in the order the serial merge would first
// update them and are inserted during the assembly. Together with the pending
// resize applied before the tasks start (as the first serial put would), this
// keeps the table layout, and so the summary row order, the same.
// ---------------------------------------------------------------------------

// Where one non-dual ATAC fragment is written.
struct AtacFragmentSink {
  FILE *evidence_fp = nullptr;        // AEV1 sidecar (sidecar-only runs)
  uint64_t *evidence_records = nullptr;
  FILE *text_fp = nullptr;            // BED or TagAlign rows
  std::vector<std::vector<macs3::FragmentRecord>> *buckets = nullptr;
  bool allow_bucket_resize = false;
  bool evidence_write_failed = false;
  bool text_write_failed = false;
  bool bucket_out_of_range = false;
};

#ifndef LEGACY_OVERFLOW
namespace lowmem_finalize_detail {

struct Counters {
  uint64_t uni = 0;
  uint64_t multi = 0;
  uint64_t passing = 0;
};

template <typename Record>
struct GroupState {
  GroupState() {
    // The mapping constructors leave most fields unset.
    last_mapping.read_id_ = 0;
    last_mapping.cell_barcode_ = 0;
    last_mapping.fragment_start_position_ = 0;
    last_mapping.fragment_length_ = 0;
    last_mapping.mapq_ = 0;
    last_mapping.direction_ = 0;
    last_mapping.is_unique_ = 0;
    last_mapping.num_dups_ = 0;
    last_mapping.positive_alignment_length_ = 0;
    last_mapping.negative_alignment_length_ = 0;
  }
  bool active = false;
  Record last_mapping;
  uint32_t num_last_mapping_dups = 0;
  std::vector<Record> bulk_dups;
};

struct Job {
  uint32_t rid = 0;
  std::vector<size_t> file_indices;  // ascending global spill-file indices
  uint64_t spill_bytes = 0;
};

template <typename Record>
struct TaskResult {
  bool done = false;
  std::string error;
  Counters counters;
  GroupState<Record> tail;
  std::string partition_path;
  uint64_t partition_records = 0;
  uint64_t records_merged = 0;
  AtacSummaryDeltaLog new_keys;
};

// Full decode: the payload string and AtacSpillRecord the serial merge uses.
struct FullDecodePolicy {
  typedef AtacSpillRecord Record;
  // 1: a record; 0: end of file; -1: error.
  static int Next(OverflowReader *reader, uint32_t expected_rid,
                  uint16_t schema, Record *record, std::string *payload,
                  std::string *error) {
    uint32_t rid = 0;
    if (!reader->ReadNext(rid, *payload)) {
      return 0;
    }
    if (rid != expected_rid) {
      *error = "Low-memory spill file " + reader->GetPath() +
               " holds a record for another reference";
      return -1;
    }
    if (!DecodeAtacKwaySpillRecord(payload->data(), payload->size(), schema,
                                   record, error)) {
      return -1;
    }
    return 1;
  }
};

// Lean decode: the fixed record header straight into the fragment fields;
// no payload string and no SAMMapping members.
struct LeanDecodePolicy {
  typedef PairedEndMappingWithBarcode Record;
  static int Next(OverflowReader *reader, uint32_t /*expected_rid*/,
                  uint16_t schema, Record *record, std::string * /*payload*/,
                  std::string *error) {
    AtacKwaySpillRecordHeaderV1 header;
    const int got = reader->ReadNextAtacRecordHeader(&header, error);
    if (got <= 0) {
      return got;
    }
    if (!DecodeAtacKwaySpillRecordLean(header, schema, record, error)) {
      return -1;
    }
    return 1;
  }
};

}  // namespace lowmem_finalize_detail
#endif  // LEGACY_OVERFLOW

class AtacLowMemFinalizer {
 public:
  explicit AtacLowMemFinalizer(MappingWriter<AtacSpillRecord> &writer)
      : w_(writer), p_(writer.mapping_parameters_) {
    summary_enabled_ = !p_.summary_metadata_file_path.empty();
    dedup_bulk_ = p_.remove_pcr_duplicates && !p_.is_bulk_data &&
                  p_.remove_pcr_duplicates_at_bulk_level;
  }

  // Output formats and modes the per-reference path covers.
  static bool ParametersEligible(const MappingParameters &p) {
    return !p.AtacDualFragmentAndBam() && !p.CreatesMergeableAtacSpill() &&
           !p.is_bulk_data &&
           (p.mapping_output_format == MAPPINGFORMAT_BED ||
            p.mapping_output_format == MAPPINGFORMAT_TAGALIGN);
  }

  // The non-dual branch of MappingWriter<AtacSpillRecord>::AppendMapping,
  // writing to `sink`.
  static void WriteNonDualFragment(MappingWriter<AtacSpillRecord> &w,
                                   AtacFragmentSink *sink, uint32_t rid,
                                   const SequenceBatch &reference,
                                   const PairedEndMappingWithBarcode &mwb);

#ifndef LEGACY_OVERFLOW
  // Returns false, before writing anything, when the spill files are not
  // ATAC k-way files of the covered schema; the caller then runs the serial
  // merge.
  bool Run(uint32_t num_reference_sequences, const SequenceBatch &reference,
           const khash_t(k64_seq) * barcode_whitelist_lookup_table,
           int requested_workers);
#endif  // LEGACY_OVERFLOW

 private:
  static void WriteText(AtacFragmentSink *sink, const std::string &line) {
    if (sink->text_fp != nullptr &&
        fwrite(line.data(), 1, line.size(), sink->text_fp) != line.size()) {
      sink->text_write_failed = true;
    }
  }

#ifndef LEGACY_OVERFLOW
  typedef lowmem_finalize_detail::Counters Counters;
  typedef lowmem_finalize_detail::Job Job;

  // direct: the serial order at assembly (tails); otherwise the task's own
  // log, aggregated per barcode and applied when the task ends.
  void Summary(bool direct, AtacSummaryDeltaLog *task_summary,
               uint64_t barcode, int type, uint64_t change) {
    if (!summary_enabled_) {
      return;
    }
    if (direct) {
      w_.summary_metadata_.UpdateCount(barcode, type, change);
      return;
    }
    task_summary->Add(barcode, type, change);
  }

  // Adds a finished task's deltas: one atomic add per barcode and field for
  // barcodes that already have a row; the others keep their first-seen order
  // and are inserted during the ordered assembly. No row is inserted here.
  void ApplyTaskSummary(const AtacSummaryDeltaLog &task_summary,
                        AtacSummaryDeltaLog *new_keys) {
    SummaryMetadata &summary = w_.summary_metadata_;
    for (const AtacSummaryDelta &delta : task_summary.keys()) {
      const khiter_t it = summary.FindExisting(delta.barcode);
      if (it == summary.End()) {
        new_keys->Append(delta);
        continue;
      }
      if (delta.duplicate != 0) {
        summary.AddExistingAtomic(it, SUMMARY_METADATA_DUP, delta.duplicate);
      }
      if (delta.lowmapq != 0) {
        summary.AddExistingAtomic(it, SUMMARY_METADATA_LOWMAPQ,
                                  delta.lowmapq);
      }
      if (delta.mapped != 0) {
        summary.AddExistingAtomic(it, SUMMARY_METADATA_MAPPED, delta.mapped);
      }
    }
  }

  AtacFragmentSink DirectSink() {
    AtacFragmentSink sink;
    if (p_.AtacSidecarOnly()) {
      sink.evidence_fp = w_.atac_evidence_fp_;
      sink.evidence_records = &w_.atac_evidence_records_written_;
    }
    sink.text_fp = w_.mapping_output_file_;
    if (p_.mapping_output_format == MAPPINGFORMAT_BED &&
        p_.macs3_frag_buffer) {
      sink.buckets = p_.macs3_frag_buffer.get();
      sink.allow_bucket_resize = true;
    }
    return sink;
  }

  template <typename Record>
  void EmitGroup(uint32_t rid, const SequenceBatch &reference,
                 const khash_t(k64_seq) * barcode_whitelist_lookup_table,
                 lowmem_finalize_detail::GroupState<Record> *state,
                 bool end_of_stream, AtacFragmentSink *sink, bool direct,
                 AtacSummaryDeltaLog *task_summary, Counters *counters);

  template <typename Record>
  void EmitTail(uint32_t rid, const SequenceBatch &reference,
                const khash_t(k64_seq) * barcode_whitelist_lookup_table,
                lowmem_finalize_detail::GroupState<Record> *state,
                bool end_of_stream, Counters *counters) {
    AtacFragmentSink sink = DirectSink();
    EmitGroup(rid, reference, barcode_whitelist_lookup_table, state,
              end_of_stream, &sink, /*direct=*/true, nullptr, counters);
    if (sink.evidence_write_failed) {
      ExitWithMessage("Failed to write ATAC evidence binary record: " +
                      w_.atac_evidence_tmp_path_);
    }
  }

  template <typename Policy>
  bool RunTask(const Job &job, uint16_t schema,
               const khash_t(k64_seq) * barcode_whitelist_lookup_table,
               const SequenceBatch &reference,
               const std::string &partition_directory,
               lowmem_finalize_detail::TaskResult<typename Policy::Record>
                   *result);

  template <typename Policy>
  void RunTasksAndAssemble(
      const std::vector<Job> &jobs, uint16_t schema, int workers,
      const SequenceBatch &reference,
      const khash_t(k64_seq) * barcode_whitelist_lookup_table,
      double start_time);

  [[noreturn]] void FailAndCleanUp(const std::vector<std::string> &partitions,
                                   const std::string &message) {
    for (const std::string &path : partitions) {
      if (!path.empty()) {
        unlink(path.c_str());
      }
    }
    std::vector<std::string> &spills =
        MappingWriter<AtacSpillRecord>::shared_overflow_file_paths_;
    for (const std::string &path : spills) {
      unlink(path.c_str());
    }
    spills.clear();
    ExitWithMessage(message);
    std::abort();
  }

#endif  // LEGACY_OVERFLOW

  MappingWriter<AtacSpillRecord> &w_;
  const MappingParameters &p_;
  bool summary_enabled_ = false;
  bool dedup_bulk_ = false;
};

void AtacLowMemFinalizer::WriteNonDualFragment(
    MappingWriter<AtacSpillRecord> &w, AtacFragmentSink *sink, uint32_t rid,
    const SequenceBatch &reference, const PairedEndMappingWithBarcode &mwb) {
  const MappingParameters &p = w.mapping_parameters_;
  if (p.mapping_output_format == MAPPINGFORMAT_BED) {
    uint32_t mapping_end_position = mwb.GetEndPosition();
    if (p.AtacSidecarOnly()) {
      // The same AEV1 record the dual BAM/CRAM branch writes, with no text
      // row formatted: start, exclusive end, collapsed duplicate count and
      // the untranslated barcode key.
      if (sink->evidence_fp != nullptr) {
        AtacEvidenceBinaryRecord rec;
        rec.chrom_id = static_cast<int32_t>(rid);
        rec.start = static_cast<int32_t>(mwb.GetStartPosition());
        rec.end = static_cast<int32_t>(mapping_end_position);
        rec.count = static_cast<uint32_t>(mwb.num_dups_);
        rec.barcode_key = mwb.cell_barcode_;
        if (fwrite(&rec, sizeof(rec), 1, sink->evidence_fp) != 1) {
          sink->evidence_write_failed = true;
        } else {
          ++*sink->evidence_records;
        }
      }
    } else if (p.is_bulk_data) {
      const std::string strand = mwb.IsPositiveStrand() ? "+" : "-";
      const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
      WriteText(sink, std::string(reference_sequence_name) + "\t" +
                          std::to_string(mwb.GetStartPosition()) + "\t" +
                          std::to_string(mapping_end_position) + "\tN\t" +
                          std::to_string(mwb.mapq_) + "\t" + strand + "\t" +
                          std::to_string(mwb.num_dups_) + "\n");
    } else {
      const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
      const std::string translated_barcode = w.barcode_translator_.Translate(
          mwb.cell_barcode_, w.cell_barcode_length_);
      WriteText(sink, std::string(reference_sequence_name) + "\t" +
                          std::to_string(mwb.GetStartPosition()) + "\t" +
                          std::to_string(mapping_end_position) + "\t" +
                          translated_barcode + "\t" +
                          std::to_string(mwb.num_dups_) + "\n");
    }
    if (sink->buckets != nullptr) {
      macs3::FragmentRecord rec;
      rec.chrom_id = static_cast<int32_t>(rid);
      rec.start = static_cast<int32_t>(mwb.GetStartPosition());
      rec.end = static_cast<int32_t>(mapping_end_position);
      rec.count = static_cast<uint32_t>(mwb.num_dups_);
      if (rec.end > rec.start && rec.count > 0) {
        auto &buckets = *sink->buckets;
        if (rid >= buckets.size()) {
          if (!sink->allow_bucket_resize) {
            sink->bucket_out_of_range = true;
            return;
          }
          buckets.resize(rid + 1);
        }
        buckets[rid].push_back(rec);
      }
    }
  } else {
    bool positive_strand = mwb.IsPositiveStrand();
    uint32_t positive_read_end =
        mwb.fragment_start_position_ + mwb.positive_alignment_length_;
    uint32_t negative_read_end =
        mwb.fragment_start_position_ + mwb.fragment_length_;
    uint32_t negative_read_start =
        negative_read_end - mwb.negative_alignment_length_;
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    if (positive_strand) {
      WriteText(sink, std::string(reference_sequence_name) + "\t" +
                          std::to_string(mwb.fragment_start_position_) + "\t" +
                          std::to_string(positive_read_end) + "\tN\t" +
                          std::to_string(mwb.mapq_) + "\t+\n" +
                          std::string(reference_sequence_name) + "\t" +
                          std::to_string(negative_read_start) + "\t" +
                          std::to_string(negative_read_end) + "\tN\t" +
                          std::to_string(mwb.mapq_) + "\t-\n");
    } else {
      WriteText(sink, std::string(reference_sequence_name) + "\t" +
                          std::to_string(negative_read_start) + "\t" +
                          std::to_string(negative_read_end) + "\tN\t" +
                          std::to_string(mwb.mapq_) + "\t-\n" +
                          std::string(reference_sequence_name) + "\t" +
                          std::to_string(mwb.fragment_start_position_) + "\t" +
                          std::to_string(positive_read_end) + "\tN\t" +
                          std::to_string(mwb.mapq_) + "\t+\n");
    }
  }
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

template <>
void MappingWriter<AtacSpillRecord>::OutputTempMappingsToOverflow(
    uint32_t num_reference_sequences,
    std::vector<std::vector<AtacSpillRecord>> &mappings_on_diff_ref_seqs) {
  if (!tls_overflow_writer_) {
    std::string base_dir = mapping_parameters_.temp_directory_path;
    tls_overflow_writer_ =
        std::unique_ptr<OverflowWriter>(new OverflowWriter(base_dir, "chromap"));
  }
  const uint16_t sm = AtacSpillSchemaMaskForParameters(
      mapping_parameters_.AtacDualFragmentAndBam() ||
          mapping_parameters_.CreatesMergeableAtacSpill(),
      mapping_parameters_.is_bulk_data,
      mapping_parameters_.CreatesMergeableAtacSpill() &&
          !mapping_parameters_.is_bulk_data);
  tls_overflow_writer_->EnableAtacSpillFileHeader(sm);
  for (uint32_t rid = 0; rid < num_reference_sequences; ++rid) {
    for (const auto &mapping : mappings_on_diff_ref_seqs[rid]) {
      tls_overflow_writer_->WriteAtac(rid, mapping);
    }
    mappings_on_diff_ref_seqs[rid].clear();
  }
}

namespace {

template <typename MappingRecord>
void LoadMappingFromOverflowPayload(MappingRecord *out, FILE *fp,
                                    uint16_t atac_spill_file_schema_mask) {
  (void)atac_spill_file_schema_mask;
  out->LoadFromFile(fp);
}

template <typename MappingRecord>
void LoadMappingFromOverflowPayloadBytes(
    MappingRecord *out, const std::string &payload,
    uint16_t atac_spill_file_schema_mask) {
  FILE *memory =
      fmemopen(const_cast<char *>(payload.data()), payload.size(), "rb");
  if (memory == nullptr) {
    ExitWithMessage("Cannot decode overflow payload");
  }
  LoadMappingFromOverflowPayload(out, memory, atac_spill_file_schema_mask);
  const long consumed = ftell(memory);
  fclose(memory);
  if (consumed < 0 || static_cast<size_t>(consumed) != payload.size()) {
    ExitWithMessage("Overflow payload has trailing or truncated bytes");
  }
}

template <>
void LoadMappingFromOverflowPayloadBytes<AtacSpillRecord>(
    AtacSpillRecord *out, const std::string &payload,
    uint16_t atac_spill_file_schema_mask) {
  std::string error;
  if (!DecodeAtacKwaySpillRecord(payload.data(), payload.size(),
                                 atac_spill_file_schema_mask, out, &error)) {
    ExitWithMessage(error);
  }
}

template <typename MappingRecord>
uint16_t OverflowSpillSchemaMaskFromReader(OverflowReader *r) {
  (void)r;
  return 0;
}

template <>
uint16_t OverflowSpillSchemaMaskFromReader<AtacSpillRecord>(
    OverflowReader *r) {
  if (!r->FileHasAtacSpillHeader()) {
    ExitWithMessage(
        "ATAC spill overflow file is missing AtacSpillFileHeader; merge "
        "refused");
  }
  return r->AtacSpillSchemaFromFileHeader();
}

}  // namespace

template <typename MappingRecord>
void MappingWriter<MappingRecord>::WriteMergeableAtacSpillFromOverflow(
    uint32_t, const SequenceBatch &, const khash_t(k64_seq) *, uint64_t) {
  ExitWithMessage(
      "Durable mergeable spill output is supported only for ATAC records");
}

template <>
void MappingWriter<AtacSpillRecord>::WriteMergeableAtacSpillFromOverflow(
    uint32_t num_reference_sequences, const SequenceBatch &reference,
    const khash_t(k64_seq) *barcode_whitelist_lookup_table,
    uint64_t local_num_sample_barcodes) {
  AtacMergeableSpillMetadata metadata;
  metadata.schema_mask = AtacSpillSchemaMaskForParameters(
      /*dual_bam_fragments=*/true, mapping_parameters_.is_bulk_data,
      /*raw_barcode_evidence=*/!mapping_parameters_.is_bulk_data);
  if (mapping_parameters_.emit_noY_stream ||
      mapping_parameters_.emit_Y_stream ||
      mapping_parameters_.emit_y_read_names ||
      mapping_parameters_.emit_y_noy_fastq) {
    metadata.schema_mask = static_cast<uint16_t>(
        metadata.schema_mask | kAtacSpillSchemaHasYHit);
  }
  if (mapping_parameters_.remove_pcr_duplicates) {
    metadata.flags = static_cast<uint16_t>(
        metadata.flags | kAtacMergeableRemovePcrDuplicates);
  }
  if (mapping_parameters_.remove_pcr_duplicates_at_bulk_level) {
    metadata.flags = static_cast<uint16_t>(
        metadata.flags | kAtacMergeableBulkLevelDedup);
  }
  if (mapping_parameters_.Tn5_shift) {
    metadata.flags =
        static_cast<uint16_t>(metadata.flags | kAtacMergeableTn5Shift);
  }
  if (mapping_parameters_.output_mappings_not_in_whitelist) {
    metadata.flags = static_cast<uint16_t>(
        metadata.flags | kAtacMergeableOutputMappingsNotInWhitelist);
  }
  if (mapping_parameters_.allocate_multi_mappings) {
    metadata.flags = static_cast<uint16_t>(
        metadata.flags | kAtacMergeableAllocateMultiMappings);
  }
  metadata.flags = static_cast<uint16_t>(
      metadata.flags | kAtacMergeableHasSummaryEvidence);
  metadata.flags = static_cast<uint16_t>(
      metadata.flags | kAtacMergeableHasHotSidecar);
  if (mapping_parameters_.output_num_uniq_cache_slots) {
    metadata.flags = static_cast<uint16_t>(
        metadata.flags | kAtacMergeableSummaryCardinality);
  }
  if (mapping_parameters_.mergeable_spill_read_range_late_bound) {
    metadata.flags = static_cast<uint16_t>(
        metadata.flags | kAtacMergeableReadRangeLateBound);
  }
  metadata.shard_ordinal = mapping_parameters_.mergeable_spill_shard_ordinal;
  metadata.shard_count = mapping_parameters_.mergeable_spill_shard_count;
  if (mapping_parameters_.mergeable_spill_read_range_late_bound) {
    metadata.first_global_read_ordinal = 0;
    metadata.input_record_count = mergeable_summary_evidence_count_;
  } else {
    metadata.first_global_read_ordinal =
        mapping_parameters_.mergeable_spill_first_global_read;
    metadata.input_record_count =
        mapping_parameters_.mergeable_spill_input_record_count;
    if (metadata.input_record_count != mergeable_summary_evidence_count_) {
      ExitWithMessage(
          "Observed synchronized input count disagrees with the explicit "
          "mergeable ATAC shard range");
    }
  }
  metadata.local_num_sample_barcodes = local_num_sample_barcodes;
  metadata.barcode_whitelist_fingerprint =
      mapping_parameters_.barcode_whitelist_fingerprint;
  metadata.barcode_correction_error_threshold = static_cast<uint8_t>(
      mapping_parameters_.barcode_correction_error_threshold);
  metadata.barcode_correction_probability_threshold =
      mapping_parameters_.barcode_correction_probability_threshold;
  metadata.multi_mapping_allocation_distance =
      mapping_parameters_.multi_mapping_allocation_distance;
  metadata.multi_mapping_allocation_seed =
      mapping_parameters_.multi_mapping_allocation_seed;
  metadata.max_num_best_mappings = static_cast<uint32_t>(
      mapping_parameters_.max_num_best_mappings);
  metadata.summary_evidence_count = mergeable_summary_evidence_count_;
  metadata.cache_size = static_cast<uint32_t>(mapping_parameters_.cache_size);
  metadata.k_for_minhash =
      static_cast<uint32_t>(mapping_parameters_.k_for_minhash);
  std::stringstream coefficient_stream(mapping_parameters_.frip_est_params);
  std::string coefficient;
  while (std::getline(coefficient_stream, coefficient, ';')) {
    try {
      metadata.frip_est_coefficients.push_back(std::stod(coefficient));
    } catch (...) {
      ExitWithMessage("Invalid mergeable ATAC FRIP summary coefficients");
    }
  }
  if (metadata.frip_est_coefficients.size() != 5) {
    ExitWithMessage("Mergeable ATAC requires five FRIP summary coefficients");
  }
  metadata.barcode_length = cell_barcode_length_;
  metadata.mapq_threshold = mapping_parameters_.mapq_threshold;
  metadata.tn5_forward_shift =
      static_cast<int8_t>(mapping_parameters_.Tn5_forward_shift);
  metadata.tn5_reverse_shift =
      static_cast<int8_t>(mapping_parameters_.Tn5_reverse_shift);
  metadata.sample_id = mapping_parameters_.mergeable_spill_sample_id;
  metadata.input_id = mapping_parameters_.mergeable_spill_input_id;
  metadata.references.reserve(num_reference_sequences);
  for (uint32_t rid = 0; rid < num_reference_sequences; ++rid) {
    AtacMergeableSpillReference entry;
    entry.name = reference.GetSequenceNameAt(rid);
    entry.length = reference.GetSequenceLengthAt(rid);
    metadata.references.push_back(std::move(entry));
  }
  if (!mapping_parameters_.is_bulk_data) {
    if (barcode_whitelist_lookup_table == nullptr) {
      ExitWithMessage(
          "barcoded ATAC mergeable spill requires a barcode whitelist");
    }
    metadata.barcode_abundance_entries.reserve(
        kh_size(barcode_whitelist_lookup_table));
    for (khiter_t it = kh_begin(barcode_whitelist_lookup_table);
         it != kh_end(barcode_whitelist_lookup_table); ++it) {
      if (kh_exist(barcode_whitelist_lookup_table, it)) {
        AtacBarcodeAbundanceEntry entry;
        entry.barcode_key = kh_key(barcode_whitelist_lookup_table, it);
        entry.count = kh_value(barcode_whitelist_lookup_table, it);
        metadata.barcode_abundance_entries.push_back(entry);
      }
    }
    std::sort(metadata.barcode_abundance_entries.begin(),
              metadata.barcode_abundance_entries.end(),
              [](const AtacBarcodeAbundanceEntry &a,
                 const AtacBarcodeAbundanceEntry &b) {
                return a.barcode_key < b.barcode_key;
              });
  }

  AtacMergeableSpillWriter writer;
  AtacHotSpillWriter hot_writer;
  std::string error;
  if (!writer.Open(mapping_parameters_.create_mergeable_spill_record_path,
                   metadata, &error)) {
    ExitWithMessage(error);
  }
  const std::string hot_spill_path = AtacHotSpillSidecarPath(
      mapping_parameters_.create_mergeable_spill_record_path);
  if (!hot_writer.Open(hot_spill_path, metadata, &error)) {
    ExitWithMessage(error);
  }
  if (mergeable_summary_evidence_file_ == nullptr ||
      fflush(mergeable_summary_evidence_file_) != 0 ||
      fseek(mergeable_summary_evidence_file_, 0, SEEK_SET) != 0) {
    ExitWithMessage("Cannot replay mergeable ATAC summary evidence");
  }
  for (uint64_t i = 0; i < mergeable_summary_evidence_count_; ++i) {
    AtacSummaryEvidence evidence;
    evidence.raw_barcode_qual.assign(
        mapping_parameters_.is_bulk_data ? 0 : cell_barcode_length_, '\0');
    if (fread(&evidence.raw_barcode_key, sizeof(evidence.raw_barcode_key), 1,
              mergeable_summary_evidence_file_) != 1 ||
        fread(&evidence.raw_barcode_n_mask,
              sizeof(evidence.raw_barcode_n_mask), 1,
              mergeable_summary_evidence_file_) != 1 ||
        fread(&evidence.cache_slot1, sizeof(evidence.cache_slot1), 1,
              mergeable_summary_evidence_file_) != 1 ||
        fread(&evidence.cache_slot2, sizeof(evidence.cache_slot2), 1,
              mergeable_summary_evidence_file_) != 1 ||
        (!evidence.raw_barcode_qual.empty() &&
         fread(&evidence.raw_barcode_qual[0], 1,
               evidence.raw_barcode_qual.size(),
               mergeable_summary_evidence_file_) !=
             evidence.raw_barcode_qual.size()) ||
        !writer.AppendSummaryEvidence(evidence, &error)) {
      ExitWithMessage(error.empty()
                          ? "Cannot read mergeable ATAC summary evidence"
                          : error);
    }
  }

  struct FileRecord {
    uint32_t rid = 0;
    AtacSpillRecord mapping;
    size_t file_index = 0;

    bool operator<(const FileRecord &other) const {
      if (rid != other.rid) {
        return rid > other.rid;
      }
      const bool a_less_b = mapping < other.mapping;
      const bool b_less_a = other.mapping < mapping;
      if (!a_less_b && !b_less_a) {
        return file_index > other.file_index;
      }
      return !a_less_b;
    }
  };

  std::unordered_map<uint32_t, std::vector<size_t>> rid_to_files;
  std::unordered_set<uint32_t> all_rids;
  uint16_t agreed_schema = 0;
  bool have_agreed_schema = false;
  for (size_t fi = 0; fi < shared_overflow_file_paths_.size(); ++fi) {
    OverflowReader scanner(shared_overflow_file_paths_[fi]);
    if (!scanner.IsValid()) {
      ExitWithMessage("Cannot open ATAC overflow run for durable spill: " +
                      shared_overflow_file_paths_[fi]);
    }
    uint32_t rid = 0;
    std::string payload;
    std::unordered_set<uint32_t> rids_in_file;
    if (scanner.FileHasAtacKwayHeader()) {
      rid = scanner.AtacKwayReferenceIdFromFileHeader();
      rids_in_file.insert(rid);
      all_rids.insert(rid);
    } else {
      while (scanner.ReadNext(rid, payload)) {
        rids_in_file.insert(rid);
        all_rids.insert(rid);
      }
    }
    const uint16_t schema =
        OverflowSpillSchemaMaskFromReader<AtacSpillRecord>(&scanner);
    if (!have_agreed_schema) {
      agreed_schema = schema;
      have_agreed_schema = true;
    } else if (schema != agreed_schema) {
      ExitWithMessage(
          "Mismatched ATAC spill schema between worker overflow runs");
    }
    if ((schema & kAtacSpillSchemaHasBamPair) == 0) {
      ExitWithMessage(
          "Mergeable ATAC spill source is missing BAM-pair payloads");
    }
    for (uint32_t file_rid : rids_in_file) {
      rid_to_files[file_rid].push_back(fi);
    }
  }

  std::vector<uint32_t> rids(all_rids.begin(), all_rids.end());
  std::sort(rids.begin(), rids.end());
  std::vector<std::unique_ptr<OverflowReader>> readers(
      shared_overflow_file_paths_.size());

  auto load_record = [&](size_t fi, uint32_t expected_rid,
                         FileRecord *output, bool *eof) {
    uint32_t rid = 0;
    std::string payload;
    if (!readers[fi]->ReadNext(rid, payload)) {
      *eof = true;
      return;
    }
    *eof = false;
    if (rid != expected_rid) {
      ExitWithMessage(
          "ATAC overflow run crossed reference buckets during durable merge");
    }
    AtacSpillRecord mapping;
    LoadMappingFromOverflowPayloadBytes(&mapping, payload, agreed_schema);
    if (!mapping.HasBamPairSection()) {
      ExitWithMessage("Invalid ATAC overflow payload during durable merge");
    }
    output->rid = rid;
    output->mapping = std::move(mapping);
    output->file_index = fi;
  };

  for (uint32_t current_rid : rids) {
    const auto &file_indices = rid_to_files[current_rid];
    std::priority_queue<FileRecord> heap;
    for (size_t fi : file_indices) {
      readers[fi].reset(new OverflowReader(shared_overflow_file_paths_[fi]));
      if (!readers[fi]->IsValid()) {
        ExitWithMessage("Cannot reopen ATAC overflow run for durable merge");
      }
      FileRecord record;
      bool eof = false;
      load_record(fi, current_rid, &record, &eof);
      if (!eof) {
        heap.push(std::move(record));
      }
    }
    while (!heap.empty()) {
      FileRecord record = heap.top();
      heap.pop();
      if (record.mapping.num_dups_ != 1 ||
          record.mapping.read_id_ >= metadata.input_record_count ||
          record.mapping.sam1.read_id_ != record.mapping.read_id_ ||
          record.mapping.sam2.read_id_ != record.mapping.read_id_) {
        ExitWithMessage(
            "ATAC mergeable spill source is not a pre-dedup local-id record");
      }
      const bool y_hit =
          reads_with_y_hit_ != nullptr &&
          reads_with_y_hit_->count(record.mapping.read_id_) != 0;
      record.mapping.SetYHit(y_hit);
      if (!writer.Append(record.rid, record.mapping, &error)) {
        ExitWithMessage(error);
      }
      if (!hot_writer.Append(record.rid, record.mapping, &error)) {
        ExitWithMessage(error);
      }

      FileRecord next;
      bool eof = false;
      load_record(record.file_index, current_rid, &next, &eof);
      if (!eof) {
        heap.push(std::move(next));
      }
    }
    for (size_t fi : file_indices) {
      readers[fi].reset();
    }
  }

  // Publish the hot companion first. The parent header advertises it only
  // after the full spill is atomically renamed, so a completed parent can
  // never point at an unpublished hot file.
  if (!hot_writer.Finalize(&error) || !writer.Finalize(&error)) {
    ExitWithMessage(error);
  }
  for (const std::string &file_path : shared_overflow_file_paths_) {
    unlink(file_path.c_str());
  }
  shared_overflow_file_paths_.clear();
  std::cerr << "Wrote mergeable pre-dedup ATAC spill shard "
            << metadata.shard_ordinal << "/" << metadata.shard_count << " ("
            << writer.record_count() << " records) to "
            << mapping_parameters_.create_mergeable_spill_record_path
            << " with fixed-width hot companion " << hot_spill_path << "\n";
}

template <typename MappingRecord>
void MappingWriter<MappingRecord>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  ProcessLowMemOverflowSerial(num_mappings_in_mem, num_reference_sequences,
                              reference, barcode_whitelist_lookup_table);
}

// The serial low-memory merge: the v1.1.0 loop with the 1.1.1 fixes, kept
// verbatim. MappingWriter<AtacSpillRecord> also has a per-reference path
// (AtacLowMemFinalizer below) that produces the same bytes.
template <typename MappingRecord>
void MappingWriter<MappingRecord>::ProcessLowMemOverflowSerial(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  
  // All thread-local writers should already be closed by now
  // Use the shared collection of overflow file paths
  
  if (shared_overflow_file_paths_.empty()) {
    return;
  }
  
  std::cerr << "Processing " << shared_overflow_file_paths_.size() << " overflow files for k-way merge\n";
  std::cerr << "Low-memory overflow: mid-batch flush count: "
            << LowMemMidBatchOverflowFlushCount() << "\n";

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
      ExitOnUnopenableOverflowFile(shared_overflow_file_paths_,
                                   shared_overflow_file_paths_[fi],
                                   scanner.OpenErrno(), 0, 0);
    }

    // Scan file to find all rids it contains
    uint32_t rid;
    std::string payload;
    std::unordered_set<uint32_t> rids_in_file;
    if (scanner.FileHasAtacKwayHeader()) {
      rid = scanner.AtacKwayReferenceIdFromFileHeader();
      rids_in_file.insert(rid);
      all_rids.insert(rid);
    } else {
      while (scanner.ReadNext(rid, payload)) {
        rids_in_file.insert(rid);
        all_rids.insert(rid);
      }
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

  uint16_t merge_spill_schema_agreed = 0;
  bool merge_spill_schema_agreed_set = false;

  // For each rid, k-way merge all files containing that rid
  for (uint32_t current_rid : rids) {
    const auto& file_indices = rid_to_files[current_rid];
    
    // Priority queue for k-way merge within this rid
    std::priority_queue<FileRecord> merge_heap;
    
    // Initialize: read first record from each file for this rid
    // Note: OverflowWriter creates one file per rid, so each file contains only records for one rid
    for (size_t fi : file_indices) {
      readers[fi].reset(new OverflowReader(shared_overflow_file_paths_[fi]));
      if (!readers[fi]->IsValid()) {
        ExitOnUnopenableOverflowFile(shared_overflow_file_paths_,
                                     shared_overflow_file_paths_[fi],
                                     readers[fi]->OpenErrno(), current_rid,
                                     file_indices.size());
      }

      uint32_t rid;
      std::string payload;
      if (readers[fi]->ReadNext(rid, payload)) {
        // Verify rid matches (should always be true, but check for safety)
        if (rid == current_rid) {
          const uint16_t spill_mask =
              OverflowSpillSchemaMaskFromReader<MappingRecord>(
                  readers[fi].get());
          if (!merge_spill_schema_agreed_set) {
            merge_spill_schema_agreed = spill_mask;
            merge_spill_schema_agreed_set = true;
          } else if (spill_mask != merge_spill_schema_agreed) {
            ExitWithMessage(
                "Mismatched ATAC spill schema_mask between overflow temp "
                "files (merge refused)");
          }
          MappingRecord mapping;
          LoadMappingFromOverflowPayloadBytes(&mapping, payload, spill_mask);

          FileRecord rec;
          rec.rid = rid;
          rec.mapping = mapping;
          rec.file_index = fi;
          merge_heap.push(rec);
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
            const uint16_t spill_mask =
                OverflowSpillSchemaMaskFromReader<MappingRecord>(
                    readers[min_rec.file_index].get());
            if (spill_mask != merge_spill_schema_agreed) {
              ExitWithMessage(
                  "Mismatched ATAC spill schema_mask within overflow merge "
                  "(merge refused)");
            }
            MappingRecord mapping;
            LoadMappingFromOverflowPayloadBytes(&mapping, payload, spill_mask);

            FileRecord rec;
            rec.rid = rid;
            rec.mapping = mapping;
            rec.file_index = min_rec.file_index;
            merge_heap.push(rec);
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

template void MappingWriter<AtacSpillRecord>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<AtacSpillRecord>>&);
// MappingWriter<AtacSpillRecord>::ProcessAndOutputMappingsInLowMemoryFromOverflow
// is an explicit specialization (after AtacLowMemFinalizer below).

// BED writers (single/paired × with/without barcode). Serialization
// added in bed_mapping.h; generic templates do all the work.
template void MappingWriter<MappingWithBarcode>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<MappingWithBarcode>>&);
template void MappingWriter<MappingWithBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t, uint32_t, const SequenceBatch&, const khash_t(k64_seq)*);

template void MappingWriter<MappingWithoutBarcode>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<MappingWithoutBarcode>>&);
template void MappingWriter<MappingWithoutBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t, uint32_t, const SequenceBatch&, const khash_t(k64_seq)*);

template void MappingWriter<PairedEndMappingWithBarcode>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<PairedEndMappingWithBarcode>>&);
template void MappingWriter<PairedEndMappingWithBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t, uint32_t, const SequenceBatch&, const khash_t(k64_seq)*);

template void MappingWriter<PairedEndMappingWithoutBarcode>::OutputTempMappingsToOverflow(
    uint32_t, std::vector<std::vector<PairedEndMappingWithoutBarcode>>&);
template void MappingWriter<PairedEndMappingWithoutBarcode>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
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

// All four BED mapping types (single/paired × with/without barcode) now
// have WriteToFile / LoadFromFile / SerializedSize in bed_mapping.h plus
// operator< / operator==, so the generic
// MappingWriter<MappingRecord>::OutputTempMappingsToOverflow and
// ProcessAndOutputMappingsInLowMemoryFromOverflow templates apply
// directly. Explicit template instantiations live below alongside the
// SAM/PAF/Pairs/Dual instantiations; no per-type stub bodies needed.

// Add explicit instantiation for CloseThreadOverflowWriter
template void MappingWriter<SAMMapping>::CloseThreadOverflowWriter();
template void MappingWriter<PAFMapping>::CloseThreadOverflowWriter();
template void MappingWriter<PairedPAFMapping>::CloseThreadOverflowWriter();
template void MappingWriter<PairsMapping>::CloseThreadOverflowWriter();
template void MappingWriter<MappingWithBarcode>::CloseThreadOverflowWriter();
template void MappingWriter<MappingWithoutBarcode>::CloseThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithBarcode>::CloseThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithoutBarcode>::CloseThreadOverflowWriter();
template void MappingWriter<AtacSpillRecord>::CloseThreadOverflowWriter();

// Add explicit instantiation for RotateThreadOverflowWriter
template void MappingWriter<SAMMapping>::RotateThreadOverflowWriter();
template void MappingWriter<PAFMapping>::RotateThreadOverflowWriter();
template void MappingWriter<PairedPAFMapping>::RotateThreadOverflowWriter();
template void MappingWriter<PairsMapping>::RotateThreadOverflowWriter();
template void MappingWriter<MappingWithBarcode>::RotateThreadOverflowWriter();
template void MappingWriter<MappingWithoutBarcode>::RotateThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithBarcode>::RotateThreadOverflowWriter();
template void MappingWriter<PairedEndMappingWithoutBarcode>::RotateThreadOverflowWriter();
template void MappingWriter<AtacSpillRecord>::RotateThreadOverflowWriter();

// The two emission blocks of the serial merge, kept statement for statement:
// end_of_stream == false is the in-loop block (best duplicate, then the MAPQ
// test); true is the block after the loop (MAPQ test on the last record, then
// best duplicate).
template <typename Record>
void AtacLowMemFinalizer::EmitGroup(
    uint32_t rid, const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table,
    lowmem_finalize_detail::GroupState<Record> *state, bool end_of_stream,
    AtacFragmentSink *sink, bool direct, AtacSummaryDeltaLog *task_summary,
    Counters *counters) {
  Record &last_mapping = state->last_mapping;
  const uint32_t num_last_mapping_dups = state->num_last_mapping_dups;
  std::vector<Record> &temp_dups_for_bulk_level_dedup = state->bulk_dups;
  if (!end_of_stream) {
    if (dedup_bulk_) {
      size_t best_mapping_index = FindBestIndexFromDuplicatesT(
          barcode_whitelist_lookup_table, temp_dups_for_bulk_level_dedup);
      last_mapping = temp_dups_for_bulk_level_dedup[best_mapping_index];
      temp_dups_for_bulk_level_dedup.clear();
    }

    if (last_mapping.mapq_ >= p_.mapq_threshold) {
      last_mapping.num_dups_ =
          std::min((uint32_t)std::numeric_limits<uint8_t>::max(),
                   num_last_mapping_dups);
      if (p_.Tn5_shift) {
        last_mapping.Tn5Shift(p_.Tn5_forward_shift, p_.Tn5_reverse_shift);
      }

      WriteNonDualFragment(w_, sink, rid, reference, last_mapping);
      ++counters->passing;
      Summary(direct, task_summary, last_mapping.GetBarcode(),
              SUMMARY_METADATA_DUP, num_last_mapping_dups - 1);
    } else {
      Summary(direct, task_summary, last_mapping.GetBarcode(),
              SUMMARY_METADATA_LOWMAPQ, num_last_mapping_dups);
    }
    Summary(direct, task_summary, last_mapping.GetBarcode(),
            SUMMARY_METADATA_MAPPED, num_last_mapping_dups);

    if (last_mapping.is_unique_ == 1) {
      ++counters->uni;
    } else {
      ++counters->multi;
    }
    return;
  }

  if (last_mapping.mapq_ >= p_.mapq_threshold) {
    if (dedup_bulk_) {
      size_t best_mapping_index = FindBestIndexFromDuplicatesT(
          barcode_whitelist_lookup_table, temp_dups_for_bulk_level_dedup);
      last_mapping = temp_dups_for_bulk_level_dedup[best_mapping_index];
      temp_dups_for_bulk_level_dedup.clear();
    }

    last_mapping.num_dups_ = std::min(
        (uint32_t)std::numeric_limits<uint8_t>::max(), num_last_mapping_dups);
    if (p_.Tn5_shift) {
      last_mapping.Tn5Shift(p_.Tn5_forward_shift, p_.Tn5_reverse_shift);
    }
    WriteNonDualFragment(w_, sink, rid, reference, last_mapping);
    ++counters->passing;

    Summary(direct, task_summary, last_mapping.GetBarcode(), SUMMARY_METADATA_DUP,
            num_last_mapping_dups - 1);
  } else {
    Summary(direct, task_summary, last_mapping.GetBarcode(),
            SUMMARY_METADATA_LOWMAPQ, num_last_mapping_dups);
  }
  Summary(direct, task_summary, last_mapping.GetBarcode(), SUMMARY_METADATA_MAPPED,
          num_last_mapping_dups);

  if (last_mapping.is_unique_ == 1) {
    ++counters->uni;
  } else {
    ++counters->multi;
  }
}

template <typename Policy>
bool AtacLowMemFinalizer::RunTask(
    const Job &job, uint16_t schema,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table,
    const SequenceBatch &reference, const std::string &partition_directory,
    lowmem_finalize_detail::TaskResult<typename Policy::Record> *result) {
  typedef typename Policy::Record Record;
  const std::vector<std::string> &paths =
      MappingWriter<AtacSpillRecord>::shared_overflow_file_paths_;

  FILE *partition = nullptr;
  if (!CreateLowMemPartitionFile(partition_directory, job.rid,
                                 &result->partition_path, &partition,
                                 &result->error)) {
    result->partition_path.clear();
    return false;
  }
  std::vector<char> partition_buffer(1u << 20);
  (void)setvbuf(partition, partition_buffer.data(), _IOFBF,
                partition_buffer.size());

  AtacFragmentSink sink;
  if (p_.AtacSidecarOnly()) {
    sink.evidence_fp = partition;
    sink.evidence_records = &result->partition_records;
  } else {
    sink.text_fp = partition;
  }
  if (p_.mapping_output_format == MAPPINGFORMAT_BED && p_.macs3_frag_buffer) {
    sink.buckets = p_.macs3_frag_buffer.get();
    sink.allow_bucket_resize = false;
  }

  // The heap order of the serial merge: the smallest record first; equal
  // records by the smaller global spill-file index.
  struct HeapEntry {
    Record mapping;
    size_t file_index;
    size_t slot;
    bool operator<(const HeapEntry &other) const {
      const bool a_less_b = mapping < other.mapping;
      const bool b_less_a = other.mapping < mapping;
      if (!a_less_b && !b_less_a) {
        return file_index > other.file_index;
      }
      return !a_less_b;
    }
  };

  std::vector<std::unique_ptr<OverflowReader>> readers(
      job.file_indices.size());
  std::priority_queue<HeapEntry> heap;
  std::string payload;
  bool ok = true;
  for (size_t k = 0; k < job.file_indices.size(); ++k) {
    const size_t fi = job.file_indices[k];
    readers[k].reset(new OverflowReader(paths[fi]));
    if (!readers[k]->IsValid()) {
      result->error = UnopenableOverflowFileMessage(
          paths[fi], readers[k]->OpenErrno(), job.rid,
          job.file_indices.size());
      ok = false;
      break;
    }
    HeapEntry entry;
    entry.file_index = fi;
    entry.slot = k;
    const int got = Policy::Next(readers[k].get(), job.rid, schema,
                                 &entry.mapping, &payload, &result->error);
    if (got < 0) {
      ok = false;
      break;
    }
    if (got > 0) {
      heap.push(entry);
    }
  }

  lowmem_finalize_detail::GroupState<Record> state;
  AtacSummaryDeltaLog task_summary;
  bool first = true;
  while (ok && !heap.empty()) {
    HeapEntry min_rec = heap.top();
    heap.pop();
    ++result->records_merged;

    const Record &current_min_mapping = min_rec.mapping;
    const bool current_mapping_is_duplicated_at_cell_level =
        !first && current_min_mapping == state.last_mapping;
    const bool current_mapping_is_duplicated_at_bulk_level =
        !first && dedup_bulk_ &&
        current_min_mapping.IsSamePosition(state.last_mapping);
    const bool current_mapping_is_duplicated =
        current_mapping_is_duplicated_at_cell_level ||
        current_mapping_is_duplicated_at_bulk_level;

    if (p_.remove_pcr_duplicates && current_mapping_is_duplicated) {
      ++state.num_last_mapping_dups;
      if (dedup_bulk_) {
        std::vector<Record> &temp = state.bulk_dups;
        if (!temp.empty() && current_min_mapping == temp.back()) {
          temp.back() = current_min_mapping;
          temp.back().num_dups_ += 1;
        } else {
          temp.push_back(current_min_mapping);
          temp.back().num_dups_ = 1;
        }
      }
      state.last_mapping = current_min_mapping;
    } else {
      if (!first) {
        EmitGroup(job.rid, reference, barcode_whitelist_lookup_table, &state,
                  /*end_of_stream=*/false, &sink, /*direct=*/false,
                  &task_summary, &result->counters);
      }
      state.last_mapping = current_min_mapping;
      state.num_last_mapping_dups = 1;
      if (dedup_bulk_) {
        state.bulk_dups.push_back(current_min_mapping);
        state.bulk_dups.back().num_dups_ = 1;
      }
      first = false;
    }

    HeapEntry next;
    next.file_index = min_rec.file_index;
    next.slot = min_rec.slot;
    const int got = Policy::Next(readers[min_rec.slot].get(), job.rid, schema,
                                 &next.mapping, &payload, &result->error);
    if (got < 0) {
      ok = false;
      break;
    }
    if (got > 0) {
      heap.push(next);
    }
  }
  readers.clear();

  if (ok && first) {
    result->error = "Low-memory spill files of reference " +
                    std::to_string(job.rid) + " hold no records";
    ok = false;
  }
  if (ok && sink.bucket_out_of_range) {
    result->error =
        "Low-memory finalization: MACS3 fragment buckets are not sized for "
        "reference " + std::to_string(job.rid);
    ok = false;
  }
  if (fclose(partition) != 0 && ok) {
    result->error = "Cannot close low-memory finalization partition " +
                    result->partition_path;
    ok = false;
  }
  if (ok && (sink.evidence_write_failed || sink.text_write_failed)) {
    result->error = "Cannot write low-memory finalization partition " +
                    result->partition_path;
    ok = false;
  }
  if (!ok) {
    unlink(result->partition_path.c_str());
    result->partition_path.clear();
    return false;
  }

  if (summary_enabled_) {
    ApplyTaskSummary(task_summary, &result->new_keys);
    task_summary.Clear();
  }
  state.active = true;
  result->tail = state;
  // This reference's spill files are no longer needed.
  for (size_t fi : job.file_indices) {
    unlink(paths[fi].c_str());
  }
  result->done = true;
  return true;
}

template <typename Policy>
void AtacLowMemFinalizer::RunTasksAndAssemble(
    const std::vector<Job> &jobs, uint16_t schema, int workers,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table,
    double start_time) {
  typedef typename Policy::Record Record;
  std::vector<std::string> &paths =
      MappingWriter<AtacSpillRecord>::shared_overflow_file_paths_;
  const std::string partition_directory = DirectoryOfPath(paths[0]);

  std::vector<lowmem_finalize_detail::TaskResult<Record>> results(
      jobs.size());
  // Longest first, by spill bytes; the output does not depend on the order.
  std::vector<size_t> order(jobs.size());
  for (size_t j = 0; j < order.size(); ++j) {
    order[j] = j;
  }
  std::sort(order.begin(), order.end(), [&jobs](size_t a, size_t b) {
    if (jobs[a].spill_bytes != jobs[b].spill_bytes) {
      return jobs[a].spill_bytes > jobs[b].spill_bytes;
    }
    return jobs[a].rid < jobs[b].rid;
  });

  std::atomic<bool> failed(false);
#pragma omp parallel for schedule(dynamic, 1) num_threads(workers)
  for (int64_t k = 0; k < static_cast<int64_t>(order.size()); ++k) {
    if (failed.load(std::memory_order_relaxed)) {
      continue;
    }
    const size_t j = order[static_cast<size_t>(k)];
    // Every task runs under a host permit when the host provides hooks.
    LowMemPermitScope permit(p_);
    const bool ok =
        RunTask<Policy>(jobs[j], schema, barcode_whitelist_lookup_table,
                        reference, partition_directory, &results[j]);
    permit.Release(results[j].records_merged, jobs[j].spill_bytes);
    if (!ok) {
      failed.store(true, std::memory_order_relaxed);
    }
  }

  std::vector<std::string> partitions(results.size());
  for (size_t j = 0; j < results.size(); ++j) {
    partitions[j] = results[j].partition_path;
  }
  if (failed.load()) {
    std::string first_error;
    for (size_t j = 0; j < results.size(); ++j) {
      if (!results[j].error.empty()) {
        first_error = results[j].error;
        break;
      }
    }
    FailAndCleanUp(partitions, first_error.empty()
                                   ? std::string("Low-memory finalization failed")
                                   : first_error);
  }

  // Ordered assembly, reference by reference.
  FILE *destination =
      p_.AtacSidecarOnly() ? w_.atac_evidence_fp_ : w_.mapping_output_file_;
  std::vector<char> copy_buffer;
  Counters totals;
  lowmem_finalize_detail::GroupState<Record> pending;
  uint32_t pending_rid = 0;
  bool has_pending = false;
  for (size_t j = 0; j < jobs.size(); ++j) {
    lowmem_finalize_detail::TaskResult<Record> &result = results[j];
    if (has_pending) {
      EmitTail(pending_rid, reference, barcode_whitelist_lookup_table,
               &pending, /*end_of_stream=*/false, &totals);
    }
    if (destination != nullptr) {
      uint64_t copied = 0;
      std::string error;
      if (!AppendFileContents(destination, result.partition_path,
                              &copy_buffer, &copied, &error)) {
        FailAndCleanUp(partitions, error);
      }
    }
    if (p_.AtacSidecarOnly()) {
      w_.atac_evidence_records_written_ += result.partition_records;
    }
    unlink(result.partition_path.c_str());
    partitions[j].clear();
    if (summary_enabled_) {
      for (const AtacSummaryDelta &key : result.new_keys.keys()) {
        // One put per field, so a put always follows the insertion.
        w_.summary_metadata_.UpdateCount(key.barcode, SUMMARY_METADATA_DUP,
                                         key.duplicate);
        w_.summary_metadata_.UpdateCount(key.barcode,
                                         SUMMARY_METADATA_LOWMAPQ,
                                         key.lowmapq);
        w_.summary_metadata_.UpdateCount(key.barcode, SUMMARY_METADATA_MAPPED,
                                         key.mapped);
      }
    }
    totals.uni += result.counters.uni;
    totals.multi += result.counters.multi;
    totals.passing += result.counters.passing;
    pending = result.tail;
    pending_rid = jobs[j].rid;
    has_pending = true;
    result.tail = lowmem_finalize_detail::GroupState<Record>();
  }
  if (has_pending) {
    EmitTail(pending_rid, reference, barcode_whitelist_lookup_table, &pending,
             /*end_of_stream=*/true, &totals);
  }

  for (const std::string &path : paths) {
    unlink(path.c_str());
  }
  paths.clear();

  if (p_.remove_pcr_duplicates) {
    std::cerr << "Sorted, deduped and outputed mappings in "
              << GetRealTime() - start_time << "s.\n";
  } else {
    std::cerr << "Sorted and outputed mappings in "
              << GetRealTime() - start_time << "s.\n";
  }
  std::cerr << "# uni-mappings: " << totals.uni
            << ", # multi-mappings: " << totals.multi
            << ", total: " << totals.uni + totals.multi << ".\n";
  std::cerr << "Number of output mappings (passed filters): "
            << totals.passing << "\n";
}

bool AtacLowMemFinalizer::Run(
    uint32_t num_reference_sequences, const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table,
    int requested_workers) {
  std::vector<std::string> &paths =
      MappingWriter<AtacSpillRecord>::shared_overflow_file_paths_;
  if (paths.empty()) {
    return true;  // As the serial merge: nothing is written.
  }

  // Probe every spill file: reference, schema and size.
  std::vector<AtacSpillProbe> probes(paths.size());
  for (size_t i = 0; i < paths.size(); ++i) {
    probes[i] = ProbeAtacKwaySpillFile(paths[i]);
    switch (probes[i].status) {
      case AtacSpillProbeStatus::kOpenFailed:
        ExitOnUnopenableOverflowFile(paths, paths[i], probes[i].open_errno, 0,
                                     0);
        break;
      case AtacSpillProbeStatus::kNotKway:
        return false;
      case AtacSpillProbeStatus::kInvalid:
        ExitWithMessage("Unsupported or invalid ATAC k-way spill file header");
        break;
      case AtacSpillProbeStatus::kKway:
        break;
    }
  }
  const uint16_t schema = probes[0].schema_mask;
  for (size_t i = 1; i < probes.size(); ++i) {
    if (probes[i].schema_mask != schema) {
      ExitWithMessage(
          "Mismatched ATAC spill schema_mask between overflow temp files "
          "(merge refused)");
    }
  }
  if ((schema & (kAtacSpillSchemaHasBamPair | kAtacSpillSchemaIsBulk |
                 kAtacSpillSchemaHasRawBarcodeEvidence)) != 0) {
    return false;
  }

  // One job per reference, in reference order.
  std::map<uint32_t, Job> by_reference;
  for (size_t i = 0; i < probes.size(); ++i) {
    Job &job = by_reference[probes[i].reference_id];
    job.rid = probes[i].reference_id;
    job.file_indices.push_back(i);
    job.spill_bytes += probes[i].bytes;
  }
  std::vector<Job> jobs;
  jobs.reserve(by_reference.size());
  size_t max_files = 0;
  uint32_t reference_with_most_files = 0;
  for (auto &entry : by_reference) {
    if (entry.second.file_indices.size() > max_files) {
      max_files = entry.second.file_indices.size();
      reference_with_most_files = entry.first;
    }
    jobs.push_back(std::move(entry.second));
  }

  std::cerr << "Processing " << paths.size()
            << " overflow files for k-way merge\n";
  std::cerr << "Low-memory overflow: mid-batch flush count: "
            << LowMemMidBatchOverflowFlushCount() << "\n";
  const double start_time = GetRealTime();

  LowMemFinalizeWorkerPlan plan;
  std::string error;
  if (!PlanLowMemFinalizeWorkers(requested_workers, jobs.size(), max_files,
                                 reference_with_most_files, &plan, &error)) {
    FailAndCleanUp(std::vector<std::string>(), error);
  }
  std::cerr << "Low-memory finalization: " << jobs.size()
            << " reference tasks on " << plan.workers << " threads";
  if (plan.limited_by_open_files) {
    std::cerr << " (limited by the open-file limit of "
              << plan.open_file_limit << ")";
  }
  std::cerr << "\n";

  if (summary_enabled_ &&
      w_.summary_metadata_.ResizePendingOnNextPut()) {
    // The serial merge's first summary update would resize the table before
    // anything else; do the same before the tasks read it.
    w_.summary_metadata_.ApplyPendingResizeLikePut();
  }
  if (p_.mapping_output_format == MAPPINGFORMAT_BED && p_.macs3_frag_buffer &&
      p_.macs3_frag_buffer->size() < num_reference_sequences) {
    p_.macs3_frag_buffer->resize(num_reference_sequences);
  }

  // The covered schema has no optional sections, so the lean decode reads
  // exactly the fields the merge uses (FullDecodePolicy gives the same bytes).
  RunTasksAndAssemble<lowmem_finalize_detail::LeanDecodePolicy>(
      jobs, schema, plan.workers, reference, barcode_whitelist_lookup_table,
      start_time);
  return true;
}

template <>
void MappingWriter<AtacSpillRecord>::ProcessAndOutputMappingsInLowMemoryFromOverflow(
    uint32_t num_mappings_in_mem, uint32_t num_reference_sequences,
    const SequenceBatch &reference,
    const khash_t(k64_seq) * barcode_whitelist_lookup_table) {
  const int requested = RequestedLowMemFinalizeThreads(mapping_parameters_);
  if (requested >= 2 &&
      AtacLowMemFinalizer::ParametersEligible(mapping_parameters_)) {
    AtacLowMemFinalizer finalizer(*this);
    if (finalizer.Run(num_reference_sequences, reference,
                      barcode_whitelist_lookup_table, requested)) {
      return;
    }
  }
  if (shared_overflow_file_paths_.empty()) {
    return;
  }
  // The serial merge, under one host permit when the host provides hooks.
  uint64_t spill_bytes = 0;
  const uint64_t spill_files = shared_overflow_file_paths_.size();
  if (mapping_parameters_.PermitHooksEnabled()) {
    for (const std::string &path : shared_overflow_file_paths_) {
      struct stat st;
      if (stat(path.c_str(), &st) == 0) {
        spill_bytes += static_cast<uint64_t>(st.st_size);
      }
    }
  }
  LowMemPermitScope permit(mapping_parameters_);
  ProcessLowMemOverflowSerial(num_mappings_in_mem, num_reference_sequences,
                              reference, barcode_whitelist_lookup_table);
  permit.Release(spill_files, spill_bytes);
}


// AtacSpillRecord has WriteToFile / LoadFromFile / SerializedSize
// (see atac_dual_mapping.h) and operator< / operator==, so the generic
// MappingWriter<MappingRecord> overflow templates apply directly. Use
// explicit instantiations below in place of the empty stubs that
// previously blocked --low-mem + --atac-fragments dual ATAC output.

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
  return DeriveReadGroupFromFilenameImpl(filename);
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
      if (!mapping_parameters_.ReadGroupSourcePath().empty()) {
        rg_id = DeriveReadGroupFromFilenameImpl(
            mapping_parameters_.ReadGroupSourcePath());
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
  return BuildBamRecordFromSamMappingFields(
      mapping_parameters_, cell_barcode_length_, barcode_translator_, rid,
      reference, mapping);
}

template <>
void MappingWriter<AtacSpillRecord>::OutputHeader(
    uint32_t num_reference_sequences, const SequenceBatch &reference) {
  if (!mapping_parameters_.AtacDualFragmentAndBam()) {
    if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED) {
      if (mapping_parameters_.macs3_frag_chrom_names) {
        auto& names = *mapping_parameters_.macs3_frag_chrom_names;
        names.clear();
        names.reserve(num_reference_sequences);
        for (uint32_t i = 0; i < num_reference_sequences; ++i) {
          names.emplace_back(reference.GetSequenceNameAt(i));
        }
      }
      if (mapping_parameters_.macs3_frag_buffer) {
        mapping_parameters_.macs3_frag_buffer->assign(
            num_reference_sequences, std::vector<macs3::FragmentRecord>());
      }
      if (mapping_parameters_.AtacSidecarOnly()) {
        OpenAtacEvidenceBinaryOutput(num_reference_sequences, reference);
      }
    }
    return;
  }
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM ||
      mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM) {
    if (!hts_out_) {
      OpenHtsOutput();
    }
    BuildHtsHeader(num_reference_sequences, reference);

    if (noY_hts_out_ && hts_hdr_) {
      if (sam_hdr_write(noY_hts_out_, hts_hdr_) < 0) {
        ExitWithMessage("Failed to write header to noY BAM/CRAM output");
      }
      if (mapping_parameters_.write_index &&
          mapping_parameters_.noY_output_path != "-" &&
          mapping_parameters_.noY_output_path != "/dev/stdout" &&
          mapping_parameters_.noY_output_path != "/dev/stderr") {
        this->noY_index_path_ = mapping_parameters_.noY_output_path;
        if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) {
          this->noY_index_path_ += ".bai";
        } else if (mapping_parameters_.mapping_output_format ==
                   MAPPINGFORMAT_CRAM) {
          this->noY_index_path_ += ".crai";
        }
        int min_shift = 0;
        if (sam_idx_init(noY_hts_out_, hts_hdr_, min_shift,
                         this->noY_index_path_.c_str()) < 0) {
          ExitWithMessage("Failed to initialize index for noY BAM/CRAM output");
        }
      }
    }
    if (Y_hts_out_ && hts_hdr_) {
      if (sam_hdr_write(Y_hts_out_, hts_hdr_) < 0) {
        ExitWithMessage("Failed to write header to Y BAM/CRAM output");
      }
      if (mapping_parameters_.write_index &&
          mapping_parameters_.Y_output_path != "-" &&
          mapping_parameters_.Y_output_path != "/dev/stdout" &&
          mapping_parameters_.Y_output_path != "/dev/stderr") {
        this->Y_index_path_ = mapping_parameters_.Y_output_path;
        if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) {
          this->Y_index_path_ += ".bai";
        } else if (mapping_parameters_.mapping_output_format ==
                   MAPPINGFORMAT_CRAM) {
          this->Y_index_path_ += ".crai";
        }
        int min_shift = 0;
        if (sam_idx_init(Y_hts_out_, hts_hdr_, min_shift,
                         this->Y_index_path_.c_str()) < 0) {
          ExitWithMessage("Failed to initialize index for Y BAM/CRAM output");
        }
      }
    }
  }
  if (mapping_parameters_.macs3_frag_memory_accumulator &&
      mapping_parameters_.macs3_frag_memory_accumulator->IsValid()) {
    std::vector<std::string> chrom_names;
    chrom_names.reserve(num_reference_sequences);
    for (uint32_t i = 0; i < num_reference_sequences; ++i) {
      chrom_names.emplace_back(reference.GetSequenceNameAt(i));
    }
    mapping_parameters_.macs3_frag_memory_accumulator->ResizeForReferenceNames(
        num_reference_sequences, chrom_names);
  }
  if (mapping_parameters_.macs3_frag_chrom_names) {
    auto& names = *mapping_parameters_.macs3_frag_chrom_names;
    names.clear();
    names.reserve(num_reference_sequences);
    for (uint32_t i = 0; i < num_reference_sequences; ++i) {
      names.emplace_back(reference.GetSequenceNameAt(i));
    }
  }
  if (mapping_parameters_.macs3_frag_buffer) {
    mapping_parameters_.macs3_frag_buffer->assign(
        num_reference_sequences, std::vector<macs3::FragmentRecord>());
  }
  OpenAtacEvidenceBinaryOutput(num_reference_sequences, reference);
}

template <>
void MappingWriter<AtacSpillRecord>::AppendMapping(
    uint32_t rid, const SequenceBatch &reference,
    const AtacSpillRecord &mapping) {
  if (!mapping_parameters_.AtacDualFragmentAndBam()) {
    // Sidecar-only AEV1 records, or BED/TagAlign rows; the same code writes
    // the per-reference partitions of the parallel low-memory merge.
    AtacFragmentSink sink;
    if (mapping_parameters_.AtacSidecarOnly()) {
      sink.evidence_fp = atac_evidence_fp_;
      sink.evidence_records = &atac_evidence_records_written_;
    }
    sink.text_fp = mapping_output_file_;
    if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BED &&
        mapping_parameters_.macs3_frag_buffer) {
      sink.buckets = mapping_parameters_.macs3_frag_buffer.get();
      sink.allow_bucket_resize = true;
    }
    AtacLowMemFinalizer::WriteNonDualFragment(*this, &sink, rid, reference,
                                              mapping);
    if (sink.evidence_write_failed) {
      ExitWithMessage("Failed to write ATAC evidence binary record: " +
                      atac_evidence_tmp_path_);
    }
    return;
  }
  if (!mapping.HasBamPairSection()) {
    ExitWithMessage(
        "ATAC dual BAM/CRAM output requires BAM mate payloads on each fragment "
        "record; refusing to emit fragments without mates (overflow schema "
        "mismatch or corrupt spill payload)");
  }
  const PairedEndMappingWithBarcode &frag = mapping;
  uint32_t mapping_end_position = frag.GetEndPosition();
  const bool bulk_no_barcode = !mapping_parameters_.HasBarcodeInput();
  // Both the barcoded 5-column and bulk 7-column fragment schemas expose the
  // collapsed duplicate count. In-memory peak calling must consume the same
  // weight as a file-source or standalone MACS3 run on the emitted fragments.
  const uint32_t peak_count = static_cast<uint32_t>(frag.num_dups_);
  if (atac_evidence_fp_) {
    AppendAtacEvidenceBinaryRecord(rid, frag.GetStartPosition(),
                                   mapping_end_position,
                                   peak_count, frag.cell_barcode_);
  }
  if (atac_fragment_gz_output_ || atac_fragment_output_file_) {
    const char *reference_sequence_name = reference.GetSequenceNameAt(rid);
    if (bulk_no_barcode) {
      const std::string strand = frag.IsPositiveStrand() ? "+" : "-";
      AppendAtacFragmentOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(frag.GetStartPosition()) + "\t" +
          std::to_string(mapping_end_position) + "\tN\t" +
          std::to_string(frag.mapq_) + "\t" + strand + "\t" +
          std::to_string(frag.num_dups_) + "\n");
    } else {
      AppendAtacFragmentOutput(
          std::string(reference_sequence_name) + "\t" +
          std::to_string(frag.GetStartPosition()) + "\t" +
          std::to_string(mapping_end_position) + "\t" +
          barcode_translator_.Translate(frag.cell_barcode_,
                                        cell_barcode_length_) +
          "\t" + std::to_string(frag.num_dups_) + "\n");
    }
  }
  if (mapping_parameters_.macs3_frag_memory_accumulator &&
      mapping_parameters_.macs3_frag_memory_accumulator->IsValid()) {
    std::string acc_err;
    if (!mapping_parameters_.macs3_frag_memory_accumulator->Add(
            rid, frag.GetStartPosition(), mapping_end_position, peak_count,
            &acc_err)) {
      ExitWithMessage("MACS3 FRAG memory accumulator: " + acc_err);
    }
  }
  if (mapping_parameters_.macs3_frag_buffer) {
    macs3::FragmentRecord rec;
    rec.chrom_id = static_cast<int32_t>(rid);
    rec.start = static_cast<int32_t>(frag.GetStartPosition());
    rec.end = static_cast<int32_t>(mapping_end_position);
    rec.count = peak_count;
    if (rec.end > rec.start && rec.count > 0) {
      auto& buckets = *mapping_parameters_.macs3_frag_buffer;
      if (rid >= buckets.size()) {
        buckets.resize(rid + 1);
      }
      buckets[rid].push_back(rec);
    }
  }

  auto write_sam = [&](const SAMMapping &m) {
    bam1_t *b =
        BuildBamRecordFromSamMappingFields(
            mapping_parameters_, cell_barcode_length_, barcode_translator_, rid,
            reference, m);
    if (mapping_parameters_.sort_bam && bam_sorter_) {
      bool hasY =
          mapping.HasYHit() ||
          (reads_with_y_hit_ && reads_with_y_hit_->count(m.read_id_) > 0);
      bam_sorter_->addRecord(b, m.read_id_, hasY);
      bam_destroy1(b);
    } else {
      if (sam_write1(hts_out_, hts_hdr_, b) < 0) {
        bam_destroy1(b);
        ExitWithMessage("Failed to write BAM/CRAM record");
      }
      if ((reads_with_y_hit_ || mapping.HasYHit()) &&
          (noY_hts_out_ || Y_hts_out_)) {
        bool is_y_hit = mapping.HasYHit() ||
                        (reads_with_y_hit_ &&
                         reads_with_y_hit_->count(m.read_id_) > 0);
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
      bam_destroy1(b);
    }
  };

  write_sam(mapping.sam1);
  write_sam(mapping.sam2);
}

template <>
void MappingWriter<AtacSpillRecord>::OutputTempMapping(
    const std::string &temp_mapping_output_file_path,
    uint32_t num_reference_sequences,
    const std::vector<std::vector<AtacSpillRecord>> &mappings) {
  FILE *temp_mapping_output_file =
      fopen(temp_mapping_output_file_path.c_str(), "wb");
  assert(temp_mapping_output_file != NULL);
  for (size_t ri = 0; ri < num_reference_sequences; ++ri) {
    size_t num_mappings = mappings[ri].size();
    fwrite(&num_mappings, sizeof(size_t), 1, temp_mapping_output_file);
    if (mappings[ri].size() > 0) {
      for (size_t mi = 0; mi < num_mappings; ++mi) {
        mappings[ri][mi].WriteToFile(temp_mapping_output_file);
      }
    }
  }
  fclose(temp_mapping_output_file);
}

template <>
void MappingWriter<AtacSpillRecord>::OpenHtsOutput() {
  if (!mapping_parameters_.AtacDualFragmentAndBam()) {
    return;
  }
  if (mapping_parameters_.mapping_output_format != MAPPINGFORMAT_BAM &&
      mapping_parameters_.mapping_output_format != MAPPINGFORMAT_CRAM) {
    return;
  }

  const char *output_path = mapping_parameters_.mapping_output_file_path.c_str();
  const char *hts_mode =
      (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) ? "wb"
                                                                       : "wc";

  hts_out_ = sam_open(output_path, hts_mode);
  if (!hts_out_) {
    ExitWithMessage("Failed to open BAM/CRAM output file: " +
                    mapping_parameters_.mapping_output_file_path);
  }

  int effective_hts_threads = mapping_parameters_.hts_threads;
  if (effective_hts_threads == 0) {
    effective_hts_threads = std::min(mapping_parameters_.num_threads, 4);
  }
  hts_set_threads(hts_out_, effective_hts_threads);

  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM &&
      !mapping_parameters_.reference_file_path.empty()) {
    int ret = hts_set_fai_filename(hts_out_,
                                   mapping_parameters_.reference_file_path.c_str());
    if (ret != 0) {
      ExitWithMessage("Failed to set CRAM reference file: " +
                      mapping_parameters_.reference_file_path);
    }
  }
}

template <>
void MappingWriter<AtacSpillRecord>::FinalizeSortedOutput() {
  if (!mapping_parameters_.AtacDualFragmentAndBam()) {
    return;
  }
  if (!bam_sorter_) return;

  bam_sorter_->finalize();

  bam1_t *b;
  bool hasY;

  while (bam_sorter_->nextRecord(&b, &hasY)) {
    if (sam_write1(hts_out_, hts_hdr_, b) < 0) {
      ExitWithMessage("Failed to write sorted BAM/CRAM record");
    }
    if (hasY && Y_hts_out_) {
      if (sam_write1(Y_hts_out_, hts_hdr_, b) < 0) {
        ExitWithMessage("Failed to write sorted BAM/CRAM record to Y stream");
      }
    }
    if (!hasY && noY_hts_out_) {
      if (sam_write1(noY_hts_out_, hts_hdr_, b) < 0) {
        ExitWithMessage(
            "Failed to write sorted BAM/CRAM record to noY stream");
      }
    }
  }

  bam_sorter_.reset();
}

template <>
void MappingWriter<AtacSpillRecord>::CloseHtsOutput() {
  if (!mapping_parameters_.AtacDualFragmentAndBam()) {
    return;
  }
  if (bam_sorter_) {
    FinalizeSortedOutput();
  }

  if (hts_out_) {
    if (mapping_parameters_.write_index &&
        mapping_parameters_.mapping_output_file_path != "-" &&
        mapping_parameters_.mapping_output_file_path != "/dev/stdout" &&
        mapping_parameters_.mapping_output_file_path != "/dev/stderr") {
      if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) {
        int ret = sam_idx_save(hts_out_);
        if (ret < 0) {
          std::cerr << "Warning: Failed to save BAM index (error code: " << ret
                    << ")\n";
        }
      }
    }
    sam_close(hts_out_);
    hts_out_ = nullptr;
  }
  if (hts_hdr_) {
    bam_hdr_destroy(hts_hdr_);
    hts_hdr_ = nullptr;
  }
  if (noY_hts_out_) {
    if (mapping_parameters_.write_index &&
        mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM &&
        mapping_parameters_.noY_output_path != "-" &&
        mapping_parameters_.noY_output_path != "/dev/stdout" &&
        mapping_parameters_.noY_output_path != "/dev/stderr") {
      int ret = sam_idx_save(noY_hts_out_);
      (void)ret;
    }
    sam_close(noY_hts_out_);
    noY_hts_out_ = nullptr;
  }
  if (Y_hts_out_) {
    if (mapping_parameters_.write_index &&
        mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM &&
        mapping_parameters_.Y_output_path != "-" &&
        mapping_parameters_.Y_output_path != "/dev/stdout" &&
        mapping_parameters_.Y_output_path != "/dev/stderr") {
      int ret = sam_idx_save(Y_hts_out_);
      (void)ret;
    }
    sam_close(Y_hts_out_);
    Y_hts_out_ = nullptr;
  }
}

template <>
void MappingWriter<AtacSpillRecord>::OpenYFilterStreams() {
  if (!mapping_parameters_.AtacDualFragmentAndBam()) {
    return;
  }
  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM ||
      mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM) {
    const char *hts_mode =
        (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) ? "wb"
                                                                         : "wc";
    int effective_hts_threads = mapping_parameters_.hts_threads;
    if (effective_hts_threads == 0) {
      effective_hts_threads = std::min(mapping_parameters_.num_threads, 4);
    }

    if (mapping_parameters_.emit_noY_stream && !noY_hts_out_) {
      noY_hts_out_ =
          sam_open(mapping_parameters_.noY_output_path.c_str(), hts_mode);
      if (!noY_hts_out_) {
        ExitWithMessage("Failed to open noY BAM/CRAM output file: " +
                        mapping_parameters_.noY_output_path);
      }
      hts_set_threads(noY_hts_out_, effective_hts_threads);
      if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM &&
          !mapping_parameters_.reference_file_path.empty()) {
        int ret = hts_set_fai_filename(
            noY_hts_out_, mapping_parameters_.reference_file_path.c_str());
        if (ret != 0) {
          ExitWithMessage(
              "Failed to set CRAM reference file for noY stream: " +
              mapping_parameters_.reference_file_path);
        }
      }
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
        int ret = hts_set_fai_filename(
            Y_hts_out_, mapping_parameters_.reference_file_path.c_str());
        if (ret != 0) {
          ExitWithMessage("Failed to set CRAM reference file for Y stream: " +
                          mapping_parameters_.reference_file_path);
        }
      }
    }
  } else {
    if (mapping_parameters_.emit_noY_stream && !noY_output_file_) {
      noY_output_file_ =
          fopen(mapping_parameters_.noY_output_path.c_str(), "w");
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
void MappingWriter<AtacSpillRecord>::CloseYFilterStreams() {
  if (!mapping_parameters_.AtacDualFragmentAndBam()) {
    return;
  }
  if (noY_hts_out_) {
    if (mapping_parameters_.write_index &&
        mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM &&
        mapping_parameters_.noY_output_path != "-" &&
        mapping_parameters_.noY_output_path != "/dev/stdout" &&
        mapping_parameters_.noY_output_path != "/dev/stderr") {
      int ret = sam_idx_save(noY_hts_out_);
      (void)ret;
    }
    sam_close(noY_hts_out_);
    noY_hts_out_ = nullptr;
  }
  if (Y_hts_out_) {
    if (mapping_parameters_.write_index &&
        mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM &&
        mapping_parameters_.Y_output_path != "-" &&
        mapping_parameters_.Y_output_path != "/dev/stdout" &&
        mapping_parameters_.Y_output_path != "/dev/stderr") {
      int ret = sam_idx_save(Y_hts_out_);
      (void)ret;
    }
    sam_close(Y_hts_out_);
    Y_hts_out_ = nullptr;
  }

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
void MappingWriter<AtacSpillRecord>::BuildHtsHeader(
    uint32_t num_ref_seqs, const SequenceBatch &reference) {
  if (!mapping_parameters_.AtacDualFragmentAndBam()) {
    return;
  }
  std::string header_text = "@HD\tVN:1.6";

  if (mapping_parameters_.sort_bam) {
    header_text += "\tSO:coordinate";
  } else {
    header_text += "\tSO:unknown";
  }
  header_text += "\n";

  for (uint32_t rid = 0; rid < num_ref_seqs; ++rid) {
    header_text += "@SQ\tSN:" + std::string(reference.GetSequenceNameAt(rid)) +
                   "\tLN:" + std::to_string(reference.GetSequenceLengthAt(rid)) +
                   "\n";
  }

  header_text +=
      "@PG\tID:chromap\tPN:chromap\tVN:" + std::string(CHROMAP_VERSION) + "\n";

  if (!mapping_parameters_.read_group_id.empty()) {
    std::string rg_id = mapping_parameters_.read_group_id;
    if (rg_id == "auto") {
      if (!mapping_parameters_.ReadGroupSourcePath().empty()) {
        rg_id = DeriveReadGroupFromFilenameImpl(
            mapping_parameters_.ReadGroupSourcePath());
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

  if (sam_hdr_write(hts_out_, hts_hdr_) < 0) {
    ExitWithMessage("Failed to write BAM/CRAM header");
  }

  if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM &&
      !mapping_parameters_.reference_file_path.empty()) {
    faidx_t *fai = fai_load3(mapping_parameters_.reference_file_path.c_str(),
                             NULL, NULL, FAI_CREATE);
    if (!fai) {
      ExitWithMessage(
          "Failed to load FASTA index for CRAM reference: " +
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
      std::string error_msg =
          "CRAM reference file is missing " +
          std::to_string(missing_contigs.size()) +
          " contig(s) present in header:\n";
      for (size_t j = 0; j < missing_contigs.size() && j < 10; ++j) {
        error_msg += "  " + missing_contigs[j] + "\n";
      }
      if (missing_contigs.size() > 10) {
        error_msg += "  ... and " +
                     std::to_string(missing_contigs.size() - 10) + " more\n";
      }
      error_msg +=
          "Ensure the reference FASTA matches the index used for mapping.";
      ExitWithMessage(error_msg);
    }
  }

  if (mapping_parameters_.write_index) {
    this->hts_index_path_ = mapping_parameters_.mapping_output_file_path;
    int min_shift = 0;
    if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_BAM) {
      this->hts_index_path_ += ".bai";
    } else if (mapping_parameters_.mapping_output_format == MAPPINGFORMAT_CRAM) {
      this->hts_index_path_ += ".crai";
    }
    int ret =
        sam_idx_init(hts_out_, hts_hdr_, min_shift, this->hts_index_path_.c_str());
    if (ret < 0) {
      ExitWithMessage("Failed to initialize BAM/CRAM index");
    }
  }

  if (mapping_parameters_.sort_bam) {
    std::string tmpDir = mapping_parameters_.temp_directory_path;
    if (tmpDir.empty()) {
      const char *env_tmp = std::getenv("TMPDIR");
      if (env_tmp) {
        tmpDir = env_tmp;
      } else {
        tmpDir = "/tmp";
      }
    }
    bam_sorter_.reset(
        new BamSorter(mapping_parameters_.sort_bam_ram_limit, tmpDir));
  }
}

}  // namespace chromap
