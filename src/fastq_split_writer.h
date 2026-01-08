#ifndef FASTQ_SPLIT_WRITER_H_
#define FASTQ_SPLIT_WRITER_H_

#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <zlib.h>

#include "sequence_batch.h"
#include "utils.h"

namespace chromap {

// Helper class to write Y/noY FASTQ files
class FastqSplitWriter {
 public:
  FastqSplitWriter(const std::vector<std::string> &y_output_paths,
                   const std::vector<std::string> &noy_output_paths,
                   const std::string &compression)
      : compression_(compression), num_written_y_(0), num_written_noy_(0) {
    // Open Y output files
    for (const auto &path : y_output_paths) {
      if (compression_ == "gz") {
        y_files_.push_back(gzopen(path.c_str(), "wb"));
        if (!y_files_.back()) {
          ExitWithMessage("Failed to open Y FASTQ output file: " + path);
        }
      } else {
        y_file_streams_.emplace_back(new std::ofstream(path));
        if (!y_file_streams_.back()->is_open()) {
          ExitWithMessage("Failed to open Y FASTQ output file: " + path);
        }
      }
    }
    
    // Open noY output files
    for (const auto &path : noy_output_paths) {
      if (compression_ == "gz") {
        noy_files_.push_back(gzopen(path.c_str(), "wb"));
        if (!noy_files_.back()) {
          ExitWithMessage("Failed to open noY FASTQ output file: " + path);
        }
      } else {
        noy_file_streams_.emplace_back(new std::ofstream(path));
        if (!noy_file_streams_.back()->is_open()) {
          ExitWithMessage("Failed to open noY FASTQ output file: " + path);
        }
      }
    }
  }

  ~FastqSplitWriter() {
    // Close gzip files
    for (gzFile f : y_files_) {
      if (f) gzclose(f);
    }
    for (gzFile f : noy_files_) {
      if (f) gzclose(f);
    }
    
    // Close regular files
    for (auto &stream : y_file_streams_) {
      if (stream && stream->is_open()) {
        stream->close();
      }
    }
    for (auto &stream : noy_file_streams_) {
      if (stream && stream->is_open()) {
        stream->close();
      }
    }
  }

  // Write a single-end read to appropriate stream
  void WriteRead(uint32_t read_index, const SequenceBatch &read_batch,
                 bool has_y_hit, uint32_t mate_index = 0) {
    const char *name = read_batch.GetSequenceNameAt(read_index);
    if (!name || name[0] == '\0') {
      return;
    }
    const char *comment = read_batch.GetSequenceCommentAt(read_index);
    uint32_t comment_len = read_batch.GetSequenceCommentLengthAt(read_index);
    const char *seq = read_batch.GetSequenceAt(read_index);
    uint32_t seq_len = read_batch.GetSequenceLengthAt(read_index);
    const char *qual = read_batch.GetSequenceQualAt(read_index);
    uint32_t qual_len = read_batch.GetSequenceQualLengthAt(read_index);
    bool has_qual = (qual != nullptr && qual_len > 0);
    
    // Select output stream based on Y-hit status
    gzFile gz_file = nullptr;
    std::ofstream *file_stream = nullptr;
    if (has_y_hit) {
      if (mate_index < y_files_.size()) {
        gz_file = y_files_[mate_index];
      } else if (mate_index < y_file_streams_.size()) {
        file_stream = y_file_streams_[mate_index].get();
      }
    } else {
      if (mate_index < noy_files_.size()) {
        gz_file = noy_files_[mate_index];
      } else if (mate_index < noy_file_streams_.size()) {
        file_stream = noy_file_streams_[mate_index].get();
      }
    }
    
    if (!gz_file && !file_stream) return;
    
    // Build header: name + comment
    std::string header;
    header += (has_qual ? "@" : ">");
    header += name;
    if (comment_len > 0) {
      header += " ";
      header.append(comment, comment_len);
    }
    header += "\n";
    
    // Write header
    if (gz_file) {
      gzwrite(gz_file, header.c_str(), header.length());
    } else {
      *file_stream << header;
    }
    
    // Write sequence
    if (gz_file) {
      gzwrite(gz_file, seq, seq_len);
      gzputc(gz_file, '\n');
    } else {
      file_stream->write(seq, seq_len);
      *file_stream << "\n";
    }
    
    // Write quality if present
    if (has_qual) {
      if (gz_file) {
        gzputc(gz_file, '+');
        gzputc(gz_file, '\n');
        gzwrite(gz_file, qual, qual_len);
        gzputc(gz_file, '\n');
      } else {
        *file_stream << "+\n";
        file_stream->write(qual, qual_len);
        *file_stream << "\n";
      }
    }
    
    if (has_y_hit) {
      ++num_written_y_;
    } else {
      ++num_written_noy_;
    }
  }

  // Write paired-end reads (both mates go to same stream if either has Y hit)
  void WritePairedReads(uint32_t pair_index,
                        const SequenceBatch &read_batch1,
                        const SequenceBatch &read_batch2,
                        bool has_y_hit) {
    WriteRead(pair_index, read_batch1, has_y_hit, 0);
    WriteRead(pair_index, read_batch2, has_y_hit, 1);
  }

  uint64_t GetNumWrittenY() const { return num_written_y_; }
  uint64_t GetNumWrittenNoY() const { return num_written_noy_; }

 private:
  std::string compression_;
  std::vector<gzFile> y_files_;
  std::vector<gzFile> noy_files_;
  std::vector<std::unique_ptr<std::ofstream>> y_file_streams_;
  std::vector<std::unique_ptr<std::ofstream>> noy_file_streams_;
  uint64_t num_written_y_;
  uint64_t num_written_noy_;
};

}  // namespace chromap

#endif  // FASTQ_SPLIT_WRITER_H_
