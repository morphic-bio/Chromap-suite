#ifndef FASTQ_SPLIT_WRITER_H_
#define FASTQ_SPLIT_WRITER_H_

#include <cstring>
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
        gzFile file = gzopen(path.c_str(), "wb");
        if (!file) {
          ExitWithMessage("Failed to open Y FASTQ output file: " + path);
        }
        if (gzbuffer(file, 1 << 20) != 0) {
          gzclose(file);
          ExitWithMessage("Failed to set gzip buffer for Y FASTQ output file: " + path);
        }
        y_files_.push_back(file);
      } else {
        y_file_streams_.emplace_back(new std::ofstream(path, std::ios::out | std::ios::binary));
        if (!y_file_streams_.back()->is_open()) {
          ExitWithMessage("Failed to open Y FASTQ output file: " + path);
        }
      }
    }
    
    // Open noY output files
    for (const auto &path : noy_output_paths) {
      if (compression_ == "gz") {
        gzFile file = gzopen(path.c_str(), "wb");
        if (!file) {
          ExitWithMessage("Failed to open noY FASTQ output file: " + path);
        }
        if (gzbuffer(file, 1 << 20) != 0) {
          gzclose(file);
          ExitWithMessage("Failed to set gzip buffer for noY FASTQ output file: " + path);
        }
        noy_files_.push_back(file);
      } else {
        noy_file_streams_.emplace_back(new std::ofstream(path, std::ios::out | std::ios::binary));
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

    record_buffer_.clear();
    record_buffer_.reserve(1 + std::strlen(name) + comment_len + seq_len +
                           qual_len + 8);
    record_buffer_ += (has_qual ? "@" : ">");
    record_buffer_ += name;
    if (comment_len > 0) {
      record_buffer_ += " ";
      record_buffer_.append(comment, comment_len);
    }
    record_buffer_ += "\n";
    record_buffer_.append(seq, seq_len);
    record_buffer_ += "\n";

    if (has_qual) {
      record_buffer_ += "+\n";
      record_buffer_.append(qual, qual_len);
      record_buffer_ += "\n";
    }

    if (gz_file) {
      const int written = gzwrite(gz_file, record_buffer_.data(),
                                  static_cast<unsigned int>(record_buffer_.size()));
      if (written != static_cast<int>(record_buffer_.size())) {
        ExitWithMessage("Failed to write Y/noY FASTQ gzip record");
      }
    } else {
      file_stream->write(record_buffer_.data(),
                         static_cast<std::streamsize>(record_buffer_.size()));
      if (!file_stream->good()) {
        ExitWithMessage("Failed to write Y/noY FASTQ record");
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
  std::string record_buffer_;
  uint64_t num_written_y_;
  uint64_t num_written_noy_;
};

}  // namespace chromap

#endif  // FASTQ_SPLIT_WRITER_H_
