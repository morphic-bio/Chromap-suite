#include "sequence_batch.h"

#include <tuple>
#include <cstdlib>
#include <cstring>

#include "utils.h"

namespace chromap {

void SequenceBatch::InitializeLoading(const std::string &sequence_file_path) {
  sequence_file_ = gzopen(sequence_file_path.c_str(), "r");
  if (sequence_file_ == NULL) {
    ExitWithMessage("Cannot find sequence file " + sequence_file_path);
  }
  sequence_kseq_ = kseq_init(sequence_file_);
}

void SequenceBatch::FinalizeLoading() {
  kseq_destroy(sequence_kseq_);
  gzclose(sequence_file_);
}

void SequenceBatch::ResetLoadedSequences() {
  num_loaded_sequences_ = 0;
  num_bases_ = 0;
}

void SequenceBatch::AssignKstring(kstring_t &dest, const char *data,
                                  size_t len) {
  if (dest.m < len + 1) {
    dest.m = len + 1;
    dest.s = static_cast<char *>(realloc(dest.s, dest.m));
    if (dest.s == NULL) {
      ExitWithMessage("Failed to allocate sequence storage");
    }
  }
  if (len > 0 && data != NULL) {
    memcpy(dest.s, data, len);
  }
  dest.s[len] = '\0';
  dest.l = len;
}

void SequenceBatch::AssignLoadedSequence(uint32_t sequence_index,
                                         const char *name, size_t name_len,
                                         const char *comment,
                                         size_t comment_len, const char *seq,
                                         size_t seq_len, const char *qual,
                                         size_t qual_len) {
  if (sequence_index >= sequence_batch_.size()) {
    ExitWithMessage("Sequence index exceeds batch capacity");
  }
  if (sequence_index > num_loaded_sequences_) {
    ExitWithMessage("Sequence batches must be filled in order");
  }

  kseq_t *sequence = sequence_batch_[sequence_index];
  AssignKstring(sequence->seq, seq, seq_len);
  ReplaceByEffectiveRange(sequence->seq, /*is_seq=*/true);
  AssignKstring(sequence->name, name, name_len);
  AssignKstring(sequence->comment, comment, comment_len);
  if (qual != NULL && qual_len > 0) {
    AssignKstring(sequence->qual, qual, qual_len);
    ReplaceByEffectiveRange(sequence->qual, /*is_seq=*/false);
  } else {
    AssignKstring(sequence->qual, NULL, 0);
  }
  sequence->id = total_num_loaded_sequences_;
  ++total_num_loaded_sequences_;
  if (sequence_index == num_loaded_sequences_) {
    ++num_loaded_sequences_;
  }
}

bool SequenceBatch::LoadOneSequenceAndSaveAt(uint32_t sequence_index) {
  if (sequence_index == 0) {
    num_loaded_sequences_ = 0;
  }

  int length = kseq_read(sequence_kseq_);
  while (length == 0) {
    length = kseq_read(sequence_kseq_);
  }

  if (length > 0) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    std::swap(sequence_kseq_->seq, sequence->seq);
    ReplaceByEffectiveRange(sequence->seq, /*is_seq=*/true);
    std::swap(sequence_kseq_->name, sequence->name);
    std::swap(sequence_kseq_->comment, sequence->comment);
    sequence->id = total_num_loaded_sequences_;
    ++total_num_loaded_sequences_;

    if (sequence_index >= num_loaded_sequences_) {
      ++num_loaded_sequences_;
    } else if (sequence_index + 1 != num_loaded_sequences_) {
      std::cerr << sequence_index << " " << num_loaded_sequences_ << "\n";
      ExitWithMessage(
          "Shouldn't override other sequences rather than the last!");
    }

    if (sequence_kseq_->qual.l != 0) {  // fastq file
      std::swap(sequence_kseq_->qual, sequence->qual);
      ReplaceByEffectiveRange(sequence->qual, /*is_seq=*/false);
    }
    return false;
  }

  // Make sure to reach the end of the file rather than meet an error.
  if (length != -1) {
    ExitWithMessage(
        "Didn't reach the end of sequence file, which might be corrupted!");
  }
  return true;
}

uint32_t SequenceBatch::LoadBatch() {
  double real_start_time = GetRealTime();
  num_loaded_sequences_ = 0;
  for (uint32_t sequence_index = 0; sequence_index < max_num_sequences_;
       ++sequence_index) {
    if (LoadOneSequenceAndSaveAt(sequence_index)) {
      break;
    }
  }

  if (num_loaded_sequences_ != 0) {
    std::cerr << "Loaded sequence batch successfully in "
              << GetRealTime() - real_start_time << "s, ";
    std::cerr << "number of sequences: " << num_loaded_sequences_ << ".\n";
  } else {
    std::cerr << "No more sequences.\n";
  }
  return num_loaded_sequences_;
}

void SequenceBatch::LoadAllSequences() {
  double real_start_time = GetRealTime();
  sequence_batch_.reserve(200);
  num_loaded_sequences_ = 0;
  num_bases_ = 0;
  int length = kseq_read(sequence_kseq_);
  while (length >= 0) {
    if (length > 0) {
      sequence_batch_.emplace_back((kseq_t *)calloc(1, sizeof(kseq_t)));
      kseq_t *sequence = sequence_batch_.back();
      std::swap(sequence_kseq_->seq, sequence->seq);
      ReplaceByEffectiveRange(sequence->seq, /*is_seq=*/true);
      std::swap(sequence_kseq_->name, sequence->name);
      std::swap(sequence_kseq_->comment, sequence->comment);
      if (sequence_kseq_->qual.l != 0) {  // fastq file
        std::swap(sequence_kseq_->qual, sequence->qual);
        ReplaceByEffectiveRange(sequence->qual, /*is_seq=*/false);
      }
      sequence->id = total_num_loaded_sequences_;
      ++total_num_loaded_sequences_;
      ++num_loaded_sequences_;
      num_bases_ += length;
    }
    length = kseq_read(sequence_kseq_);
  }

  // Make sure to reach the end of the file rather than meet an error.
  if (length != -1) {
    ExitWithMessage(
        "Didn't reach the end of sequence file, which might be corrupted!");
  }

  std::cerr << "Loaded all sequences successfully in "
            << GetRealTime() - real_start_time << "s, ";
  std::cerr << "number of sequences: " << num_loaded_sequences_ << ", ";
  std::cerr << "number of bases: " << num_bases_ << ".\n";
}

void SequenceBatch::ReplaceByEffectiveRange(kstring_t &seq, bool is_seq) {
  seq.l = effective_range_.Replace(seq.s, seq.l, is_seq);
}

}  // namespace chromap
