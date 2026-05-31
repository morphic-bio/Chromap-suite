#ifndef SEQUENCEBATCH_H_
#define SEQUENCEBATCH_H_

#include <unistd.h>
#include <zlib.h>

#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "kseq.h"
#include "sequence_effective_range.h"
#include "utils.h"

namespace chromap {

class SequenceBatch {
 public:
  KSEQ_INIT(gzFile, gzread);

  // When 'max_num_sequences' is not specified. This batch can be used to load
  // any number of sequences with a positive full effective range.
  SequenceBatch() : effective_range_(SequenceEffectiveRange()) {}

  // Construct once and use update sequences when loading each batch.
  SequenceBatch(uint32_t max_num_sequences,
                const SequenceEffectiveRange &effective_range)
      : max_num_sequences_(max_num_sequences),
        effective_range_(effective_range) {
    sequence_batch_.reserve(max_num_sequences_);
    for (uint32_t i = 0; i < max_num_sequences_; ++i) {
      sequence_batch_.emplace_back((kseq_t *)calloc(1, sizeof(kseq_t)));
      sequence_batch_.back()->f = NULL;
    }
    negative_sequence_batch_.assign(max_num_sequences_, "");
    negative_sequence_prepared_.assign(max_num_sequences_, 0);
  }

  ~SequenceBatch() {
    if (sequence_batch_.size() > 0) {
      for (uint32_t i = 0; i < sequence_batch_.size(); ++i) {
        kseq_destroy(sequence_batch_[i]);
      }
    }
  }

  inline uint64_t GetNumSequences() const { return num_loaded_sequences_; }

  inline uint32_t GetMaxBatchSize() const { return max_num_sequences_; }

  inline uint64_t GetNumBases() const { return num_bases_; }

  inline std::vector<kseq_t *> &GetSequenceBatch() { return sequence_batch_; }

  inline std::vector<std::string> &GetNegativeSequenceBatch() {
    return negative_sequence_batch_;
  }

  inline std::vector<uint8_t> &GetNegativeSequencePreparedBatch() {
    return negative_sequence_prepared_;
  }

  inline const char *GetSequenceAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->seq.s;
  }

  inline uint32_t GetSequenceLengthAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->seq.l;
  }

  inline const char *GetSequenceNameAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->name.s;
  }

  inline uint32_t GetSequenceNameLengthAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->name.l;
  }

  inline const char *GetSequenceCommentAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->comment.s;
  }

  inline uint32_t GetSequenceCommentLengthAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->comment.l;
  }

  inline const char *GetSequenceQualAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->qual.s;
  }

  inline uint32_t GetSequenceQualLengthAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->qual.l;
  }
  inline uint32_t GetSequenceIdAt(uint32_t sequence_index) const {
    return sequence_batch_[sequence_index]->id;
  }

  inline const std::string &GetNegativeSequenceAt(
      uint32_t sequence_index) const {
    return negative_sequence_batch_[sequence_index];
  }

  inline bool IsNegativeSequencePreparedAt(uint32_t sequence_index) const {
    return sequence_index < negative_sequence_prepared_.size() &&
           negative_sequence_prepared_[sequence_index] != 0;
  }

  inline bool HasFullPositiveEffectiveRange() const {
    return effective_range_.IsFullRangeAndPositiveStrand();
  }

  // big_endian: N_pos is in the order of sequence
  // little_endian: N_pos is in the order from the sequence right side to left,
  //                this is the order of the GenerateSeed
  // e.g: If the sequence is "ACN", big endian returns N at 2,
  //      little endian returns N at 0.
  inline void GetSequenceNsAt(uint32_t sequence_index, bool little_endian,
                              std::vector<int> &N_pos) {
    const int l = sequence_batch_[sequence_index]->seq.l;
    const char *s = sequence_batch_[sequence_index]->seq.s;
    N_pos.clear();
    if (little_endian) {
      for (int i = l - 1; i >= 0; --i) {
        if (s[i] == 'N') N_pos.push_back(l - 1 - i);
      }
    } else {
      for (int i = 0; i < l; ++i) {
        if (s[i] == 'N') N_pos.push_back(i);
      }
    }
  }

  inline bool IsNInSequenceAt(uint32_t sequence_index) {
    const int l = sequence_batch_[sequence_index]->seq.l;
    const char *s = sequence_batch_[sequence_index]->seq.s;
    for (int i = 0 ; i < l ; ++i)
      if (s[i] == 'N')
        return true;
    return false;
  }

  //  inline char GetReverseComplementBaseOfSequenceAt(uint32_t sequence_index,
  //  uint32_t position) {
  //    kseq_t *sequence = sequence_batch_[sequence_index];
  //    return Uint8ToChar(((uint8_t)3) ^
  //    (CharToUint8((sequence->seq.s)[sequence->seq.l - position - 1])));
  //  }

  inline void PrepareNegativeSequenceAt(uint32_t sequence_index) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    uint32_t sequence_length = sequence->seq.l;
    std::string &negative_sequence = negative_sequence_batch_[sequence_index];
    static const std::array<char, 256> complement_lookup = [] {
      std::array<char, 256> table = {};
      table.fill('N');
      table[static_cast<unsigned char>('A')] = 'T';
      table[static_cast<unsigned char>('C')] = 'G';
      table[static_cast<unsigned char>('G')] = 'C';
      table[static_cast<unsigned char>('T')] = 'A';
      table[static_cast<unsigned char>('N')] = 'N';
      table[static_cast<unsigned char>('a')] = 'T';
      table[static_cast<unsigned char>('c')] = 'G';
      table[static_cast<unsigned char>('g')] = 'C';
      table[static_cast<unsigned char>('t')] = 'A';
      table[static_cast<unsigned char>('n')] = 'N';
      return table;
    }();
    negative_sequence.resize(sequence_length);
    for (uint32_t i = 0; i < sequence_length; ++i) {
      negative_sequence[i] = complement_lookup[static_cast<unsigned char>(
          sequence->seq.s[sequence_length - i - 1])];
    }
    if (sequence_index < negative_sequence_prepared_.size()) {
      negative_sequence_prepared_[sequence_index] = 1;
    }
  }

  inline char *PrepareLoadedNegativeSequenceBuffer(uint32_t sequence_index,
                                                   size_t seq_len) {
    if (sequence_index >= negative_sequence_batch_.size()) {
      ExitWithMessage("Sequence index exceeds batch capacity");
    }
    std::string &negative_sequence = negative_sequence_batch_[sequence_index];
    negative_sequence.resize(seq_len);
    if (sequence_index < negative_sequence_prepared_.size()) {
      negative_sequence_prepared_[sequence_index] = 0;
    }
    return seq_len == 0 ? nullptr : &negative_sequence[0];
  }

  inline void CommitLoadedNegativeSequenceBuffer(uint32_t sequence_index) {
    if (sequence_index >= negative_sequence_batch_.size()) {
      ExitWithMessage("Sequence index exceeds batch capacity");
    }
    if (sequence_index < negative_sequence_prepared_.size()) {
      negative_sequence_prepared_[sequence_index] = 1;
    }
  }

  inline void TrimSequenceAt(uint32_t sequence_index, int length_after_trim) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    if (length_after_trim >= (int)sequence->seq.l) {
      return;
    }

    if (IsNegativeSequencePreparedAt(sequence_index)) {
      negative_sequence_batch_[sequence_index].erase(
          negative_sequence_batch_[sequence_index].begin(),
          negative_sequence_batch_[sequence_index].begin() + sequence->seq.l -
              length_after_trim);
    }

    sequence->seq.l = length_after_trim;
    sequence->seq.s[sequence->seq.l] = '\0';
    sequence->qual.l = length_after_trim;
    sequence->qual.s[sequence->qual.l] = '\0';
  }

  inline void SwapSequenceBatch(SequenceBatch &batch) {
    sequence_batch_.swap(batch.GetSequenceBatch());
    negative_sequence_batch_.swap(batch.GetNegativeSequenceBatch());
    negative_sequence_prepared_.swap(batch.GetNegativeSequencePreparedBatch());
  }

  void InitializeLoading(const std::string &sequence_file_path);

  void FinalizeLoading();

  void ResetLoadedSequences();

  void AssignLoadedSequence(uint32_t sequence_index, const char *name,
                            size_t name_len, const char *comment,
                            size_t comment_len, const char *seq,
                            size_t seq_len, const char *qual,
                            size_t qual_len);

  char *PrepareLoadedSequenceBuffer(uint32_t sequence_index, size_t seq_len);

  void CommitLoadedSequenceBuffer(uint32_t sequence_index, const char *name,
                                  size_t name_len, const char *comment,
                                  size_t comment_len, size_t seq_len,
                                  const char *qual, size_t qual_len);

  void CommitLoadedSequenceBufferWithId(
      uint32_t sequence_index, uint32_t sequence_id, const char *name,
      size_t name_len, const char *comment, size_t comment_len, size_t seq_len,
      const char *qual, size_t qual_len);

  // The func should never override other sequences rather than the last, which
  // means 'sequence_index' cannot be smaller than 'num_loaded_sequences_' - 1.
  // Return true when reaching the end of the file.
  bool LoadOneSequenceAndSaveAt(uint32_t sequence_index);

  // Return the number of sequences loaded into the batch and return 0 if there
  // is no more sequences. This func now is only used to load barcodes.
  uint32_t LoadBatch();

  // Load all sequences in a file. This function should only be used to load
  // reference. And once the reference is loaded, the batch should never be
  // updated. This func is slow when there are large number of sequences.
  void LoadAllSequences();

  inline void CorrectBaseAt(uint32_t sequence_index, uint32_t base_position,
                            char correct_base) {
    kseq_t *sequence = sequence_batch_[sequence_index];
    sequence->seq.s[base_position] = correct_base;
    if (sequence_index < negative_sequence_prepared_.size()) {
      negative_sequence_prepared_[sequence_index] = 0;
    }
  }

  inline uint64_t GenerateSeedFromSequenceAt(uint32_t sequence_index,
                                             uint32_t start_position,
                                             uint32_t seed_length) const {
    const char *sequence = GetSequenceAt(sequence_index);
    const uint32_t sequence_length = GetSequenceLengthAt(sequence_index);
    return GenerateSeedFromSequence(sequence, sequence_length, start_position,
                                    seed_length);
  }

  inline void ReorderSequences(const std::vector<int> &rid_rank) {
    std::vector<kseq_t *> tmp_sequence_batch_ = sequence_batch_;
    std::vector<std::string> tmp_negative_sequence_batch_ =
        negative_sequence_batch_;
    std::vector<uint8_t> tmp_negative_sequence_prepared =
        negative_sequence_prepared_;
    for (size_t i = 0; i < sequence_batch_.size(); ++i) {
      sequence_batch_[rid_rank[i]] = tmp_sequence_batch_[i];
    }

    if (negative_sequence_batch_.size() > 0) {
      for (size_t i = 0; i < sequence_batch_.size(); ++i) {
        negative_sequence_batch_[rid_rank[i]] = tmp_negative_sequence_batch_[i];
      }
    }
    if (negative_sequence_prepared_.size() > 0) {
      for (size_t i = 0; i < sequence_batch_.size(); ++i) {
        negative_sequence_prepared_[rid_rank[i]] =
            tmp_negative_sequence_prepared[i];
      }
    }
  }

 protected:
  // When 'is_seq' is set to true, this func will complement the base when
  // necessary. Otherwise, it will just reverse the sequence.
  void ReplaceByEffectiveRange(kstring_t &seq, bool is_seq);

  void AssignKstring(kstring_t &dest, const char *data, size_t len);

  // This is the accumulated number of sequences that have ever been loaded into
  // the batch. It is useful for tracking read ids.
  uint32_t total_num_loaded_sequences_ = 0;

  // This is the number of sequences loaded into the current batch.
  uint32_t num_loaded_sequences_ = 0;

  // This is the number of bases loaded into the current batch. It is only
  // populated for the reference.
  uint64_t num_bases_ = 0;

  // This is the max number of sequences that can be loaded into the batch. It
  // is set to 0 when there is no such restriction.
  uint32_t max_num_sequences_ = 0;

  gzFile sequence_file_;
  kseq_t *sequence_kseq_ = nullptr;
  std::vector<kseq_t *> sequence_batch_;

  // TODO: avoid constructing the negative sequence batch.
  std::vector<std::string> negative_sequence_batch_;
  std::vector<uint8_t> negative_sequence_prepared_;

  // Actual range within each sequence.
  const SequenceEffectiveRange effective_range_;
};

}  // namespace chromap

#endif  // SEQUENCEBATCH_H_
