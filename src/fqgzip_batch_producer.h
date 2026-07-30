#ifndef FQGZIP_BATCH_PRODUCER_H_
#define FQGZIP_BATCH_PRODUCER_H_

#include <condition_variable>
#include <cstdint>
#include <deque>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "sequence_batch.h"
#include "sequence_effective_range.h"

namespace chromap {

struct FqgzipPairedEndBatch {
  FqgzipPairedEndBatch(uint32_t batch_size,
                       const SequenceEffectiveRange &read1_effective_range,
                       const SequenceEffectiveRange &read2_effective_range);

  SequenceBatch read_batch1;
  SequenceBatch read_batch2;
  uint32_t num_loaded_pairs = 0;
};

// Ordered Chromap batch producer backed by fqgzip's progressive paired plan.
// fqgzip independently decodes/indexes R1 and R2; this adapter joins them only
// at fqgzip's validated paired-record boundary and leaves mapping unchanged.
class FqgzipPairedEndBatchProducer {
 public:
  FqgzipPairedEndBatchProducer(
      const std::string &read1_path, const std::string &read2_path,
      uint32_t batch_size,
      const SequenceEffectiveRange &read1_effective_range,
      const SequenceEffectiveRange &read2_effective_range,
      uint64_t global_record_offset, uint32_t requested_shards,
      uint32_t worker_threads, size_t queue_depth);
  FqgzipPairedEndBatchProducer(const FqgzipPairedEndBatchProducer &) = delete;
  FqgzipPairedEndBatchProducer &operator=(
      const FqgzipPairedEndBatchProducer &) = delete;
  ~FqgzipPairedEndBatchProducer();

  void Start();
  bool Pop(FqgzipPairedEndBatch **batch, std::string *error);
  void Release(FqgzipPairedEndBatch *batch);
  void Stop();

 private:
  struct ReaderState;

  FqgzipPairedEndBatch *AcquireFreeSlot();
  uint32_t LoadBatch(FqgzipPairedEndBatch *batch);
  void MarkDone();
  void MarkFailed(const std::string &message);
  void Run();

  std::string read1_path_;
  std::string read2_path_;
  uint32_t batch_size_;
  uint64_t global_record_offset_;
  uint32_t requested_shards_;
  uint32_t worker_threads_;
  std::unique_ptr<ReaderState> reader_;
  std::vector<std::unique_ptr<FqgzipPairedEndBatch>> slots_;
  std::deque<FqgzipPairedEndBatch *> free_slots_;
  std::deque<FqgzipPairedEndBatch *> ready_slots_;
  std::mutex mutex_;
  std::condition_variable free_cv_;
  std::condition_variable ready_cv_;
  std::thread worker_;
  bool stop_requested_ = false;
  bool done_ = false;
  bool failed_ = false;
  std::string error_;
};

uint32_t PopFqgzipBatchIntoSequenceBatches(
    FqgzipPairedEndBatchProducer *producer, SequenceBatch &read_batch1,
    SequenceBatch &read_batch2);

}  // namespace chromap

#endif  // FQGZIP_BATCH_PRODUCER_H_
