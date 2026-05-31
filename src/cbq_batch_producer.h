#ifndef CBQ_BATCH_PRODUCER_H_
#define CBQ_BATCH_PRODUCER_H_

#include <condition_variable>
#include <cstdint>
#include <deque>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "sequence_batch.h"
#include "sequence_effective_range.h"

namespace chromap {

struct CbqPairedEndBatch {
  CbqPairedEndBatch(uint32_t batch_size,
                    const SequenceEffectiveRange &read1_effective_range,
                    const SequenceEffectiveRange &read2_effective_range,
                    const SequenceEffectiveRange &barcode_effective_range);

  SequenceBatch read_batch1;
  SequenceBatch read_batch2;
  SequenceBatch barcode_batch;
  uint32_t num_loaded_pairs = 0;
  uint64_t batch_index = 0;
  uint64_t first_record = 0;
};

class CbqPairedEndBatchProducer {
 public:
  using LoadFn = std::function<uint32_t(CbqPairedEndBatch *)>;

  CbqPairedEndBatchProducer(
      uint32_t batch_size, const SequenceEffectiveRange &read1_effective_range,
      const SequenceEffectiveRange &read2_effective_range,
      const SequenceEffectiveRange &barcode_effective_range, size_t queue_depth,
      LoadFn load_fn);
  CbqPairedEndBatchProducer(const CbqPairedEndBatchProducer &) = delete;
  CbqPairedEndBatchProducer &operator=(const CbqPairedEndBatchProducer &) =
      delete;
  ~CbqPairedEndBatchProducer();

  void Start();
  bool Pop(CbqPairedEndBatch **batch, std::string *error);
  void Release(CbqPairedEndBatch *batch);
  void Stop();

 private:
  CbqPairedEndBatch *AcquireFreeSlot();
  void MarkDone();
  void MarkFailed(const std::string &message);
  void Run();

  LoadFn load_fn_;
  std::vector<std::unique_ptr<CbqPairedEndBatch>> slots_;
  std::deque<CbqPairedEndBatch *> free_slots_;
  std::deque<CbqPairedEndBatch *> ready_slots_;
  std::mutex mutex_;
  std::condition_variable free_cv_;
  std::condition_variable ready_cv_;
  std::thread worker_;
  bool stop_requested_ = false;
  bool done_ = false;
  bool failed_ = false;
  std::string error_;
};

uint32_t PopCbqBatchIntoSequenceBatches(CbqPairedEndBatchProducer *producer,
                                        SequenceBatch &read_batch1,
                                        SequenceBatch &read_batch2,
                                        SequenceBatch &barcode_batch);

class CbqPairedEndRangeBatchProducer {
 public:
  using RangeLoadFn =
      std::function<uint32_t(uint64_t first_record, uint32_t record_count,
                             CbqPairedEndBatch *batch)>;

  CbqPairedEndRangeBatchProducer(
      uint32_t batch_size, uint64_t total_records,
      const SequenceEffectiveRange &read1_effective_range,
      const SequenceEffectiveRange &read2_effective_range,
      const SequenceEffectiveRange &barcode_effective_range,
      size_t worker_count, size_t queue_depth, RangeLoadFn load_fn);
  CbqPairedEndRangeBatchProducer(
      const CbqPairedEndRangeBatchProducer &) = delete;
  CbqPairedEndRangeBatchProducer &operator=(
      const CbqPairedEndRangeBatchProducer &) = delete;
  ~CbqPairedEndRangeBatchProducer();

  void Start();
  bool Pop(CbqPairedEndBatch **batch, std::string *error);
  void Release(CbqPairedEndBatch *batch);
  void Stop();

  uint64_t total_records() const { return total_records_; }
  size_t worker_count() const { return worker_count_; }

 private:
  CbqPairedEndBatch *AcquireFreeSlot();
  void MarkWorkerDone();
  void MarkFailed(const std::string &message);
  void RunWorker();

  uint32_t batch_size_ = 0;
  uint64_t total_records_ = 0;
  uint64_t total_batches_ = 0;
  size_t worker_count_ = 1;
  RangeLoadFn load_fn_;
  std::vector<std::unique_ptr<CbqPairedEndBatch>> slots_;
  std::deque<CbqPairedEndBatch *> free_slots_;
  std::map<uint64_t, CbqPairedEndBatch *> ready_slots_;
  std::vector<std::thread> workers_;
  std::mutex mutex_;
  std::condition_variable free_cv_;
  std::condition_variable ready_cv_;
  uint64_t next_batch_to_assign_ = 0;
  uint64_t next_batch_to_pop_ = 0;
  size_t running_workers_ = 0;
  bool started_ = false;
  bool stop_requested_ = false;
  bool done_ = false;
  bool failed_ = false;
  std::string error_;
};

uint32_t PopCbqBatchIntoSequenceBatches(
    CbqPairedEndRangeBatchProducer *producer, SequenceBatch &read_batch1,
    SequenceBatch &read_batch2, SequenceBatch &barcode_batch);

}  // namespace chromap

#endif  // CBQ_BATCH_PRODUCER_H_
