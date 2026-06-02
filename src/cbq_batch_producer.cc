#include "cbq_batch_producer.h"

#include <exception>
#include <limits>
#include <utility>

#include "utils.h"

namespace chromap {

CbqPairedEndBatch::CbqPairedEndBatch(
    uint32_t batch_size, const SequenceEffectiveRange &read1_effective_range,
    const SequenceEffectiveRange &read2_effective_range,
    const SequenceEffectiveRange &barcode_effective_range)
    : read_batch1(batch_size, read1_effective_range),
      read_batch2(batch_size, read2_effective_range),
      barcode_batch(batch_size, barcode_effective_range) {}

CbqPairedEndBatchProducer::CbqPairedEndBatchProducer(
    uint32_t batch_size, const SequenceEffectiveRange &read1_effective_range,
    const SequenceEffectiveRange &read2_effective_range,
    const SequenceEffectiveRange &barcode_effective_range, size_t queue_depth,
    LoadFn load_fn)
    : load_fn_(std::move(load_fn)) {
  if (queue_depth == 0) {
    queue_depth = 1;
  }
  slots_.reserve(queue_depth);
  for (size_t i = 0; i < queue_depth; ++i) {
    slots_.emplace_back(new CbqPairedEndBatch(batch_size,
                                              read1_effective_range,
                                              read2_effective_range,
                                              barcode_effective_range));
    free_slots_.push_back(slots_.back().get());
  }
}

CbqPairedEndBatchProducer::~CbqPairedEndBatchProducer() { Stop(); }

void CbqPairedEndBatchProducer::Start() {
  worker_ = std::thread(&CbqPairedEndBatchProducer::Run, this);
}

bool CbqPairedEndBatchProducer::Pop(CbqPairedEndBatch **batch,
                                    std::string *error) {
  std::unique_lock<std::mutex> lock(mutex_);
  ready_cv_.wait(lock, [this] {
    return failed_ || done_ || !ready_slots_.empty();
  });
  if (failed_) {
    if (error != nullptr) {
      *error = error_;
    }
    return false;
  }
  if (ready_slots_.empty()) {
    *batch = nullptr;
    return true;
  }
  *batch = ready_slots_.front();
  ready_slots_.pop_front();
  return true;
}

void CbqPairedEndBatchProducer::Release(CbqPairedEndBatch *batch) {
  if (batch == nullptr) {
    return;
  }
  {
    std::lock_guard<std::mutex> lock(mutex_);
    free_slots_.push_back(batch);
  }
  free_cv_.notify_one();
}

void CbqPairedEndBatchProducer::Stop() {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    stop_requested_ = true;
  }
  free_cv_.notify_all();
  ready_cv_.notify_all();
  if (worker_.joinable()) {
    worker_.join();
  }
}

CbqPairedEndBatch *CbqPairedEndBatchProducer::AcquireFreeSlot() {
  std::unique_lock<std::mutex> lock(mutex_);
  free_cv_.wait(lock, [this] {
    return stop_requested_ || !free_slots_.empty();
  });
  if (stop_requested_) {
    return nullptr;
  }
  CbqPairedEndBatch *slot = free_slots_.front();
  free_slots_.pop_front();
  return slot;
}

void CbqPairedEndBatchProducer::MarkDone() {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    done_ = true;
  }
  ready_cv_.notify_all();
}

void CbqPairedEndBatchProducer::MarkFailed(const std::string &message) {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    failed_ = true;
    done_ = true;
    error_ = message;
  }
  ready_cv_.notify_all();
  free_cv_.notify_all();
}

void CbqPairedEndBatchProducer::Run() {
  try {
    for (;;) {
      CbqPairedEndBatch *slot = AcquireFreeSlot();
      if (slot == nullptr) {
        return;
      }
      slot->num_loaded_pairs = load_fn_(slot);
      if (slot->num_loaded_pairs == 0) {
        Release(slot);
        MarkDone();
        return;
      }
      {
        std::lock_guard<std::mutex> lock(mutex_);
        ready_slots_.push_back(slot);
      }
      ready_cv_.notify_one();
    }
  } catch (const std::exception &ex) {
    MarkFailed(ex.what());
  } catch (...) {
    MarkFailed("unknown CBQ producer failure");
  }
}

uint32_t PopCbqBatchIntoSequenceBatches(CbqPairedEndBatchProducer *producer,
                                        SequenceBatch &read_batch1,
                                        SequenceBatch &read_batch2,
                                        SequenceBatch &barcode_batch) {
  CbqPairedEndBatch *batch = nullptr;
  std::string error;
  if (!producer->Pop(&batch, &error)) {
    ExitWithMessage("CBQ producer failed: " + error);
  }
  if (batch == nullptr) {
    return 0;
  }

  const uint32_t num_loaded_pairs = batch->num_loaded_pairs;
  batch->read_batch1.SwapSequenceBatch(read_batch1);
  batch->read_batch2.SwapSequenceBatch(read_batch2);
  batch->barcode_batch.SwapSequenceBatch(barcode_batch);
  producer->Release(batch);
  return num_loaded_pairs;
}

CbqPairedEndRangeBatchProducer::CbqPairedEndRangeBatchProducer(
    uint32_t batch_size, uint64_t total_records,
    const SequenceEffectiveRange &read1_effective_range,
    const SequenceEffectiveRange &read2_effective_range,
    const SequenceEffectiveRange &barcode_effective_range, size_t worker_count,
    size_t queue_depth, RangeLoadFn load_fn)
    : batch_size_(batch_size),
      total_records_(total_records),
      worker_count_(worker_count == 0 ? 1 : worker_count),
      load_fn_(std::move(load_fn)) {
  if (batch_size_ == 0) {
    batch_size_ = 1;
  }
  if (total_records_ >
      std::numeric_limits<uint64_t>::max() -
          static_cast<uint64_t>(batch_size_) + 1U) {
    ExitWithMessage("CBQ range producer record count overflows batch count");
  }
  total_batches_ =
      (total_records_ + static_cast<uint64_t>(batch_size_) - 1U) /
      static_cast<uint64_t>(batch_size_);
  if (queue_depth == 0) {
    queue_depth = 1;
  }
  if (queue_depth < worker_count_) {
    queue_depth = worker_count_;
  }
  slots_.reserve(queue_depth);
  for (size_t i = 0; i < queue_depth; ++i) {
    slots_.emplace_back(new CbqPairedEndBatch(batch_size_,
                                              read1_effective_range,
                                              read2_effective_range,
                                              barcode_effective_range));
    free_slots_.push_back(slots_.back().get());
  }
}

CbqPairedEndRangeBatchProducer::~CbqPairedEndRangeBatchProducer() { Stop(); }

void CbqPairedEndRangeBatchProducer::Start() {
  std::lock_guard<std::mutex> lock(mutex_);
  if (started_) {
    return;
  }
  started_ = true;
  if (total_batches_ == 0) {
    done_ = true;
    ready_cv_.notify_all();
    return;
  }
  running_workers_ = worker_count_;
  workers_.reserve(worker_count_);
  for (size_t i = 0; i < worker_count_; ++i) {
    workers_.push_back(
        std::thread(&CbqPairedEndRangeBatchProducer::RunWorker, this));
  }
}

bool CbqPairedEndRangeBatchProducer::Pop(CbqPairedEndBatch **batch,
                                         std::string *error) {
  std::unique_lock<std::mutex> lock(mutex_);
  ready_cv_.wait(lock, [this] {
    return failed_ || ready_slots_.count(next_batch_to_pop_) != 0 ||
           (done_ && next_batch_to_pop_ >= total_batches_);
  });
  if (failed_) {
    if (error != nullptr) {
      *error = error_;
    }
    return false;
  }
  std::map<uint64_t, CbqPairedEndBatch *>::iterator it =
      ready_slots_.find(next_batch_to_pop_);
  if (it == ready_slots_.end()) {
    *batch = nullptr;
    return true;
  }
  *batch = it->second;
  ready_slots_.erase(it);
  ++next_batch_to_pop_;
  return true;
}

void CbqPairedEndRangeBatchProducer::Release(CbqPairedEndBatch *batch) {
  if (batch == nullptr) {
    return;
  }
  {
    std::lock_guard<std::mutex> lock(mutex_);
    free_slots_.push_back(batch);
  }
  free_cv_.notify_one();
}

void CbqPairedEndRangeBatchProducer::Stop() {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    stop_requested_ = true;
  }
  free_cv_.notify_all();
  ready_cv_.notify_all();
  for (size_t i = 0; i < workers_.size(); ++i) {
    if (workers_[i].joinable()) {
      workers_[i].join();
    }
  }
  workers_.clear();
}

CbqPairedEndBatch *CbqPairedEndRangeBatchProducer::AcquireFreeSlot() {
  std::unique_lock<std::mutex> lock(mutex_);
  free_cv_.wait(lock, [this] {
    return stop_requested_ || failed_ || !free_slots_.empty();
  });
  if (stop_requested_ || failed_) {
    return nullptr;
  }
  CbqPairedEndBatch *slot = free_slots_.front();
  free_slots_.pop_front();
  return slot;
}

void CbqPairedEndRangeBatchProducer::MarkWorkerDone() {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    if (running_workers_ > 0) {
      --running_workers_;
    }
    if (running_workers_ == 0) {
      done_ = true;
    }
  }
  ready_cv_.notify_all();
  free_cv_.notify_all();
}

void CbqPairedEndRangeBatchProducer::MarkFailed(
    const std::string &message) {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    failed_ = true;
    done_ = true;
    error_ = message;
  }
  ready_cv_.notify_all();
  free_cv_.notify_all();
}

void CbqPairedEndRangeBatchProducer::RunWorker() {
  try {
    for (;;) {
      uint64_t batch_index = 0;
      {
        std::lock_guard<std::mutex> lock(mutex_);
        if (stop_requested_ || failed_ ||
            next_batch_to_assign_ >= total_batches_) {
          break;
        }
        batch_index = next_batch_to_assign_++;
      }

      CbqPairedEndBatch *slot = AcquireFreeSlot();
      if (slot == nullptr) {
        break;
      }

      const uint64_t first_record =
          batch_index * static_cast<uint64_t>(batch_size_);
      const uint64_t remaining = total_records_ - first_record;
      const uint32_t record_count =
          remaining < static_cast<uint64_t>(batch_size_)
              ? static_cast<uint32_t>(remaining)
              : batch_size_;
      slot->batch_index = batch_index;
      slot->first_record = first_record;
      slot->num_loaded_pairs = load_fn_(first_record, record_count, slot);
      if (slot->num_loaded_pairs != record_count) {
        Release(slot);
        MarkFailed("CBQ range producer loaded an unexpected record count");
        break;
      }

      {
        std::lock_guard<std::mutex> lock(mutex_);
        ready_slots_[batch_index] = slot;
      }
      ready_cv_.notify_all();
    }
  } catch (const std::exception &ex) {
    MarkFailed(ex.what());
  } catch (...) {
    MarkFailed("unknown CBQ range producer failure");
  }
  MarkWorkerDone();
}

uint32_t PopCbqBatchIntoSequenceBatches(
    CbqPairedEndRangeBatchProducer *producer, SequenceBatch &read_batch1,
    SequenceBatch &read_batch2, SequenceBatch &barcode_batch) {
  CbqPairedEndBatch *batch = nullptr;
  std::string error;
  if (!producer->Pop(&batch, &error)) {
    ExitWithMessage("CBQ range producer failed: " + error);
  }
  if (batch == nullptr) {
    return 0;
  }

  const uint32_t num_loaded_pairs = batch->num_loaded_pairs;
  batch->read_batch1.SwapSequenceBatch(read_batch1);
  batch->read_batch2.SwapSequenceBatch(read_batch2);
  batch->barcode_batch.SwapSequenceBatch(barcode_batch);
  producer->Release(batch);
  return num_loaded_pairs;
}

}  // namespace chromap
