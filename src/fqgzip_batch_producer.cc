#include "fqgzip_batch_producer.h"

#include <fqgzip/fqgzip.h>

#include <cstring>
#include <exception>
#include <limits>
#include <stdexcept>
#include <utility>

#include "utils.h"

namespace chromap {
namespace {

struct HeaderFields {
  const char *name = nullptr;
  size_t name_size = 0;
  const char *comment = nullptr;
  size_t comment_size = 0;
};

HeaderFields ParseHeader(const fqgz_bytes_view_t &header) {
  if (header.size < 2 || header.data == nullptr || header.data[0] != '@') {
    throw std::runtime_error("fqgzip returned a FASTQ header without @");
  }
  size_t name_begin = 1;
  size_t name_end = name_begin;
  while (name_end < header.size && header.data[name_end] != ' ' &&
         header.data[name_end] != '\t') {
    ++name_end;
  }
  if (name_end == name_begin) {
    throw std::runtime_error("fqgzip returned an empty FASTQ read name");
  }
  size_t comment_begin = name_end;
  while (comment_begin < header.size &&
         (header.data[comment_begin] == ' ' ||
          header.data[comment_begin] == '\t')) {
    ++comment_begin;
  }
  HeaderFields result;
  result.name = header.data + name_begin;
  result.name_size = name_end - name_begin;
  if (comment_begin < header.size) {
    result.comment = header.data + comment_begin;
    result.comment_size = header.size - comment_begin;
  }
  return result;
}

void AssignRecord(const fqgz_fastq_record_view_t &record,
                  uint32_t sequence_index, uint32_t sequence_id,
                  SequenceBatch &batch) {
  if (record.sequence.size > std::numeric_limits<uint32_t>::max() ||
      record.quality.size > std::numeric_limits<uint32_t>::max()) {
    throw std::runtime_error("fqgzip FASTQ record exceeds Chromap limits");
  }
  const HeaderFields header = ParseHeader(record.header);
  char *sequence =
      batch.PrepareLoadedSequenceBuffer(sequence_index, record.sequence.size);
  if (record.sequence.size != 0) {
    std::memcpy(sequence, record.sequence.data, record.sequence.size);
  }
  batch.CommitLoadedSequenceBufferWithId(
      sequence_index, sequence_id, header.name, header.name_size,
      header.comment, header.comment_size, record.sequence.size,
      record.quality.data, record.quality.size);
}

}  // namespace

struct FqgzipPairedEndBatchProducer::ReaderState {
  fqgz_progressive_plan_t *plan = nullptr;
  fqgz_shard_t *shard = nullptr;
  uint64_t next_pair_ordinal = 0;

  ~ReaderState() {
    fqgz_shard_close(shard);
    fqgz_progressive_close(plan);
  }
};

FqgzipPairedEndBatch::FqgzipPairedEndBatch(
    uint32_t batch_size,
    const SequenceEffectiveRange &read1_effective_range,
    const SequenceEffectiveRange &read2_effective_range)
    : read_batch1(batch_size, read1_effective_range),
      read_batch2(batch_size, read2_effective_range) {}

FqgzipPairedEndBatchProducer::FqgzipPairedEndBatchProducer(
    const std::string &read1_path, const std::string &read2_path,
    uint32_t batch_size,
    const SequenceEffectiveRange &read1_effective_range,
    const SequenceEffectiveRange &read2_effective_range,
    uint64_t global_record_offset, uint32_t requested_shards,
    uint32_t worker_threads, size_t queue_depth)
    : read1_path_(read1_path),
      read2_path_(read2_path),
      batch_size_(batch_size == 0 ? 1 : batch_size),
      global_record_offset_(global_record_offset),
      requested_shards_(requested_shards == 0 ? 1 : requested_shards),
      worker_threads_(worker_threads < 2 ? 2 : worker_threads),
      reader_(new ReaderState()) {
  if (queue_depth == 0) {
    queue_depth = 1;
  }
  slots_.reserve(queue_depth);
  for (size_t i = 0; i < queue_depth; ++i) {
    slots_.emplace_back(new FqgzipPairedEndBatch(
        batch_size_, read1_effective_range, read2_effective_range));
    free_slots_.push_back(slots_.back().get());
  }
}

FqgzipPairedEndBatchProducer::~FqgzipPairedEndBatchProducer() { Stop(); }

void FqgzipPairedEndBatchProducer::Start() {
  worker_ = std::thread(&FqgzipPairedEndBatchProducer::Run, this);
}

bool FqgzipPairedEndBatchProducer::Pop(FqgzipPairedEndBatch **batch,
                                       std::string *error) {
  std::unique_lock<std::mutex> lock(mutex_);
  ready_cv_.wait(lock,
                 [this] { return failed_ || done_ || !ready_slots_.empty(); });
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

void FqgzipPairedEndBatchProducer::Release(FqgzipPairedEndBatch *batch) {
  if (batch == nullptr) {
    return;
  }
  {
    std::lock_guard<std::mutex> lock(mutex_);
    free_slots_.push_back(batch);
  }
  free_cv_.notify_one();
}

void FqgzipPairedEndBatchProducer::Stop() {
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

FqgzipPairedEndBatch *FqgzipPairedEndBatchProducer::AcquireFreeSlot() {
  std::unique_lock<std::mutex> lock(mutex_);
  free_cv_.wait(lock,
                [this] { return stop_requested_ || !free_slots_.empty(); });
  if (stop_requested_) {
    return nullptr;
  }
  FqgzipPairedEndBatch *slot = free_slots_.front();
  free_slots_.pop_front();
  return slot;
}

uint32_t FqgzipPairedEndBatchProducer::LoadBatch(
    FqgzipPairedEndBatch *batch) {
  batch->read_batch1.ResetLoadedSequences();
  batch->read_batch2.ResetLoadedSequences();
  uint32_t count = 0;
  while (count < batch_size_) {
    if (reader_->shard == nullptr) {
      fqgz_shard_info_t info{};
      const fqgz_status_t status =
          fqgz_progressive_next_shard(reader_->plan, &reader_->shard, &info);
      if (status == FQGZ_EOF) {
        break;
      }
      if (status != FQGZ_OK) {
        throw std::runtime_error(
            std::string("fqgzip shard discovery failed: ") +
            fqgz_progressive_error(reader_->plan));
      }
    }

    fqgz_pair_record_view_t pair{};
    const fqgz_status_t status =
        fqgz_shard_next_pair(reader_->shard, &pair);
    if (status == FQGZ_EOF) {
      fqgz_shard_close(reader_->shard);
      reader_->shard = nullptr;
      continue;
    }
    if (status != FQGZ_OK) {
      throw std::runtime_error(std::string("fqgzip paired shard failed: ") +
                               fqgz_shard_error(reader_->shard));
    }
    if (pair.pair_ordinal != reader_->next_pair_ordinal) {
      throw std::runtime_error("fqgzip delivered a non-contiguous pair ordinal");
    }
    const uint64_t global_id =
        global_record_offset_ + reader_->next_pair_ordinal;
    if (global_id > std::numeric_limits<uint32_t>::max()) {
      throw std::runtime_error("fqgzip record ordinal exceeds Chromap's uint32 read-id limit");
    }
    AssignRecord(pair.r1, count, static_cast<uint32_t>(global_id),
                 batch->read_batch1);
    AssignRecord(pair.r2, count, static_cast<uint32_t>(global_id),
                 batch->read_batch2);
    ++reader_->next_pair_ordinal;
    ++count;
  }
  batch->num_loaded_pairs = count;
  return count;
}

void FqgzipPairedEndBatchProducer::MarkDone() {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    done_ = true;
  }
  ready_cv_.notify_all();
}

void FqgzipPairedEndBatchProducer::MarkFailed(const std::string &message) {
  {
    std::lock_guard<std::mutex> lock(mutex_);
    failed_ = true;
    done_ = true;
    error_ = message;
  }
  ready_cv_.notify_all();
  free_cv_.notify_all();
}

void FqgzipPairedEndBatchProducer::Run() {
  try {
    fqgz_pair_options_t options;
    fqgz_pair_options_init(&options);
    options.requested_shards = requested_shards_;
    options.worker_threads = worker_threads_;
    const fqgz_status_t open_status = fqgz_progressive_open(
        read1_path_.c_str(), read2_path_.c_str(), &options, &reader_->plan);
    if (open_status != FQGZ_OK) {
      const std::string detail =
          reader_->plan == nullptr ? fqgz_status_string(open_status)
                                   : fqgz_progressive_error(reader_->plan);
      throw std::runtime_error("cannot open fqgzip input: " + detail);
    }

    for (;;) {
      FqgzipPairedEndBatch *slot = AcquireFreeSlot();
      if (slot == nullptr) {
        return;
      }
      if (LoadBatch(slot) == 0) {
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
    MarkFailed("unknown fqgzip producer failure");
  }
}

uint32_t PopFqgzipBatchIntoSequenceBatches(
    FqgzipPairedEndBatchProducer *producer, SequenceBatch &read_batch1,
    SequenceBatch &read_batch2) {
  FqgzipPairedEndBatch *batch = nullptr;
  std::string error;
  if (!producer->Pop(&batch, &error)) {
    ExitWithMessage("fqgzip producer failed: " + error);
  }
  if (batch == nullptr) {
    return 0;
  }
  const uint32_t count = batch->num_loaded_pairs;
  batch->read_batch1.SwapSequenceBatch(read_batch1);
  batch->read_batch2.SwapSequenceBatch(read_batch2);
  producer->Release(batch);
  return count;
}

}  // namespace chromap
