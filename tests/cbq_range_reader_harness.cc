#include "cbq_reader.h"

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <thread>
#include <vector>

namespace {

const uint64_t kFnvOffset = 1469598103934665603ULL;
const uint64_t kFnvPrime = 1099511628211ULL;

struct Args {
  std::string cbq_path;
  uint32_t mate_count = 2;
  uint32_t batch_size = 131072;
  int threads = 4;
  bool materialize_sequence = false;
};

struct Stats {
  uint64_t records = 0;
  uint64_t batches = 0;
  uint64_t segments = 0;
  uint64_t bases = 0;
  uint64_t read_ordinal_sum = 0;
  uint64_t lane_ordinal_sum = 0;
  uint64_t record_hash_sum = 0;
  uint64_t record_hash_xor = 0;
};

void Usage(const char *prog) {
  std::cerr
      << "Usage: " << prog
      << " --cbq PATH [--mate-count 1|2] [--threads N]\n"
      << "       [--batch-size N] [--materialize none|sequence]\n";
}

bool ParseUint32(const char *text, uint32_t *out) {
  char *end = nullptr;
  const unsigned long value = std::strtoul(text, &end, 10);
  if (end == text || *end != '\0' ||
      value > static_cast<unsigned long>(std::numeric_limits<uint32_t>::max())) {
    return false;
  }
  *out = static_cast<uint32_t>(value);
  return true;
}

bool ParseInt(const char *text, int *out) {
  char *end = nullptr;
  const long value = std::strtol(text, &end, 10);
  if (end == text || *end != '\0' ||
      value < static_cast<long>(std::numeric_limits<int>::min()) ||
      value > static_cast<long>(std::numeric_limits<int>::max())) {
    return false;
  }
  *out = static_cast<int>(value);
  return true;
}

bool ParseArgs(int argc, char **argv, Args *args) {
  for (int i = 1; i < argc; ++i) {
    const std::string key = argv[i];
    auto need_value = [&](const char *name) -> const char * {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << name << "\n";
        return nullptr;
      }
      return argv[++i];
    };

    if (key == "--cbq") {
      const char *value = need_value("--cbq");
      if (value == nullptr) return false;
      args->cbq_path = value;
    } else if (key == "--mate-count") {
      const char *value = need_value("--mate-count");
      if (value == nullptr || !ParseUint32(value, &args->mate_count)) {
        std::cerr << "Invalid --mate-count\n";
        return false;
      }
    } else if (key == "--threads") {
      const char *value = need_value("--threads");
      if (value == nullptr || !ParseInt(value, &args->threads)) {
        std::cerr << "Invalid --threads\n";
        return false;
      }
    } else if (key == "--batch-size") {
      const char *value = need_value("--batch-size");
      if (value == nullptr || !ParseUint32(value, &args->batch_size)) {
        std::cerr << "Invalid --batch-size\n";
        return false;
      }
    } else if (key == "--materialize") {
      const char *value = need_value("--materialize");
      if (value == nullptr) return false;
      const std::string mode = value;
      if (mode == "none") {
        args->materialize_sequence = false;
      } else if (mode == "sequence") {
        args->materialize_sequence = true;
      } else {
        std::cerr << "--materialize must be none or sequence\n";
        return false;
      }
    } else if (key == "--help" || key == "-h") {
      return false;
    } else {
      std::cerr << "Unknown argument: " << key << "\n";
      return false;
    }
  }

  if (args->cbq_path.empty()) {
    std::cerr << "--cbq is required\n";
    return false;
  }
  if (args->mate_count != 1 && args->mate_count != 2) {
    std::cerr << "--mate-count must be 1 or 2\n";
    return false;
  }
  if (args->batch_size == 0) {
    std::cerr << "--batch-size must be positive\n";
    return false;
  }
  if (args->threads <= 0) {
    std::cerr << "--threads must be positive\n";
    return false;
  }
  return true;
}

void HashByte(uint64_t *hash, unsigned char value) {
  *hash ^= static_cast<uint64_t>(value);
  *hash *= kFnvPrime;
}

void HashU64(uint64_t *hash, uint64_t value) {
  for (unsigned i = 0; i < 8; ++i) {
    HashByte(hash, static_cast<unsigned char>((value >> (i * 8)) & 0xffU));
  }
}

void HashBytes(uint64_t *hash, const char *data, size_t size) {
  if (data == nullptr) {
    HashU64(hash, 0);
    HashU64(hash, size);
    return;
  }
  for (size_t i = 0; i < size; ++i) {
    HashByte(hash, static_cast<unsigned char>(data[i]));
  }
}

void MergeStats(const Stats &src, Stats *dest) {
  dest->records += src.records;
  dest->batches += src.batches;
  dest->segments += src.segments;
  dest->bases += src.bases;
  dest->read_ordinal_sum += src.read_ordinal_sum;
  dest->lane_ordinal_sum += src.lane_ordinal_sum;
  dest->record_hash_sum += src.record_hash_sum;
  dest->record_hash_xor ^= src.record_hash_xor;
}

bool AccumulateRecord(const chromap::CbqReadView &record,
                      bool materialize_sequence, Stats *stats,
                      std::vector<char> *scratch, std::string *error) {
  uint64_t record_hash = kFnvOffset;
  ++stats->records;
  stats->read_ordinal_sum += record.read_ordinal;
  stats->lane_ordinal_sum += record.lane_read_ordinal;

  HashU64(&record_hash, record.read_ordinal);
  HashU64(&record_hash, record.lane_read_ordinal);
  HashU64(&record_hash, record.segment_count);
  HashByte(&record_hash, static_cast<unsigned char>(record.read_filter));
  HashU64(&record_hash, record.read_name.size);
  HashBytes(&record_hash, record.read_name.data, record.read_name.size);
  HashU64(&record_hash, record.read_name_extra.size);
  HashBytes(&record_hash, record.read_name_extra.data,
            record.read_name_extra.size);

  if (record.segment_count != 0 && record.segments == nullptr) {
    *error = "CBQ record has segments but segment pointer is null";
    return false;
  }

  for (uint32_t isegment = 0; isegment < record.segment_count; ++isegment) {
    const chromap::CbqSegmentView &segment = record.segments[isegment];
    const size_t length = chromap::CbqSegmentSequenceLength(segment);
    ++stats->segments;
    stats->bases += static_cast<uint64_t>(length);
    HashU64(&record_hash, isegment);
    HashU64(&record_hash, length);
    HashU64(&record_hash, segment.quality.size);
    HashBytes(&record_hash, segment.quality.data, segment.quality.size);

    if (materialize_sequence) {
      if (scratch->size() < length + 1) {
        scratch->resize(length + 1);
      }
      size_t written = 0;
      if (!chromap::MaterializeCbqSegmentSequenceToBuffer(
              segment, scratch->data(), scratch->size(), &written, error)) {
        return false;
      }
      if (written != length) {
        *error = "CBQ materialized sequence length mismatch";
        return false;
      }
      HashBytes(&record_hash, scratch->data(), written);
    }
  }

  stats->record_hash_sum += record_hash;
  stats->record_hash_xor ^= record_hash;
  return true;
}

bool RunSequentialRange(const Args &args, uint64_t first_record,
                        uint64_t record_count, Stats *stats,
                        std::string *error) {
  chromap::CbqLaneReader reader(args.cbq_path, args.mate_count);
  if (!reader.Open(error)) {
    return false;
  }

  std::vector<char> scratch;
  const uint64_t range_end =
      record_count == std::numeric_limits<uint64_t>::max()
          ? std::numeric_limits<uint64_t>::max()
          : first_record + record_count;
  uint64_t seen = 0;
  for (;;) {
    chromap::CbqReadBatchView batch;
    const chromap::CbqReadStatus status =
        reader.NextBatch(args.batch_size, &batch, error);
    if (status == chromap::CbqReadStatus::kError) {
      reader.Close();
      return false;
    }
    if (status == chromap::CbqReadStatus::kEnd) {
      break;
    }
    ++stats->batches;
    for (uint32_t i = 0; i < batch.record_count; ++i, ++seen) {
      if (seen < first_record) {
        continue;
      }
      if (seen >= range_end) {
        reader.Close();
        return true;
      }
      if (!AccumulateRecord(batch.records[i], args.materialize_sequence, stats,
                            &scratch, error)) {
        reader.Close();
        return false;
      }
    }
  }
  reader.Close();
  return true;
}

bool RunIndexedRange(const Args &args, uint64_t first_record,
                     uint64_t record_count, Stats *stats,
                     std::string *error) {
  chromap::CbqLaneReader reader(args.cbq_path, args.mate_count);
  if (!reader.OpenRange(first_record, record_count, error)) {
    return false;
  }

  std::vector<char> scratch;
  for (;;) {
    chromap::CbqReadBatchView batch;
    const chromap::CbqReadStatus status =
        reader.NextBatch(args.batch_size, &batch, error);
    if (status == chromap::CbqReadStatus::kError) {
      reader.Close();
      return false;
    }
    if (status == chromap::CbqReadStatus::kEnd) {
      break;
    }
    ++stats->batches;
    for (uint32_t i = 0; i < batch.record_count; ++i) {
      if (!AccumulateRecord(batch.records[i], args.materialize_sequence, stats,
                            &scratch, error)) {
        reader.Close();
        return false;
      }
    }
  }
  reader.Close();
  return true;
}

bool RunCachedIndexedRange(const Args &args, const chromap::CbqLaneIndex &index,
                           uint64_t first_record, uint64_t record_count,
                           Stats *stats, std::string *error) {
  chromap::CbqLaneReader reader(args.cbq_path, args.mate_count);
  if (!reader.OpenRangeWithIndex(index, first_record, record_count, error)) {
    return false;
  }

  std::vector<char> scratch;
  for (;;) {
    chromap::CbqReadBatchView batch;
    const chromap::CbqReadStatus status =
        reader.NextBatch(args.batch_size, &batch, error);
    if (status == chromap::CbqReadStatus::kError) {
      reader.Close();
      return false;
    }
    if (status == chromap::CbqReadStatus::kEnd) {
      break;
    }
    ++stats->batches;
    for (uint32_t i = 0; i < batch.record_count; ++i) {
      if (!AccumulateRecord(batch.records[i], args.materialize_sequence, stats,
                            &scratch, error)) {
        reader.Close();
        return false;
      }
    }
  }
  reader.Close();
  return true;
}

bool SameStats(const Stats &a, const Stats &b) {
  return a.records == b.records && a.segments == b.segments &&
         a.bases == b.bases && a.read_ordinal_sum == b.read_ordinal_sum &&
         a.lane_ordinal_sum == b.lane_ordinal_sum &&
         a.record_hash_sum == b.record_hash_sum &&
         a.record_hash_xor == b.record_hash_xor;
}

void PrintStats(const char *label, const Stats &stats) {
  std::cout << label << "\t"
            << "records=" << stats.records << "\t"
            << "batches=" << stats.batches << "\t"
            << "segments=" << stats.segments << "\t"
            << "bases=" << stats.bases << "\t"
            << "read_ordinal_sum=" << stats.read_ordinal_sum << "\t"
            << "lane_ordinal_sum=" << stats.lane_ordinal_sum << "\t"
            << "hash_sum=" << stats.record_hash_sum << "\t"
            << "hash_xor=" << stats.record_hash_xor << "\n";
}

bool CheckRange(const Args &args, const chromap::CbqLaneIndex &index,
                uint64_t first_record, uint64_t record_count,
                const std::string &label) {
  Stats sequential;
  Stats uncached_indexed;
  Stats indexed;
  std::string error;
  if (!RunSequentialRange(args, first_record, record_count, &sequential,
                          &error)) {
    std::cerr << "Sequential range failed for " << label << ": " << error
              << "\n";
    return false;
  }
  if (!RunIndexedRange(args, first_record, record_count, &uncached_indexed,
                       &error)) {
    std::cerr << "Indexed range failed for " << label << ": " << error << "\n";
    return false;
  }
  if (!RunCachedIndexedRange(args, index, first_record, record_count,
                             &indexed, &error)) {
    std::cerr << "Cached indexed range failed for " << label << ": " << error
              << "\n";
    return false;
  }
  if (!SameStats(sequential, uncached_indexed) ||
      !SameStats(sequential, indexed)) {
    std::cerr << "Range stats mismatch for " << label << "\n";
    PrintStats("sequential", sequential);
    PrintStats("indexed", uncached_indexed);
    PrintStats("cached_indexed", indexed);
    return false;
  }
  PrintStats(label.c_str(), indexed);
  return true;
}

}  // namespace

int main(int argc, char **argv) {
  Args args;
  if (!ParseArgs(argc, argv, &args)) {
    Usage(argv[0]);
    return 2;
  }

  chromap::CbqLaneReader metadata_reader(args.cbq_path, args.mate_count);
  chromap::CbqLaneIndex lane_index;
  std::string error;
  if (!metadata_reader.LoadIndex(&lane_index, &error)) {
    std::cerr << "Failed to load CBQ range metadata: " << error << "\n";
    return 1;
  }
  const uint64_t total_records = lane_index.total_records;
  std::cout << "metadata\ttotal_records=" << total_records << "\tthreads="
            << args.threads << "\tbatch_size=" << args.batch_size
            << "\tmaterialize="
            << (args.materialize_sequence ? "sequence" : "none") << "\n";

  Stats sequential_all;
  if (!RunSequentialRange(args, 0, std::numeric_limits<uint64_t>::max(),
                          &sequential_all, &error)) {
    std::cerr << "Sequential full decode failed: " << error << "\n";
    return 1;
  }

  int nthreads = args.threads;
  if (total_records == 0) {
    nthreads = 1;
  } else if (total_records < static_cast<uint64_t>(nthreads)) {
    nthreads = static_cast<int>(total_records);
  }
  const uint64_t chunk_size =
      total_records == 0
          ? 0
          : (total_records + static_cast<uint64_t>(nthreads) - 1U) /
                static_cast<uint64_t>(nthreads);

  std::vector<Stats> thread_stats(static_cast<size_t>(nthreads));
  std::vector<std::string> thread_errors(static_cast<size_t>(nthreads));
  std::vector<int> thread_ok(static_cast<size_t>(nthreads), 1);
  std::vector<std::thread> workers;
  workers.reserve(static_cast<size_t>(nthreads));
  for (int ithread = 0; ithread < nthreads; ++ithread) {
    const uint64_t first = static_cast<uint64_t>(ithread) * chunk_size;
    const uint64_t count =
        first >= total_records ? 0 : std::min(chunk_size, total_records - first);
    workers.push_back(std::thread([&, ithread, first, count]() {
      if (!RunCachedIndexedRange(
              args, lane_index, first, count,
              &thread_stats[static_cast<size_t>(ithread)],
              &thread_errors[static_cast<size_t>(ithread)])) {
        thread_ok[static_cast<size_t>(ithread)] = 0;
      }
    }));
  }
  for (size_t i = 0; i < workers.size(); ++i) {
    workers[i].join();
  }

  Stats indexed_parallel;
  for (int ithread = 0; ithread < nthreads; ++ithread) {
    if (!thread_ok[static_cast<size_t>(ithread)]) {
      std::cerr << "Indexed worker " << ithread
                << " failed: " << thread_errors[static_cast<size_t>(ithread)]
                << "\n";
      return 1;
    }
    MergeStats(thread_stats[static_cast<size_t>(ithread)], &indexed_parallel);
  }

  PrintStats("sequential_full", sequential_all);
  PrintStats("indexed_parallel", indexed_parallel);
  if (!SameStats(sequential_all, indexed_parallel)) {
    std::cerr << "Parallel indexed range totals do not match sequential decode\n";
    return 1;
  }

  if (!CheckRange(args, lane_index, 0, 0, "range_empty_start")) return 1;
  if (!CheckRange(args, lane_index, total_records, 0, "range_empty_end")) {
    return 1;
  }
  if (!CheckRange(args, lane_index, 0, std::min<uint64_t>(13, total_records),
                  "range_prefix")) {
    return 1;
  }
  if (total_records > 1 &&
      !CheckRange(args, lane_index, 1,
                  std::min<uint64_t>(17, total_records - 1),
                  "range_offset_one")) {
    return 1;
  }
  if (total_records > 0 &&
      !CheckRange(args, lane_index, total_records / 2,
                  std::min<uint64_t>(31, total_records - total_records / 2),
                  "range_middle")) {
    return 1;
  }
  if (total_records > 17 &&
      !CheckRange(args, lane_index, total_records - 17, 17, "range_suffix")) {
    return 1;
  }

  std::cout << "PASS\tCBQ range reader totals match sequential decode\n";
  return 0;
}
