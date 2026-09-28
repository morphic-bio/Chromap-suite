#include "fastq_bgzf_input.h"

#include <sys/stat.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <limits>
#include <sstream>
#include <thread>

#include "star_input/BgzfBlockReader.h"
#include "star_input/BgzfRangeReader.h"

namespace chromap {
namespace {

const char *const kStreamRoles[3] = {"read 1", "read 2", "barcode"};

// kseq ends the read name at the first whitespace character and takes the
// rest of the header line, after that one character, as the comment.
void SplitHeader(const char *header, size_t header_length, size_t *name_length,
                 const char **comment, size_t *comment_length) {
  size_t length = 0;
  while (length < header_length &&
         !std::isspace(static_cast<unsigned char>(header[length]))) {
    ++length;
  }
  *name_length = length;
  if (length < header_length) {
    *comment = header + length + 1;
    *comment_length = header_length - length - 1;
  } else {
    *comment = nullptr;
    *comment_length = 0;
  }
}

// The read-name stem STAR Suite compares across synchronized streams.
size_t ReadNameStemLength(const char *name, size_t length) {
  if (length >= 2 && name[length - 2] == '/' && name[length - 1] >= '1' &&
      name[length - 1] <= '3') {
    return length - 2;
  }
  return length;
}

// Parses the first record of a candidate file with the reader that would read
// it, inflating on the calling thread. An empty *why means the file qualifies.
void ProbeBgzfFastq(const std::string &path, std::string *first_name,
                    std::string *why) {
  why->clear();
  first_name->clear();
  struct stat info;
  if (::stat(path.c_str(), &info) != 0) {
    *why = "cannot stat the file";
    return;
  }
  if (!S_ISREG(info.st_mode)) {
    *why = "not a regular file";
    return;
  }
  star_input::BgzfDetection detection;
  std::string error;
  if (!star_input::detect_bgzf(path, &detection, &error)) {
    *why = error;
    return;
  }
  if (!detection.isBgzf) {
    *why = "not BGZF";
    return;
  }
  star_input::BgzfRangeReader reader;
  if (!reader.open(path, 0, std::numeric_limits<uint64_t>::max(),
                   /*worker_threads=*/0, /*check_crc=*/true, &error)) {
    *why = error;
    return;
  }
  star_input::BgzfFastqRecord record;
  if (!reader.next(&record, &error)) {
    *why = error.empty() ? "no FASTQ records" : "first record: " + error;
    return;
  }
  size_t name_length = 0;
  const char *comment = nullptr;
  size_t comment_length = 0;
  SplitHeader(record.name_data(), record.nameLength, &name_length, &comment,
              &comment_length);
  first_name->assign(record.name_data(), name_length);
}

}  // namespace

bool SelectBgzfFastqInput(FastqBgzfMode mode,
                          const std::vector<std::string> &paths,
                          bool *use_bgzf, std::string *message) {
  *use_bgzf = false;
  message->clear();
  if (mode == FastqBgzfMode::kOff) {
    *message = "--input-bgzf-mode off";
    return true;
  }
  const bool required = mode == FastqBgzfMode::kOn;
  std::vector<std::string> first_names;
  for (const std::string &path : paths) {
    std::string first_name;
    std::string why;
    ProbeBgzfFastq(path, &first_name, &why);
    if (!why.empty()) {
      *message = path + ": " + why;
      if (required) {
        *message = "--input-bgzf-mode on requires regular BGZF FASTQ input; " +
                   *message;
        return false;
      }
      return true;
    }
    first_names.push_back(first_name);
  }
  for (size_t i = 1; i < first_names.size(); ++i) {
    const std::string &a = first_names[0];
    const std::string &b = first_names[i];
    const size_t stem_a = ReadNameStemLength(a.data(), a.size());
    const size_t stem_b = ReadNameStemLength(b.data(), b.size());
    if (stem_a != stem_b || a.compare(0, stem_a, b, 0, stem_b) != 0) {
      *message = "first read names do not pair ('" + a + "' in " + paths[0] +
                 ", '" + b + "' in " + paths[i] + ")";
      if (required) {
        *message = "--input-bgzf-mode on: " + *message;
        return false;
      }
      return true;
    }
  }
  *use_bgzf = true;
  return true;
}

std::vector<uint32_t> BgzfInflateWorkersPerStream(int requested,
                                                  int num_threads,
                                                  size_t num_streams) {
  std::vector<uint32_t> workers(num_streams, 0);
  if (num_streams == 0) {
    return workers;
  }
  const int streams = static_cast<int>(num_streams);
  int total = 0;
  if (requested > 0) {
    total = requested;
  } else if (num_threads >= 3) {
    total = std::max(num_threads - streams, streams);
  }
  for (int i = 0; i < streams; ++i) {
    workers[i] = static_cast<uint32_t>(total / streams + (i < total % streams));
  }
  return workers;
}

struct BgzfFastqStream::Impl {
  std::string path;
  star_input::BgzfRangeReader reader;
  star_input::BgzfFastqRecord record;
  // Keeps the inflated window behind the record's zero-copy views alive until
  // the record has been copied into the batch.
  star_input::BgzfBatchLease lease;
  bool open = false;
};

BgzfFastqStream::BgzfFastqStream() : impl_(new Impl()) {}

BgzfFastqStream::~BgzfFastqStream() = default;

const std::string &BgzfFastqStream::path() const { return impl_->path; }

bool BgzfFastqStream::Open(const std::string &path, uint32_t inflate_workers,
                           std::string *error) {
  impl_->path = path;
  impl_->open = false;
  impl_->lease.clear();
  std::string reader_error;
  if (!impl_->reader.open(path, 0, std::numeric_limits<uint64_t>::max(),
                          inflate_workers, /*check_crc=*/true, &reader_error)) {
    *error = path + ": " + reader_error;
    return false;
  }
  impl_->open = true;
  error->clear();
  return true;
}

bool BgzfFastqStream::LoadNext(SequenceBatch &batch, uint32_t index,
                               uint64_t *source_ordinal, std::string *error) {
  error->clear();
  if (!impl_->open) {
    *error = "BGZF FASTQ stream is not open";
    return false;
  }
  star_input::BgzfFastqRecord &record = impl_->record;
  std::string reader_error;
  do {
    // The previous record has been copied; release its window.
    impl_->lease.clear();
    if (!impl_->reader.next(&record, &reader_error, &impl_->lease)) {
      if (!reader_error.empty()) {
        *error = impl_->path + ": " + reader_error;
      }
      return false;
    }
  } while (record.sequenceLength == 0);  // kseq skips empty records

  size_t name_length = 0;
  const char *comment = nullptr;
  size_t comment_length = 0;
  SplitHeader(record.name_data(), record.nameLength, &name_length, &comment,
              &comment_length);
  batch.AssignLoadedSequence(index, record.name_data(), name_length, comment,
                             comment_length, record.sequence_data(),
                             record.sequenceLength, record.quality_data(),
                             record.qualityLength);
  *source_ordinal = record.ordinal;
  return true;
}

bool BgzfFastqStream::LoadBatch(SequenceBatch &batch, uint32_t max_records,
                                uint32_t *loaded, std::string *error) {
  batch.ResetLoadedSequences();
  *loaded = 0;
  uint64_t ordinal = 0;
  while (*loaded < max_records) {
    if (!LoadNext(batch, *loaded, &ordinal, error)) {
      return error->empty();
    }
    ++*loaded;
  }
  return true;
}

BgzfPairedEndReadProvider::BgzfPairedEndReadProvider(
    const std::string &read1_path, const std::string &read2_path,
    const std::string &barcode_path,
    const std::vector<uint32_t> &inflate_workers)
    : inflate_workers_(inflate_workers) {
  paths_.push_back(read1_path);
  paths_.push_back(read2_path);
  if (!barcode_path.empty()) {
    paths_.push_back(barcode_path);
  }
  inflate_workers_.resize(paths_.size(), 0);
}

bool BgzfPairedEndReadProvider::Open(std::string *error) {
  streams_.clear();
  num_pairs_ = 0;
  for (size_t i = 0; i < paths_.size(); ++i) {
    streams_.emplace_back(new BgzfFastqStream());
    if (!streams_.back()->Open(paths_[i], inflate_workers_[i], error)) {
      *error = std::string("BGZF ") + kStreamRoles[i] + " input " + *error;
      streams_.clear();
      return false;
    }
  }
  return true;
}

bool BgzfPairedEndReadProvider::HasBarcode() const {
  return paths_.size() == 3;
}

bool BgzfPairedEndReadProvider::LoadBatch(uint32_t max_pairs,
                                          SequenceBatch &read_batch1,
                                          SequenceBatch &read_batch2,
                                          SequenceBatch &barcode_batch,
                                          uint32_t &num_loaded_pairs,
                                          std::string &error) {
  num_loaded_pairs = 0;
  error.clear();
  if (streams_.size() != paths_.size()) {
    error = "BGZF paired-end read provider is not open";
    return false;
  }
  SequenceBatch *const batches[3] = {&read_batch1, &read_batch2,
                                     &barcode_batch};
  const size_t num_streams = streams_.size();
  ordinals_.resize(num_streams);
  uint32_t loaded[3] = {0, 0, 0};
  std::string stream_errors[3];
  // Each file is consumed in order on its own thread, as STAR Suite feeds
  // each mate through its own producer; the batch's records are paired below.
  auto fill = [&](size_t i) {
    SequenceBatch &batch = *batches[i];
    std::vector<uint64_t> &ordinals = ordinals_[i];
    if (ordinals.size() < max_pairs) {
      ordinals.resize(max_pairs);
    }
    batch.ResetLoadedSequences();
    uint32_t count = 0;
    while (count < max_pairs &&
           streams_[i]->LoadNext(batch, count, &ordinals[count],
                                 &stream_errors[i])) {
      ++count;
    }
    loaded[i] = count;
  };
  std::vector<std::thread> threads;
  for (size_t i = 1; i < num_streams; ++i) {
    threads.emplace_back(fill, i);
  }
  fill(0);
  for (std::thread &thread : threads) {
    thread.join();
  }
  for (size_t i = 0; i < num_streams; ++i) {
    if (!stream_errors[i].empty()) {
      error = std::string("BGZF ") + kStreamRoles[i] + ": " + stream_errors[i];
      return false;
    }
  }

  // Pair the records as STAR Suite pairs synchronized BGZF streams: record k
  // of every file has the same ordinal and the same read-name stem, and no
  // file ends before the others.
  uint32_t paired = loaded[0];
  for (size_t i = 1; i < num_streams; ++i) {
    paired = std::min(paired, loaded[i]);
  }
  for (uint32_t k = 0; k < paired; ++k) {
    const char *name0 = batches[0]->GetSequenceNameAt(k);
    const size_t stem0 =
        ReadNameStemLength(name0, batches[0]->GetSequenceNameLengthAt(k));
    for (size_t i = 1; i < num_streams; ++i) {
      const char *name = batches[i]->GetSequenceNameAt(k);
      const size_t length = batches[i]->GetSequenceNameLengthAt(k);
      const size_t stem = ReadNameStemLength(name, length);
      if (ordinals_[i][k] != ordinals_[0][k] || stem != stem0 ||
          std::memcmp(name0, name, stem0) != 0) {
        std::ostringstream message;
        message << "BGZF FASTQ records do not pair at pair " << num_pairs_ + k
                << ": " << kStreamRoles[0] << " '" << std::string(name0, stem0)
                << "' (record " << ordinals_[0][k] << " of " << paths_[0]
                << "), " << kStreamRoles[i] << " '" << std::string(name, length)
                << "' (record " << ordinals_[i][k] << " of " << paths_[i]
                << "). If these files pair by position despite their names, "
                   "--input-bgzf-mode off reads them without this check";
        error = message.str();
        return false;
      }
    }
  }
  for (size_t i = 1; i < num_streams; ++i) {
    if (loaded[i] != loaded[0]) {
      std::ostringstream message;
      message << "BGZF FASTQ record counts differ: ";
      for (size_t j = 0; j < num_streams; ++j) {
        message << (j == 0 ? "" : ", ") << kStreamRoles[j] << " file "
                << paths_[j];
      }
      message << " do not all have record " << num_pairs_ + paired;
      error = message.str();
      return false;
    }
  }
  num_loaded_pairs = loaded[0];
  num_pairs_ += num_loaded_pairs;
  return true;
}

}  // namespace chromap
