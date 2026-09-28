// FASTQ intake harness: reads one paired-end lane (read 1, read 2 and an
// optional barcode file) through a chosen Chromap reader, optionally writes
// every record Chromap would receive, and reports the loading rate.
//
//   fastq_intake_harness --reader READER [--threads N] [--bgzf-threads N]
//                        [--batch N] [--dump FILE|-] R1 R2 [BC]
//
// READER is one of
//   kseq-serial   zlib (kseq), one record per file in turn on one thread
//                 (Chromap's loader below three threads);
//   kseq-threads  zlib (kseq), each file on its own thread for every batch
//                 (Chromap's loader from three threads up);
//   bgzf          the BGZF paired-end read provider;
//   auto          what --input-bgzf-mode auto selects for the lane.
//
// --dump writes, per pair and per file, "name<TAB>comment<TAB>seq<TAB>qual".

#include <sys/stat.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <string>
#include <thread>
#include <vector>

#include "fastq_bgzf_input.h"
#include "sequence_batch.h"
#include "sequence_effective_range.h"
#include "utils.h"

using namespace chromap;

namespace {

int Usage() {
  std::cerr << "usage: fastq_intake_harness --reader "
               "kseq-serial|kseq-threads|bgzf|auto [--threads N] "
               "[--bgzf-threads N] [--batch N] [--dump FILE|-] R1 R2 [BC]\n";
  return 2;
}

void DumpBatches(FILE *out, std::vector<SequenceBatch *> &batches,
                 uint32_t num_pairs) {
  for (uint32_t i = 0; i < num_pairs; ++i) {
    for (SequenceBatch *batch : batches) {
      fwrite(batch->GetSequenceNameAt(i), 1,
             batch->GetSequenceNameLengthAt(i), out);
      fputc('\t', out);
      if (batch->GetSequenceCommentLengthAt(i) != 0) {
        fwrite(batch->GetSequenceCommentAt(i), 1,
               batch->GetSequenceCommentLengthAt(i), out);
      }
      fputc('\t', out);
      fwrite(batch->GetSequenceAt(i), 1, batch->GetSequenceLengthAt(i), out);
      fputc('\t', out);
      if (batch->GetSequenceQualLengthAt(i) != 0) {
        fwrite(batch->GetSequenceQualAt(i), 1,
               batch->GetSequenceQualLengthAt(i), out);
      }
      fputc('\n', out);
    }
  }
}

}  // namespace

int main(int argc, char **argv) {
  std::string reader;
  std::string dump_path;
  int num_threads = 3;
  int bgzf_threads = 0;
  uint32_t batch_size = 500000;
  std::vector<std::string> paths;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--reader" && i + 1 < argc) {
      reader = argv[++i];
    } else if (arg == "--threads" && i + 1 < argc) {
      num_threads = std::atoi(argv[++i]);
    } else if (arg == "--bgzf-threads" && i + 1 < argc) {
      bgzf_threads = std::atoi(argv[++i]);
    } else if (arg == "--batch" && i + 1 < argc) {
      batch_size = static_cast<uint32_t>(std::strtoul(argv[++i], nullptr, 10));
    } else if (arg == "--dump" && i + 1 < argc) {
      dump_path = argv[++i];
    } else if (!arg.empty() && arg[0] == '-') {
      return Usage();
    } else {
      paths.push_back(arg);
    }
  }
  if (paths.size() < 2 || paths.size() > 3 || batch_size == 0) {
    return Usage();
  }
  if (reader == "auto") {
    bool use_bgzf = false;
    std::string message;
    if (!SelectBgzfFastqInput(FastqBgzfMode::kAuto, paths, &use_bgzf,
                              &message)) {
      std::cerr << message << "\n";
      return 1;
    }
    reader = use_bgzf ? "bgzf" : "kseq-threads";
    std::cerr << "auto selected " << reader
              << (message.empty() ? "" : " (" + message + ")") << "\n";
  }
  if (reader != "kseq-serial" && reader != "kseq-threads" &&
      reader != "bgzf") {
    return Usage();
  }

  FILE *out = nullptr;
  if (!dump_path.empty()) {
    out = dump_path == "-" ? stdout : fopen(dump_path.c_str(), "wb");
    if (out == nullptr) {
      std::cerr << "cannot open " << dump_path << "\n";
      return 1;
    }
  }

  const SequenceEffectiveRange full_range;
  std::vector<std::unique_ptr<SequenceBatch>> owned;
  std::vector<SequenceBatch *> batches;
  for (size_t i = 0; i < 3; ++i) {
    owned.emplace_back(new SequenceBatch(batch_size, full_range));
  }
  for (size_t i = 0; i < paths.size(); ++i) {
    batches.push_back(owned[i].get());
  }

  std::unique_ptr<BgzfPairedEndReadProvider> provider;
  if (reader == "bgzf") {
    const std::vector<uint32_t> workers =
        BgzfInflateWorkersPerStream(bgzf_threads, num_threads, paths.size());
    provider.reset(new BgzfPairedEndReadProvider(
        paths[0], paths[1], paths.size() == 3 ? paths[2] : std::string(),
        workers));
    std::string error;
    if (!provider->Open(&error)) {
      std::cerr << error << "\n";
      return 1;
    }
    std::cerr << "bgzf inflate workers";
    for (uint32_t w : workers) std::cerr << " " << w;
    std::cerr << "\n";
  } else {
    for (size_t i = 0; i < paths.size(); ++i) {
      batches[i]->InitializeLoading(paths[i]);
    }
  }

  uint64_t compressed_bytes = 0;
  for (const std::string &path : paths) {
    struct stat info;
    if (::stat(path.c_str(), &info) == 0 && S_ISREG(info.st_mode)) {
      compressed_bytes += static_cast<uint64_t>(info.st_size);
    }
  }

  const double start = GetRealTime();
  uint64_t total_pairs = 0;
  while (true) {
    uint32_t num_pairs = 0;
    if (reader == "bgzf") {
      std::string error;
      if (!provider->LoadBatch(batch_size, *owned[0], *owned[1], *owned[2],
                               num_pairs, error)) {
        std::cerr << error << "\n";
        return 1;
      }
    } else if (reader == "kseq-serial") {
      while (num_pairs < batch_size) {
        size_t ended = 0;
        for (SequenceBatch *batch : batches) {
          ended += batch->LoadOneSequenceAndSaveAt(num_pairs) ? 1 : 0;
        }
        if (ended == batches.size()) break;
        if (ended != 0) {
          std::cerr << "Numbers of reads and barcodes don't match!\n";
          return 1;
        }
        ++num_pairs;
      }
    } else {
      std::vector<uint32_t> counts(batches.size(), 0);
      std::vector<std::thread> threads;
      for (size_t s = 0; s < batches.size(); ++s) {
        threads.emplace_back([&, s]() {
          uint32_t i = 0;
          for (; i < batch_size; ++i) {
            if (batches[s]->LoadOneSequenceAndSaveAt(i)) break;
          }
          counts[s] = i;
        });
      }
      for (std::thread &thread : threads) thread.join();
      for (uint32_t count : counts) {
        if (count != counts[0]) {
          std::cerr << "Numbers of reads and barcodes don't match!\n";
          return 1;
        }
      }
      num_pairs = counts[0];
    }
    if (num_pairs == 0) break;
    if (out != nullptr) DumpBatches(out, batches, num_pairs);
    total_pairs += num_pairs;
  }
  const double seconds = GetRealTime() - start;
  if (reader != "bgzf") {
    for (SequenceBatch *batch : batches) batch->FinalizeLoading();
  }
  if (out != nullptr && out != stdout) fclose(out);
  std::fprintf(stderr,
               "reader=%s pairs=%llu seconds=%.3f pairs_per_second=%.0f "
               "compressed_MB_per_second=%.1f\n",
               reader.c_str(), static_cast<unsigned long long>(total_pairs),
               seconds, seconds > 0 ? total_pairs / seconds : 0.0,
               seconds > 0 ? compressed_bytes / 1e6 / seconds : 0.0);
  return 0;
}
