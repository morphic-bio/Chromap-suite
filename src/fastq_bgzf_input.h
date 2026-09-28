#ifndef FASTQ_BGZF_INPUT_H_
#define FASTQ_BGZF_INPUT_H_

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "mapping_parameters.h"
#include "paired_end_read_provider.h"
#include "sequence_batch.h"

namespace chromap {

// BGZF FASTQ intake. This file is the only Chromap code that uses the BGZF
// reader mirrored from STAR Suite under src/star_input/ (see MIRROR.md there),
// so replacing that copy with a shared FASTQ library changes the implementation
// in fastq_bgzf_input.cc and the build, and nothing else.

// Decides whether the FASTQ files of one lane use the BGZF reader. A file
// qualifies when it is a regular BGZF file whose first record the reader
// parses; the first read names of all files must also pair. kOff never uses
// the reader. kAuto uses it when every file qualifies and otherwise returns
// true with *use_bgzf false and *message saying why. kOn returns false with
// *message describing the file that does not qualify.
bool SelectBgzfFastqInput(FastqBgzfMode mode,
                          const std::vector<std::string> &paths,
                          bool *use_bgzf, std::string *message);

// Inflate workers for each of num_streams files. requested > 0 is the total for
// the lane; 0 derives it from num_threads as num_threads - num_streams, with at
// least one worker per stream from three threads up (below three threads the
// files inflate on the loading thread). The total is split evenly, the
// remainder going to the first streams, as STAR Suite splits it across mates.
std::vector<uint32_t> BgzfInflateWorkersPerStream(int requested,
                                                  int num_threads,
                                                  size_t num_streams);

// One ordered FASTQ file read with the BGZF reader: its members are inflated in
// parallel by `inflate_workers` threads and consumed in file order. Records are
// stored as kseq stores them for Chromap: the read name ends at the first
// whitespace, the rest of the header line is the comment, and records with an
// empty sequence are skipped.
class BgzfFastqStream {
 public:
  BgzfFastqStream();
  ~BgzfFastqStream();

  BgzfFastqStream(const BgzfFastqStream &) = delete;
  BgzfFastqStream &operator=(const BgzfFastqStream &) = delete;

  bool Open(const std::string &path, uint32_t inflate_workers,
            std::string *error);

  // Stores the next record in slot `index` of `batch`. Returns false at a
  // clean end (error empty) or on failure (error set). *source_ordinal is the
  // record's ordinal in the file, counting skipped empty records.
  bool LoadNext(SequenceBatch &batch, uint32_t index, uint64_t *source_ordinal,
                std::string *error);

  // Resets `batch` and fills it with up to max_records records.
  bool LoadBatch(SequenceBatch &batch, uint32_t max_records, uint32_t *loaded,
                 std::string *error);

  const std::string &path() const;

 private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

// Read 1, read 2 and an optional barcode read of one lane from BGZF FASTQ
// files. Each file is inflated in parallel by its own workers and consumed in
// order on its own thread; the records of a batch are then paired as STAR
// Suite pairs its synchronized BGZF streams: the same ordinal in every file,
// the same read-name stem (a trailing /1, /2 or /3 removed), and no file
// ending before the others.
class BgzfPairedEndReadProvider : public PairedEndReadProvider {
 public:
  // An empty barcode_path means no barcode read. inflate_workers has one entry
  // per file in the order read 1, read 2, barcode.
  BgzfPairedEndReadProvider(const std::string &read1_path,
                            const std::string &read2_path,
                            const std::string &barcode_path,
                            const std::vector<uint32_t> &inflate_workers);

  bool Open(std::string *error);

  bool HasBarcode() const override;

  bool LoadBatch(uint32_t max_pairs, SequenceBatch &read_batch1,
                 SequenceBatch &read_batch2, SequenceBatch &barcode_batch,
                 uint32_t &num_loaded_pairs, std::string &error) override;

 private:
  std::vector<std::string> paths_;
  std::vector<uint32_t> inflate_workers_;
  std::vector<std::unique_ptr<BgzfFastqStream>> streams_;
  // Source ordinals of the current batch, per file.
  std::vector<std::vector<uint64_t>> ordinals_;
  uint64_t num_pairs_ = 0;
};

}  // namespace chromap

#endif  // FASTQ_BGZF_INPUT_H_
