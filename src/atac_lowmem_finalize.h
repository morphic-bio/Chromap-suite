#ifndef ATAC_LOWMEM_FINALIZE_H_
#define ATAC_LOWMEM_FINALIZE_H_

// Record-independent helpers for the per-reference low-memory finalisation of
// paired-end ATAC spills (MappingWriter<AtacSpillRecord>): spill-file probes,
// the worker plan (threads, host permits, open-file budget), temporary
// partition files and the ordered log of summary barcodes that are new to the
// summary table.

#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <vector>

#include "mapping_parameters.h"

namespace chromap {

// Result of reading the fixed header of one low-memory spill file.
enum class AtacSpillProbeStatus {
  kKway,         // an ATAC k-way spill file with a valid header
  kNotKway,      // any other spill format (the serial merge handles it)
  kInvalid,      // k-way magic with an invalid header
  kOpenFailed,   // the file could not be opened; open_errno is set
};

struct AtacSpillProbe {
  AtacSpillProbeStatus status = AtacSpillProbeStatus::kNotKway;
  uint32_t reference_id = 0;
  uint16_t schema_mask = 0;
  uint64_t bytes = 0;
  int open_errno = 0;
};

AtacSpillProbe ProbeAtacKwaySpillFile(const std::string &path);

// Workers requested by the parameters: 0 means num_threads; never below 1.
int RequestedLowMemFinalizeThreads(const MappingParameters &parameters);

// How many per-reference tasks may run at once. Every task keeps all spill
// files of its reference open, plus its partition file.
struct LowMemFinalizeWorkerPlan {
  int workers = 1;
  bool limited_by_open_files = false;
  uint64_t open_file_limit = 0;  // soft RLIMIT_NOFILE; 0 when unlimited
  uint64_t open_files_now = 0;
};

// Returns false (with a message) when even one task cannot open all spill
// files of the reference with the most files.
bool PlanLowMemFinalizeWorkers(int requested_workers, size_t num_tasks,
                               size_t max_files_per_reference,
                               uint32_t reference_with_most_files,
                               LowMemFinalizeWorkerPlan *plan,
                               std::string *error);

// Raises the soft open-file limit to the hard limit. Used by the command-line
// front ends only; the library never changes process limits.
void RaiseOpenFileSoftLimitToHard();

// Holds one host permit (MappingParameters permit hooks) for its lifetime and
// reports the work done when it is released. Does nothing without hooks.
class LowMemPermitScope {
 public:
  explicit LowMemPermitScope(const MappingParameters &parameters);
  ~LowMemPermitScope();
  void Release(uint64_t work_units, uint64_t work_bytes);

 private:
  LowMemPermitScope(const LowMemPermitScope &) = delete;
  LowMemPermitScope &operator=(const LowMemPermitScope &) = delete;

  void (*release_hook_)(void *, uint64_t, uint64_t, uint64_t,
                        uint64_t) = nullptr;
  void *hook_ctx_ = nullptr;
  uint64_t wait_ns_ = 0;
  double start_sec_ = 0.0;
  bool held_ = false;
};

// A temporary per-reference output file, written by one task and appended to
// the final output in reference order.
bool CreateLowMemPartitionFile(const std::string &directory, uint32_t rid,
                               std::string *path, FILE **file,
                               std::string *error);

// Appends the whole content of `path` to `destination` with a buffered copy.
bool AppendFileContents(FILE *destination, const std::string &path,
                        std::vector<char> *buffer, uint64_t *bytes_copied,
                        std::string *error);

// Directory part of a path ("." when it has none).
std::string DirectoryOfPath(const std::string &path);

// Summary counts for barcodes that had no row in the summary table when the
// parallel merge started, in the order the serial merge would first update
// them, aggregated per barcode.
struct AtacSummaryNewKey {
  uint64_t barcode = 0;
  uint64_t duplicate = 0;
  uint64_t lowmapq = 0;
  uint64_t mapped = 0;
};

class AtacSummaryNewKeyLog {
 public:
  // type is SUMMARY_METADATA_DUP, SUMMARY_METADATA_LOWMAPQ or
  // SUMMARY_METADATA_MAPPED.
  void Add(uint64_t barcode, int type, uint64_t change);
  const std::vector<AtacSummaryNewKey> &keys() const { return keys_; }

 private:
  std::vector<AtacSummaryNewKey> keys_;
  std::unordered_map<uint64_t, size_t> index_;
};

}  // namespace chromap

#endif  // ATAC_LOWMEM_FINALIZE_H_
