#include "atac_lowmem_finalize.h"

#include <dirent.h>
#include <fcntl.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <limits>
#include <sstream>

#include "atac_kway_spill.h"
#include "atac_spill_record.h"
#include "summary_metadata.h"
#include "utils.h"

namespace chromap {

AtacSpillProbe ProbeAtacKwaySpillFile(const std::string &path) {
  AtacSpillProbe probe;
  const int fd = open(path.c_str(), O_RDONLY);
  if (fd < 0) {
    probe.status = AtacSpillProbeStatus::kOpenFailed;
    probe.open_errno = errno;
    return probe;
  }
  struct stat st;
  if (fstat(fd, &st) == 0) {
    probe.bytes = static_cast<uint64_t>(st.st_size);
  }
  AtacKwaySpillFileHeaderV1 header;
  memset(&header, 0, sizeof(header));
  const ssize_t got = pread(fd, &header, sizeof(header), 0);
  close(fd);
  if (got < static_cast<ssize_t>(sizeof(header.magic)) ||
      memcmp(header.magic, kAtacKwaySpillMagicV1, sizeof(header.magic)) != 0) {
    probe.status = AtacSpillProbeStatus::kNotKway;
    return probe;
  }
  // The same checks as OverflowReader::ConsumeAtacSpillFilePrefixIfPresent.
  const uint16_t known_schema = static_cast<uint16_t>(
      kAtacSpillSchemaHasBamPair | kAtacSpillSchemaHasYHit |
      kAtacSpillSchemaIsBulk | kAtacSpillSchemaHasRawBarcodeEvidence);
  if (got != static_cast<ssize_t>(sizeof(header)) ||
      header.format_version != kAtacKwaySpillFormatVersion ||
      header.fixed_header_bytes != sizeof(header) ||
      header.record_codec_version != kAtacKwaySpillRecordCodecVersion ||
      header.endian_marker != kAtacKwaySpillEndianMarker ||
      (header.schema_mask & ~known_schema) != 0 || header.flags != 0 ||
      header.reserved != 0) {
    probe.status = AtacSpillProbeStatus::kInvalid;
    return probe;
  }
  probe.status = AtacSpillProbeStatus::kKway;
  probe.reference_id = header.reference_id;
  probe.schema_mask = header.schema_mask;
  return probe;
}

int RequestedLowMemFinalizeThreads(const MappingParameters &parameters) {
  int requested = parameters.low_mem_finalize_threads;
  if (requested == 0) {
    requested = parameters.num_threads;
  }
  return std::max(1, requested);
}

namespace {

uint64_t CountOpenFileDescriptors() {
  DIR *d = opendir("/proc/self/fd");
  if (d == nullptr) {
    return 64;
  }
  uint64_t count = 0;
  while (struct dirent *entry = readdir(d)) {
    if (entry->d_name[0] != '.') {
      ++count;
    }
  }
  closedir(d);
  // The directory stream itself holds one descriptor.
  return count > 0 ? count - 1 : 0;
}

}  // namespace

bool PlanLowMemFinalizeWorkers(int requested_workers, size_t num_tasks,
                               size_t max_files_per_reference,
                               uint32_t reference_with_most_files,
                               LowMemFinalizeWorkerPlan *plan,
                               std::string *error) {
  LowMemFinalizeWorkerPlan result;
  const uint64_t wanted = static_cast<uint64_t>(std::max<size_t>(
      1, std::min<size_t>(static_cast<size_t>(std::max(1, requested_workers)),
                          std::max<size_t>(1, num_tasks))));
  result.open_files_now = CountOpenFileDescriptors();
  struct rlimit limit;
  uint64_t budget = std::numeric_limits<uint64_t>::max();
  if (getrlimit(RLIMIT_NOFILE, &limit) == 0 &&
      limit.rlim_cur != RLIM_INFINITY) {
    result.open_file_limit = static_cast<uint64_t>(limit.rlim_cur);
    const uint64_t reserve =
        result.open_files_now + std::max<uint64_t>(64, result.open_files_now);
    budget = result.open_file_limit > reserve
                 ? result.open_file_limit - reserve
                 : 0;
  }
  // Each task: every spill file of its reference, its partition, and one
  // spare for the temporary buffers of the standard library.
  const uint64_t per_task = static_cast<uint64_t>(max_files_per_reference) + 2;
  const uint64_t fit = budget / per_task;
  if (fit == 0) {
    std::ostringstream message;
    message << "Low-memory finalization needs " << per_task
            << " open files at once for reference " << reference_with_most_files
            << " (" << max_files_per_reference
            << " spill files), but the open-file limit is "
            << result.open_file_limit << " with " << result.open_files_now
            << " files already open. Raise the limit (ulimit -n) or use a "
               "larger --low-mem-ram so that fewer spill files are written";
    *error = message.str();
    return false;
  }
  result.workers = static_cast<int>(std::min(wanted, fit));
  result.limited_by_open_files = fit < wanted;
  *plan = result;
  return true;
}

void RaiseOpenFileSoftLimitToHard() {
  struct rlimit limit;
  if (getrlimit(RLIMIT_NOFILE, &limit) != 0) {
    return;
  }
  if (limit.rlim_cur != limit.rlim_max) {
    limit.rlim_cur = limit.rlim_max;
    (void)setrlimit(RLIMIT_NOFILE, &limit);
  }
}

LowMemPermitScope::LowMemPermitScope(const MappingParameters &parameters) {
  if (!parameters.PermitHooksEnabled()) {
    return;
  }
  release_hook_ = parameters.permit_release_hook;
  hook_ctx_ = parameters.permit_hook_ctx;
  wait_ns_ = parameters.permit_acquire_hook(hook_ctx_);
  start_sec_ = GetRealTime();
  held_ = true;
}

LowMemPermitScope::~LowMemPermitScope() { Release(0, 0); }

void LowMemPermitScope::Release(uint64_t work_units, uint64_t work_bytes) {
  if (!held_) {
    return;
  }
  held_ = false;
  const double end_sec = GetRealTime();
  const uint64_t work_ns =
      end_sec > start_sec_
          ? static_cast<uint64_t>((end_sec - start_sec_) * 1e9)
          : 0ULL;
  release_hook_(hook_ctx_, wait_ns_, work_units, work_bytes, work_ns);
}

bool CreateLowMemPartitionFile(const std::string &directory, uint32_t rid,
                               std::string *path, FILE **file,
                               std::string *error) {
  std::ostringstream pattern;
  pattern << directory;
  if (!directory.empty() && directory.back() != '/') {
    pattern << '/';
  }
  pattern << "chromap_lowmem_part_" << getpid() << "_" << rid << "_XXXXXX";
  std::string name = pattern.str();
  std::vector<char> writable(name.begin(), name.end());
  writable.push_back('\0');
  const int fd = mkstemp(writable.data());
  if (fd < 0) {
    *error = "Cannot create low-memory finalization partition in " +
             directory + ": " + std::strerror(errno);
    return false;
  }
  *path = writable.data();
  *file = fdopen(fd, "wb");
  if (*file == nullptr) {
    const int saved = errno;
    close(fd);
    unlink(path->c_str());
    *error = "Cannot open low-memory finalization partition " + *path + ": " +
             std::strerror(saved);
    return false;
  }
  return true;
}

bool AppendFileContents(FILE *destination, const std::string &path,
                        std::vector<char> *buffer, uint64_t *bytes_copied,
                        std::string *error) {
  *bytes_copied = 0;
  FILE *source = fopen(path.c_str(), "rb");
  if (source == nullptr) {
    *error = "Cannot reopen low-memory finalization partition " + path + ": " +
             std::strerror(errno);
    return false;
  }
  if (buffer->empty()) {
    buffer->resize(4u << 20);
  }
  bool ok = true;
  while (true) {
    const size_t got = fread(buffer->data(), 1, buffer->size(), source);
    if (got > 0 && fwrite(buffer->data(), 1, got, destination) != got) {
      *error = "Cannot append low-memory finalization partition " + path;
      ok = false;
      break;
    }
    *bytes_copied += got;
    if (got < buffer->size()) {
      if (ferror(source)) {
        *error = "Cannot read low-memory finalization partition " + path;
        ok = false;
      }
      break;
    }
  }
  fclose(source);
  return ok;
}

std::string DirectoryOfPath(const std::string &path) {
  const size_t slash = path.find_last_of('/');
  if (slash == std::string::npos) {
    return ".";
  }
  if (slash == 0) {
    return "/";
  }
  return path.substr(0, slash);
}

void AtacSummaryDeltaLog::Add(uint64_t barcode, int type, uint64_t change) {
  auto found = index_.find(barcode);
  size_t slot;
  if (found == index_.end()) {
    slot = keys_.size();
    index_.emplace(barcode, slot);
    AtacSummaryDelta key;
    key.barcode = barcode;
    keys_.push_back(key);
  } else {
    slot = found->second;
  }
  AtacSummaryDelta &key = keys_[slot];
  if (type == SUMMARY_METADATA_DUP) {
    key.duplicate += change;
  } else if (type == SUMMARY_METADATA_LOWMAPQ) {
    key.lowmapq += change;
  } else if (type == SUMMARY_METADATA_MAPPED) {
    key.mapped += change;
  } else {
    ExitWithMessage("Unexpected summary field in the low-memory merge");
  }
}

void AtacSummaryDeltaLog::Append(const AtacSummaryDelta &delta) {
  index_.emplace(delta.barcode, keys_.size());
  keys_.push_back(delta);
}

void AtacSummaryDeltaLog::Clear() {
  std::vector<AtacSummaryDelta>().swap(keys_);
  std::unordered_map<uint64_t, size_t>().swap(index_);
}

}  // namespace chromap
