#include "sequence_batch.h"

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <zlib.h>

#include "utils.h"

namespace chromap {
namespace {

const char kReferenceSidecarMagic[8] = {'C', 'H', 'R', 'R',
                                        'E', 'F', 'S', '\0'};
const uint32_t kReferenceSidecarVersion = 1;
const uint32_t kEndianMarker = 0x01020304U;
const uint64_t kFnvOffset = 1469598103934665603ULL;
const uint64_t kFnvPrime = 1099511628211ULL;
const size_t kReferenceDirectAlignment = 4096;
const size_t kReferenceReadBlockBytes = 256U * 1024U * 1024U;

#pragma pack(push, 1)
struct ReferenceSidecarHeaderV1 {
  char magic[8];
  uint32_t version;
  uint32_t header_bytes;
  uint32_t endian_marker;
  uint32_t sequence_count;
  uint64_t total_bases;
  uint64_t directory_offset;
  uint64_t names_offset;
  uint64_t bases_offset;
  uint64_t file_bytes;
  uint64_t reference_fingerprint;
  uint64_t reserved[4];
};

struct ReferenceSidecarEntryV1 {
  uint64_t base_offset;
  uint64_t base_length;
  uint64_t name_offset;
  uint32_t name_bytes;
  uint32_t sequence_crc32;
};
#pragma pack(pop)

static_assert(sizeof(ReferenceSidecarHeaderV1) == 104,
              "unexpected materialized-reference header layout");
static_assert(sizeof(ReferenceSidecarEntryV1) == 32,
              "unexpected materialized-reference entry layout");
static_assert(kReferenceReadBlockBytes % kReferenceDirectAlignment == 0,
              "reference read block must be direct-I/O aligned");

struct ReferenceBaseSpan {
  uint64_t file_offset;
  uint64_t bytes;
  char *destination;
};

struct ReferenceReadTask {
  uint64_t file_offset;
  size_t bytes;
  char *destination;
};

void SetError(std::string *error, const std::string &message) {
  if (error != nullptr) *error = message;
}

std::string ErrnoMessage(const std::string &prefix) {
  return prefix + ": " + std::strerror(errno);
}

bool AddU64(uint64_t a, uint64_t b, uint64_t *result) {
  if (b > std::numeric_limits<uint64_t>::max() - a) return false;
  *result = a + b;
  return true;
}

bool MultiplyU64(uint64_t a, uint64_t b, uint64_t *result) {
  if (a != 0 && b > std::numeric_limits<uint64_t>::max() / a) return false;
  *result = a * b;
  return true;
}

void HashByte(uint8_t value, uint64_t *hash) {
  *hash ^= value;
  *hash *= kFnvPrime;
}

void HashU32(uint32_t value, uint64_t *hash) {
  for (unsigned int i = 0; i < 4; ++i) {
    HashByte(static_cast<uint8_t>(value >> (8 * i)), hash);
  }
}

void HashU64(uint64_t value, uint64_t *hash) {
  for (unsigned int i = 0; i < 8; ++i) {
    HashByte(static_cast<uint8_t>(value >> (8 * i)), hash);
  }
}

uint32_t SequenceCrc32(const char *data, uint64_t length) {
  uLong crc = crc32(0L, Z_NULL, 0);
  while (length > 0) {
    const uInt chunk = static_cast<uInt>(std::min<uint64_t>(
        length, std::numeric_limits<uInt>::max()));
    crc = crc32(crc, reinterpret_cast<const Bytef *>(data), chunk);
    data += chunk;
    length -= chunk;
  }
  return static_cast<uint32_t>(crc);
}

uint64_t ReferenceFingerprint(
    const SequenceBatch &reference,
    const std::vector<ReferenceSidecarEntryV1> &entries) {
  uint64_t hash = kFnvOffset;
  HashU32(kReferenceSidecarVersion, &hash);
  HashU64(reference.GetNumSequences(), &hash);
  HashU64(reference.GetNumBases(), &hash);
  for (uint32_t rid = 0; rid < reference.GetNumSequences(); ++rid) {
    const uint32_t name_length = reference.GetSequenceNameLengthAt(rid);
    HashU32(name_length, &hash);
    const char *name = reference.GetSequenceNameAt(rid);
    for (uint32_t i = 0; i < name_length; ++i) {
      HashByte(static_cast<uint8_t>(name[i]), &hash);
    }
    HashU64(reference.GetSequenceLengthAt(rid), &hash);
    HashU32(entries[rid].sequence_crc32, &hash);
  }
  return hash == 0 ? 1 : hash;
}

bool WriteExact(FILE *file, const void *data, size_t bytes) {
  return bytes == 0 || fwrite(data, 1, bytes, file) == bytes;
}

bool ReadExact(FILE *file, void *data, size_t bytes) {
  return bytes == 0 || fread(data, 1, bytes, file) == bytes;
}

bool PreadExact(int fd, void *data, size_t bytes, uint64_t offset,
                size_t required_alignment = 1) {
  char *cursor = static_cast<char *>(data);
  size_t remaining = bytes;
  while (remaining > 0) {
    const ssize_t got = pread(fd, cursor, remaining,
                              static_cast<off_t>(offset));
    if (got < 0 && errno == EINTR) continue;
    if (got <= 0) return false;
    if (required_alignment > 1 &&
        static_cast<size_t>(got) < remaining &&
        static_cast<size_t>(got) % required_alignment != 0) {
      errno = EIO;
      return false;
    }
    cursor += got;
    remaining -= static_cast<size_t>(got);
    offset += static_cast<uint64_t>(got);
  }
  return true;
}

void CopyReferenceBlock(const char *block, uint64_t block_offset,
                        size_t block_bytes,
                        const std::vector<ReferenceBaseSpan> &spans) {
  const uint64_t block_end = block_offset + block_bytes;
  for (size_t i = 0; i < spans.size(); ++i) {
    const ReferenceBaseSpan &span = spans[i];
    const uint64_t span_end = span.file_offset + span.bytes;
    const uint64_t copy_begin = std::max(block_offset, span.file_offset);
    const uint64_t copy_end = std::min(block_end, span_end);
    if (copy_begin < copy_end) {
      memcpy(span.destination + (copy_begin - span.file_offset),
             block + (copy_begin - block_offset), copy_end - copy_begin);
    }
  }
}

bool LoadReferenceBasesDirectParallel(
    const std::string &path, int buffered_fd, uint64_t payload_begin,
    uint64_t file_bytes, const std::vector<ReferenceBaseSpan> &spans,
    int requested_threads, int *threads_used, std::string *error) {
#ifdef O_DIRECT
  const int direct_fd = open(path.c_str(), O_RDONLY | O_DIRECT);
  if (direct_fd < 0) {
    SetError(error, ErrnoMessage("O_DIRECT open failed"));
    return false;
  }

  const uint64_t aligned_begin =
      payload_begin - payload_begin % kReferenceDirectAlignment;
  const uint64_t aligned_end =
      file_bytes - file_bytes % kReferenceDirectAlignment;
  if (aligned_begin >= aligned_end) {
    close(direct_fd);
    SetError(error, "materialized reference payload is smaller than one "
                    "direct-I/O alignment block");
    return false;
  }
  const uint64_t direct_bytes = aligned_end - aligned_begin;
  const uint64_t num_blocks =
      (direct_bytes + kReferenceReadBlockBytes - 1) /
      kReferenceReadBlockBytes;
  const int loader_threads = std::max<int>(
      1, std::min<uint64_t>(std::max(1, requested_threads), num_blocks));
  if (threads_used != nullptr) *threads_used = loader_threads;

  int read_failed = 0;
  int allocation_failed = 0;
#pragma omp parallel num_threads(loader_threads) reduction(| : read_failed) reduction(| : allocation_failed)
  {
    void *raw_buffer = nullptr;
    const int allocation_result =
        posix_memalign(&raw_buffer, kReferenceDirectAlignment,
                       kReferenceReadBlockBytes);
    if (allocation_result != 0) {
      allocation_failed = 1;
      read_failed = 1;
    }

#pragma omp for schedule(static)
    for (int64_t block_id = 0;
         block_id < static_cast<int64_t>(num_blocks); ++block_id) {
      if (raw_buffer == nullptr) continue;
      const uint64_t block_offset =
          aligned_begin + static_cast<uint64_t>(block_id) *
                              kReferenceReadBlockBytes;
      const size_t block_bytes = static_cast<size_t>(std::min<uint64_t>(
          kReferenceReadBlockBytes, aligned_end - block_offset));
      if (!PreadExact(direct_fd, raw_buffer, block_bytes, block_offset,
                      kReferenceDirectAlignment)) {
        read_failed = 1;
        continue;
      }
      CopyReferenceBlock(static_cast<const char *>(raw_buffer), block_offset,
                         block_bytes, spans);
    }
    free(raw_buffer);
  }
  close(direct_fd);

  if (read_failed) {
    if (allocation_failed) {
      SetError(error,
               "cannot allocate aligned materialized-reference buffer");
    } else {
      SetError(error, "parallel O_DIRECT pread failed");
    }
    return false;
  }

  if (aligned_end < file_bytes) {
    const size_t tail_bytes = static_cast<size_t>(file_bytes - aligned_end);
    std::vector<char> tail(tail_bytes);
    if (!PreadExact(buffered_fd, tail.data(), tail.size(), aligned_end)) {
      SetError(error, "cannot read materialized-reference tail");
      return false;
    }
    CopyReferenceBlock(tail.data(), aligned_end, tail.size(), spans);
  }
  return true;
#else
  (void)path;
  (void)buffered_fd;
  (void)payload_begin;
  (void)file_bytes;
  (void)spans;
  (void)requested_threads;
  (void)threads_used;
  SetError(error, "O_DIRECT is unavailable on this platform");
  return false;
#endif
}

bool LoadReferenceBasesBufferedParallel(
    int fd, const std::vector<ReferenceBaseSpan> &spans,
    int requested_threads, int *threads_used) {
  std::vector<ReferenceReadTask> tasks;
  for (size_t i = 0; i < spans.size(); ++i) {
    uint64_t offset = 0;
    while (offset < spans[i].bytes) {
      const size_t bytes = static_cast<size_t>(std::min<uint64_t>(
          kReferenceReadBlockBytes, spans[i].bytes - offset));
      tasks.push_back({spans[i].file_offset + offset, bytes,
                       spans[i].destination + offset});
      offset += bytes;
    }
  }
  const int loader_threads = std::max<int>(
      1, std::min<size_t>(std::max(1, requested_threads), tasks.size()));
  if (threads_used != nullptr) *threads_used = loader_threads;
  int read_failed = 0;
#pragma omp parallel for schedule(dynamic, 1) num_threads(loader_threads) reduction(| : read_failed)
  for (int64_t task_id = 0; task_id < static_cast<int64_t>(tasks.size());
       ++task_id) {
    const ReferenceReadTask &task = tasks[task_id];
    if (!PreadExact(fd, task.destination, task.bytes, task.file_offset)) {
      read_failed = 1;
    }
  }
  return read_failed == 0;
}

}  // namespace

bool SequenceBatch::SaveMaterializedReference(
    const std::string &path, MaterializedReferenceInfo *info,
    std::string *error) const {
  if (path.empty()) {
    SetError(error, "materialized reference path is empty");
    return false;
  }
  if (num_loaded_sequences_ == 0 ||
      num_loaded_sequences_ != sequence_batch_.size()) {
    SetError(error, "materialized reference requires a complete reference");
    return false;
  }
  if (num_loaded_sequences_ > std::numeric_limits<uint32_t>::max()) {
    SetError(error, "materialized reference has too many sequences");
    return false;
  }

  std::vector<ReferenceSidecarEntryV1> entries(num_loaded_sequences_);
  uint64_t names_bytes = 0;
  uint64_t bases_bytes = 0;
  for (uint32_t rid = 0; rid < num_loaded_sequences_; ++rid) {
    const uint32_t name_bytes = GetSequenceNameLengthAt(rid);
    const uint64_t base_bytes = GetSequenceLengthAt(rid);
    if (!AddU64(names_bytes, name_bytes, &names_bytes) ||
        !AddU64(bases_bytes, base_bytes, &bases_bytes)) {
      SetError(error, "materialized reference size overflows uint64");
      return false;
    }
    entries[rid].name_bytes = name_bytes;
    entries[rid].base_length = base_bytes;
    entries[rid].sequence_crc32 =
        SequenceCrc32(GetSequenceAt(rid), base_bytes);
  }

  uint64_t directory_bytes = 0;
  if (!MultiplyU64(entries.size(), sizeof(ReferenceSidecarEntryV1),
                   &directory_bytes)) {
    SetError(error, "materialized reference directory size overflows uint64");
    return false;
  }

  ReferenceSidecarHeaderV1 header = {};
  memcpy(header.magic, kReferenceSidecarMagic, sizeof(header.magic));
  header.version = kReferenceSidecarVersion;
  header.header_bytes = sizeof(header);
  header.endian_marker = kEndianMarker;
  header.sequence_count = static_cast<uint32_t>(num_loaded_sequences_);
  header.total_bases = bases_bytes;
  header.directory_offset = sizeof(header);
  if (!AddU64(header.directory_offset, directory_bytes,
              &header.names_offset) ||
      !AddU64(header.names_offset, names_bytes, &header.bases_offset) ||
      !AddU64(header.bases_offset, bases_bytes, &header.file_bytes)) {
    SetError(error, "materialized reference file size overflows uint64");
    return false;
  }

  uint64_t next_name = header.names_offset;
  uint64_t next_base = header.bases_offset;
  for (uint32_t rid = 0; rid < num_loaded_sequences_; ++rid) {
    entries[rid].name_offset = next_name;
    entries[rid].base_offset = next_base;
    next_name += entries[rid].name_bytes;
    next_base += entries[rid].base_length;
  }
  header.reference_fingerprint = ReferenceFingerprint(*this, entries);

  const std::string temporary_path =
      path + ".tmp." + std::to_string(static_cast<unsigned long long>(getpid()));
  FILE *file = fopen(temporary_path.c_str(), "wb");
  if (file == nullptr) {
    SetError(error, ErrnoMessage("cannot create " + temporary_path));
    return false;
  }

  bool ok = WriteExact(file, &header, sizeof(header)) &&
            WriteExact(file, entries.data(),
                       entries.size() * sizeof(entries[0]));
  for (uint32_t rid = 0; ok && rid < num_loaded_sequences_; ++rid) {
    ok = WriteExact(file, GetSequenceNameAt(rid), entries[rid].name_bytes);
  }
  for (uint32_t rid = 0; ok && rid < num_loaded_sequences_; ++rid) {
    ok = WriteExact(file, GetSequenceAt(rid), entries[rid].base_length);
  }
  if (ok && fflush(file) != 0) ok = false;
  if (fclose(file) != 0) ok = false;
  if (!ok) {
    const int saved_errno = errno;
    std::remove(temporary_path.c_str());
    errno = saved_errno;
    SetError(error, ErrnoMessage("cannot write " + temporary_path));
    return false;
  }
  if (std::rename(temporary_path.c_str(), path.c_str()) != 0) {
    const int saved_errno = errno;
    std::remove(temporary_path.c_str());
    errno = saved_errno;
    SetError(error, ErrnoMessage("cannot publish " + path));
    return false;
  }

  if (info != nullptr) {
    info->fingerprint = header.reference_fingerprint;
    info->num_bases = header.total_bases;
    info->num_sequences = header.sequence_count;
  }
  return true;
}

bool SequenceBatch::LoadMaterializedReference(
    const std::string &path, int num_threads, MaterializedReferenceInfo *info,
    std::string *error) {
  const double real_start_time = GetRealTime();
  if (!sequence_batch_.empty() || num_loaded_sequences_ != 0) {
    SetError(error, "materialized reference destination is not empty");
    return false;
  }

  FILE *file = fopen(path.c_str(), "rb");
  if (file == nullptr) {
    SetError(error, ErrnoMessage("cannot open " + path));
    return false;
  }
  const int fd = fileno(file);
  struct stat file_stat = {};
  ReferenceSidecarHeaderV1 header = {};
  bool ok = fd >= 0 && fstat(fd, &file_stat) == 0 &&
            ReadExact(file, &header, sizeof(header));
  if (!ok) {
    fclose(file);
    SetError(error, "cannot read materialized reference header from " + path);
    return false;
  }
  if (memcmp(header.magic, kReferenceSidecarMagic, sizeof(header.magic)) != 0 ||
      header.version != kReferenceSidecarVersion ||
      header.header_bytes != sizeof(header) ||
      header.endian_marker != kEndianMarker) {
    fclose(file);
    SetError(error, "unsupported materialized reference header in " + path);
    return false;
  }
  if (header.sequence_count == 0 || header.total_bases == 0 ||
      header.directory_offset != sizeof(header) ||
      header.file_bytes != static_cast<uint64_t>(file_stat.st_size)) {
    fclose(file);
    SetError(error, "invalid materialized reference dimensions in " + path);
    return false;
  }

  uint64_t directory_bytes = 0;
  uint64_t expected_names_offset = 0;
  if (!MultiplyU64(header.sequence_count, sizeof(ReferenceSidecarEntryV1),
                   &directory_bytes) ||
      !AddU64(header.directory_offset, directory_bytes,
              &expected_names_offset) ||
      header.names_offset != expected_names_offset ||
      header.names_offset > header.bases_offset ||
      header.bases_offset > header.file_bytes) {
    fclose(file);
    SetError(error, "invalid materialized reference offsets in " + path);
    return false;
  }

  std::vector<ReferenceSidecarEntryV1> entries(header.sequence_count);
  if (!ReadExact(file, entries.data(), entries.size() * sizeof(entries[0]))) {
    fclose(file);
    SetError(error, "truncated materialized reference directory in " + path);
    return false;
  }

  uint64_t next_name = header.names_offset;
  uint64_t next_base = header.bases_offset;
  uint64_t observed_bases = 0;
  for (uint32_t rid = 0; rid < header.sequence_count; ++rid) {
    uint64_t name_end = 0;
    uint64_t base_end = 0;
    if (entries[rid].name_bytes == 0 ||
        entries[rid].base_length == 0 ||
        entries[rid].base_length > std::numeric_limits<uint32_t>::max() ||
        entries[rid].name_offset != next_name ||
        entries[rid].base_offset != next_base ||
        !AddU64(next_name, entries[rid].name_bytes, &name_end) ||
        !AddU64(next_base, entries[rid].base_length, &base_end) ||
        name_end > header.bases_offset || base_end > header.file_bytes ||
        !AddU64(observed_bases, entries[rid].base_length, &observed_bases)) {
      fclose(file);
      SetError(error, "invalid materialized reference entry in " + path);
      return false;
    }
    next_name = name_end;
    next_base = base_end;
  }
  if (next_name != header.bases_offset || next_base != header.file_bytes ||
      observed_bases != header.total_bases) {
    fclose(file);
    SetError(error, "inconsistent materialized reference payload in " + path);
    return false;
  }

  const uint64_t names_bytes_u64 = header.bases_offset - header.names_offset;
  if (names_bytes_u64 > std::numeric_limits<size_t>::max()) {
    fclose(file);
    SetError(error, "materialized reference names exceed addressable memory");
    return false;
  }
  std::vector<char> names(static_cast<size_t>(names_bytes_u64));
  if (!PreadExact(fd, names.data(), names.size(), header.names_offset)) {
    fclose(file);
    SetError(error, "cannot read materialized reference names from " + path);
    return false;
  }

  sequence_batch_.reserve(header.sequence_count);
  negative_sequence_batch_.assign(header.sequence_count, "");
  negative_sequence_prepared_.assign(header.sequence_count, 0);
  std::vector<ReferenceBaseSpan> base_spans;
  base_spans.reserve(header.sequence_count);
  for (uint32_t rid = 0; rid < header.sequence_count; ++rid) {
    sequence_batch_.emplace_back(static_cast<kseq_t *>(calloc(1, sizeof(kseq_t))));
    if (sequence_batch_.back() == nullptr) {
      fclose(file);
      SetError(error, "cannot allocate materialized reference sequence");
      return false;
    }
    kseq_t *sequence = sequence_batch_.back();
    const size_t name_offset = static_cast<size_t>(
        entries[rid].name_offset - header.names_offset);
    AssignKstring(sequence->name, names.data() + name_offset,
                  entries[rid].name_bytes);
    AssignKstring(sequence->comment, nullptr, 0);
    AssignKstring(sequence->qual, nullptr, 0);
    sequence->seq.l = entries[rid].base_length;
    sequence->seq.m = entries[rid].base_length + 1;
    sequence->seq.s = static_cast<char *>(malloc(sequence->seq.m));
    if (sequence->seq.s == nullptr) {
      fclose(file);
      SetError(error, "cannot allocate materialized reference bases");
      return false;
    }
    sequence->seq.s[sequence->seq.l] = '\0';
    sequence->id = rid;
    base_spans.push_back({entries[rid].base_offset,
                          entries[rid].base_length, sequence->seq.s});
  }

  int loader_threads = 1;
  std::string direct_error;
  const bool used_direct = LoadReferenceBasesDirectParallel(
      path, fd, header.bases_offset, header.file_bytes, base_spans,
      num_threads, &loader_threads, &direct_error);
  if (!used_direct) {
    if (!LoadReferenceBasesBufferedParallel(fd, base_spans, num_threads,
                                            &loader_threads)) {
      fclose(file);
      SetError(error, "cannot read materialized reference bases from " + path);
      return false;
    }
    std::cerr << "Materialized reference input path: parallel buffered "
                 "positioned-reader fallback ("
              << direct_error << "), block size "
              << kReferenceReadBlockBytes << " bytes, " << loader_threads
              << " thread(s).\n";
  } else {
    std::cerr << "Materialized reference input path: parallel direct I/O "
                 "positioned reader, block size "
              << kReferenceReadBlockBytes << " bytes, " << loader_threads
              << " thread(s).\n";
  }
  fclose(file);
  max_num_sequences_ = header.sequence_count;
  num_loaded_sequences_ = header.sequence_count;
  total_num_loaded_sequences_ = header.sequence_count;
  num_bases_ = header.total_bases;
  const uint64_t fingerprint = ReferenceFingerprint(*this, entries);
  if (fingerprint != header.reference_fingerprint) {
    SetError(error, "materialized reference fingerprint mismatch in " + path);
    return false;
  }
  if (info != nullptr) {
    info->fingerprint = fingerprint;
    info->num_bases = header.total_bases;
    info->num_sequences = header.sequence_count;
  }
  std::cerr << "Loaded materialized reference successfully in "
            << GetRealTime() - real_start_time << "s using " << loader_threads
            << " thread(s), number of sequences: " << header.sequence_count
            << ", number of bases: " << header.total_bases << ".\n";
  return true;
}

}  // namespace chromap
