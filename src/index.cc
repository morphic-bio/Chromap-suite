#include "index.h"

#include <assert.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "minimizer_generator.h"

namespace chromap {
namespace {

const char kIndexReferenceFooterMagic[8] = {'C', 'H', 'R', 'I',
                                           'D', 'X', 'R', 'F'};
const uint32_t kIndexReferenceFooterVersion = 1;
const size_t kIndexDirectAlignment = 4096;
const size_t kIndexDirectBlockBytes = 256U * 1024U * 1024U;

#pragma pack(push, 1)
struct IndexPrefixV1 {
  int32_t kmer_size;
  int32_t window_size;
  uint32_t declared_lookup_size;
  uint32_t n_buckets;
  uint32_t size;
  uint32_t n_occupied;
  uint32_t upper_bound;
};

struct IndexReferenceFooterV1 {
  char magic[8];
  uint32_t version;
  uint32_t footer_bytes;
  uint64_t reference_fingerprint;
  uint64_t total_bases;
  uint32_t sequence_count;
  uint32_t reserved;
};
#pragma pack(pop)

static_assert(sizeof(int) == sizeof(int32_t),
              "Chromap index requires a 32-bit int");
static_assert(sizeof(IndexPrefixV1) == 28,
              "unexpected Chromap index prefix layout");
static_assert(sizeof(IndexReferenceFooterV1) == 40,
              "unexpected Chromap index reference footer layout");
static_assert(kIndexDirectBlockBytes % kIndexDirectAlignment == 0,
              "direct-I/O block must be alignment-sized");

struct IndexDiskLayout {
  IndexPrefixV1 prefix = {};
  uint32_t occurrence_count = 0;
  uint64_t flags_offset = 0;
  uint64_t flags_bytes = 0;
  uint64_t keys_offset = 0;
  uint64_t keys_bytes = 0;
  uint64_t values_offset = 0;
  uint64_t values_bytes = 0;
  uint64_t occurrence_count_offset = 0;
  uint64_t occurrences_offset = 0;
  uint64_t occurrences_bytes = 0;
  uint64_t payload_end = 0;
  uint64_t file_bytes = 0;
  bool has_reference_footer = false;
  IndexReferenceFooterV1 reference_footer = {};
};

struct IndexLoadSpan {
  uint64_t file_offset;
  uint64_t bytes;
  char *destination;
};

void SetIndexLoadError(std::string *error, const std::string &message) {
  if (error != nullptr) *error = message;
}

bool AddIndexU64(uint64_t left, uint64_t right, uint64_t *result) {
  if (right > std::numeric_limits<uint64_t>::max() - left) return false;
  *result = left + right;
  return true;
}

bool MultiplyIndexU64(uint64_t left, uint64_t right, uint64_t *result) {
  if (left != 0 && right > std::numeric_limits<uint64_t>::max() / left) {
    return false;
  }
  *result = left * right;
  return true;
}

bool PreadIndexExact(int fd, void *destination, size_t bytes,
                     uint64_t offset, std::string *error) {
  char *cursor = static_cast<char *>(destination);
  size_t remaining = bytes;
  while (remaining > 0) {
    const ssize_t got =
        pread(fd, cursor, remaining, static_cast<off_t>(offset));
    if (got < 0 && errno == EINTR) continue;
    if (got < 0) {
      SetIndexLoadError(error, std::strerror(errno));
      return false;
    }
    if (got == 0) {
      SetIndexLoadError(error, "unexpected end of file");
      return false;
    }
    cursor += got;
    remaining -= static_cast<size_t>(got);
    offset += static_cast<uint64_t>(got);
  }
  return true;
}

bool InspectIndexLayout(int fd, uint64_t file_bytes, IndexDiskLayout *layout,
                        std::string *error) {
  layout->file_bytes = file_bytes;
  if (!PreadIndexExact(fd, &layout->prefix, sizeof(layout->prefix), 0,
                       error)) {
    SetIndexLoadError(error, "cannot read Chromap index prefix: " + *error);
    return false;
  }
  const IndexPrefixV1 &prefix = layout->prefix;
  if (prefix.kmer_size <= 0 || prefix.window_size <= 0 ||
      prefix.n_buckets == 0 ||
      (prefix.n_buckets & (prefix.n_buckets - 1)) != 0 ||
      prefix.declared_lookup_size != prefix.size ||
      prefix.size > prefix.n_occupied ||
      prefix.n_occupied > prefix.n_buckets ||
      prefix.upper_bound > prefix.n_buckets) {
    SetIndexLoadError(error, "invalid Chromap index dimensions");
    return false;
  }

  const uint64_t flag_words =
      prefix.n_buckets < 16 ? 1 : prefix.n_buckets >> 4;
  layout->flags_offset = sizeof(IndexPrefixV1);
  if (!MultiplyIndexU64(flag_words, sizeof(uint32_t),
                        &layout->flags_bytes) ||
      !AddIndexU64(layout->flags_offset, layout->flags_bytes,
                   &layout->keys_offset) ||
      !MultiplyIndexU64(prefix.n_buckets, sizeof(uint64_t),
                        &layout->keys_bytes) ||
      !AddIndexU64(layout->keys_offset, layout->keys_bytes,
                   &layout->values_offset)) {
    SetIndexLoadError(error, "Chromap index khash layout overflows");
    return false;
  }
  layout->values_bytes = layout->keys_bytes;
  if (!AddIndexU64(layout->values_offset, layout->values_bytes,
                   &layout->occurrence_count_offset) ||
      !AddIndexU64(layout->occurrence_count_offset, sizeof(uint32_t),
                   &layout->occurrences_offset) ||
      layout->occurrence_count_offset + sizeof(uint32_t) > file_bytes ||
      !PreadIndexExact(fd, &layout->occurrence_count,
                       sizeof(layout->occurrence_count),
                       layout->occurrence_count_offset, error) ||
      !MultiplyIndexU64(layout->occurrence_count, sizeof(uint64_t),
                        &layout->occurrences_bytes) ||
      !AddIndexU64(layout->occurrences_offset, layout->occurrences_bytes,
                   &layout->payload_end)) {
    SetIndexLoadError(error, "invalid Chromap occurrence-table layout");
    return false;
  }

  if (file_bytes == layout->payload_end) return true;
  uint64_t footer_end = 0;
  if (!AddIndexU64(layout->payload_end, sizeof(IndexReferenceFooterV1),
                   &footer_end) ||
      file_bytes != footer_end ||
      !PreadIndexExact(fd, &layout->reference_footer,
                       sizeof(layout->reference_footer), layout->payload_end,
                       error)) {
    SetIndexLoadError(error, "invalid trailing Chromap index data");
    return false;
  }
  const IndexReferenceFooterV1 &footer = layout->reference_footer;
  if (memcmp(footer.magic, kIndexReferenceFooterMagic,
             sizeof(footer.magic)) != 0 ||
      footer.version != kIndexReferenceFooterVersion ||
      footer.footer_bytes != sizeof(footer) ||
      footer.reference_fingerprint == 0 || footer.total_bases == 0 ||
      footer.sequence_count == 0) {
    SetIndexLoadError(error, "invalid Chromap index reference footer");
    return false;
  }
  layout->has_reference_footer = true;
  return true;
}

void CopyIndexBlock(const char *block, uint64_t block_offset,
                    size_t block_bytes,
                    const std::vector<IndexLoadSpan> &spans) {
  const uint64_t block_end = block_offset + block_bytes;
  for (size_t i = 0; i < spans.size(); ++i) {
    const IndexLoadSpan &span = spans[i];
    const uint64_t span_end = span.file_offset + span.bytes;
    const uint64_t copy_begin = std::max(block_offset, span.file_offset);
    const uint64_t copy_end = std::min(block_end, span_end);
    if (copy_begin < copy_end) {
      memcpy(span.destination + (copy_begin - span.file_offset),
             block + (copy_begin - block_offset), copy_end - copy_begin);
    }
  }
}

bool LoadIndexSpansDirect(const std::string &path, int buffered_fd,
                          uint64_t file_bytes,
                          const std::vector<IndexLoadSpan> &spans,
                          std::string *error) {
#ifdef O_DIRECT
  const int direct_fd = open(path.c_str(), O_RDONLY | O_DIRECT);
  if (direct_fd < 0) {
    SetIndexLoadError(error, "O_DIRECT open failed: " +
                                 std::string(std::strerror(errno)));
    return false;
  }
  void *raw_buffer = nullptr;
  if (posix_memalign(&raw_buffer, kIndexDirectAlignment,
                     kIndexDirectBlockBytes) != 0) {
    close(direct_fd);
    SetIndexLoadError(error, "cannot allocate aligned direct-I/O buffer");
    return false;
  }
  char *buffer = static_cast<char *>(raw_buffer);
  const uint64_t aligned_bytes =
      file_bytes - file_bytes % kIndexDirectAlignment;
  bool ok = true;
  uint64_t offset = 0;
  while (offset < aligned_bytes) {
    const size_t bytes = static_cast<size_t>(std::min<uint64_t>(
        kIndexDirectBlockBytes, aligned_bytes - offset));
    if (!PreadIndexExact(direct_fd, buffer, bytes, offset, error)) {
      ok = false;
      break;
    }
    CopyIndexBlock(buffer, offset, bytes, spans);
    offset += bytes;
  }
  if (ok && aligned_bytes < file_bytes) {
    const size_t tail_bytes = static_cast<size_t>(file_bytes - aligned_bytes);
    if (!PreadIndexExact(buffered_fd, buffer, tail_bytes, aligned_bytes,
                         error)) {
      ok = false;
    } else {
      CopyIndexBlock(buffer, aligned_bytes, tail_bytes, spans);
    }
  }
  free(raw_buffer);
  close(direct_fd);
  if (!ok) {
    SetIndexLoadError(error, "direct-I/O read failed: " + *error);
  }
  return ok;
#else
  (void)path;
  (void)buffered_fd;
  (void)file_bytes;
  (void)spans;
  SetIndexLoadError(error, "O_DIRECT is unavailable on this platform");
  return false;
#endif
}

bool LoadIndexSpansBuffered(int fd, const std::vector<IndexLoadSpan> &spans,
                            std::string *error) {
  for (size_t i = 0; i < spans.size(); ++i) {
    if (!PreadIndexExact(fd, spans[i].destination,
                         static_cast<size_t>(spans[i].bytes),
                         spans[i].file_offset, error)) {
      return false;
    }
  }
  return true;
}

}  // namespace

void Index::Construct(uint32_t num_sequences, const SequenceBatch &reference) {
  const double real_start_time = GetRealTime();

  std::vector<Minimizer> minimizers;
  minimizers.reserve(reference.GetNumBases() / window_size_ * 2);
  std::cerr << "Collecting minimizers.\n";
  MinimizerGenerator minimizer_generator(kmer_size_, window_size_);
  for (uint32_t sequence_index = 0; sequence_index < num_sequences;
       ++sequence_index) {
    minimizer_generator.GenerateMinimizers(reference, sequence_index,
                                           minimizers);
  }
  std::cerr << "Collected " << minimizers.size() << " minimizers.\n";
  std::cerr << "Sorting minimizers.\n";
  std::stable_sort(minimizers.begin(), minimizers.end());
  std::cerr << "Sorted all minimizers.\n";
  const size_t num_minimizers = minimizers.size();
  assert(num_minimizers > 0);
  // TODO: check this assert!
  // Here I make sure the # minimizers is less than the limit of signed int32,
  // so that I can use int to store position later.
  assert(num_minimizers <= static_cast<size_t>(INT_MAX));

  occurrence_table_.reserve(num_minimizers);
  uint64_t previous_lookup_hash =
      GenerateHashInLookupTable(minimizers[0].GetHash());
  uint32_t num_previous_minimizer_occurrences = 0;
  uint64_t num_nonsingletons = 0;
  uint32_t num_singletons = 0;
  for (size_t mi = 0; mi <= num_minimizers; ++mi) {
    const bool is_last_iteration = mi == num_minimizers;
    const uint64_t current_lookup_hash =
        is_last_iteration ? previous_lookup_hash + 1
                          : GenerateHashInLookupTable(minimizers[mi].GetHash());

    if (current_lookup_hash != previous_lookup_hash) {
      int khash_return_code = 0;
      khiter_t khash_iterator =
          kh_put(k64, lookup_table_, previous_lookup_hash, &khash_return_code);
      assert(khash_return_code != -1 && khash_return_code != 0);

      if (num_previous_minimizer_occurrences == 1) {
        // We set the lowest bit of the key value to 1 if the minimizer only
        // occurs once. And the occurrence is directly saved in the lookup
        // table.
        kh_key(lookup_table_, khash_iterator) |= 1;
        kh_value(lookup_table_, khash_iterator) = occurrence_table_.back();
        occurrence_table_.pop_back();
        ++num_singletons;
      } else {
        kh_value(lookup_table_, khash_iterator) =
            GenerateEntryValueInLookupTable(num_nonsingletons,
                                            num_previous_minimizer_occurrences);
        num_nonsingletons += num_previous_minimizer_occurrences;
      }
      num_previous_minimizer_occurrences = 1;
    } else {
      num_previous_minimizer_occurrences++;
    }

    if (is_last_iteration) {
      break;
    }

    occurrence_table_.push_back(minimizers[mi].GetHit());
    previous_lookup_hash = current_lookup_hash;
  }
  assert(num_nonsingletons + num_singletons == num_minimizers);

  std::cerr << "Kmer size: " << kmer_size_ << ", window size: " << window_size_
            << ".\n";
  std::cerr << "Lookup table size: " << kh_size(lookup_table_)
            << ", # buckets: " << kh_n_buckets(lookup_table_)
            << ", occurrence table size: " << occurrence_table_.size()
            << ", # singletons: " << num_singletons << ".\n";
  std::cerr << "Built index successfully in " << GetRealTime() - real_start_time
            << "s.\n";
}

void Index::Save(const MaterializedReferenceInfo *reference_info) const {
  const double real_start_time = GetRealTime();
  FILE *index_file = fopen(index_file_path_.c_str(), "wb");
  assert(index_file != nullptr);

  uint64_t num_bytes = 0;
  int err = 0;

  err = fwrite(&kmer_size_, sizeof(int), 1, index_file);
  num_bytes += sizeof(int);
  assert(err != 0);

  err = fwrite(&window_size_, sizeof(int), 1, index_file);
  num_bytes += sizeof(int);
  assert(err != 0);

  const uint32_t lookup_table_size = kh_size(lookup_table_);
  err = fwrite(&lookup_table_size, sizeof(uint32_t), 1, index_file);
  num_bytes += sizeof(uint32_t);
  assert(err != 0);

  kh_save(k64, lookup_table_, index_file);
  num_bytes += sizeof(uint64_t) * 2 * lookup_table_size;

  const uint32_t occurrence_table_size = occurrence_table_.size();
  err = fwrite(&occurrence_table_size, sizeof(uint32_t), 1, index_file);
  num_bytes += sizeof(uint32_t);
  assert(err != 0);

  if (occurrence_table_size > 0) {
    err = fwrite(occurrence_table_.data(), sizeof(uint64_t),
                 occurrence_table_size, index_file);
    num_bytes += sizeof(uint64_t) * occurrence_table_size;
    assert(err != 0);
  }

  if (reference_info != nullptr) {
    IndexReferenceFooterV1 footer = {};
    memcpy(footer.magic, kIndexReferenceFooterMagic, sizeof(footer.magic));
    footer.version = kIndexReferenceFooterVersion;
    footer.footer_bytes = sizeof(footer);
    footer.reference_fingerprint = reference_info->fingerprint;
    footer.total_bases = reference_info->num_bases;
    footer.sequence_count = reference_info->num_sequences;
    err = fwrite(&footer, sizeof(footer), 1, index_file);
    assert(err != 0);
  }

  fclose(index_file);
  // std::cerr << "Index size: " << num_bytes / (1024.0 * 1024 * 1024) << "GB,
  std::cerr << "Saved in " << GetRealTime() - real_start_time << "s.\n";
}

void Index::Load() {
  const double real_start_time = GetRealTime();
  const int buffered_fd = open(index_file_path_.c_str(), O_RDONLY);
  if (buffered_fd < 0) {
    ExitWithMessage("Cannot open Chromap index " + index_file_path_ + ": " +
                    std::strerror(errno));
  }
  struct stat file_stat = {};
  if (fstat(buffered_fd, &file_stat) != 0 || file_stat.st_size < 0) {
    close(buffered_fd);
    ExitWithMessage("Cannot stat Chromap index " + index_file_path_);
  }
  IndexDiskLayout layout;
  std::string load_error;
  if (!InspectIndexLayout(buffered_fd,
                          static_cast<uint64_t>(file_stat.st_size), &layout,
                          &load_error)) {
    close(buffered_fd);
    ExitWithMessage("Cannot inspect Chromap index " + index_file_path_ +
                    ": " + load_error);
  }

  kmer_size_ = layout.prefix.kmer_size;
  window_size_ = layout.prefix.window_size;
  lookup_table_->n_buckets = layout.prefix.n_buckets;
  lookup_table_->size = layout.prefix.size;
  lookup_table_->n_occupied = layout.prefix.n_occupied;
  lookup_table_->upper_bound = layout.prefix.upper_bound;
  lookup_table_->flags =
      static_cast<khint32_t *>(malloc(layout.flags_bytes));
  lookup_table_->keys =
      static_cast<uint64_t *>(malloc(layout.keys_bytes));
  lookup_table_->vals =
      static_cast<uint64_t *>(malloc(layout.values_bytes));
  if (lookup_table_->flags == nullptr || lookup_table_->keys == nullptr ||
      lookup_table_->vals == nullptr) {
    close(buffered_fd);
    ExitWithMessage("Cannot allocate Chromap index khash arrays");
  }
  occurrence_table_.resize(layout.occurrence_count);

  const std::vector<IndexLoadSpan> spans = {
      {layout.flags_offset, layout.flags_bytes,
       reinterpret_cast<char *>(lookup_table_->flags)},
      {layout.keys_offset, layout.keys_bytes,
       reinterpret_cast<char *>(lookup_table_->keys)},
      {layout.values_offset, layout.values_bytes,
       reinterpret_cast<char *>(lookup_table_->vals)},
      {layout.occurrences_offset, layout.occurrences_bytes,
       reinterpret_cast<char *>(occurrence_table_.data())}};

  const bool used_direct = LoadIndexSpansDirect(
      index_file_path_, buffered_fd, layout.file_bytes, spans, &load_error);
  if (!used_direct) {
    const std::string direct_error = load_error;
    if (!LoadIndexSpansBuffered(buffered_fd, spans, &load_error)) {
      close(buffered_fd);
      ExitWithMessage("Cannot load Chromap index " + index_file_path_ +
                      ": " + load_error);
    }
    std::cerr << "Index input path: buffered positioned-read fallback ("
              << direct_error << ").\n";
  } else {
    std::cerr << "Index input path: direct I/O block reader, block size "
              << kIndexDirectBlockBytes << " bytes.\n";
  }
  close(buffered_fd);

  has_materialized_reference_binding_ = false;
  materialized_reference_binding_ = MaterializedReferenceInfo();
  if (layout.has_reference_footer) {
    const IndexReferenceFooterV1 &footer = layout.reference_footer;
    has_materialized_reference_binding_ = true;
    materialized_reference_binding_.fingerprint =
        footer.reference_fingerprint;
    materialized_reference_binding_.num_bases = footer.total_bases;
    materialized_reference_binding_.num_sequences = footer.sequence_count;
  }

  std::cerr << "Kmer size: " << kmer_size_ << ", window size: " << window_size_
            << ".\n";
  std::cerr << "Lookup table size: " << kh_size(lookup_table_)
            << ", occurrence table size: " << occurrence_table_.size() << ".\n";
  std::cerr << "Loaded index successfully in "
            << GetRealTime() - real_start_time << "s.\n";
}

bool Index::ValidateMaterializedReference(
    const MaterializedReferenceInfo &reference_info,
    std::string *error) const {
  if (!has_materialized_reference_binding_) {
    if (error != nullptr) {
      *error = "Chromap index does not declare a materialized reference; "
               "rebuild it with --reference-sidecar";
    }
    return false;
  }
  if (reference_info.fingerprint !=
          materialized_reference_binding_.fingerprint ||
      reference_info.num_bases != materialized_reference_binding_.num_bases ||
      reference_info.num_sequences !=
          materialized_reference_binding_.num_sequences) {
    if (error != nullptr) {
      *error = "materialized reference does not match the Chromap index";
    }
    return false;
  }
  return true;
}

void Index::Statistics(uint32_t num_sequences,
                       const SequenceBatch &reference) const {
  double real_start_time = GetRealTime();
  int n = 0, n1 = 0;
  uint32_t i;
  uint64_t sum = 0, len = 0;
  fprintf(stderr, "[M::%s] kmer size: %d; skip: %d; #seq: %d\n", __func__,
          kmer_size_, window_size_, num_sequences);
  for (i = 0; i < num_sequences; ++i) {
    len += reference.GetSequenceLengthAt(i);
  }
  assert(len == reference.GetNumBases());
  if (lookup_table_) {
    n += kh_size(lookup_table_);
  }
  for (khint_t k = 0; k < kh_end(lookup_table_); ++k) {
    if (kh_exist(lookup_table_, k)) {
      sum +=
          kh_key(lookup_table_, k) & 1 ? 1 : (uint32_t)kh_val(lookup_table_, k);
      if (kh_key(lookup_table_, k) & 1) ++n1;
    }
  }
  fprintf(stderr,
          "[M::%s::%.3f] distinct minimizers: %d (%.2f%% are singletons); "
          "average occurrences: %.3lf; average spacing: %.3lf\n",
          __func__, GetRealTime() - real_start_time, n, 100.0 * n1 / n,
          (double)sum / n, (double)len / sum);
}

void Index::CheckIndex(uint32_t num_sequences,
                       const SequenceBatch &reference) const {
  std::vector<Minimizer> minimizers;
  minimizers.reserve(reference.GetNumBases() / window_size_ * 2);
  MinimizerGenerator minimizer_generator(kmer_size_, window_size_);
  for (uint32_t sequence_index = 0; sequence_index < num_sequences;
       ++sequence_index) {
    minimizer_generator.GenerateMinimizers(reference, sequence_index,
                                           minimizers);
  }
  std::cerr << "Collected " << minimizers.size() << " minimizers.\n";
  std::stable_sort(minimizers.begin(), minimizers.end());
  std::cerr << "Sorted minimizers.\n";

  uint32_t count = 0;
  for (uint32_t i = 0; i < minimizers.size(); ++i) {
    khiter_t khash_iterator = kh_get(
        k64, lookup_table_, GenerateHashInLookupTable(minimizers[i].GetHash()));
    assert(khash_iterator != kh_end(lookup_table_));
    uint64_t key = kh_key(lookup_table_, khash_iterator);
    uint64_t value = kh_value(lookup_table_, khash_iterator);
    if (key & 1) {  // singleton
      assert(minimizers[i].GetHit() == value);
      count = 0;
    } else {
      uint32_t offset = GenerateOffsetInOccurrenceTable(value);
      uint32_t num_occ = GenerateNumOccurrenceInOccurrenceTable(value);
      uint64_t value_in_index = occurrence_table_[offset + count];
      assert(value_in_index == minimizers[i].GetHit());
      ++count;
      if (count == num_occ) {
        count = 0;
      }
    }
  }
}

int Index::GenerateCandidatePositions(
    const CandidatePositionGeneratingConfig &generating_config,
    MappingMetadata &mapping_metadata) const {
  const uint32_t num_minimizers = mapping_metadata.GetNumMinimizers();
  const std::vector<Minimizer> &minimizers = mapping_metadata.minimizers_;

  std::vector<std::vector<uint64_t>> positive_candidate_position_lists;
  std::vector<std::vector<uint64_t>> negative_candidate_position_lists;
  if (generating_config.UseHeapMerge()) {
    for (uint32_t i = 0; i < num_minimizers; ++i) {
      positive_candidate_position_lists.emplace_back(std::vector<uint64_t>());
      negative_candidate_position_lists.emplace_back(std::vector<uint64_t>());
    }
  }
  bool is_candidate_position_list_sorted = true;

  mapping_metadata.positive_hits_.reserve(
      generating_config.GetMaxSeedFrequency() * 2);
  mapping_metadata.negative_hits_.reserve(
      generating_config.GetMaxSeedFrequency() * 2);

  RepetitiveSeedStats repetitive_seed_stats;
  for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
    khiter_t khash_iterator =
        kh_get(k64, lookup_table_,
               GenerateHashInLookupTable(minimizers[mi].GetHash()));
    if (khash_iterator == kh_end(lookup_table_)) {
      // std::cerr << "The minimizer is not in reference!\n";
      continue;
    }

    std::vector<uint64_t> &positive_candidate_positions =
        generating_config.UseHeapMerge() ? positive_candidate_position_lists[mi]
                                         : mapping_metadata.positive_hits_;
    std::vector<uint64_t> &negative_candidate_positions =
        generating_config.UseHeapMerge() ? negative_candidate_position_lists[mi]
                                         : mapping_metadata.negative_hits_;

    const uint64_t lookup_key = kh_key(lookup_table_, khash_iterator);
    const uint64_t lookup_value = kh_value(lookup_table_, khash_iterator);
    const uint64_t read_hit = minimizers[mi].GetHit();
    if (IsSingletonLookupKey(lookup_key)) {
      const uint64_t candidate_position = GenerateCandidatePositionFromHits(
          /*reference_hit=*/lookup_value, read_hit);
      if (AreTwoHitsOnTheSameStrand(/*reference_hit=*/lookup_value, read_hit)) {
        positive_candidate_positions.push_back(candidate_position);
      } else {
        negative_candidate_positions.push_back(candidate_position);
      }
      continue;
    }

    const uint32_t num_occurrences =
        GenerateNumOccurrenceInOccurrenceTable(lookup_value);
    if (!generating_config.IsFrequentSeed(num_occurrences)) {
      const uint32_t read_position = HitToSequencePosition(read_hit);
      const uint32_t occ_offset = GenerateOffsetInOccurrenceTable(lookup_value);
      for (uint32_t oi = 0; oi < num_occurrences; ++oi) {
        const uint64_t reference_hit = occurrence_table_[occ_offset + oi];
        const uint64_t candidate_position =
            GenerateCandidatePositionFromHits(reference_hit, read_hit);
        if (AreTwoHitsOnTheSameStrand(reference_hit, read_hit)) {
          const uint32_t reference_position =
              HitToSequencePosition(reference_hit);
          if (reference_position < read_position) {
            is_candidate_position_list_sorted = false;
          }
          positive_candidate_positions.push_back(candidate_position);
        } else {
          negative_candidate_positions.push_back(candidate_position);
        }
      }
    }

    if (generating_config.IsRepetitiveSeed(num_occurrences)) {
      const uint32_t read_position = HitToSequencePosition(read_hit);
      UpdateRepetitiveSeedStats(read_position, repetitive_seed_stats);
    }
  }

  if (generating_config.UseHeapMerge()) {
    // TODO: try to remove this sorting.
    if (!is_candidate_position_list_sorted) {
      for (uint32_t mi = 0; mi < num_minimizers; ++mi) {
        std::sort(positive_candidate_position_lists[mi].begin(),
                  positive_candidate_position_lists[mi].end());
      }
    }
    HeapMergeCandidatePositionLists(positive_candidate_position_lists,
                                    mapping_metadata.positive_hits_);
    HeapMergeCandidatePositionLists(negative_candidate_position_lists,
                                    mapping_metadata.negative_hits_);
  } else {
    std::sort(mapping_metadata.positive_hits_.begin(),
              mapping_metadata.positive_hits_.end());
    std::sort(mapping_metadata.negative_hits_.begin(),
              mapping_metadata.negative_hits_.end());
  }

#ifdef LI_DEBUG
  for (uint32_t mi = 0; mi < mapping_metadata.positive_hits_.size(); ++mi)
    printf("+ %llu %d %d\n", mapping_metadata.positive_hits_[mi],
           (int)(mapping_metadata.positive_hits_[mi] >> 32), (int)(mapping_metadata.positive_hits_[mi]));

  for (uint32_t mi = 0; mi < mapping_metadata.negative_hits_.size(); ++mi)
    printf("- %llu %d %d\n", mapping_metadata.negative_hits_[mi],
           (int)(mapping_metadata.negative_hits_[mi] >> 32), (int)(mapping_metadata.negative_hits_[mi]));
#endif

  mapping_metadata.repetitive_seed_length_ =
      repetitive_seed_stats.repetitive_seed_length;
  return repetitive_seed_stats.repetitive_seed_count;
}

int Index::GenerateCandidatePositionsFromRepetitiveReadWithMateInfoOnOneStrand(
    const Strand strand, uint32_t search_range,
    int min_num_seeds_required_for_mapping, int max_seed_frequency0,
    int error_threshold, const std::vector<Minimizer> &minimizers,
    const std::vector<Candidate> &mate_candidates,
    uint32_t &repetitive_seed_length,
    std::vector<uint64_t> &candidate_positions) const {
  const uint32_t mate_candidates_size = mate_candidates.size();
  int max_minimizer_count = 0;
  int best_candidate_num = 0;
  for (uint32_t i = 0; i < mate_candidates_size; ++i) {
    int count = mate_candidates[i].count;
    if (count > max_minimizer_count) {
      max_minimizer_count = count;
      best_candidate_num = 1;
    } else if (count == max_minimizer_count) {
      ++best_candidate_num;
    }
  }

  const bool mate_has_too_many_candidates =
      best_candidate_num >= 300 ||
      mate_candidates_size > static_cast<uint32_t>(max_seed_frequency0);
  const bool mate_has_too_many_low_support_candidates =
      max_minimizer_count <= min_num_seeds_required_for_mapping &&
      best_candidate_num >= 200;
  if (mate_has_too_many_candidates ||
      mate_has_too_many_low_support_candidates) {
    return -max_minimizer_count;
  }

  // TODO: reduce the search range based on the strand.
  std::vector<std::pair<uint64_t, uint64_t>> boundaries;
  boundaries.reserve(best_candidate_num);
  for (uint32_t ci = 0; ci < mate_candidates_size; ++ci) {
    if (mate_candidates[ci].count == max_minimizer_count) {
      const uint64_t boundary_start =
          (mate_candidates[ci].position < search_range)
              ? 0
              : (mate_candidates[ci].position - search_range);
      const uint64_t boundary_end = mate_candidates[ci].position + search_range;
      boundaries.emplace_back(boundary_start, boundary_end);
    }
  }

  const uint32_t raw_boundary_size = boundaries.size();
  if (raw_boundary_size == 0) {
    return max_minimizer_count;
  }

  // Merge adjacent boundary point. Assume the candidates are sorted by
  // coordinate, and thus boundaries are also sorted.
  uint32_t boundary_size = 1;
  for (uint32_t bi = 1; bi < raw_boundary_size; ++bi) {
    if (boundaries[boundary_size - 1].second < boundaries[bi].first) {
      boundaries[boundary_size] = boundaries[bi];
      ++boundary_size;
    } else {
      boundaries[boundary_size - 1].second = boundaries[bi].second;
    }
  }
  boundaries.resize(boundary_size);

  RepetitiveSeedStats repetitive_seed_stats;
  for (uint32_t mi = 0; mi < minimizers.size(); ++mi) {
    khiter_t khash_iterator =
        kh_get(k64, lookup_table_,
               GenerateHashInLookupTable(minimizers[mi].GetHash()));
    if (khash_iterator == kh_end(lookup_table_)) {
      // std::cerr << "The minimizer is not in reference!\n";
      continue;
    }

    const uint64_t lookup_key = kh_key(lookup_table_, khash_iterator);
    const uint64_t lookup_value = kh_value(lookup_table_, khash_iterator);
    const uint64_t read_hit = minimizers[mi].GetHit();
    const uint32_t read_position = HitToSequencePosition(read_hit);
    if (IsSingletonLookupKey(lookup_key)) {
      const uint64_t candidate_position =
          GenerateCandidatePositionFromHits(lookup_value, read_hit);
      const bool on_same_strand =
          AreTwoHitsOnTheSameStrand(lookup_value, read_hit);
      if ((on_same_strand && strand == kPositive) ||
          (!on_same_strand && strand == kNegative)) {
        candidate_positions.push_back(candidate_position);
      }
      continue;
    }

    const uint32_t offset = GenerateOffsetInOccurrenceTable(lookup_value);
    const uint32_t num_occurrences =
        GenerateNumOccurrenceInOccurrenceTable(lookup_value);
    int32_t prev_l = 0;
    for (uint32_t bi = 0; bi < boundary_size; ++bi) {
      // Use binary search to locate the coordinate near mate position.
      int32_t l = prev_l, m = 0, r = num_occurrences - 1;
      uint64_t boundary = boundaries[bi].first;
      while (l <= r) {
        m = (l + r) / 2;
        uint64_t candidate_position =
            GenerateCandidatePositionFromOccurrenceTableEntry(
                occurrence_table_[offset + m]);
        if (candidate_position < boundary) {
          l = m + 1;
        } else if (candidate_position > boundary) {
          r = m - 1;
        } else {
          break;
        }
      }
      // For next boundary, we don't have to start from l=0.
      prev_l = m;

      for (uint32_t oi = m; oi < num_occurrences; ++oi) {
        const uint64_t reference_hit = occurrence_table_[offset + oi];
        if ((GenerateCandidatePositionFromOccurrenceTableEntry(reference_hit)) >
            boundaries[bi].second) {
          break;
        }
        const uint64_t candidate_position =
            GenerateCandidatePositionFromHits(reference_hit, read_hit);
        const bool on_same_strand =
            AreTwoHitsOnTheSameStrand(reference_hit, read_hit);
        if ((on_same_strand && strand == kPositive) ||
            (!on_same_strand && strand == kNegative)) {
          candidate_positions.push_back(candidate_position);
        }
      }
    }

    if (num_occurrences >= (uint32_t)max_seed_frequency0) {
      UpdateRepetitiveSeedStats(read_position, repetitive_seed_stats);
    }
  }

  std::sort(candidate_positions.begin(), candidate_positions.end());
  repetitive_seed_length = repetitive_seed_stats.repetitive_seed_length;
  return max_minimizer_count;
}

uint64_t Index::GenerateCandidatePositionFromHits(uint64_t reference_hit,
                                                  uint64_t read_hit) const {
  const uint32_t reference_position = HitToSequencePosition(reference_hit);
  const uint32_t read_position = HitToSequencePosition(read_hit);
  // For now we can't see the reference here. So let us don't validate this
  // candidate position. Instead, we do it later some time when we check the
  // candidates.
  const uint32_t reference_start_position =
      AreTwoHitsOnTheSameStrand(reference_hit, read_hit)
          ? reference_position - read_position
          : reference_position + read_position - kmer_size_ + 1;
  const uint64_t reference_id = HitToSequenceIndex(reference_hit);
  return SequenceIndexAndPositionToCandidatePosition(reference_id,
                                                     reference_start_position);
}

void Index::UpdateRepetitiveSeedStats(uint32_t read_position,
                                      RepetitiveSeedStats &stats) const {
  if (stats.previous_repetitive_seed_position > read_position) {
    // First minimizer.
    stats.repetitive_seed_length += kmer_size_;
  } else {
    if (read_position < stats.previous_repetitive_seed_position + kmer_size_ +
                            window_size_ - 1) {
      stats.repetitive_seed_length +=
          read_position - stats.previous_repetitive_seed_position;
    } else {
      stats.repetitive_seed_length += kmer_size_;
    }
  }
  stats.previous_repetitive_seed_position = read_position;
  ++stats.repetitive_seed_count;
}

}  // namespace chromap
