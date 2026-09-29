// Parallel low-memory finalisation of paired-end ATAC spills: every variant
// (serial merge, per-reference tasks on 1-64 threads, host permits, open-file
// limits) must write the same bytes as the serial merge.
//
// The harness drives MappingWriter<AtacSpillRecord> through its public
// low-memory API with synthetic spill runs; no FASTQ, index or fixture is
// needed. Each (case, variant) runs in a forked child that writes its outputs
// (sidecar or text, summary CSV, MACS3 bucket dump, stderr counter lines) to
// its own directory; the parent compares them with the serial variant.
//
// Environment:
//   CHROMAP_ARTIFACT_ROOT       where the per-run directories are created
//   LOWMEM_TEST_GOLDENS_DIR     compare serial outputs with saved ones
//   LOWMEM_TEST_WRITE_GOLDENS   with the above: save them instead
//
// Built with -DLOWMEM_FINALIZE_HAS_THREADS against a library that has
// MappingParameters::low_mem_finalize_threads; without it (for saving
// goldens from an earlier library) only the serial variant runs.
#include <dirent.h>
#include <fcntl.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <condition_variable>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "atac_kway_spill.h"
#include "atac_spill_record.h"
#include "khash.h"
#include "mapping_parameters.h"
#include "mapping_writer.h"
#include "sequence_batch.h"
#include "sequence_effective_range.h"
#include "summary_metadata.h"
#include "utils.h"

using namespace chromap;

namespace {

const uint32_t kBarcodeLength = 16;

// ------------------------------------------------------------------ files

std::string ReadFile(const std::string &path) {
  std::ifstream in(path.c_str(), std::ios::binary);
  std::ostringstream out;
  out << in.rdbuf();
  return out.str();
}

bool FileExists(const std::string &path) {
  struct stat st;
  return stat(path.c_str(), &st) == 0;
}

void WriteFile(const std::string &path, const std::string &data) {
  std::ofstream out(path.c_str(), std::ios::binary);
  out << data;
}

std::vector<std::string> ListDir(const std::string &dir) {
  std::vector<std::string> names;
  DIR *d = opendir(dir.c_str());
  if (d == nullptr) return names;
  while (struct dirent *entry = readdir(d)) {
    const std::string name = entry->d_name;
    if (name != "." && name != "..") names.push_back(name);
  }
  closedir(d);
  std::sort(names.begin(), names.end());
  return names;
}

void RemoveTree(const std::string &dir) {
  for (const std::string &name : ListDir(dir)) {
    const std::string path = dir + "/" + name;
    struct stat st;
    if (lstat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
      chmod(path.c_str(), 0755);
      RemoveTree(path);
    } else {
      unlink(path.c_str());
    }
  }
  rmdir(dir.c_str());
}

void MakeDirs(const std::string &path) {
  std::string partial;
  std::istringstream parts(path);
  std::string part;
  if (!path.empty() && path[0] == '/') partial = "/";
  while (std::getline(parts, part, '/')) {
    if (part.empty()) continue;
    partial += part;
    mkdir(partial.c_str(), 0755);
    partial += "/";
  }
}

std::string SeedToSequence(uint64_t seed, uint32_t length) {
  static const char kBases[] = "ACGT";
  std::string seq;
  for (uint32_t i = 0; i < length; ++i) {
    seq.push_back(kBases[(seed >> ((length - 1 - i) * 2)) & 3]);
  }
  return seq;
}

// ------------------------------------------------------------------ cases

enum class OutputMode { kSidecar, kBed, kTagAlign };

struct CaseSpec {
  std::string name;
  uint32_t num_references = 1;
  OutputMode output = OutputMode::kSidecar;
  bool remove_pcr_duplicates = true;
  bool bulk_level_dedup = false;
  bool tn5_shift = true;
  bool summary = true;
  bool macs3_buffer = true;
  bool translate_table = false;
  std::vector<std::vector<std::vector<AtacSpillRecord>>> flushes;
  std::vector<std::pair<uint64_t, uint64_t>> whitelist;
  std::vector<uint64_t> summary_preseed;  // barcodes credited TOTAL first
  // Applied in the child after the spills are written.
  std::function<void(const std::string &tmp_dir)> before_merge;
};

struct Variant {
  std::string name;
  int finalize_threads = 1;   // 1 = the serial merge
  int permits = 0;            // >0: fake host permit hooks with this pool
  int rlimit_nofile = 0;      // >0: soft and hard open-file limit
  bool expect_failure = false;
  std::string expected_message;
  bool compare_with_serial = true;
};

AtacSpillRecord Rec(uint64_t read_id, uint64_t barcode, uint32_t start,
                    uint16_t length, uint8_t mapq, uint8_t direction = 1,
                    uint8_t is_unique = 1) {
  PairedEndMappingWithBarcode bed(read_id, barcode, start, length, mapq,
                                  direction, is_unique, /*num_dups=*/1,
                                  /*positive_alignment_length=*/50,
                                  /*negative_alignment_length=*/50);
  return AtacSpillRecord(bed);
}

uint64_t Draw(std::mt19937_64 &rng, uint64_t n) { return rng() % n; }

// Random records spread over flushes: small barcode and position pools make
// duplicates (cell and bulk level) and cross-file ties common.
void AddRandomFlushes(CaseSpec *spec, uint64_t seed, int flushes,
                      int records_per_flush, int barcodes, int positions,
                      uint64_t *read_id) {
  std::mt19937_64 rng(seed);
  for (int f = 0; f < flushes; ++f) {
    std::vector<std::vector<AtacSpillRecord>> flush(spec->num_references);
    for (int i = 0; i < records_per_flush; ++i) {
      const uint32_t rid =
          static_cast<uint32_t>(Draw(rng, spec->num_references));
      const uint64_t barcode = 0x1000 + Draw(rng, barcodes);
      const uint32_t start = 1000 + 37 * static_cast<uint32_t>(Draw(rng, positions));
      const uint16_t length = static_cast<uint16_t>(120 + 10 * Draw(rng, 3));
      const uint8_t mapq = static_cast<uint8_t>(20 + Draw(rng, 25));
      const uint8_t direction = static_cast<uint8_t>(Draw(rng, 2));
      const uint8_t unique = static_cast<uint8_t>(Draw(rng, 4) != 0);
      flush[rid].push_back(
          Rec((*read_id)++, barcode, start, length, mapq, direction, unique));
    }
    spec->flushes.push_back(flush);
  }
}

// Barcodes credited TOTAL before the merge, as mapping does for whitelisted
// reads: every barcode of the pool except every fifth, so both existing and
// new summary rows occur.
std::vector<uint64_t> MostBarcodes(int barcodes) {
  std::vector<uint64_t> list;
  for (int b = 0; b < barcodes; ++b) {
    if (b % 5 != 4) list.push_back(0x1000 + b);
  }
  return list;
}

std::vector<std::pair<uint64_t, uint64_t>> WhitelistFor(int barcodes,
                                                        uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::vector<std::pair<uint64_t, uint64_t>> list;
  for (int b = 0; b < barcodes; ++b) {
    // Every third barcode is absent from the table (1.1.1 semantics);
    // abundances repeat so that ties occur.
    if (b % 3 == 2) continue;
    list.push_back(std::make_pair(0x1000 + b, Draw(rng, 4)));
  }
  return list;
}

// ------------------------------------------------------------------ permits

struct FakePermits {
  std::mutex mutex;
  std::condition_variable cv;
  int available = 0;
  int in_use = 0;
  int max_in_use = 0;
  uint64_t acquires = 0;
  uint64_t releases = 0;
  uint64_t work_units = 0;
};

FakePermits g_permits;

uint64_t FakeAcquire(void *) {
  std::unique_lock<std::mutex> lock(g_permits.mutex);
  g_permits.cv.wait(lock, [] { return g_permits.available > 0; });
  --g_permits.available;
  ++g_permits.in_use;
  ++g_permits.acquires;
  g_permits.max_in_use = std::max(g_permits.max_in_use, g_permits.in_use);
  return 0;
}

void FakeRelease(void *, uint64_t, uint64_t units, uint64_t, uint64_t) {
  {
    std::lock_guard<std::mutex> lock(g_permits.mutex);
    ++g_permits.available;
    --g_permits.in_use;
    ++g_permits.releases;
    g_permits.work_units += units;
  }
  g_permits.cv.notify_one();
}

// ------------------------------------------------------------------ child

int RunChild(const CaseSpec &spec, const Variant &variant,
             const std::string &dir) {
  const std::string tmp_dir = dir + "/tmp";
  MappingParameters params;
  params.low_memory_mode = true;
  params.mapping_output_format = spec.output == OutputMode::kTagAlign
                                     ? MAPPINGFORMAT_TAGALIGN
                                     : MAPPINGFORMAT_BED;
  if (spec.output == OutputMode::kSidecar) {
    params.atac_sidecar_only = true;
    params.atac_fragment_binary_output_file_path = dir + "/out.bin";
  } else {
    params.mapping_output_file_path = dir + "/out.txt";
  }
  params.temp_directory_path = tmp_dir;
  params.remove_pcr_duplicates = spec.remove_pcr_duplicates;
  params.remove_pcr_duplicates_at_bulk_level = spec.bulk_level_dedup;
  params.is_bulk_data = false;
  params.mapq_threshold = 30;
  params.Tn5_shift = spec.tn5_shift;
  params.num_threads = 4;
  if (spec.summary) params.summary_metadata_file_path = dir + "/summary.csv";
  if (spec.macs3_buffer) {
    params.macs3_frag_buffer =
        std::make_shared<std::vector<std::vector<macs3::FragmentRecord>>>();
    params.macs3_frag_chrom_names = std::make_shared<std::vector<std::string>>();
  }
  if (spec.translate_table) {
    std::ostringstream table;
    for (uint64_t b = 0x1000; b < 0x1100; ++b) {
      const std::string seq = SeedToSequence(b, kBarcodeLength);
      std::string translated(seq.rbegin(), seq.rend());
      table << seq << "\t" << translated << "\n";
    }
    WriteFile(dir + "/translate.tsv", table.str());
    params.barcode_translate_table_file_path = dir + "/translate.tsv";
    params.barcode_translate_from_first_column = true;
  }
#ifdef LOWMEM_FINALIZE_HAS_THREADS
  params.low_mem_finalize_threads = variant.finalize_threads;
  if (variant.permits > 0) {
    g_permits.available = variant.permits;
    params.permit_acquire_hook = &FakeAcquire;
    params.permit_release_hook = &FakeRelease;
    params.permit_hook_ctx = &g_permits;
  }
#endif

  SequenceBatch reference(spec.num_references, SequenceEffectiveRange());
  for (uint32_t rid = 0; rid < spec.num_references; ++rid) {
    reference.AssignLoadedReferenceMetadata(rid, "chr" + std::to_string(rid),
                                            10000000);
  }
  khash_t(k64_seq) *whitelist = kh_init(k64_seq);
  for (const auto &entry : spec.whitelist) {
    int ret = 0;
    khiter_t it = kh_put(k64_seq, whitelist, entry.first, &ret);
    kh_value(whitelist, it) = entry.second;
  }

  uint64_t spilled_records = 0;
  {
    MappingWriter<AtacSpillRecord> writer(params, kBarcodeLength,
                                          std::vector<int>());
    writer.OutputHeader(spec.num_references, reference);
    for (uint64_t barcode : spec.summary_preseed) {
      writer.UpdateSummaryMetadata(barcode, SUMMARY_METADATA_TOTAL, 1);
    }
    for (const auto &flush : spec.flushes) {
      std::vector<std::vector<AtacSpillRecord>> by_reference = flush;
      by_reference.resize(spec.num_references);
      for (auto &records : by_reference) {
        std::sort(records.begin(), records.end());
        spilled_records += records.size();
      }
      writer.OutputTempMappingsToOverflow(spec.num_references, by_reference);
      MappingWriter<AtacSpillRecord>::RotateThreadOverflowWriter();
    }
    if (spec.before_merge) spec.before_merge(tmp_dir);
    if (variant.rlimit_nofile > 0) {
      struct rlimit limit;
      limit.rlim_cur = static_cast<rlim_t>(variant.rlimit_nofile);
      limit.rlim_max = static_cast<rlim_t>(variant.rlimit_nofile);
      setrlimit(RLIMIT_NOFILE, &limit);
    }
    writer.ProcessAndOutputMappingsInLowMemoryFromOverflow(
        0, spec.num_references, reference, whitelist);
    writer.OutputSummaryMetadata();
  }
  kh_destroy(k64_seq, whitelist);

  if (spec.macs3_buffer) {
    std::ostringstream dump;
    const auto &buckets = *params.macs3_frag_buffer;
    for (size_t rid = 0; rid < buckets.size(); ++rid) {
      dump << "rid " << rid << " " << buckets[rid].size() << "\n";
      for (const auto &rec : buckets[rid]) {
        dump << rec.chrom_id << " " << rec.start << " " << rec.end << " "
             << rec.count << "\n";
      }
    }
    WriteFile(dir + "/buckets.txt", dump.str());
  }
  if (variant.permits > 0) {
    std::ostringstream stats;
    stats << "max_in_use " << g_permits.max_in_use << "\nacquires "
          << g_permits.acquires << "\nreleases " << g_permits.releases
          << "\nwork_units " << g_permits.work_units << "\nspilled "
          << spilled_records << "\n";
    WriteFile(dir + "/permits.txt", stats.str());
    const bool any_work = spilled_records > 0;
    if (g_permits.max_in_use > variant.permits ||
        g_permits.acquires != g_permits.releases ||
        (any_work ? g_permits.acquires == 0 : g_permits.acquires != 0)) {
      std::cerr << "[child] permit accounting failed\n";
      return 5;
    }
    if (variant.finalize_threads != 1 && g_permits.work_units != spilled_records) {
      std::cerr << "[child] permit work units " << g_permits.work_units
                << " != spilled records " << spilled_records << "\n";
      return 5;
    }
  }
  return 0;
}

// ------------------------------------------------------------------ parent

const char *kStderrPatterns[] = {
    "Processing ", "Low-memory overflow: mid-batch flush count:",
    "# uni-mappings:", "Number of output mappings (passed filters):"};

std::string CounterLines(const std::string &stderr_text) {
  std::istringstream lines(stderr_text);
  std::string line, out;
  while (std::getline(lines, line)) {
    for (const char *pattern : kStderrPatterns) {
      if (line.compare(0, strlen(pattern), pattern) == 0) {
        out += line + "\n";
        break;
      }
    }
  }
  return out;
}

const char *kOutputs[] = {"out.bin", "out.bin.chroms.tsv", "out.txt",
                          "summary.csv", "buckets.txt"};

bool RunVariant(const CaseSpec &spec, const Variant &variant,
                const std::string &root, std::string *dir_out) {
  const std::string dir = root + "/" + spec.name + "/" + variant.name;
  RemoveTree(dir);
  MakeDirs(dir + "/tmp");
  *dir_out = dir;
  std::cout.flush();
  std::cerr.flush();
  const pid_t pid = fork();
  if (pid == 0) {
    const int fd = open((dir + "/stderr.txt").c_str(),
                        O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd >= 0) {
      dup2(fd, STDERR_FILENO);
      close(fd);
    }
    const int code = RunChild(spec, variant, dir);
    std::cerr.flush();
    _exit(code);
  }
  int status = 0;
  waitpid(pid, &status, 0);
  const std::string child_stderr = ReadFile(dir + "/stderr.txt");
  const std::vector<std::string> leftovers = ListDir(dir + "/tmp");
  std::string problem;
  if (variant.expect_failure) {
    if (!WIFEXITED(status) || WEXITSTATUS(status) == 0) {
      problem = "expected an error exit";
    } else if (child_stderr.find(variant.expected_message) ==
               std::string::npos) {
      problem = "error without '" + variant.expected_message + "'";
    }
  } else if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
    problem = WIFSIGNALED(status)
                  ? "killed by signal " + std::to_string(WTERMSIG(status))
                  : "exit " + std::to_string(WEXITSTATUS(status));
  }
  if (problem.empty() && !leftovers.empty()) {
    problem = "left " + std::to_string(leftovers.size()) +
              " files in the temp directory (first: " + leftovers[0] + ")";
  }
  if (!problem.empty()) {
    std::cerr << "[lowmem-parallel] FAIL " << spec.name << "/" << variant.name
              << ": " << problem << " (" << dir << ")\n";
    std::istringstream lines(child_stderr);
    std::string line;
    std::vector<std::string> tail;
    while (std::getline(lines, line)) tail.push_back(line);
    for (size_t i = tail.size() > 6 ? tail.size() - 6 : 0; i < tail.size(); ++i)
      std::cerr << "    " << tail[i] << "\n";
    return false;
  }
  WriteFile(dir + "/counters.txt", CounterLines(child_stderr));
  return true;
}

bool SameOutputs(const std::string &a, const std::string &b,
                 std::string *difference) {
  std::vector<std::string> names(std::begin(kOutputs), std::end(kOutputs));
  names.push_back("counters.txt");
  for (const std::string &name : names) {
    const bool ea = FileExists(a + "/" + name);
    const bool eb = FileExists(b + "/" + name);
    if (ea != eb) {
      *difference = name + " exists in only one run";
      return false;
    }
    if (ea && ReadFile(a + "/" + name) != ReadFile(b + "/" + name)) {
      *difference = name + " differs";
      return false;
    }
  }
  return true;
}

// ------------------------------------------------------------------ case set

std::vector<CaseSpec> BuildCases() {
  std::vector<CaseSpec> cases;
  uint64_t read_id = 1;

  {  // U01
    CaseSpec c;
    c.name = "U01_one_reference_one_file";
    std::vector<std::vector<AtacSpillRecord>> flush(1);
    for (int i = 0; i < 50; ++i)
      flush[0].push_back(Rec(read_id++, 0x1000 + i, 5000 + 100 * i, 150, 40));
    c.flushes.push_back(flush);
    cases.push_back(c);
  }
  {  // U02
    CaseSpec c;
    c.name = "U02_three_references_400_flushes";
    c.num_references = 3;
    AddRandomFlushes(&c, 11, 400, 12, 8, 40, &read_id);
    c.whitelist = WhitelistFor(8, 12);
    c.summary_preseed = MostBarcodes(8);
    cases.push_back(c);
  }
  {  // U03
    CaseSpec c;
    c.name = "U03_duplicates_across_reference_boundary";
    c.num_references = 3;
    for (int f = 0; f < 3; ++f) {
      std::vector<std::vector<AtacSpillRecord>> flush(3);
      flush[0].push_back(Rec(read_id++, 0x1001, 9000, 150, 40));
      flush[1].push_back(Rec(read_id++, 0x1001, 9000, 150, 40));
      flush[1].push_back(Rec(read_id++, 0x1002, 100, 150, 40));
      flush[2].push_back(Rec(read_id++, 0x1002, 100, 150, 40));
      c.flushes.push_back(flush);
    }
    cases.push_back(c);
  }
  {  // U04
    CaseSpec c;
    c.name = "U04_identical_records_in_several_files";
    c.num_references = 2;
    for (int f = 0; f < 5; ++f) {
      std::vector<std::vector<AtacSpillRecord>> flush(2);
      flush[0].push_back(Rec(77, 0x1003, 4000, 150, 40));
      flush[1].push_back(Rec(78, 0x1004, 4000, 150, 35, 0, 0));
      c.flushes.push_back(flush);
    }
    cases.push_back(c);
  }
  {  // U05
    CaseSpec c;
    c.name = "U05_bulk_level_dedup";
    c.num_references = 4;
    c.bulk_level_dedup = true;
    AddRandomFlushes(&c, 21, 60, 30, 12, 15, &read_id);
    c.whitelist = WhitelistFor(12, 22);
    c.summary_preseed = MostBarcodes(12);
    cases.push_back(c);
  }
  // U06-U08: end-of-stream groups (bulk-level dedup). At each position the
  // last record in sort order and the bulk-selected best differ in MAPQ.
  auto end_group = [&read_id](std::vector<AtacSpillRecord> *records,
                              uint32_t start, bool last_passes) {
    // 0x1010 has two duplicates (the best); 0x1011 sorts after it.
    const uint8_t best_mapq = last_passes ? 25 : 40;
    const uint8_t last_mapq = last_passes ? 40 : 25;
    records->push_back(Rec(read_id++, 0x1010, start, 150, best_mapq));
    records->push_back(Rec(read_id++, 0x1010, start, 150, best_mapq));
    records->push_back(Rec(read_id++, 0x1011, start, 150, last_mapq));
  };
  for (int pattern = 0; pattern < 2; ++pattern) {  // U06 (A), U07 (B)
    CaseSpec c;
    c.name = pattern == 0 ? "U06_end_of_stream_last_passes"
                          : "U07_end_of_stream_best_passes";
    c.num_references = 2;
    c.bulk_level_dedup = true;
    c.whitelist = {{0x1010, 5}, {0x1011, 1}};
    std::vector<std::vector<AtacSpillRecord>> flush(2);
    flush[0].push_back(Rec(read_id++, 0x1012, 100, 150, 40));
    end_group(&flush[0], 500, pattern == 0);
    flush[1].push_back(Rec(read_id++, 0x1012, 100, 150, 40));
    end_group(&flush[1], 700, pattern == 0);
    c.flushes.push_back(flush);
    cases.push_back(c);
  }
  {  // U08
    CaseSpec c;
    c.name = "U08_end_groups_before_empty_references";
    c.num_references = 5;
    c.bulk_level_dedup = true;
    c.whitelist = {{0x1010, 5}, {0x1011, 1}};
    std::vector<std::vector<AtacSpillRecord>> flush(5);
    end_group(&flush[0], 300, true);
    end_group(&flush[1], 300, false);
    end_group(&flush[2], 300, true);  // references 3 and 4 stay empty
    c.flushes.push_back(flush);
    cases.push_back(c);
  }
  {  // U09
    CaseSpec c;
    c.name = "U09_empty_references";
    c.num_references = 6;
    for (int f = 0; f < 4; ++f) {
      std::vector<std::vector<AtacSpillRecord>> flush(6);
      flush[1].push_back(Rec(read_id++, 0x1001 + f, 200 + f, 150, 40));
      flush[3].push_back(Rec(read_id++, 0x1001, 200, 150, 40));
      c.flushes.push_back(flush);
    }
    cases.push_back(c);
  }
  {  // U09b
    CaseSpec c;
    c.name = "U09b_no_records";
    c.num_references = 3;
    cases.push_back(c);
  }
  {  // U10
    CaseSpec c;
    c.name = "U10_300_duplicates";
    for (int f = 0; f < 3; ++f) {
      std::vector<std::vector<AtacSpillRecord>> flush(1);
      for (int i = 0; i < 100; ++i)
        flush[0].push_back(Rec(read_id++, 0x1005, 3000, 150, 40));
      c.flushes.push_back(flush);
    }
    cases.push_back(c);
  }
  {  // U11: pending resize at the start (788 keys fill 1,024 buckets)
    CaseSpec c;
    c.name = "U11_summary_pending_resize";
    c.num_references = 3;
    AddRandomFlushes(&c, 31, 20, 40, 60, 50, &read_id);
    for (uint64_t b = 0; b < 788; ++b) c.summary_preseed.push_back(0x20000 + b);
    cases.push_back(c);
  }
  {  // U11b: new keys cross a threshold as the last insertions
    CaseSpec c;
    c.name = "U11b_summary_new_keys_cross_threshold";
    c.num_references = 3;
    AddRandomFlushes(&c, 32, 20, 40, 16, 50, &read_id);
    // 780 pre-seeded keys plus the 16 barcodes of the flushes: the 8th new
    // key reaches 788 (1,024 buckets), later ones follow the resize.
    for (uint64_t b = 0; b < 780; ++b) c.summary_preseed.push_back(0x20000 + b);
    cases.push_back(c);
  }
  {  // U11c: some barcodes seeded, others new, new ones in several references
    CaseSpec c;
    c.name = "U11c_summary_mixed_keys";
    c.num_references = 4;
    AddRandomFlushes(&c, 33, 30, 30, 40, 30, &read_id);
    for (uint64_t b = 0; b < 40; b += 2) c.summary_preseed.push_back(0x1000 + b);
    cases.push_back(c);
  }
  {  // U12
    CaseSpec c;
    c.name = "U12_no_tn5_no_dedup";
    c.num_references = 2;
    c.tn5_shift = false;
    c.remove_pcr_duplicates = false;
    AddRandomFlushes(&c, 41, 10, 30, 6, 10, &read_id);
    c.summary_preseed = MostBarcodes(6);
    cases.push_back(c);
  }
  {  // U13 BED text with a translate table
    CaseSpec c;
    c.name = "U13_bed_text_translated";
    c.num_references = 3;
    c.output = OutputMode::kBed;
    c.translate_table = true;
    AddRandomFlushes(&c, 51, 30, 20, 10, 20, &read_id);
    c.summary_preseed = MostBarcodes(10);
    cases.push_back(c);
  }
  {  // U13b BED text, bulk-level dedup
    CaseSpec c;
    c.name = "U13b_bed_text_bulk_dedup";
    c.num_references = 3;
    c.output = OutputMode::kBed;
    c.bulk_level_dedup = true;
    AddRandomFlushes(&c, 52, 30, 20, 10, 20, &read_id);
    c.whitelist = WhitelistFor(10, 53);
    c.summary_preseed = MostBarcodes(10);
    cases.push_back(c);
  }
  {  // U13c TagAlign
    CaseSpec c;
    c.name = "U13c_tagalign";
    c.num_references = 3;
    c.output = OutputMode::kTagAlign;
    c.macs3_buffer = false;
    AddRandomFlushes(&c, 54, 30, 20, 10, 20, &read_id);
    c.summary_preseed = MostBarcodes(10);
    cases.push_back(c);
  }
  {  // U14 open-file limits (40 files for reference 0)
    CaseSpec c;
    c.name = "U14_open_file_limits";
    c.num_references = 2;
    AddRandomFlushes(&c, 61, 40, 10, 6, 20, &read_id);
    c.summary_preseed = MostBarcodes(6);
    cases.push_back(c);
  }
  {  // U15 a task fails: corrupt the first record of one spill file
    CaseSpec c;
    c.name = "U15_task_failure_cleanup";
    c.num_references = 3;
    AddRandomFlushes(&c, 71, 5, 20, 6, 20, &read_id);
    c.before_merge = [](const std::string &tmp_dir) {
      const std::vector<std::string> names = ListDir(tmp_dir);
      if (names.empty()) return;
      const std::string path = tmp_dir + "/" + names[names.size() / 2];
      const int fd = open(path.c_str(), O_WRONLY);
      if (fd >= 0) {
        const uint32_t zero = 0;
        // File header (32) + block header (16) + record length (4).
        (void)pwrite(fd, &zero, sizeof(zero), 52);
        close(fd);
      }
    };
    cases.push_back(c);
  }
  return cases;
}

std::vector<Variant> VariantsFor(const CaseSpec &spec) {
  std::vector<Variant> variants;
  // U15 exercises the per-reference path's cleanup; the serial merge stops
  // on the corrupt record without removing its spill files, as in v1.1.x.
  if (spec.name == "U15_task_failure_cleanup") {
#ifdef LOWMEM_FINALIZE_HAS_THREADS
    Variant v;
    v.name = "t7_failure";
    v.finalize_threads = 7;
    v.expect_failure = true;
    v.expected_message = "invalid ATAC k-way record header";
    v.compare_with_serial = false;
    variants.push_back(v);
#endif
    return variants;
  }
  Variant serial;
  serial.name = "serial";
  variants.push_back(serial);
#ifdef LOWMEM_FINALIZE_HAS_THREADS
  const int thread_counts[] = {2, 7, 32, 64};
  for (int t : thread_counts) {
    Variant v;
    v.name = "t" + std::to_string(t);
    v.finalize_threads = t;
    variants.push_back(v);
  }
  Variant automatic;
  automatic.name = "auto";
  automatic.finalize_threads = 0;
  variants.push_back(automatic);
  Variant one_permit;
  one_permit.name = "t7_permits1";
  one_permit.finalize_threads = 7;
  one_permit.permits = 1;
  variants.push_back(one_permit);
  if (spec.name == "U02_three_references_400_flushes" ||
      spec.name == "U05_bulk_level_dedup" ||
      spec.name == "U11c_summary_mixed_keys") {
    Variant three;
    three.name = "t7_permits3";
    three.finalize_threads = 7;
    three.permits = 3;
    variants.push_back(three);
    Variant serial_permit;
    serial_permit.name = "serial_permits1";
    serial_permit.finalize_threads = 1;
    serial_permit.permits = 1;
    variants.push_back(serial_permit);
  }
  if (spec.name == "U14_open_file_limits") {
    Variant tight;
    tight.name = "t7_nofile_150";  // room for one or two tasks
    tight.finalize_threads = 7;
    tight.rlimit_nofile = 150;
    variants.push_back(tight);
    Variant too_tight;
    too_tight.name = "t7_nofile_40";  // below one reference's 40 files
    too_tight.finalize_threads = 7;
    too_tight.rlimit_nofile = 40;
    too_tight.expect_failure = true;
    too_tight.expected_message = "Low-memory finalization needs";
    too_tight.compare_with_serial = false;
    variants.push_back(too_tight);
  }
#else
  (void)spec;
#endif
  return variants;
}

#ifdef LOWMEM_FINALIZE_HAS_THREADS
// U16: the lean decode yields the same fragment fields as the full decode and
// rejects the same invalid record headers.
bool SameFragment(const PairedEndMappingWithBarcode &a,
                  const PairedEndMappingWithBarcode &b) {
  return a.read_id_ == b.read_id_ && a.cell_barcode_ == b.cell_barcode_ &&
         a.fragment_start_position_ == b.fragment_start_position_ &&
         a.fragment_length_ == b.fragment_length_ && a.mapq_ == b.mapq_ &&
         a.direction_ == b.direction_ && a.is_unique_ == b.is_unique_ &&
         a.num_dups_ == b.num_dups_ &&
         a.positive_alignment_length_ == b.positive_alignment_length_ &&
         a.negative_alignment_length_ == b.negative_alignment_length_;
}

bool RunLeanDecodeCheck() {
  std::mt19937_64 rng(91);
  std::vector<uint8_t> encoded;
  std::string error;
  for (int i = 0; i < 20000; ++i) {
    AtacSpillRecord rec =
        Rec(rng(), rng() & 0xffffffffULL, static_cast<uint32_t>(rng() % 250000000),
            static_cast<uint16_t>(1 + rng() % 60000),
            static_cast<uint8_t>(rng() % 64), static_cast<uint8_t>(rng() % 2),
            static_cast<uint8_t>(rng() % 2));
    rec.num_dups_ = static_cast<uint8_t>(rng() % 256);
    rec.positive_alignment_length_ = static_cast<uint16_t>(rng() % 400);
    rec.negative_alignment_length_ = static_cast<uint16_t>(rng() % 400);
    rec.SetYHit((rng() % 2) != 0);
    if (!EncodeAtacKwaySpillRecord(rec, 0, &encoded, &error) ||
        encoded.size() != sizeof(AtacKwaySpillRecordHeaderV1)) {
      std::cerr << "[lowmem-parallel] U16 encode failed: " << error << "\n";
      return false;
    }
    AtacSpillRecord full;
    PairedEndMappingWithBarcode lean;
    AtacKwaySpillRecordHeaderV1 header;
    memcpy(&header, encoded.data(), sizeof(header));
    if (!DecodeAtacKwaySpillRecord(encoded.data(), encoded.size(), 0, &full,
                                   &error) ||
        !DecodeAtacKwaySpillRecordLean(header, 0, &lean, &error) ||
        !SameFragment(full, lean)) {
      std::cerr << "[lowmem-parallel] U16 record " << i
                << " decodes differently (" << error << ")\n";
      return false;
    }
  }
  // Each invalid header must be rejected by both decoders.
  AtacSpillRecord base = Rec(5, 0x1000, 1000, 150, 40);
  if (!EncodeAtacKwaySpillRecord(base, 0, &encoded, &error)) return false;
  AtacKwaySpillRecordHeaderV1 good;
  memcpy(&good, encoded.data(), sizeof(good));
  std::vector<std::pair<std::string, std::function<void(AtacKwaySpillRecordHeaderV1 *)>>> bad = {
      {"magic", [](AtacKwaySpillRecordHeaderV1 *h) { h->magic ^= 1; }},
      {"codec version", [](AtacKwaySpillRecordHeaderV1 *h) { h->codec_version = 9; }},
      {"fixed header bytes", [](AtacKwaySpillRecordHeaderV1 *h) { h->fixed_header_bytes = 40; }},
      {"zero length", [](AtacKwaySpillRecordHeaderV1 *h) { h->fragment_length = 0; }},
      {"row flags", [](AtacKwaySpillRecordHeaderV1 *h) { h->row_flags = 2; }},
      {"quality bytes", [](AtacKwaySpillRecordHeaderV1 *h) { h->barcode_quality_bytes = 1; }},
      {"n mask", [](AtacKwaySpillRecordHeaderV1 *h) { h->raw_barcode_n_mask = 1; }},
      {"bam pair bytes", [](AtacKwaySpillRecordHeaderV1 *h) { h->bam_pair_bytes = 4; }},
  };
  for (const auto &entry : bad) {
    AtacKwaySpillRecordHeaderV1 h = good;
    entry.second(&h);
    AtacSpillRecord full;
    PairedEndMappingWithBarcode lean;
    const bool full_ok = DecodeAtacKwaySpillRecord(&h, sizeof(h), 0, &full, &error);
    const bool lean_ok = DecodeAtacKwaySpillRecordLean(h, 0, &lean, &error);
    if (full_ok || lean_ok) {
      std::cerr << "[lowmem-parallel] U16 invalid " << entry.first
                << " accepted (full " << full_ok << ", lean " << lean_ok
                << ")\n";
      return false;
    }
  }
  return true;
}
#endif

}  // namespace

int main() {
  const char *artifact_root = getenv("CHROMAP_ARTIFACT_ROOT");
  const std::string base = (artifact_root != nullptr && artifact_root[0])
                               ? std::string(artifact_root) +
                                     "/lowmem_parallel_finalize"
                               : std::string("/tmp");
  MakeDirs(base);
  std::string pattern = base + "/run.XXXXXX";
  std::vector<char> buffer(pattern.begin(), pattern.end());
  buffer.push_back('\0');
  if (mkdtemp(buffer.data()) == nullptr) {
    std::cerr << "[lowmem-parallel] cannot create a work directory under "
              << base << ": " << strerror(errno) << "\n";
    return 2;
  }
  const std::string root = buffer.data();
  const char *goldens_env = getenv("LOWMEM_TEST_GOLDENS_DIR");
  const std::string goldens = goldens_env != nullptr ? goldens_env : "";
  const bool write_goldens = getenv("LOWMEM_TEST_WRITE_GOLDENS") != nullptr;

  int failures = 0;
  int runs = 0;
  for (const CaseSpec &spec : BuildCases()) {
    std::string serial_dir;
    bool case_ok = true;
    for (const Variant &variant : VariantsFor(spec)) {
      ++runs;
      std::string dir;
      if (!RunVariant(spec, variant, root, &dir)) {
        ++failures;
        case_ok = false;
        continue;
      }
      if (variant.name == "serial") {
        serial_dir = dir;
        if (!goldens.empty()) {
          const std::string golden_dir = goldens + "/" + spec.name;
          if (write_goldens) {
            MakeDirs(golden_dir);
            for (const char *name : kOutputs) {
              const std::string path = dir + "/" + name;
              if (FileExists(path)) WriteFile(golden_dir + "/" + name, ReadFile(path));
            }
            WriteFile(golden_dir + "/counters.txt",
                      ReadFile(dir + "/counters.txt"));
          } else {
            std::string difference;
            if (!SameOutputs(golden_dir, dir, &difference)) {
              std::cerr << "[lowmem-parallel] FAIL " << spec.name
                        << "/serial against the saved goldens: " << difference
                        << "\n";
              ++failures;
              case_ok = false;
            }
          }
        }
        continue;
      }
      if (!variant.compare_with_serial) continue;
      std::string difference;
      if (serial_dir.empty() || !SameOutputs(serial_dir, dir, &difference)) {
        std::cerr << "[lowmem-parallel] FAIL " << spec.name << "/"
                  << variant.name << ": "
                  << (serial_dir.empty() ? "no serial run" : difference)
                  << " (" << dir << ")\n";
        ++failures;
        case_ok = false;
      }
    }
    std::cerr << "[lowmem-parallel] " << (case_ok ? "PASS " : "FAIL ")
              << spec.name << "\n";
    if (case_ok) RemoveTree(root + "/" + spec.name);
  }
#ifdef LOWMEM_FINALIZE_HAS_THREADS
  ++runs;
  if (RunLeanDecodeCheck()) {
    std::cerr << "[lowmem-parallel] PASS U16_lean_decode\n";
  } else {
    std::cerr << "[lowmem-parallel] FAIL U16_lean_decode\n";
    ++failures;
  }
#endif
  if (failures == 0) {
    rmdir(root.c_str());
    std::cerr << "[lowmem-parallel] PASS: " << runs << " runs\n";
    return 0;
  }
  std::cerr << "[lowmem-parallel] FAIL: " << failures << " failures in "
            << runs << " runs (artifacts in " << root << ")\n";
  return 1;
}
