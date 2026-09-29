// Low-memory overflow merge edge cases (Chromap Suite 1.1.1).
//
// Drives MappingWriter<AtacSpillRecord> in sidecar-only mode through its
// public low-memory API with synthetic spill runs (no FASTQ, index or
// fixture). Each case runs in a forked child so that its exit status and
// stderr can be checked.
//
//   missing_spill_file   A spill file is removed before the merge. The run
//                        must stop with an error instead of dropping the
//                        records in that file.
//   fd_limit             RLIMIT_NOFILE is below the number of spill files of
//                        one reference. The run must stop with an error
//                        instead of merging only the files it could open.
//   fd_limit_control     The same spills with a normal limit: every record
//                        is written.
//   bulk_dedup_empty_whitelist
//                        Bulk-level dedup with an empty whitelist table (a
//                        barcoded run without --barcode-whitelist). No
//                        abundance is read from the empty table; the record
//                        with the most duplicates is kept.
//   bulk_dedup_absent_barcode
//                        Bulk-level dedup with barcodes that are not in the
//                        whitelist table (--output-mappings-not-in-whitelist).
//                        An absent barcode has abundance 0.
#include <dirent.h>
#include <fcntl.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "atac_spill_record.h"
#include "khash.h"
#include "mapping_parameters.h"
#include "mapping_writer.h"
#include "sequence_batch.h"
#include "sequence_effective_range.h"
#include "utils.h"

// The khash helpers for k64_seq are declared in namespace chromap.
using namespace chromap;

namespace {

const uint32_t kBarcodeLength = 16;
const int kExitContentMismatch = 3;
const int kExitNoErrorRaised = 4;

struct OutputRecord {
  int32_t chrom_id;
  int32_t start;
  int32_t end;
  uint32_t count;
  uint64_t barcode;
};

struct CaseSpec {
  std::string name;
  uint32_t num_references = 1;
  bool bulk_level_dedup = false;
  // One entry per flush; each flush holds records per reference.
  std::vector<std::vector<std::vector<AtacSpillRecord>>> flushes;
  // Whitelist table content: (barcode, abundance). Empty means an empty table.
  std::vector<std::pair<uint64_t, uint64_t>> whitelist;
  // Applied in the child after the spills are written and before the merge.
  std::function<void(const std::string &tmp_dir)> before_merge;
  bool expect_error = false;
  std::string expected_message;
  std::vector<OutputRecord> expected_records;
};

AtacSpillRecord MakeRecord(uint64_t read_id, uint64_t barcode, uint32_t start,
                           uint16_t length, uint8_t mapq) {
  PairedEndMappingWithBarcode bed(read_id, barcode, start, length, mapq,
                                  /*direction=*/1, /*is_unique=*/1,
                                  /*num_dups=*/1,
                                  /*positive_alignment_length=*/50,
                                  /*negative_alignment_length=*/50);
  return AtacSpillRecord(bed);
}

std::vector<std::string> ListSpillFiles(const std::string &dir) {
  std::vector<std::string> paths;
  DIR *d = opendir(dir.c_str());
  if (d == nullptr) {
    return paths;
  }
  while (struct dirent *entry = readdir(d)) {
    const std::string name = entry->d_name;
    if (name.compare(0, 8, "chromap_") == 0) {
      paths.push_back(dir + "/" + name);
    }
  }
  closedir(d);
  std::sort(paths.begin(), paths.end());
  return paths;
}

void RemoveTree(const std::string &dir) {
  DIR *d = opendir(dir.c_str());
  if (d == nullptr) {
    return;
  }
  while (struct dirent *entry = readdir(d)) {
    const std::string name = entry->d_name;
    if (name == "." || name == "..") {
      continue;
    }
    const std::string path = dir + "/" + name;
    struct stat st;
    if (lstat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode)) {
      RemoveTree(path);
    } else {
      unlink(path.c_str());
    }
  }
  closedir(d);
  rmdir(dir.c_str());
}

bool ReadSidecar(const std::string &path, std::vector<OutputRecord> *records,
                 std::string *error) {
  std::ifstream in(path.c_str(), std::ios::binary);
  if (!in) {
    *error = "cannot open sidecar " + path;
    return false;
  }
  char header[32];
  if (!in.read(header, sizeof(header)) || memcmp(header, "AEV1", 4) != 0) {
    *error = "bad AEV1 header in " + path;
    return false;
  }
  uint64_t num_records = 0;
  memcpy(&num_records, header + 24, sizeof(num_records));
  records->clear();
  for (uint64_t i = 0; i < num_records; ++i) {
    char raw[24];
    if (!in.read(raw, sizeof(raw))) {
      *error = "truncated AEV1 body in " + path;
      return false;
    }
    OutputRecord r;
    memcpy(&r.chrom_id, raw, 4);
    memcpy(&r.start, raw + 4, 4);
    memcpy(&r.end, raw + 8, 4);
    memcpy(&r.count, raw + 12, 4);
    memcpy(&r.barcode, raw + 16, 8);
    records->push_back(r);
  }
  char extra;
  if (in.read(&extra, 1)) {
    *error = "trailing bytes after the declared AEV1 records";
    return false;
  }
  return true;
}

std::string Describe(const OutputRecord &r) {
  std::ostringstream out;
  out << "(chrom " << r.chrom_id << ", " << r.start << "-" << r.end
      << ", count " << r.count << ", barcode " << r.barcode << ")";
  return out.str();
}

// Runs in the child process. Returns the exit code.
int RunCaseInChild(const CaseSpec &spec, const std::string &case_dir) {
  const std::string tmp_dir = case_dir + "/tmp";
  const std::string sidecar = case_dir + "/out.bin";

  MappingParameters params;
  params.low_memory_mode = true;
  params.atac_sidecar_only = true;
  params.atac_fragment_binary_output_file_path = sidecar;
  params.mapping_output_format = chromap::MAPPINGFORMAT_BED;
  params.temp_directory_path = tmp_dir;
  params.remove_pcr_duplicates = true;
  params.remove_pcr_duplicates_at_bulk_level = spec.bulk_level_dedup;
  params.is_bulk_data = false;
  params.mapq_threshold = 30;
  params.Tn5_shift = false;
  params.num_threads = 1;

  SequenceBatch reference(spec.num_references, SequenceEffectiveRange());
  for (uint32_t rid = 0; rid < spec.num_references; ++rid) {
    reference.AssignLoadedReferenceMetadata(
        rid, "ref" + std::to_string(rid), 1000000);
  }

  khash_t(k64_seq) *whitelist = kh_init(k64_seq);
  for (const auto &entry : spec.whitelist) {
    int ret = 0;
    khiter_t it = kh_put(k64_seq, whitelist, entry.first, &ret);
    kh_value(whitelist, it) = entry.second;
  }

  {
    MappingWriter<AtacSpillRecord> writer(params, kBarcodeLength,
                                          std::vector<int>());
    writer.OutputHeader(spec.num_references, reference);
    for (const auto &flush : spec.flushes) {
      std::vector<std::vector<AtacSpillRecord>> by_reference = flush;
      by_reference.resize(spec.num_references);
      for (auto &records : by_reference) {
        std::sort(records.begin(), records.end());
      }
      writer.OutputTempMappingsToOverflow(spec.num_references, by_reference);
      MappingWriter<AtacSpillRecord>::RotateThreadOverflowWriter();
    }
    if (spec.before_merge) {
      spec.before_merge(tmp_dir);
    }
    writer.ProcessAndOutputMappingsInLowMemoryFromOverflow(
        0, spec.num_references, reference, whitelist);
  }
  kh_destroy(k64_seq, whitelist);

  std::vector<OutputRecord> records;
  std::string error;
  const bool read_ok = ReadSidecar(sidecar, &records, &error);
  if (spec.expect_error) {
    std::cerr << "[child] the merge completed without an error; "
              << (read_ok ? std::to_string(records.size()) +
                                " records were written"
                          : error)
              << "\n";
    return kExitNoErrorRaised;
  }
  if (!read_ok) {
    std::cerr << "[child] " << error << "\n";
    return kExitContentMismatch;
  }
  if (records.size() != spec.expected_records.size()) {
    std::cerr << "[child] " << records.size() << " records written, expected "
              << spec.expected_records.size() << "\n";
    return kExitContentMismatch;
  }
  for (size_t i = 0; i < records.size(); ++i) {
    const OutputRecord &got = records[i];
    const OutputRecord &want = spec.expected_records[i];
    if (got.chrom_id != want.chrom_id || got.start != want.start ||
        got.end != want.end || got.count != want.count ||
        got.barcode != want.barcode) {
      std::cerr << "[child] record " << i << " is " << Describe(got)
                << ", expected " << Describe(want) << "\n";
      return kExitContentMismatch;
    }
  }
  if (!ListSpillFiles(tmp_dir).empty()) {
    std::cerr << "[child] spill files were left in " << tmp_dir << "\n";
    return kExitContentMismatch;
  }
  return 0;
}

std::string ReadFile(const std::string &path) {
  std::ifstream in(path.c_str());
  std::ostringstream out;
  out << in.rdbuf();
  return out.str();
}

bool RunCase(const CaseSpec &spec, const std::string &root) {
  const std::string case_dir = root + "/" + spec.name;
  RemoveTree(case_dir);
  if (mkdir(case_dir.c_str(), 0755) != 0 ||
      mkdir((case_dir + "/tmp").c_str(), 0755) != 0) {
    std::cerr << "[lowmem-edge] FAIL " << spec.name << ": cannot create "
              << case_dir << "\n";
    return false;
  }
  const std::string stderr_path = case_dir + "/stderr.txt";
  std::cout.flush();
  std::cerr.flush();
  const pid_t pid = fork();
  if (pid < 0) {
    std::cerr << "[lowmem-edge] FAIL " << spec.name << ": fork failed\n";
    return false;
  }
  if (pid == 0) {
    const int fd = open(stderr_path.c_str(), O_WRONLY | O_CREAT | O_TRUNC,
                        0644);
    if (fd >= 0) {
      dup2(fd, STDERR_FILENO);
      close(fd);
    }
    const int code = RunCaseInChild(spec, case_dir);
    std::cerr.flush();
    _exit(code);
  }
  int status = 0;
  if (waitpid(pid, &status, 0) != pid) {
    std::cerr << "[lowmem-edge] FAIL " << spec.name << ": waitpid failed\n";
    return false;
  }
  const std::string child_stderr = ReadFile(stderr_path);
  std::string verdict;
  bool pass = false;
  if (WIFSIGNALED(status)) {
    verdict = "child killed by signal " + std::to_string(WTERMSIG(status));
  } else if (!WIFEXITED(status)) {
    verdict = "child did not exit normally";
  } else if (spec.expect_error) {
    const int code = WEXITSTATUS(status);
    if (code == 0 || code == kExitNoErrorRaised ||
        code == kExitContentMismatch) {
      verdict = "expected an error exit, got exit " + std::to_string(code);
    } else if (child_stderr.find(spec.expected_message) == std::string::npos) {
      verdict = "error exit without the expected message '" +
                spec.expected_message + "'";
    } else {
      pass = true;
    }
  } else if (WEXITSTATUS(status) == 0) {
    pass = true;
  } else {
    verdict = "exit " + std::to_string(WEXITSTATUS(status));
  }
  if (pass) {
    std::cerr << "[lowmem-edge] PASS " << spec.name << "\n";
    RemoveTree(case_dir);
    return true;
  }
  std::cerr << "[lowmem-edge] FAIL " << spec.name << ": " << verdict
            << " (artifacts in " << case_dir << ")\n";
  if (!child_stderr.empty()) {
    std::cerr << "  child stderr (last lines):\n";
    std::istringstream lines(child_stderr);
    std::vector<std::string> tail;
    std::string line;
    while (std::getline(lines, line)) {
      tail.push_back(line);
    }
    const size_t first = tail.size() > 6 ? tail.size() - 6 : 0;
    for (size_t i = first; i < tail.size(); ++i) {
      std::cerr << "    " << tail[i] << "\n";
    }
  }
  return false;
}

// 1 reference, `num_flushes` flushes of one distinct record each, so the
// reference has one spill file per flush.
CaseSpec ManyFilesOneReference(const std::string &name, int num_flushes) {
  CaseSpec spec;
  spec.name = name;
  spec.num_references = 1;
  for (int f = 0; f < num_flushes; ++f) {
    std::vector<std::vector<AtacSpillRecord>> flush(1);
    const uint32_t start = 1000 + 100 * static_cast<uint32_t>(f);
    flush[0].push_back(MakeRecord(static_cast<uint64_t>(f), 0x1111u + f,
                                  start, 200, 40));
    spec.flushes.push_back(flush);
    OutputRecord out = {0, static_cast<int32_t>(start),
                        static_cast<int32_t>(start + 200), 1,
                        static_cast<uint64_t>(0x1111u + f)};
    spec.expected_records.push_back(out);
  }
  return spec;
}

}  // namespace

int main() {
  std::string base;
  const char *artifact_root = getenv("CHROMAP_ARTIFACT_ROOT");
  if (artifact_root != nullptr && artifact_root[0] != '\0') {
    base = std::string(artifact_root) + "/lowmem_overflow_edge_cases";
    std::string partial;
    std::istringstream parts(base);
    std::string part;
    while (std::getline(parts, part, '/')) {
      if (part.empty()) {
        partial += "/";
        continue;
      }
      partial += part;
      mkdir(partial.c_str(), 0755);
      partial += "/";
    }
  } else {
    base = "/tmp";
  }
  std::string root_template = base + "/run.XXXXXX";
  std::vector<char> root_buffer(root_template.begin(), root_template.end());
  root_buffer.push_back('\0');
  if (mkdtemp(root_buffer.data()) == nullptr) {
    std::cerr << "[lowmem-edge] cannot create a work directory under " << base
              << ": " << strerror(errno) << "\n";
    return 2;
  }
  const std::string root = root_buffer.data();

  std::vector<CaseSpec> cases;

  // A spill file removed before the merge.
  {
    CaseSpec spec;
    spec.name = "missing_spill_file";
    spec.num_references = 2;
    for (int f = 0; f < 3; ++f) {
      std::vector<std::vector<AtacSpillRecord>> flush(2);
      for (uint32_t rid = 0; rid < 2; ++rid) {
        flush[rid].push_back(MakeRecord(10 * f + rid, 0x2000u + f,
                                        5000 + 300 * f, 150, 40));
      }
      spec.flushes.push_back(flush);
    }
    spec.before_merge = [](const std::string &tmp_dir) {
      const std::vector<std::string> spills = ListSpillFiles(tmp_dir);
      if (!spills.empty()) {
        unlink(spills[spills.size() / 2].c_str());
      }
    };
    spec.expect_error = true;
    spec.expected_message = "Cannot open low-memory overflow file";
    cases.push_back(spec);
  }

  // More spill files for one reference than open files allowed.
  {
    CaseSpec spec = ManyFilesOneReference("fd_limit", 40);
    spec.before_merge = [](const std::string &) {
      struct rlimit limit;
      limit.rlim_cur = 24;
      limit.rlim_max = 24;
      setrlimit(RLIMIT_NOFILE, &limit);
    };
    spec.expect_error = true;
    spec.expected_message = "Cannot open low-memory overflow file";
    spec.expected_records.clear();
    cases.push_back(spec);
  }

  // Control: the same spills with the normal limit.
  cases.push_back(ManyFilesOneReference("fd_limit_control", 40));

  // Bulk-level dedup, empty whitelist table. At 7000/180 barcode 0x30 has one
  // record and 0x31 has two; 0x31 has more duplicates and is kept. The count
  // is the whole position group (3). 9000/200 is a single record.
  {
    CaseSpec spec;
    spec.name = "bulk_dedup_empty_whitelist";
    spec.num_references = 1;
    spec.bulk_level_dedup = true;
    std::vector<std::vector<AtacSpillRecord>> flush1(1), flush2(1);
    flush1[0].push_back(MakeRecord(1, 0x30, 7000, 180, 40));
    flush1[0].push_back(MakeRecord(2, 0x31, 7000, 180, 40));
    flush2[0].push_back(MakeRecord(3, 0x31, 7000, 180, 41));
    flush2[0].push_back(MakeRecord(4, 0x32, 9000, 200, 40));
    spec.flushes.push_back(flush1);
    spec.flushes.push_back(flush2);
    spec.expected_records.push_back({0, 7000, 7180, 3, 0x31});
    spec.expected_records.push_back({0, 9000, 9200, 1, 0x32});
    cases.push_back(spec);
  }

  // Bulk-level dedup with barcodes absent from a non-empty whitelist table.
  // At 7000/180: 0x40 (absent) sorts before 0x41 (abundance 7); equal
  // duplicates, so the higher abundance (0x41) is kept. At 8000/180: 0x50
  // (present, abundance 0) sorts before 0x51 (absent); equal duplicates and
  // equal abundance (0), so the first (0x50) is kept. At 9000/200: 0x60
  // (absent) is alone and kept.
  {
    CaseSpec spec;
    spec.name = "bulk_dedup_absent_barcode";
    spec.num_references = 1;
    spec.bulk_level_dedup = true;
    for (uint64_t filler = 0x1000; filler < 0x1010; ++filler) {
      spec.whitelist.push_back(std::make_pair(filler, filler));
    }
    spec.whitelist.push_back(std::make_pair(0x41u, 7u));
    spec.whitelist.push_back(std::make_pair(0x50u, 0u));
    std::vector<std::vector<AtacSpillRecord>> flush1(1), flush2(1);
    flush1[0].push_back(MakeRecord(1, 0x40, 7000, 180, 40));
    flush2[0].push_back(MakeRecord(2, 0x41, 7000, 180, 40));
    flush1[0].push_back(MakeRecord(3, 0x50, 8000, 180, 40));
    flush2[0].push_back(MakeRecord(4, 0x51, 8000, 180, 40));
    flush1[0].push_back(MakeRecord(5, 0x60, 9000, 200, 40));
    spec.flushes.push_back(flush1);
    spec.flushes.push_back(flush2);
    spec.expected_records.push_back({0, 7000, 7180, 2, 0x41});
    spec.expected_records.push_back({0, 8000, 8180, 2, 0x50});
    spec.expected_records.push_back({0, 9000, 9200, 1, 0x60});
    cases.push_back(spec);
  }

  int failures = 0;
  for (const CaseSpec &spec : cases) {
    if (!RunCase(spec, root)) {
      ++failures;
    }
  }
  if (failures == 0) {
    rmdir(root.c_str());
    std::cerr << "[lowmem-edge] PASS: " << cases.size() << " cases\n";
    return 0;
  }
  std::cerr << "[lowmem-edge] FAIL: " << failures << " of " << cases.size()
            << " cases\n";
  return 1;
}
