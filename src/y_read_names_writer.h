#ifndef Y_READ_NAMES_WRITER_H_
#define Y_READ_NAMES_WRITER_H_

#include <cstdint>
#include <fstream>
#include <string>
#include <unordered_set>
#include <vector>

#include "utils.h"

namespace chromap {

// Helper class to write normalized Y read names
class YReadNamesWriter {
 public:
  YReadNamesWriter(const std::string &output_path)
      : output_path_(output_path), num_written_(0) {
    output_file_.open(output_path_);
    if (!output_file_.is_open()) {
      ExitWithMessage("Failed to open Y read names output file: " + output_path_);
    }
  }

  ~YReadNamesWriter() {
    if (output_file_.is_open()) {
      output_file_.close();
    }
  }

  // Normalize read name: strip leading @, stop at whitespace, strip trailing /1 or /2
  static std::string NormalizeReadName(const char *read_name) {
    if (!read_name) return "";
    
    std::string name(read_name);
    
    // Strip leading @ if present
    if (!name.empty() && name[0] == '@') {
      name = name.substr(1);
    }
    
    // Stop at first whitespace
    size_t space_pos = name.find_first_of(" \t\n\r");
    if (space_pos != std::string::npos) {
      name = name.substr(0, space_pos);
    }
    
    // Strip trailing /1 or /2
    if (name.length() >= 2) {
      if (name.substr(name.length() - 2) == "/1" || name.substr(name.length() - 2) == "/2") {
        name = name.substr(0, name.length() - 2);
      }
    }
    
    return name;
  }

  // Write a normalized read name (deduplicates automatically)
  void WriteReadName(uint32_t read_id, const char *read_name) {
    if (written_read_ids_.count(read_id) > 0) {
      return;  // Already written (read-level deduplication)
    }
    
    std::string normalized = NormalizeReadName(read_name);
    if (!normalized.empty()) {
      output_file_ << normalized << "\n";
      written_read_ids_.insert(read_id);
      ++num_written_;
    }
  }

  // Write multiple read names from a set of read IDs
  void WriteReadNames(const std::unordered_set<uint32_t> &read_ids,
                      const std::vector<const char*> &read_names) {
    for (uint32_t read_id : read_ids) {
      if (read_id < read_names.size() && read_names[read_id]) {
        WriteReadName(read_id, read_names[read_id]);
      }
    }
  }

  uint64_t GetNumWritten() const { return num_written_; }

 private:
  std::string output_path_;
  std::ofstream output_file_;
  std::unordered_set<uint32_t> written_read_ids_;
  uint64_t num_written_;
};

}  // namespace chromap

#endif  // Y_READ_NAMES_WRITER_H_
