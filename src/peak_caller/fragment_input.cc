#include "fragment_input.h"

#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <zlib.h>

namespace chromap {
namespace peaks {

namespace {

void TrimInPlace(std::string* s) {
  while (!s->empty() && (s->back() == '\r' || s->back() == '\n')) {
    s->pop_back();
  }
  size_t a = 0;
  while (a < s->size() && std::isspace(static_cast<unsigned char>((*s)[a]))) {
    ++a;
  }
  if (a > 0) {
    s->erase(0, a);
  }
  while (!s->empty() && std::isspace(static_cast<unsigned char>(s->back()))) {
    s->pop_back();
  }
}

bool ParseInt32(const char* p, int32_t* out) {
  char* endp = nullptr;
  long v = std::strtol(p, &endp, 10);
  if (endp == p || (endp && *endp != '\0')) {
    return false;
  }
  *out = static_cast<int32_t>(v);
  return true;
}

}  // namespace

bool LoadFragmentsFromTsv(const std::string& path,
                          std::vector<ChromFragments>* by_chrom) {
  by_chrom->clear();
  std::unordered_map<std::string, size_t> index_of;
  const bool use_gz = path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
  gzFile gz = nullptr;
  std::ifstream fin;
  if (use_gz) {
    gz = gzopen(path.c_str(), "rb");
    if (gz == nullptr) {
      return false;
    }
  } else {
    fin.open(path);
    if (!fin) {
      return false;
    }
  }
  auto process_line = [&](const std::string& raw) {
    std::string s = raw;
    TrimInPlace(&s);
    if (s.empty() || s[0] == '#') {
      return;
    }
    std::vector<std::string> fields;
    std::istringstream is(s);
    std::string col;
    while (std::getline(is, col, '\t')) {
      fields.push_back(std::move(col));
    }
    if (fields.size() < 5) {
      return;
    }
    int32_t start = 0, end = 0, dups = 1;
    if (!ParseInt32(fields[1].c_str(), &start) ||
        !ParseInt32(fields[2].c_str(), &end) ||
        !ParseInt32(fields[4].c_str(), &dups)) {
      return;
    }
    if (end <= start || dups < 0) {
      return;
    }
    if (dups == 0) {
      return;
    }
    const std::string& chrom = fields[0];
    size_t c = 0;
    auto it = index_of.find(chrom);
    if (it == index_of.end()) {
      c = by_chrom->size();
      index_of.emplace(chrom, c);
      ChromFragments block;
      block.name = chrom;
      by_chrom->push_back(std::move(block));
    } else {
      c = it->second;
    }
    Fragment f;
    f.start = start;
    f.end = end;
    f.count = dups;
    ChromFragments& ch = (*by_chrom)[c];
    ch.frags.push_back(f);
    if (end > ch.max_end) {
      ch.max_end = end;
    }
  };
  if (use_gz) {
    char buf[1 << 16];
    for (;;) {
      if (gzgets(gz, buf, sizeof(buf)) == nullptr) {
        break;
      }
      process_line(std::string(buf));
    }
    if (gzclose(gz) != Z_OK) {
      return false;
    }
  } else {
    std::string line;
    while (std::getline(fin, line)) {
      process_line(line);
    }
  }
  if (by_chrom->empty()) {
    return false;
  }
  for (auto& ch : *by_chrom) {
    std::sort(
        ch.frags.begin(), ch.frags.end(),
        [](const Fragment& a, const Fragment& b) {
          if (a.start != b.start) {
            return a.start < b.start;
          }
          return a.end < b.end;
        });
  }
  std::sort(by_chrom->begin(), by_chrom->end(),
            [](const ChromFragments& a, const ChromFragments& b) {
              return a.name < b.name;
            });
  return true;
}

}  // namespace peaks
}  // namespace chromap
