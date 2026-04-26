#ifndef BARCODETRANSLATOR_H_
#define BARCODETRANSLATOR_H_

#include <cinttypes>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <zlib.h>

#include "khash.h"
#include "utils.h"

namespace chromap {

KHASH_INIT(k64_str, uint64_t, char *, 1, kh_int64_hash_func,
           kh_int64_hash_equal);

// The class for handling barcode convertion.
class BarcodeTranslator {
 public:
  BarcodeTranslator() {
    barcode_translate_table_ = NULL;
    from_bc_length_ = -1;
  }

  ~BarcodeTranslator() {
    if (barcode_translate_table_ != NULL) {
      khiter_t k;
      for (k = kh_begin(barcode_translate_table_);
           k != kh_end(barcode_translate_table_); ++k) {
        if (kh_exist(barcode_translate_table_, k))
          free(kh_value(barcode_translate_table_, k));
      }
      kh_destroy(k64_str, barcode_translate_table_);
    }
  }

  // `from_first_column` selects the table column convention. The historical
  // Chromap default is `to_bc<TAB>from_bc` (col1 is the destination /
  // value, col2 is the source / hash key); pass `from_first_column=true`
  // for the natural `from_bc<TAB>to_bc` order (e.g. ARC-style atac2gex
  // with ATAC barcode in column 1 and GEX barcode in column 2).
  void SetTranslateTable(const std::string &file,
                         bool from_first_column = false) {
    barcode_translate_table_ = kh_init(k64_str);
    table_file_path_ = file;
    from_first_column_ = from_first_column;

    gzFile barcode_translate_file = gzopen(file.c_str(), "r");
    if (barcode_translate_file == NULL) {
      std::cerr << "BarcodeTranslator: failed to open translate table file: "
                << file << std::endl;
      exit(-1);
    }
    const uint32_t line_buffer_size = 512;
    char file_line[line_buffer_size];
    while (gzgets(barcode_translate_file, file_line, line_buffer_size) !=
           NULL) {
      int line_len = strlen(file_line);
      if (line_len > 0 && file_line[line_len - 1] == '\n') {
        file_line[line_len - 1] = '\0';
      }
      std::string tmp_string(file_line);
      ProcessTranslateFileLine(tmp_string, from_first_column);
    }
    gzclose(barcode_translate_file);

    mask_ = (1ull << (2 * from_bc_length_)) - 1;
  }

  std::string Translate(uint64_t bc, uint32_t bc_length) {
    if (barcode_translate_table_ == NULL) {
      return Seed2Sequence(bc, bc_length);
    }

    std::string ret;
    uint64_t i;
    for (i = 0; i < bc_length / from_bc_length_; ++i) {
      uint64_t seed = (bc << (2 * i * from_bc_length_)) >>
                      (2 * (bc_length / from_bc_length_ - 1) * from_bc_length_);
      seed &= mask_;
      khiter_t barcode_translate_table_iter =
          kh_get(k64_str, barcode_translate_table_, seed);
      if (barcode_translate_table_iter == kh_end(barcode_translate_table_)) {
        std::cerr
            << "BarcodeTranslator: barcode '"
            << Seed2Sequence(seed, from_bc_length_)
            << "' (segment " << i
            << ", from_bc_length=" << from_bc_length_
            << ") not found in translation table.\n"
            << "  table file       : " << table_file_path_ << "\n"
            << "  column convention: "
            << (from_first_column_ ? "from_bc<TAB>to_bc (col1=hash key)"
                                   : "to_bc<TAB>from_bc (col2=hash key)")
            << "\n"
            << "  Hint: if your file is in <source><TAB><dest> order, pass\n"
            << "        --barcode-translate-from-first (chromap CLI) or set\n"
            << "        barcode_translate_from_first_column = true\n"
            << "        in libchromap MappingParameters.\n";
        exit(-1);
      }
      std::string bc_to(
          kh_value(barcode_translate_table_, barcode_translate_table_iter));
      if (i == 0) {
        ret = bc_to;
      } else {
        ret += "-" + bc_to;
      }
    }
    return ret;
  }

 private:
  khash_t(k64_str) * barcode_translate_table_;
  int from_bc_length_;
  uint64_t mask_;
  std::string table_file_path_;
  bool from_first_column_ = false;

  std::string Seed2Sequence(uint64_t seed, uint32_t seed_length) const {
    std::string sequence;
    sequence.reserve(seed_length);
    uint64_t mask_ = 3;
    for (uint32_t i = 0; i < seed_length; ++i) {
      sequence.push_back(
          Uint8ToChar((seed >> ((seed_length - 1 - i) * 2)) & mask_));
    }
    return sequence;
  }

  void ProcessTranslateFileLine(std::string &line, bool from_first_column) {
    int i;
    int len = line.length();
    if (len == 0) {
      return;
    }
    for (i = 0; i < len; ++i) {
      if (line[i] == ',' || line[i] == '\t') break;
    }
    if (i >= len) {
      return;  // single-column line; ignore
    }

    int from_offset, from_length;
    int to_offset, to_length;
    if (from_first_column) {
      // Natural "source<TAB>destination" order.
      from_offset = 0;
      from_length = i;
      to_offset = i + 1;
      to_length = len - i - 1;
    } else {
      // Historical Chromap "destination<TAB>source" order.
      to_offset = 0;
      to_length = i;
      from_offset = i + 1;
      from_length = len - i - 1;
    }

    std::string to = line.substr(to_offset, to_length);
    from_bc_length_ = from_length;
    uint64_t from_seed = GenerateSeedFromSequence(line.c_str(), len, from_offset,
                                                  from_bc_length_);

    int khash_return_code;
    khiter_t barcode_translate_table_iter = kh_put(
        k64_str, barcode_translate_table_, from_seed, &khash_return_code);
    kh_value(barcode_translate_table_, barcode_translate_table_iter) =
        strdup(to.c_str());
  }
};

}  // namespace chromap
#endif
