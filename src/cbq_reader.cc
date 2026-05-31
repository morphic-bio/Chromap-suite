#include "cbq_reader.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cstring>
#include <dlfcn.h>
#include <limits>
#include <sstream>

namespace chromap {
namespace {

const uint64_t PRESENCE_PAIRED = 1ULL << 0;
const uint64_t PRESENCE_QUALITIES = 1ULL << 1;
const uint64_t PRESENCE_HEADERS = 1ULL << 2;
const uint64_t PRESENCE_FLAGS = 1ULL << 3;

// Defensive ceilings on per-block structural dimensions. Real CBQ blocks are
// orders of magnitude below these; the limits exist only so a malformed or
// hostile header turns into a clean error instead of an astronomical
// allocation. They are well above any legitimate encoder block and below the
// thresholds where the size arithmetic below could overflow.
const uint64_t kMaxBasesPerBlock = 1ULL << 33;       // 8 GiB of bases (nuclen)
const uint64_t kMaxSequencesPerBlock = 1ULL << 30;   // ~1.07e9 segments
const uint64_t kMaxRecordsPerBlock = 1ULL << 30;     // ~1.07e9 records
const uint64_t kMaxColumnBytes = 1ULL << 33;         // 8 GiB per decoded column
const uint64_t kMaxIndexBytes = 1ULL << 30;          // 67M block entries

bool SetError(std::string *error, const std::string &message) {
  if (error != nullptr) {
    *error = message;
  }
  return false;
}

uint64_t ReadLe64(const uint8_t *bytes) {
  uint64_t value = 0;
  for (unsigned i = 0; i < 8; ++i) {
    value |= static_cast<uint64_t>(bytes[i]) << (8 * i);
  }
  return value;
}

uint16_t ReadLe16(const uint8_t *bytes) {
  return static_cast<uint16_t>(bytes[0]) |
         static_cast<uint16_t>(static_cast<uint16_t>(bytes[1]) << 8);
}

bool CheckedSize(uint64_t value, size_t *out) {
  if (value > static_cast<uint64_t>(std::numeric_limits<size_t>::max())) {
    return false;
  }
  *out = static_cast<size_t>(value);
  return true;
}

bool IsSpaceChar(char value) {
  return std::isspace(static_cast<unsigned char>(value)) != 0;
}

CbqByteSpan MakeSpan(const std::string &value, size_t begin, size_t size) {
  CbqByteSpan span;
  span.data = value.data() + begin;
  span.size = size;
  return span;
}

CbqByteSpan MakeEmptySpan() {
  CbqByteSpan span;
  span.data = "";
  span.size = 0;
  return span;
}

char ParseReadFilterSpan(char record_type, CbqByteSpan read_name_extra) {
  if (record_type != '@') {
    return 'N';
  }
  size_t field_size = 0;
  while (field_size < read_name_extra.size &&
         !IsSpaceChar(read_name_extra.data[field_size])) {
    ++field_size;
  }
  if (field_size > 3 && read_name_extra.data[1] == ':' &&
      read_name_extra.data[2] == 'Y' && read_name_extra.data[3] == ':') {
    return 'Y';
  }
  return 'N';
}

struct ParsedHeaderSpan {
  CbqByteSpan read_name;
  CbqByteSpan read_name_extra;
  char read_filter = 'N';
};

ParsedHeaderSpan ParseHeaderPayloadSpan(const std::string &backing,
                                        size_t payload_begin,
                                        size_t payload_size,
                                        bool has_qualities) {
  static const char separators[] = {'/', ' '};
  char record_type = has_qualities ? '@' : '>';
  size_t cursor = payload_begin;
  const size_t payload_end = payload_begin + payload_size;
  if (cursor < payload_end &&
      (backing[cursor] == '@' || backing[cursor] == '>')) {
    record_type = backing[cursor];
    ++cursor;
  }

  while (cursor < payload_end && IsSpaceChar(backing[cursor])) {
    ++cursor;
  }

  const size_t name_begin = cursor;
  while (cursor < payload_end && !IsSpaceChar(backing[cursor])) {
    ++cursor;
  }
  size_t name_end = cursor;
  for (char separator : separators) {
    for (size_t pos = name_begin; pos < name_end; ++pos) {
      if (backing[pos] == separator) {
        name_end = pos;
        break;
      }
    }
  }

  while (cursor < payload_end && IsSpaceChar(backing[cursor])) {
    ++cursor;
  }

  ParsedHeaderSpan parsed;
  parsed.read_name = MakeSpan(backing, name_begin, name_end - name_begin);
  parsed.read_name_extra = MakeSpan(backing, cursor, payload_end - cursor);
  parsed.read_filter = ParseReadFilterSpan(record_type, parsed.read_name_extra);
  return parsed;
}

struct FileHeaderFields {
  uint8_t version = 0;
  uint64_t presence_flags = 0;
  uint64_t compression_level = 0;
  uint64_t block_size = 0;

  bool IsPaired() const { return (presence_flags & PRESENCE_PAIRED) != 0; }
  bool HasQualities() const {
    return (presence_flags & PRESENCE_QUALITIES) != 0;
  }
  bool HasHeaders() const {
    return (presence_flags & PRESENCE_HEADERS) != 0;
  }
  bool HasFlags() const { return (presence_flags & PRESENCE_FLAGS) != 0; }
};

bool ParseFileHeader(const std::array<uint8_t, 64> &bytes,
                     FileHeaderFields *header, std::string *error) {
  if (std::memcmp(bytes.data(), "CBQFILE", 7) != 0) {
    return SetError(error, "CBQ file header magic is not CBQFILE");
  }
  header->version = bytes[7];
  if (header->version != 1) {
    std::ostringstream msg;
    msg << "unsupported CBQ file version: "
        << static_cast<unsigned>(header->version);
    return SetError(error, msg.str());
  }
  header->presence_flags = ReadLe64(bytes.data() + 8);
  header->compression_level = ReadLe64(bytes.data() + 16);
  header->block_size = ReadLe64(bytes.data() + 24);
  return true;
}

struct BlockHeaderFields {
  uint8_t version = 0;
  uint64_t len_z_seq_len = 0;
  uint64_t len_z_header_len = 0;
  uint64_t len_z_npos = 0;
  uint64_t len_z_seq = 0;
  uint64_t len_z_flags = 0;
  uint64_t len_z_headers = 0;
  uint64_t len_z_qual = 0;
  uint64_t nuclen = 0;
  uint64_t len_nef = 0;
  uint64_t num_records = 0;
  uint64_t num_sequences = 0;
};

bool ParseBlockHeader(const std::array<uint8_t, 96> &bytes,
                      BlockHeaderFields *header, std::string *error) {
  if (std::memcmp(bytes.data(), "BLK", 3) != 0) {
    return SetError(error, "CBQ block header magic is not BLK");
  }
  header->version = bytes[3];
  if (header->version != 1) {
    std::ostringstream msg;
    msg << "unsupported CBQ block version: "
        << static_cast<unsigned>(header->version);
    return SetError(error, msg.str());
  }
  header->len_z_seq_len = ReadLe64(bytes.data() + 8);
  header->len_z_header_len = ReadLe64(bytes.data() + 16);
  header->len_z_npos = ReadLe64(bytes.data() + 24);
  header->len_z_seq = ReadLe64(bytes.data() + 32);
  header->len_z_flags = ReadLe64(bytes.data() + 40);
  header->len_z_headers = ReadLe64(bytes.data() + 48);
  header->len_z_qual = ReadLe64(bytes.data() + 56);
  header->nuclen = ReadLe64(bytes.data() + 64);
  header->len_nef = ReadLe64(bytes.data() + 72);
  header->num_records = ReadLe64(bytes.data() + 80);
  header->num_sequences = ReadLe64(bytes.data() + 88);
  return true;
}

class ZstdRuntime {
 public:
  ZstdRuntime() = default;
  ~ZstdRuntime() {
    if (handle_ != nullptr) {
      dlclose(handle_);
    }
  }

  bool Load(std::string *error) {
    if (decompress_ != nullptr) {
      return true;
    }
    handle_ = dlopen("libzstd.so.1", RTLD_LAZY | RTLD_LOCAL);
    if (handle_ == nullptr) {
      handle_ = dlopen("libzstd.so", RTLD_LAZY | RTLD_LOCAL);
    }
    if (handle_ == nullptr) {
      return SetError(error, "could not load libzstd.so.1 or libzstd.so");
    }
    decompress_ = reinterpret_cast<DecompressFn>(
        dlsym(handle_, "ZSTD_decompress"));
    is_error_ = reinterpret_cast<IsErrorFn>(dlsym(handle_, "ZSTD_isError"));
    get_error_name_ = reinterpret_cast<GetErrorNameFn>(
        dlsym(handle_, "ZSTD_getErrorName"));
    if (decompress_ == nullptr || is_error_ == nullptr ||
        get_error_name_ == nullptr) {
      return SetError(error, "libzstd is missing required decompress symbols");
    }
    return true;
  }

  bool Decompress(const std::vector<uint8_t> &compressed,
                  size_t expected_size, std::vector<uint8_t> *out,
                  const std::string &column_name, std::string *error) {
    if (!Load(error)) {
      return false;
    }
    out->assign(expected_size, 0);
    if (expected_size == 0 && compressed.empty()) {
      return true;
    }
    const size_t nbytes = decompress_(out->data(), out->size(),
                                      compressed.data(), compressed.size());
    if (is_error_(nbytes)) {
      std::ostringstream msg;
      msg << "zstd decompression failed for CBQ column " << column_name
          << ": " << get_error_name_(nbytes);
      return SetError(error, msg.str());
    }
    if (nbytes != expected_size) {
      std::ostringstream msg;
      msg << "zstd decompressed CBQ column " << column_name << " to "
          << nbytes << " bytes, expected " << expected_size;
      return SetError(error, msg.str());
    }
    return true;
  }

 private:
  using DecompressFn = size_t (*)(void *, size_t, const void *, size_t);
  using IsErrorFn = unsigned (*)(size_t);
  using GetErrorNameFn = const char *(*)(size_t);

  void *handle_ = nullptr;
  DecompressFn decompress_ = nullptr;
  IsErrorFn is_error_ = nullptr;
  GetErrorNameFn get_error_name_ = nullptr;
};

ZstdRuntime &GetZstdRuntime() {
  static ZstdRuntime runtime;
  return runtime;
}

struct BitVectorLite {
  std::vector<uint64_t> words;
  uint64_t len = 0;

  bool Bit(uint64_t pos) const {
    if (pos >= len) {
      return false;
    }
    const size_t word_index = static_cast<size_t>(pos / 64);
    if (word_index >= words.size()) {
      return false;
    }
    return ((words[word_index] >> (pos % 64)) & 1ULL) != 0;
  }

  uint64_t Bits(uint64_t pos, uint64_t nbits) const {
    uint64_t value = 0;
    for (uint64_t i = 0; i < nbits; ++i) {
      if (Bit(pos + i)) {
        value |= 1ULL << i;
      }
    }
    return value;
  }
};

class SliceReader {
 public:
  explicit SliceReader(const std::vector<uint8_t> &bytes) : bytes_(bytes) {}

  bool ReadU8(uint8_t *value, std::string *error) {
    if (offset_ + 1 > bytes_.size()) {
      return SetError(error, "truncated sucds payload while reading u8");
    }
    *value = bytes_[offset_++];
    return true;
  }

  bool ReadBool(bool *value, std::string *error) {
    uint8_t byte = 0;
    if (!ReadU8(&byte, error)) {
      return false;
    }
    *value = byte != 0;
    return true;
  }

  bool ReadU16(uint16_t *value, std::string *error) {
    if (offset_ + 2 > bytes_.size()) {
      return SetError(error, "truncated sucds payload while reading u16");
    }
    *value = ReadLe16(bytes_.data() + offset_);
    offset_ += 2;
    return true;
  }

  bool ReadU64(uint64_t *value, std::string *error) {
    if (offset_ + 8 > bytes_.size()) {
      return SetError(error, "truncated sucds payload while reading u64");
    }
    *value = ReadLe64(bytes_.data() + offset_);
    offset_ += 8;
    return true;
  }

  bool SkipBytes(uint64_t nbytes64, std::string *error) {
    size_t nbytes = 0;
    // Compare as (nbytes > remaining) rather than (offset_ + nbytes > size) so
    // a huge nbytes cannot wrap the addition and pass. offset_ <= size always.
    if (!CheckedSize(nbytes64, &nbytes) ||
        nbytes > bytes_.size() - offset_) {
      return SetError(error, "truncated sucds payload while skipping bytes");
    }
    offset_ += nbytes;
    return true;
  }

  bool ReadBitVector(BitVectorLite *bv, std::string *error) {
    if (!ReadVecU64(&bv->words, error) || !ReadU64(&bv->len, error)) {
      return false;
    }
    // Overflow-safe ceil-div: (len + 63) / 64 wraps when len is within 63 of
    // 2^64, which would let a tiny word vector satisfy the check below and
    // leave Bit() indexing past the vector.
    const uint64_t expected_words =
        bv->len / 64 + (bv->len % 64 != 0 ? 1 : 0);
    if (expected_words > bv->words.size()) {
      return SetError(error,
                      "sucds BitVector has fewer words than bit length needs");
    }
    return true;
  }

  bool ReadDarrayBitVector(BitVectorLite *high_bits, std::string *error) {
    return ReadBitVector(high_bits, error) && SkipDarrayIndex(error) &&
           SkipOptionalDarrayIndex(error) &&
           SkipOptionalRank9SelIndex(error);
  }

  bool AtEnd() const { return offset_ == bytes_.size(); }

 private:
  bool ReadVecU64(std::vector<uint64_t> *values, std::string *error) {
    uint64_t len64 = 0;
    if (!ReadU64(&len64, error)) {
      return false;
    }
    size_t len = 0;
    if (!CheckedSize(len64, &len)) {
      return SetError(error, "sucds Vec length exceeds platform size");
    }
    values->clear();
    // Cap the reservation to what the remaining payload could actually hold so
    // a bogus length field cannot trigger a huge allocation before the read
    // loop detects the truncation.
    const size_t max_possible = (bytes_.size() - offset_) / 8;
    values->reserve(std::min(len, max_possible));
    for (size_t i = 0; i < len; ++i) {
      uint64_t value = 0;
      if (!ReadU64(&value, error)) {
        return false;
      }
      values->push_back(value);
    }
    return true;
  }

  // Skip a sucds Vec of `len` elements of `element_size` bytes each, guarding
  // the len * element_size multiply against 64-bit overflow (a wrapped product
  // would skip too few bytes and desynchronize the parser).
  bool SkipVecElements(uint64_t element_size, std::string *error) {
    uint64_t len = 0;
    if (!ReadU64(&len, error)) {
      return false;
    }
    if (element_size != 0 &&
        len > std::numeric_limits<uint64_t>::max() / element_size) {
      return SetError(error, "sucds Vec byte length overflows 64 bits");
    }
    return SkipBytes(len * element_size, error);
  }

  bool SkipVecI64(std::string *error) { return SkipVecElements(8, error); }

  bool SkipVecU16(std::string *error) { return SkipVecElements(2, error); }

  bool SkipVecU64(std::string *error) { return SkipVecElements(8, error); }

  bool SkipOptionalVecU64(std::string *error) {
    bool present = false;
    if (!ReadBool(&present, error)) {
      return false;
    }
    return present ? SkipVecU64(error) : true;
  }

  bool SkipDarrayIndex(std::string *error) {
    uint64_t unused = 0;
    bool unused_bool = false;
    return SkipVecI64(error) && SkipVecU16(error) && SkipVecU64(error) &&
           ReadU64(&unused, error) && ReadBool(&unused_bool, error);
  }

  bool SkipOptionalDarrayIndex(std::string *error) {
    bool present = false;
    if (!ReadBool(&present, error)) {
      return false;
    }
    return present ? SkipDarrayIndex(error) : true;
  }

  bool SkipRank9SelIndex(std::string *error) {
    uint64_t unused = 0;
    return ReadU64(&unused, error) && SkipVecU64(error) &&
           SkipOptionalVecU64(error) && SkipOptionalVecU64(error);
  }

  bool SkipOptionalRank9SelIndex(std::string *error) {
    bool present = false;
    if (!ReadBool(&present, error)) {
      return false;
    }
    return present ? SkipRank9SelIndex(error) : true;
  }

  const std::vector<uint8_t> &bytes_;
  size_t offset_ = 0;
};

bool DecodeEliasFanoPositions(const std::vector<uint8_t> &bytes,
                              std::vector<uint64_t> *positions,
                              std::string *error) {
  positions->clear();
  if (bytes.empty()) {
    return true;
  }

  SliceReader reader(bytes);
  BitVectorLite high_bits;
  BitVectorLite low_bits;
  uint64_t low_len = 0;
  uint64_t universe = 0;

  if (!reader.ReadDarrayBitVector(&high_bits, error) ||
      !reader.ReadBitVector(&low_bits, error) ||
      !reader.ReadU64(&low_len, error) ||
      !reader.ReadU64(&universe, error)) {
    return false;
  }
  if (low_len > 63) {
    return SetError(error,
                    "unsupported sucds Elias-Fano low-bit length > 63");
  }

  uint64_t ordinal = 0;
  for (uint64_t pos = 0; pos < high_bits.len; ++pos) {
    if (!high_bits.Bit(pos)) {
      continue;
    }
    if (low_len != 0 && ordinal * low_len + low_len > low_bits.len) {
      return SetError(error,
                      "sucds Elias-Fano low-bit vector is too short");
    }
    const uint64_t high = pos - ordinal;
    const uint64_t low = low_bits.Bits(ordinal * low_len, low_len);
    const uint64_t value = (high << low_len) | low;
    if (value >= universe) {
      return SetError(error, "decoded Elias-Fano value exceeds universe");
    }
    positions->push_back(value);
    ++ordinal;
  }

  if (!reader.AtEnd()) {
    return SetError(error, "trailing bytes after sucds Elias-Fano payload");
  }
  return true;
}

bool ReadExact(std::ifstream &stream, uint8_t *dest, size_t nbytes,
               const std::string &description, std::string *error) {
  stream.read(reinterpret_cast<char *>(dest),
              static_cast<std::streamsize>(nbytes));
  if (stream.gcount() != static_cast<std::streamsize>(nbytes)) {
    return SetError(error,
                    "truncated CBQ stream while reading " + description);
  }
  return true;
}

bool ReadVector(std::ifstream &stream, uint64_t nbytes64,
                std::vector<uint8_t> *dest, const std::string &description,
                std::string *error) {
  size_t nbytes = 0;
  if (!CheckedSize(nbytes64, &nbytes)) {
    return SetError(error,
                    "CBQ column length exceeds platform size for " +
                        description);
  }
  if (nbytes64 > kMaxColumnBytes) {
    return SetError(error,
                    "CBQ compressed column " + description +
                        " declares an implausibly large size");
  }
  // Reject a column that claims more bytes than the file actually has left, so
  // a malformed length cannot force a huge allocation before the read detects
  // the truncation.
  const std::streampos cur = stream.tellg();
  if (cur >= 0) {
    stream.seekg(0, std::ios::end);
    const std::streampos end = stream.tellg();
    stream.seekg(cur);
    if (!stream.good()) {
      return SetError(error, "CBQ stream seek failed for " + description);
    }
    if (end >= cur && nbytes64 > static_cast<uint64_t>(end - cur)) {
      return SetError(error, "CBQ column " + description +
                                 " declares more bytes than remain in file");
    }
  }
  dest->assign(nbytes, 0);
  if (nbytes == 0) {
    return true;
  }
  return ReadExact(stream, dest->data(), nbytes, description, error);
}

bool BytesToU64Vector(const std::vector<uint8_t> &bytes,
                      std::vector<uint64_t> *values,
                      const std::string &description, std::string *error) {
  if (bytes.size() % 8 != 0) {
    return SetError(error, "CBQ column " + description +
                               " byte length is not a multiple of 8");
  }
  values->clear();
  values->reserve(bytes.size() / 8);
  for (size_t offset = 0; offset < bytes.size(); offset += 8) {
    values->push_back(ReadLe64(bytes.data() + offset));
  }
  return true;
}

struct DecodedBlock {
  uint64_t num_records = 0;
  uint64_t num_sequences = 0;
  std::vector<uint64_t> seq_lengths;
  std::vector<uint64_t> seq_offsets;
  std::vector<uint64_t> header_lengths;
  std::vector<uint64_t> header_offsets;
  std::vector<uint8_t> seq_word_bytes;
  std::vector<uint64_t> n_positions;
  std::string headers;
  std::string qualities;
  size_t record_index = 0;
};

bool MakeOffsets(const std::vector<uint64_t> &lengths, uint64_t expected_total,
                 std::vector<uint64_t> *offsets,
                 const std::string &description, std::string *error) {
  offsets->assign(lengths.size() + 1, 0);
  uint64_t total = 0;
  for (size_t i = 0; i < lengths.size(); ++i) {
    // Detect 64-bit wrap so the prefix sums stay strictly monotonic. Without
    // this an attacker could pick lengths whose true sum overflows back to
    // expected_total, defeating the sum check below and yielding offsets that
    // point past the backing buffer.
    if (total > std::numeric_limits<uint64_t>::max() - lengths[i]) {
      return SetError(error, "CBQ " + description + " lengths overflow 64 bits");
    }
    total += lengths[i];
    (*offsets)[i + 1] = total;
  }
  if (total != expected_total) {
    std::ostringstream msg;
    msg << "CBQ " << description << " lengths sum to " << total
        << ", expected " << expected_total;
    return SetError(error, msg.str());
  }
  return true;
}

enum class BlockLoadStatus { kBlock, kEnd, kError };

struct CbqBlockBacking {
  DecodedBlock block;
  std::vector<std::string> synthetic_read_names;
  std::vector<CbqSegmentView> segment_views;
  std::vector<CbqReadView> read_views;
};

bool PackedBaseIsN(const CbqPackedSequenceView &packed,
                   uint64_t global_offset) {
  if (packed.n_positions == nullptr || packed.n_positions_count == 0) {
    return false;
  }
  const uint64_t *begin = packed.n_positions;
  const uint64_t *end = begin + packed.n_positions_count;
  return std::binary_search(begin, end, global_offset);
}

unsigned char PackedBaseNumber(const CbqPackedSequenceView &packed,
                               size_t index) {
  const uint64_t global_offset = packed.base_offset + index;
  if (PackedBaseIsN(packed, global_offset)) {
    return 4;
  }
  const uint64_t word_index = global_offset / 32;
  const uint64_t offset_in_word = global_offset % 32;
  const size_t byte_offset = static_cast<size_t>(word_index * 8);
  if (packed.words == nullptr || byte_offset + 8 > packed.word_bytes) {
    return 4;
  }
  const uint64_t word = ReadLe64(packed.words + byte_offset);
  return static_cast<unsigned char>((word >> (offset_in_word * 2)) & 0x3ULL);
}

char NumberToAscii(unsigned char base) {
  static const char lookup[5] = {'A', 'C', 'G', 'T', 'N'};
  return lookup[base < 5 ? base : 4];
}

struct CbqPackedByteDecode {
  char ascii[4];
};

struct CbqPackedBytePairDecode {
  char ascii[8];
};

struct CbqPackedByteReverseComplementDecode {
  char reverse_complement[4];
};

struct CbqPackedBytePairReverseComplementDecode {
  char reverse_complement[8];
};

const std::array<char, 256> &AsciiComplementLookup() {
  static const std::array<char, 256> lookup = [] {
    std::array<char, 256> table = {};
    table.fill('N');
    table[static_cast<unsigned char>('A')] = 'T';
    table[static_cast<unsigned char>('C')] = 'G';
    table[static_cast<unsigned char>('G')] = 'C';
    table[static_cast<unsigned char>('T')] = 'A';
    table[static_cast<unsigned char>('N')] = 'N';
    table[static_cast<unsigned char>('a')] = 'T';
    table[static_cast<unsigned char>('c')] = 'G';
    table[static_cast<unsigned char>('g')] = 'C';
    table[static_cast<unsigned char>('t')] = 'A';
    table[static_cast<unsigned char>('n')] = 'N';
    return table;
  }();
  return lookup;
}

const std::array<CbqPackedByteDecode, 256> &CbqPackedAsciiLookup() {
  static const std::array<CbqPackedByteDecode, 256> lookup = [] {
    std::array<CbqPackedByteDecode, 256> table = {};
    static const char ascii_lookup[4] = {'A', 'C', 'G', 'T'};
    for (unsigned byte = 0; byte < table.size(); ++byte) {
      for (unsigned base = 0; base < 4; ++base) {
        const unsigned code = (byte >> (base * 2)) & 0x3U;
        table[byte].ascii[base] = ascii_lookup[code];
      }
    }
    return table;
  }();
  return lookup;
}

const std::array<CbqPackedByteReverseComplementDecode, 256> &
CbqPackedAsciiReverseComplementLookup() {
  static const std::array<CbqPackedByteReverseComplementDecode, 256> lookup =
      [] {
        std::array<CbqPackedByteReverseComplementDecode, 256> table = {};
        static const char complement_lookup[4] = {'T', 'G', 'C', 'A'};
        for (unsigned byte = 0; byte < table.size(); ++byte) {
          for (unsigned base = 0; base < 4; ++base) {
            const unsigned reverse_code =
                (byte >> ((3 - base) * 2)) & 0x3U;
            table[byte].reverse_complement[base] =
                complement_lookup[reverse_code];
          }
        }
        return table;
      }();
  return lookup;
}

const std::array<CbqPackedBytePairDecode, 65536> &CbqPackedAscii16Lookup() {
  static const std::array<CbqPackedBytePairDecode, 65536> lookup = [] {
    std::array<CbqPackedBytePairDecode, 65536> table = {};
    static const char ascii_lookup[4] = {'A', 'C', 'G', 'T'};
    for (uint32_t word = 0; word < table.size(); ++word) {
      for (unsigned base = 0; base < 8; ++base) {
        const unsigned code = (word >> (base * 2)) & 0x3U;
        table[word].ascii[base] = ascii_lookup[code];
      }
    }
    return table;
  }();
  return lookup;
}

const std::array<CbqPackedBytePairReverseComplementDecode, 65536> &
CbqPackedAscii16ReverseComplementLookup() {
  static const std::array<CbqPackedBytePairReverseComplementDecode, 65536>
      lookup = [] {
        std::array<CbqPackedBytePairReverseComplementDecode, 65536> table = {};
        static const char complement_lookup[4] = {'T', 'G', 'C', 'A'};
        for (uint32_t word = 0; word < table.size(); ++word) {
          for (unsigned base = 0; base < 8; ++base) {
            const unsigned reverse_code =
                (word >> ((7 - base) * 2)) & 0x3U;
            table[word].reverse_complement[base] =
                complement_lookup[reverse_code];
          }
        }
        return table;
      }();
  return lookup;
}

}  // namespace

size_t CbqSegmentSequenceLength(const CbqSegmentView &segment) {
  if (segment.packed_sequence.available) {
    return segment.packed_sequence.length;
  }
  return segment.sequence.size;
}

char CbqSegmentBaseAscii(const CbqSegmentView &segment, size_t index) {
  if (segment.sequence.data != nullptr && index < segment.sequence.size) {
    return segment.sequence.data[index];
  }
  if (segment.packed_sequence.available &&
      index < segment.packed_sequence.length) {
    return NumberToAscii(PackedBaseNumber(segment.packed_sequence, index));
  }
  return 'N';
}

void MaterializeCbqSegmentSequence(const CbqSegmentView &segment,
                                   std::string *sequence) {
  if (sequence == nullptr) {
    return;
  }
  const size_t length = CbqSegmentSequenceLength(segment);
  sequence->resize(length);
  for (size_t i = 0; i < length; ++i) {
    (*sequence)[i] = CbqSegmentBaseAscii(segment, i);
  }
}

// Mirrors STAR-suite's CbqInputModule direct materializer so Chromap can fill
// its existing SequenceBatch buffers without an intermediate sequence string.
bool MaterializeCbqSegmentSequenceToBuffer(const CbqSegmentView &segment,
                                           char *dest, size_t capacity,
                                           size_t *length_out,
                                           std::string *error) {
  if (dest == nullptr) {
    return SetError(error, "CBQ sequence materialization destination is null");
  }
  const size_t length = CbqSegmentSequenceLength(segment);
  if (length >= capacity) {
    std::ostringstream msg;
    msg << "CBQ sequence length " << length
        << " exceeds destination capacity " << capacity;
    return SetError(error, msg.str());
  }

  if (segment.sequence.data != nullptr && segment.sequence.size == length) {
    if (length != 0) {
      std::memcpy(dest, segment.sequence.data, length);
    }
    dest[length] = '\0';
    if (length_out != nullptr) {
      *length_out = length;
    }
    return true;
  }

  if (!segment.packed_sequence.available) {
    std::fill(dest, dest + length, 'N');
    dest[length] = '\0';
    if (length_out != nullptr) {
      *length_out = length;
    }
    return true;
  }

  const CbqPackedSequenceView &packed = segment.packed_sequence;
  const auto &lookup8 = CbqPackedAsciiLookup();
  const auto &lookup16 = CbqPackedAscii16Lookup();
  size_t i = 0;
  while (i < length) {
    const uint64_t global_offset = packed.base_offset + static_cast<uint64_t>(i);
    const size_t byte_offset = static_cast<size_t>(global_offset >> 2);
    const unsigned base_in_byte = static_cast<unsigned>(global_offset & 0x3ULL);
    const size_t take =
        std::min<size_t>(static_cast<size_t>(4U - base_in_byte), length - i);

    if (packed.words == nullptr || byte_offset >= packed.word_bytes) {
      std::fill(dest + i, dest + i + take, 'N');
    } else if (base_in_byte == 0 && take == 4) {
      size_t local_byte_offset = byte_offset;
      while (i + 8 <= length && local_byte_offset + 1 < packed.word_bytes) {
        const uint16_t packed_word =
            static_cast<uint16_t>(packed.words[local_byte_offset]) |
            static_cast<uint16_t>(
                static_cast<uint16_t>(packed.words[local_byte_offset + 1])
                << 8);
        const CbqPackedBytePairDecode &decoded = lookup16[packed_word];
        std::memcpy(dest + i, decoded.ascii, 8);
        i += 8;
        local_byte_offset += 2;
      }
      if (i >= length) {
        continue;
      }
      if (local_byte_offset >= packed.word_bytes) {
        std::fill(dest + i, dest + length, 'N');
        i = length;
        continue;
      }
      const CbqPackedByteDecode &decoded =
          lookup8[packed.words[local_byte_offset]];
      const size_t tail = std::min<size_t>(4, length - i);
      if (tail == 4) {
        std::memcpy(dest + i, decoded.ascii, 4);
      } else {
        for (size_t j = 0; j < tail; ++j) {
          dest[i + j] = decoded.ascii[j];
        }
      }
      i += tail;
      continue;
    } else {
      const CbqPackedByteDecode &decoded = lookup8[packed.words[byte_offset]];
      for (size_t j = 0; j < take; ++j) {
        dest[i + j] = decoded.ascii[base_in_byte + j];
      }
    }

    i += take;
  }

  if (packed.n_positions != nullptr) {
    for (size_t inpos = 0; inpos < packed.n_positions_count; ++inpos) {
      const uint64_t npos = packed.n_positions[inpos];
      if (npos < packed.base_offset) {
        continue;
      }
      const uint64_t local = npos - packed.base_offset;
      if (local >= length) {
        break;
      }
      dest[local] = 'N';
    }
  }

  dest[length] = '\0';
  if (length_out != nullptr) {
    *length_out = length;
  }
  return true;
}

bool MaterializeCbqSegmentSequenceAndReverseComplementToBuffers(
    const CbqSegmentView &segment, char *dest, size_t capacity,
    char *reverse_complement_dest, size_t reverse_complement_capacity,
    size_t *length_out, std::string *error) {
  if (dest == nullptr) {
    return SetError(error, "CBQ sequence materialization destination is null");
  }
  const size_t length = CbqSegmentSequenceLength(segment);
  if (reverse_complement_dest == nullptr && length != 0) {
    return SetError(error,
                    "CBQ reverse-complement materialization destination is null");
  }
  if (length >= capacity) {
    std::ostringstream msg;
    msg << "CBQ sequence length " << length
        << " exceeds destination capacity " << capacity;
    return SetError(error, msg.str());
  }
  if (length > reverse_complement_capacity) {
    std::ostringstream msg;
    msg << "CBQ sequence length " << length
        << " exceeds reverse-complement destination capacity "
        << reverse_complement_capacity;
    return SetError(error, msg.str());
  }

  const auto &complement = AsciiComplementLookup();
  if (segment.sequence.data != nullptr && segment.sequence.size == length) {
    if (length != 0) {
      std::memcpy(dest, segment.sequence.data, length);
    }
    for (size_t i = 0; i < length; ++i) {
      const unsigned char base =
          static_cast<unsigned char>(segment.sequence.data[i]);
      reverse_complement_dest[length - 1 - i] = complement[base];
    }
    dest[length] = '\0';
    if (length_out != nullptr) {
      *length_out = length;
    }
    return true;
  }

  if (!segment.packed_sequence.available) {
    std::fill(dest, dest + length, 'N');
    std::fill(reverse_complement_dest, reverse_complement_dest + length, 'N');
    dest[length] = '\0';
    if (length_out != nullptr) {
      *length_out = length;
    }
    return true;
  }

  const CbqPackedSequenceView &packed = segment.packed_sequence;
  const auto &lookup8 = CbqPackedAsciiLookup();
  const auto &lookup8_rc = CbqPackedAsciiReverseComplementLookup();
  const auto &lookup16 = CbqPackedAscii16Lookup();
  const auto &lookup16_rc = CbqPackedAscii16ReverseComplementLookup();
  size_t i = 0;
  while (i < length) {
    const uint64_t global_offset = packed.base_offset + static_cast<uint64_t>(i);
    const size_t byte_offset = static_cast<size_t>(global_offset >> 2);
    const unsigned base_in_byte = static_cast<unsigned>(global_offset & 0x3ULL);
    const size_t take =
        std::min<size_t>(static_cast<size_t>(4U - base_in_byte), length - i);

    if (packed.words == nullptr || byte_offset >= packed.word_bytes) {
      std::fill(dest + i, dest + i + take, 'N');
      std::fill(reverse_complement_dest + (length - i - take),
                reverse_complement_dest + (length - i), 'N');
    } else if (base_in_byte == 0 && take == 4) {
      size_t local_byte_offset = byte_offset;
      while (i + 8 <= length && local_byte_offset + 1 < packed.word_bytes) {
        const uint16_t packed_word =
            static_cast<uint16_t>(packed.words[local_byte_offset]) |
            static_cast<uint16_t>(
                static_cast<uint16_t>(packed.words[local_byte_offset + 1])
                << 8);
        const CbqPackedBytePairDecode &decoded = lookup16[packed_word];
        const CbqPackedBytePairReverseComplementDecode &decoded_rc =
            lookup16_rc[packed_word];
        std::memcpy(dest + i, decoded.ascii, 8);
        std::memcpy(reverse_complement_dest + (length - i - 8),
                    decoded_rc.reverse_complement, 8);
        i += 8;
        local_byte_offset += 2;
      }
      if (i >= length) {
        continue;
      }
      if (local_byte_offset >= packed.word_bytes) {
        std::fill(dest + i, dest + length, 'N');
        std::fill(reverse_complement_dest,
                  reverse_complement_dest + (length - i), 'N');
        i = length;
        continue;
      }
      const CbqPackedByteDecode &decoded =
          lookup8[packed.words[local_byte_offset]];
      const CbqPackedByteReverseComplementDecode &decoded_rc =
          lookup8_rc[packed.words[local_byte_offset]];
      const size_t tail = std::min<size_t>(4, length - i);
      if (tail == 4) {
        std::memcpy(dest + i, decoded.ascii, 4);
        std::memcpy(reverse_complement_dest + (length - i - 4),
                    decoded_rc.reverse_complement, 4);
      } else {
        for (size_t j = 0; j < tail; ++j) {
          dest[i + j] = decoded.ascii[j];
          reverse_complement_dest[length - 1 - i - j] =
              complement[static_cast<unsigned char>(decoded.ascii[j])];
        }
      }
      i += tail;
      continue;
    } else {
      const CbqPackedByteDecode &decoded = lookup8[packed.words[byte_offset]];
      for (size_t j = 0; j < take; ++j) {
        const char base = decoded.ascii[base_in_byte + j];
        dest[i + j] = base;
        reverse_complement_dest[length - 1 - i - j] =
            complement[static_cast<unsigned char>(base)];
      }
    }

    i += take;
  }

  if (packed.n_positions != nullptr) {
    for (size_t inpos = 0; inpos < packed.n_positions_count; ++inpos) {
      const uint64_t npos = packed.n_positions[inpos];
      if (npos < packed.base_offset) {
        continue;
      }
      const uint64_t local = npos - packed.base_offset;
      if (local >= length) {
        break;
      }
      dest[local] = 'N';
      reverse_complement_dest[length - 1 - local] = 'N';
    }
  }

  dest[length] = '\0';
  if (length_out != nullptr) {
    *length_out = length;
  }
  return true;
}

struct CbqLaneReader::Impl {
  std::ifstream stream;
  uint64_t read_ordinal = 0;
  uint64_t lane_record_index = 0;
  uint64_t batch_first_record = 0;
  uint64_t current_lane_records = 0;
  uint64_t range_first_record = 0;
  uint64_t range_end_record = 0;
  bool opened = false;
  bool exhausted = false;
  bool range_mode = false;
  FileHeaderFields file_header;
  std::vector<CbqBlockIndexEntry> block_index;
  std::shared_ptr<CbqBlockBacking> batch;

  bool OpenStreamAndReadHeader(const std::string &path, uint32_t mate_count,
                               std::string *error) {
    stream.open(path.c_str(), std::ios::binary);
    if (!stream.good()) {
      return SetError(error, "could not open CBQ file: " + path);
    }

    std::array<uint8_t, 64> header_bytes{};
    if (!ReadExact(stream, header_bytes.data(), header_bytes.size(),
                   "file header", error) ||
        !ParseFileHeader(header_bytes, &file_header, error)) {
      return false;
    }

    const bool file_paired = file_header.IsPaired();
    if ((mate_count == 2) != file_paired) {
      std::ostringstream msg;
      msg << "CBQ mate-count mismatch for " << path << ": file paired="
          << (file_paired ? "true" : "false")
          << " but requested mate_count=" << mate_count;
      return SetError(error, msg.str());
    }
    return true;
  }

  bool ReadCurrentLaneIndex(const std::string &path, std::string *error) {
    stream.clear();
    stream.seekg(0, std::ios::end);
    const std::streampos end_pos = stream.tellg();
    if (end_pos < 0) {
      return SetError(error, "could not seek CBQ file for index: " + path);
    }
    const uint64_t file_size = static_cast<uint64_t>(end_pos);
    const uint64_t index_header_size = 24U;
    const uint64_t index_footer_size = 16U;
    if (file_size < 64U + index_header_size + index_footer_size) {
      return SetError(error,
                      "CBQ range mode requires a CBQINDEX footer: " + path);
    }

    std::array<uint8_t, 16> footer{};
    stream.seekg(static_cast<std::streamoff>(file_size - index_footer_size),
                 std::ios::beg);
    if (!ReadExact(stream, footer.data(), footer.size(), "index footer",
                   error)) {
      return false;
    }
    if (std::memcmp(footer.data() + 8, "CBQINDEX", 8) != 0) {
      return SetError(error,
                      "CBQ range mode requires a CBQINDEX footer: " + path);
    }

    const uint64_t z_index_size = ReadLe64(footer.data());
    if (file_size < 64U + index_header_size + index_footer_size +
                        z_index_size) {
      return SetError(error,
                      "CBQ index footer points before file header: " + path);
    }
    const uint64_t index_header_offset =
        file_size - index_footer_size - z_index_size - index_header_size;

    std::array<uint8_t, 24> index_header{};
    stream.seekg(static_cast<std::streamoff>(index_header_offset),
                 std::ios::beg);
    if (!ReadExact(stream, index_header.data(), index_header.size(),
                   "index header", error)) {
      return false;
    }
    if (std::memcmp(index_header.data(), "CBQINDEX", 8) != 0) {
      return SetError(error, "CBQ index header magic mismatch: " + path);
    }

    const uint64_t index_size = ReadLe64(index_header.data() + 8);
    const uint64_t compressed_index_size =
        ReadLe64(index_header.data() + 16);
    if (compressed_index_size != z_index_size) {
      return SetError(
          error, "CBQ index footer/header compressed-size mismatch: " + path);
    }
    if (index_size % 16U != 0) {
      return SetError(error,
                      "CBQ index size is not a multiple of 16 bytes: " + path);
    }
    if (index_size > kMaxIndexBytes) {
      return SetError(error, "CBQ index declares an implausibly large size: " +
                                 path);
    }
    const uint64_t index_entry_count = index_size / 16U;
    const uint64_t max_physical_blocks = (index_header_offset - 64U) / 96U;
    if (index_entry_count > max_physical_blocks) {
      return SetError(error,
                      "CBQ index has more entries than possible blocks: " +
                          path);
    }
    size_t index_size_size_t = 0;
    if (!CheckedSize(index_size, &index_size_size_t)) {
      return SetError(error, "CBQ index size exceeds platform size: " + path);
    }

    std::vector<uint8_t> z_index;
    if (!ReadVector(stream, compressed_index_size, &z_index, "z_index",
                    error)) {
      return false;
    }

    std::vector<uint8_t> index_bytes;
    if (!GetZstdRuntime().Decompress(z_index, index_size_size_t, &index_bytes,
                                     "index", error)) {
      return false;
    }
    if (index_bytes.size() != index_size_size_t) {
      return SetError(error,
                      "CBQ index decompressed to an unexpected size: " + path);
    }

    block_index.clear();
    block_index.reserve(index_bytes.size() / 16U);
    uint64_t previous_cumulative_records = 0;
    for (size_t offset = 0; offset < index_bytes.size(); offset += 16U) {
      CbqBlockIndexEntry entry;
      entry.offset = ReadLe64(index_bytes.data() + offset);
      entry.cumulative_records = ReadLe64(index_bytes.data() + offset + 8U);
      if (entry.offset < 64U || entry.offset >= index_header_offset) {
        return SetError(error,
                        "CBQ index contains an invalid block offset: " + path);
      }
      if (entry.cumulative_records < previous_cumulative_records) {
        return SetError(error,
                        "CBQ index cumulative records are not monotonic: " +
                            path);
      }
      previous_cumulative_records = entry.cumulative_records;
      block_index.push_back(entry);
    }

    current_lane_records =
        block_index.empty() ? 0U : block_index.back().cumulative_records;
    stream.clear();
    return true;
  }

  bool SeekToRangeStart(uint64_t first_record, std::string *error) {
    if (first_record >= current_lane_records || block_index.empty()) {
      lane_record_index = current_lane_records;
      read_ordinal = current_lane_records;
      batch.reset();
      return true;
    }

    std::vector<CbqBlockIndexEntry>::const_iterator it =
        std::upper_bound(block_index.begin(), block_index.end(), first_record,
                         [](uint64_t value,
                            const CbqBlockIndexEntry &entry) {
                           return value < entry.cumulative_records;
                         });
    if (it == block_index.end()) {
      lane_record_index = current_lane_records;
      read_ordinal = current_lane_records;
      batch.reset();
      return true;
    }

    const size_t block_index_position = static_cast<size_t>(
        std::distance(block_index.cbegin(), it));
    const uint64_t block_first_record =
        block_index_position == 0
            ? 0U
            : block_index[block_index_position - 1U].cumulative_records;

    stream.clear();
    stream.seekg(static_cast<std::streamoff>(it->offset), std::ios::beg);
    if (!stream.good()) {
      return SetError(error, "could not seek CBQ stream to indexed block");
    }
    lane_record_index = block_first_record;
    read_ordinal = block_first_record;
    batch.reset();
    return true;
  }

  bool BuildBatchViews(uint32_t mate_count, std::string *error) {
    if (!batch) {
      return SetError(error, "CBQ batch storage is not initialized");
    }
    DecodedBlock &block = batch->block;
    std::vector<std::string> &synthetic_read_names =
        batch->synthetic_read_names;
    std::vector<CbqSegmentView> &segment_views = batch->segment_views;
    std::vector<CbqReadView> &read_views = batch->read_views;

    size_t num_records = 0;
    if (!CheckedSize(block.num_records, &num_records)) {
      return SetError(error, "CBQ block record count exceeds platform size");
    }
    if (mate_count == 0 ||
        num_records > std::numeric_limits<size_t>::max() / mate_count) {
      return SetError(error, "CBQ block segment count exceeds platform size");
    }

    const size_t segment_count = num_records * mate_count;
    read_views.assign(num_records, CbqReadView());
    segment_views.assign(segment_count, CbqSegmentView());
    synthetic_read_names.clear();
    if (!file_header.HasHeaders() || block.headers.empty()) {
      synthetic_read_names.reserve(num_records);
    }

    for (size_t irecord = 0; irecord < num_records; ++irecord) {
      const uint64_t first_sequence_index =
          static_cast<uint64_t>(irecord) * mate_count;

      ParsedHeaderSpan parsed;
      if (file_header.HasHeaders() && !block.headers.empty()) {
        const size_t begin =
            static_cast<size_t>(block.header_offsets[first_sequence_index]);
        const size_t end =
            static_cast<size_t>(block.header_offsets[first_sequence_index + 1]);
        if (begin > end || end > block.headers.size()) {
          return SetError(error, "CBQ header offset is out of range");
        }
        parsed = ParseHeaderPayloadSpan(block.headers, begin, end - begin,
                                        file_header.HasQualities());
      } else {
        synthetic_read_names.push_back(
            std::to_string(lane_record_index + irecord + 1));
        parsed.read_name = MakeSpan(synthetic_read_names.back(), 0,
                                    synthetic_read_names.back().size());
        parsed.read_name_extra = MakeEmptySpan();
        parsed.read_filter = 'N';
      }

      CbqReadView &record = read_views[irecord];
      record.read_name = parsed.read_name;
      record.read_name_extra = parsed.read_name_extra;
      record.read_ordinal = read_ordinal + irecord + 1;
      record.lane_read_ordinal = lane_record_index + irecord + 1;
      record.read_filter = parsed.read_filter;
      record.segment_count = mate_count;
      record.segments = segment_views.data() + (irecord * mate_count);

      for (uint32_t isegment = 0; isegment < mate_count; ++isegment) {
        const uint64_t seq_index = first_sequence_index + isegment;
        const size_t begin =
            static_cast<size_t>(block.seq_offsets[seq_index]);
        const size_t end =
            static_cast<size_t>(block.seq_offsets[seq_index + 1]);
        if (begin > end || end > block.qualities.size()) {
          return SetError(error, "CBQ sequence offset is out of range");
        }
        const size_t length = end - begin;
        if (length > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
          return SetError(error, "CBQ sequence length exceeds uint32_t");
        }

        CbqSegmentView &segment =
            segment_views[irecord * mate_count + isegment];
        segment.source_index = isegment;
        segment.sequence = MakeEmptySpan();
        segment.quality = MakeSpan(block.qualities, begin, length);
        segment.packed_sequence.words =
            block.seq_word_bytes.empty() ? nullptr : block.seq_word_bytes.data();
        segment.packed_sequence.word_bytes = block.seq_word_bytes.size();
        segment.packed_sequence.base_offset = begin;
        segment.packed_sequence.length = static_cast<uint32_t>(length);
        const auto n_begin =
            std::lower_bound(block.n_positions.begin(),
                             block.n_positions.end(), begin);
        const auto n_end =
            std::lower_bound(block.n_positions.begin(),
                             block.n_positions.end(), end);
        segment.packed_sequence.n_positions_count =
            static_cast<size_t>(std::distance(n_begin, n_end));
        segment.packed_sequence.n_positions =
            segment.packed_sequence.n_positions_count == 0 ? nullptr
                                                           : &(*n_begin);
        segment.packed_sequence.available = !block.seq_word_bytes.empty();
        segment.original_length = static_cast<uint32_t>(length);
        segment.has_quality = file_header.HasQualities();
      }
    }

    read_ordinal += static_cast<uint64_t>(num_records);
    lane_record_index += static_cast<uint64_t>(num_records);
    return true;
  }

  BlockLoadStatus LoadNextBlock(uint32_t mate_count, std::string *error) {
    batch.reset();

    std::array<uint8_t, 8> first{};
    stream.read(reinterpret_cast<char *>(first.data()),
                static_cast<std::streamsize>(first.size()));
    if (stream.gcount() == 0 && stream.eof()) {
      return BlockLoadStatus::kEnd;
    }
    if (stream.gcount() != static_cast<std::streamsize>(first.size())) {
      SetError(error, "truncated CBQ stream while reading block/index magic");
      return BlockLoadStatus::kError;
    }
    if (std::memcmp(first.data(), "CBQINDEX", 8) == 0) {
      return BlockLoadStatus::kEnd;
    }

    std::array<uint8_t, 96> header_bytes{};
    std::copy(first.begin(), first.end(), header_bytes.begin());
    if (!ReadExact(stream, header_bytes.data() + first.size(),
                   header_bytes.size() - first.size(), "block header",
                   error)) {
      return BlockLoadStatus::kError;
    }

    BlockHeaderFields header;
    if (!ParseBlockHeader(header_bytes, &header, error)) {
      return BlockLoadStatus::kError;
    }
    if (header.nuclen > kMaxBasesPerBlock ||
        header.num_sequences > kMaxSequencesPerBlock ||
        header.num_records > kMaxRecordsPerBlock ||
        header.len_nef > kMaxColumnBytes) {
      SetError(error, "CBQ block declares implausibly large dimensions");
      return BlockLoadStatus::kError;
    }
    // Bounded above, so num_records * mate_count (mate_count <= 2) cannot wrap.
    const uint64_t expected_sequences = header.num_records * mate_count;
    if (header.num_sequences != expected_sequences) {
      std::ostringstream msg;
      msg << "CBQ block has num_sequences=" << header.num_sequences
          << " but expected " << expected_sequences
          << " for num_records=" << header.num_records
          << " mate_count=" << mate_count;
      SetError(error, msg.str());
      return BlockLoadStatus::kError;
    }

    batch = std::make_shared<CbqBlockBacking>();
    DecodedBlock &block = batch->block;

    std::vector<uint8_t> z_seq_len;
    std::vector<uint8_t> z_header_len;
    std::vector<uint8_t> z_npos;
    std::vector<uint8_t> z_seq;
    std::vector<uint8_t> z_flags;
    std::vector<uint8_t> z_headers;
    std::vector<uint8_t> z_qual;
    if (!ReadVector(stream, header.len_z_seq_len, &z_seq_len, "z_seq_len",
                    error) ||
        !ReadVector(stream, header.len_z_header_len, &z_header_len,
                    "z_header_len", error) ||
        !ReadVector(stream, header.len_z_npos, &z_npos, "z_npos", error) ||
        !ReadVector(stream, header.len_z_seq, &z_seq, "z_seq", error) ||
        !ReadVector(stream, header.len_z_flags, &z_flags, "z_flags", error) ||
        !ReadVector(stream, header.len_z_headers, &z_headers, "z_headers",
                    error) ||
        !ReadVector(stream, header.len_z_qual, &z_qual, "z_qual", error)) {
      return BlockLoadStatus::kError;
    }

    size_t num_sequences = 0;
    size_t nuclen = 0;
    size_t ef_len = 0;
    if (!CheckedSize(header.num_sequences, &num_sequences) ||
        !CheckedSize(header.nuclen, &nuclen) ||
        !CheckedSize(header.len_nef, &ef_len)) {
      SetError(error, "CBQ block metadata exceeds platform size");
      return BlockLoadStatus::kError;
    }

    std::vector<uint8_t> seq_len_bytes;
    std::vector<uint8_t> header_len_bytes;
    std::vector<uint8_t> npos_bytes;
    std::vector<uint8_t> flags_bytes;
    std::vector<uint8_t> header_bytes_column;
    std::vector<uint8_t> qual_bytes;
    const size_t seq_len_bytes_expected = num_sequences * 8;
    const size_t seq_word_bytes_expected = ((nuclen + 31) / 32) * 8;

    if (!GetZstdRuntime().Decompress(z_seq_len, seq_len_bytes_expected,
                                     &seq_len_bytes, "seq_len", error) ||
        !BytesToU64Vector(seq_len_bytes, &block.seq_lengths, "seq_len",
                          error) ||
        !MakeOffsets(block.seq_lengths, header.nuclen, &block.seq_offsets,
                     "sequence", error)) {
      return BlockLoadStatus::kError;
    }

    if (file_header.HasHeaders()) {
      if (!GetZstdRuntime().Decompress(z_header_len, seq_len_bytes_expected,
                                       &header_len_bytes, "header_len",
                                       error) ||
          !BytesToU64Vector(header_len_bytes, &block.header_lengths,
                            "header_len", error)) {
        return BlockLoadStatus::kError;
      }
    } else {
      block.header_lengths.assign(num_sequences, 0);
    }

    if (!z_npos.empty() &&
        !GetZstdRuntime().Decompress(z_npos, ef_len, &npos_bytes, "npos",
                                     error)) {
      return BlockLoadStatus::kError;
    }
    if (!DecodeEliasFanoPositions(npos_bytes, &block.n_positions, error)) {
      return BlockLoadStatus::kError;
    }
    if (!GetZstdRuntime().Decompress(z_seq, seq_word_bytes_expected,
                                     &block.seq_word_bytes, "seq", error)) {
      return BlockLoadStatus::kError;
    }
    for (uint64_t npos : block.n_positions) {
      if (npos >= header.nuclen) {
        SetError(error, "CBQ N-position exceeds sequence length");
        return BlockLoadStatus::kError;
      }
    }

    if (file_header.HasFlags()) {
      const size_t flags_expected = static_cast<size_t>(header.num_records) * 8;
      if (!GetZstdRuntime().Decompress(z_flags, flags_expected, &flags_bytes,
                                       "flags", error)) {
        return BlockLoadStatus::kError;
      }
    }

    uint64_t total_header_len = 0;
    for (uint64_t length : block.header_lengths) {
      if (total_header_len > std::numeric_limits<uint64_t>::max() - length) {
        SetError(error, "CBQ header lengths overflow 64 bits");
        return BlockLoadStatus::kError;
      }
      total_header_len += length;
    }
    if (total_header_len > kMaxColumnBytes) {
      SetError(error, "CBQ header payload declares implausibly large size");
      return BlockLoadStatus::kError;
    }
    if (!MakeOffsets(block.header_lengths, total_header_len,
                     &block.header_offsets, "header", error)) {
      return BlockLoadStatus::kError;
    }
    if (file_header.HasHeaders()) {
      size_t total_header_size = 0;
      if (!CheckedSize(total_header_len, &total_header_size) ||
          !GetZstdRuntime().Decompress(z_headers, total_header_size,
                                       &header_bytes_column, "headers",
                                       error)) {
        SetError(error, "CBQ header payload exceeds platform size");
        return BlockLoadStatus::kError;
      }
      block.headers.assign(
          reinterpret_cast<const char *>(header_bytes_column.data()),
          header_bytes_column.size());
    }

    if (file_header.HasQualities()) {
      if (!GetZstdRuntime().Decompress(z_qual, nuclen, &qual_bytes, "qual",
                                       error)) {
        return BlockLoadStatus::kError;
      }
      block.qualities.assign(reinterpret_cast<const char *>(qual_bytes.data()),
                             qual_bytes.size());
    } else {
      block.qualities.assign(nuclen, 'A');
    }

    block.num_records = header.num_records;
    block.num_sequences = header.num_sequences;
    block.record_index = 0;
    batch_first_record = lane_record_index;
    if (!BuildBatchViews(mate_count, error)) {
      return BlockLoadStatus::kError;
    }
    return BlockLoadStatus::kBlock;
  }

  BlockLoadStatus LoadNextAvailableBlock(uint32_t mate_count,
                                         std::string *error) {
    if (range_mode && range_first_record >= range_end_record) {
      opened = false;
      exhausted = true;
      return BlockLoadStatus::kEnd;
    }
    if (batch &&
        batch->block.record_index < static_cast<size_t>(batch->block.num_records)) {
      return BlockLoadStatus::kBlock;
    }
    if (range_mode && lane_record_index >= range_end_record) {
      opened = false;
      exhausted = true;
      return BlockLoadStatus::kEnd;
    }
    for (;;) {
      const BlockLoadStatus status = LoadNextBlock(mate_count, error);
      if (status == BlockLoadStatus::kError ||
          status == BlockLoadStatus::kEnd) {
        if (status == BlockLoadStatus::kEnd) {
          opened = false;
          exhausted = true;
        }
        return status;
      }
      if (batch && batch->block.num_records > 0) {
        return BlockLoadStatus::kBlock;
      }
      if (range_mode && lane_record_index >= range_end_record) {
        opened = false;
        exhausted = true;
        return BlockLoadStatus::kEnd;
      }
    }
  }
};

CbqLaneReader::CbqLaneReader(const std::string &path, uint32_t mate_count)
    : path_(path), mate_count_(mate_count), impl_(new Impl()) {}

CbqLaneReader::~CbqLaneReader() = default;

bool CbqLaneReader::Open(std::string *error) {
  if (mate_count_ != 1 && mate_count_ != 2) {
    return SetError(error, "CBQ reader supports mate_count 1 or 2");
  }
  Close();
  if (!impl_->OpenStreamAndReadHeader(path_, mate_count_, error)) {
    Close();
    return false;
  }

  impl_->read_ordinal = 0;
  impl_->lane_record_index = 0;
  impl_->batch_first_record = 0;
  impl_->current_lane_records = 0;
  impl_->range_first_record = 0;
  impl_->range_end_record = 0;
  impl_->range_mode = false;
  impl_->block_index.clear();
  impl_->batch.reset();
  impl_->opened = true;
  impl_->exhausted = false;
  return true;
}

bool CbqLaneReader::OpenRange(uint64_t first_record, uint64_t record_count,
                              std::string *error) {
  if (mate_count_ != 1 && mate_count_ != 2) {
    return SetError(error, "CBQ reader supports mate_count 1 or 2");
  }
  Close();
  if (!impl_->OpenStreamAndReadHeader(path_, mate_count_, error)) {
    Close();
    return false;
  }
  if (!impl_->ReadCurrentLaneIndex(path_, error)) {
    Close();
    return false;
  }
  if (first_record > impl_->current_lane_records) {
    std::ostringstream msg;
    msg << "CBQ range starts at record " << first_record
        << " but lane has only " << impl_->current_lane_records
        << " records";
    Close();
    return SetError(error, msg.str());
  }

  uint64_t end_record = impl_->current_lane_records;
  if (record_count != std::numeric_limits<uint64_t>::max()) {
    const uint64_t remaining = impl_->current_lane_records - first_record;
    end_record =
        record_count <= remaining ? first_record + record_count
                                  : impl_->current_lane_records;
  }

  impl_->read_ordinal = 0;
  impl_->lane_record_index = 0;
  impl_->batch_first_record = 0;
  impl_->range_mode = true;
  impl_->range_first_record = first_record;
  impl_->range_end_record = end_record;
  impl_->batch.reset();
  if (!impl_->SeekToRangeStart(first_record, error)) {
    Close();
    return false;
  }
  impl_->opened = true;
  impl_->exhausted = false;
  return true;
}

bool CbqLaneReader::LoadIndex(CbqLaneIndex *index, std::string *error) {
  if (index == nullptr) {
    return SetError(error, "CBQ LoadIndex requires a non-null index");
  }
  if (mate_count_ != 1 && mate_count_ != 2) {
    return SetError(error, "CBQ reader supports mate_count 1 or 2");
  }
  Close();
  if (!impl_->OpenStreamAndReadHeader(path_, mate_count_, error)) {
    Close();
    return false;
  }
  if (!impl_->ReadCurrentLaneIndex(path_, error)) {
    Close();
    return false;
  }

  index->has_headers = impl_->file_header.HasHeaders();
  index->total_records = impl_->current_lane_records;
  index->blocks = impl_->block_index;
  Close();
  return true;
}

bool CbqLaneReader::OpenRangeWithIndex(const CbqLaneIndex &index,
                                       uint64_t first_record,
                                       uint64_t record_count,
                                       std::string *error) {
  if (mate_count_ != 1 && mate_count_ != 2) {
    return SetError(error, "CBQ reader supports mate_count 1 or 2");
  }
  if (index.total_records != 0 && index.blocks.empty()) {
    return SetError(error, "CBQ cached index has records but no blocks");
  }
  if (!index.blocks.empty() &&
      index.blocks.back().cumulative_records != index.total_records) {
    return SetError(error,
                    "CBQ cached index total does not match last block total");
  }
  uint64_t previous_cumulative_records = 0;
  for (size_t i = 0; i < index.blocks.size(); ++i) {
    if (index.blocks[i].offset < 64U) {
      return SetError(error, "CBQ cached index contains an invalid offset");
    }
    if (index.blocks[i].cumulative_records < previous_cumulative_records) {
      return SetError(error,
                      "CBQ cached index cumulative records are not monotonic");
    }
    previous_cumulative_records = index.blocks[i].cumulative_records;
  }

  Close();
  if (!impl_->OpenStreamAndReadHeader(path_, mate_count_, error)) {
    Close();
    return false;
  }
  if (impl_->file_header.HasHeaders() != index.has_headers) {
    Close();
    return SetError(error,
                    "CBQ cached index header metadata does not match file");
  }
  impl_->block_index = index.blocks;
  impl_->current_lane_records = index.total_records;
  if (first_record > impl_->current_lane_records) {
    std::ostringstream msg;
    msg << "CBQ range starts at record " << first_record
        << " but lane has only " << impl_->current_lane_records
        << " records";
    Close();
    return SetError(error, msg.str());
  }

  uint64_t end_record = impl_->current_lane_records;
  if (record_count != std::numeric_limits<uint64_t>::max()) {
    const uint64_t remaining = impl_->current_lane_records - first_record;
    end_record =
        record_count <= remaining ? first_record + record_count
                                  : impl_->current_lane_records;
  }

  impl_->read_ordinal = 0;
  impl_->lane_record_index = 0;
  impl_->batch_first_record = 0;
  impl_->range_mode = true;
  impl_->range_first_record = first_record;
  impl_->range_end_record = end_record;
  impl_->batch.reset();
  if (!impl_->SeekToRangeStart(first_record, error)) {
    Close();
    return false;
  }
  impl_->opened = true;
  impl_->exhausted = false;
  return true;
}

CbqReadStatus CbqLaneReader::Next(CbqReadView *record, std::string *error) {
  if (record == nullptr) {
    SetError(error, "CBQ Next requires a non-null record");
    return CbqReadStatus::kError;
  }
  CbqReadBatchView batch;
  const CbqReadStatus status = NextBatch(1, &batch, error);
  if (status != CbqReadStatus::kRecord) {
    return status;
  }
  if (batch.records == nullptr || batch.record_count == 0) {
    SetError(error, "CBQ Next yielded an empty batch");
    return CbqReadStatus::kError;
  }
  *record = batch.records[0];
  return CbqReadStatus::kRecord;
}

CbqReadStatus CbqLaneReader::NextBatch(uint32_t max_records,
                                       CbqReadBatchView *batch,
                                       std::string *error) {
  if (batch == nullptr) {
    SetError(error, "CBQ NextBatch requires a non-null batch");
    return CbqReadStatus::kError;
  }
  if (max_records == 0) {
    SetError(error, "CBQ NextBatch requires a positive max_records");
    return CbqReadStatus::kError;
  }
  *batch = CbqReadBatchView();
  if (!impl_->opened && impl_->exhausted) {
    return CbqReadStatus::kEnd;
  }
  if (!impl_->opened) {
    SetError(error, "CBQ reader is not open");
    return CbqReadStatus::kError;
  }

  for (;;) {
    const BlockLoadStatus status =
        impl_->LoadNextAvailableBlock(mate_count_, error);
    if (status == BlockLoadStatus::kError) {
      return CbqReadStatus::kError;
    }
    if (status == BlockLoadStatus::kEnd) {
      return CbqReadStatus::kEnd;
    }
    if (!impl_->batch) {
      SetError(error, "CBQ batch storage is not available");
      return CbqReadStatus::kError;
    }

    DecodedBlock &block = impl_->batch->block;
    size_t start = block.record_index;
    const size_t total = static_cast<size_t>(block.num_records);
    if (start > total) {
      SetError(error, "CBQ record cursor is out of range");
      return CbqReadStatus::kError;
    }
    size_t remaining = total - start;

    if (impl_->range_mode) {
      const uint64_t block_first_record = impl_->batch_first_record;
      const uint64_t current_record =
          block_first_record + static_cast<uint64_t>(start);
      if (current_record < impl_->range_first_record) {
        const uint64_t local_start =
            impl_->range_first_record - block_first_record;
        if (local_start >= block.num_records) {
          block.record_index = total;
          continue;
        }
        start = static_cast<size_t>(local_start);
        block.record_index = start;
        remaining = total - start;
      }

      const uint64_t absolute_start =
          block_first_record + static_cast<uint64_t>(start);
      if (absolute_start >= impl_->range_end_record) {
        block.record_index = total;
        continue;
      }
      const uint64_t range_remaining =
          impl_->range_end_record - absolute_start;
      if (range_remaining < static_cast<uint64_t>(remaining)) {
        remaining = static_cast<size_t>(range_remaining);
      }
    }

    if (remaining == 0) {
      block.record_index = total;
      continue;
    }
    const size_t take =
        std::min<size_t>(remaining, static_cast<size_t>(max_records));
    if (take > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
      SetError(error, "CBQ unread batch size exceeds uint32_t");
      return CbqReadStatus::kError;
    }

    batch->records = impl_->batch->read_views.data() + start;
    batch->record_count = static_cast<uint32_t>(take);
    batch->preserves_source_order = true;
    batch->backing_storage_owned_by_reader = true;
    batch->backing = impl_->batch;
    block.record_index = start + take;
    return CbqReadStatus::kRecord;
  }
}

bool CbqLaneReader::HasHeaders() const {
  return impl_ && impl_->file_header.HasHeaders();
}

uint64_t CbqLaneReader::CurrentLaneRecordCount() const {
  return impl_ ? impl_->current_lane_records : 0U;
}

void CbqLaneReader::Close() {
  if (impl_) {
    if (impl_->stream.is_open()) {
      impl_->stream.close();
    }
    impl_->batch.reset();
    impl_->block_index.clear();
    impl_->read_ordinal = 0;
    impl_->lane_record_index = 0;
    impl_->batch_first_record = 0;
    impl_->current_lane_records = 0;
    impl_->range_first_record = 0;
    impl_->range_end_record = 0;
    impl_->range_mode = false;
    impl_->opened = false;
    impl_->exhausted = false;
  }
}

}  // namespace chromap
