#ifndef CBQ_READER_H_
#define CBQ_READER_H_

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace chromap {

enum class CbqReadStatus { kRecord, kEnd, kError };

struct CbqByteSpan {
  const char *data = nullptr;
  size_t size = 0;
};

struct CbqPackedSequenceView {
  const uint8_t *words = nullptr;
  size_t word_bytes = 0;
  uint64_t base_offset = 0;
  uint32_t length = 0;
  const uint64_t *n_positions = nullptr;
  size_t n_positions_count = 0;
  bool available = false;
};

struct CbqSegmentView {
  uint32_t source_index = 0;
  CbqByteSpan sequence;
  CbqByteSpan quality;
  CbqPackedSequenceView packed_sequence;
  uint32_t original_length = 0;
  bool has_quality = true;
};

struct CbqReadView {
  CbqByteSpan read_name;
  CbqByteSpan read_name_extra;
  uint64_t read_ordinal = 0;
  uint64_t lane_read_ordinal = 0;
  char read_filter = 'Y';
  uint32_t segment_count = 0;
  const CbqSegmentView *segments = nullptr;
};

// Borrowed view over decoded CBQ records. The backing shared_ptr keeps the
// decoded block alive while callers consume records after the reader advances.
struct CbqReadBatchView {
  const CbqReadView *records = nullptr;
  uint32_t record_count = 0;
  bool preserves_source_order = true;
  bool backing_storage_owned_by_reader = true;
  std::shared_ptr<const void> backing;
};

size_t CbqSegmentSequenceLength(const CbqSegmentView &segment);
char CbqSegmentBaseAscii(const CbqSegmentView &segment, size_t index);
void MaterializeCbqSegmentSequence(const CbqSegmentView &segment,
                                   std::string *sequence);
bool MaterializeCbqSegmentSequenceToBuffer(const CbqSegmentView &segment,
                                           char *dest, size_t capacity,
                                           size_t *length_out,
                                           std::string *error);
bool MaterializeCbqSegmentSequenceAndReverseComplementToBuffers(
    const CbqSegmentView &segment, char *dest, size_t capacity,
    char *reverse_complement_dest, size_t reverse_complement_capacity,
    size_t *length_out, std::string *error);

class CbqLaneReader {
 public:
  CbqLaneReader(const std::string &path, uint32_t mate_count);
  CbqLaneReader(const CbqLaneReader &) = delete;
  CbqLaneReader &operator=(const CbqLaneReader &) = delete;
  ~CbqLaneReader();

  bool Open(std::string *error);
  bool OpenRange(uint64_t first_record, uint64_t record_count,
                 std::string *error);
  CbqReadStatus Next(CbqReadView *record, std::string *error);
  CbqReadStatus NextBatch(uint32_t max_records, CbqReadBatchView *batch,
                          std::string *error);
  void Close();

  const std::string &path() const { return path_; }
  uint32_t mate_count() const { return mate_count_; }

  // Whether the CBQ file carries per-record names (headers). Valid only after a
  // successful Open(). Barcoded ATAC requires this on both lanes so the
  // read/barcode record-alignment check has names to compare.
  bool HasHeaders() const;

  // Total records in the current lane, populated after OpenRange() has parsed
  // the CBQINDEX footer. Sequential Open() does not require or populate it.
  uint64_t CurrentLaneRecordCount() const;

 private:
  struct Impl;

  std::string path_;
  uint32_t mate_count_ = 0;
  std::unique_ptr<Impl> impl_;
};

}  // namespace chromap

#endif  // CBQ_READER_H_
