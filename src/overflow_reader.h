#ifndef OVERFLOW_READER_H_
#define OVERFLOW_READER_H_

#include <string>
#include <cstdio>
#include <cstdint>

namespace chromap {
struct AtacKwaySpillRecordHeaderV1;
}

// Simple reader for overflow files
// Reads length-prefixed records sequentially
class OverflowReader {
public:
    explicit OverflowReader(const std::string& path);
    ~OverflowReader();

    // Read next record header and payload
    // Returns false on EOF, true on success
    bool ReadNext(uint32_t& out_rid, std::string& out_payload);

    // ATAC k-way spill files without optional sections: reads the next
    // fixed-size record header straight into `header`, with the block
    // framing checks of ReadNext. Returns 1 for a record, 0 at end of file
    // and -1 with a message on any error (it never exits).
    int ReadNextAtacRecordHeader(chromap::AtacKwaySpillRecordHeaderV1* header,
                                 std::string* error);

    // When the file begins with AtacSpillFileHeader, the prefix is consumed
    // on the first ReadNext and these reflect the header contents.
    bool FileHasAtacSpillHeader() const { return file_has_atac_spill_header_; }
    uint16_t AtacSpillSchemaFromFileHeader() const {
        return atac_spill_schema_from_file_header_;
    }
    bool FileHasAtacKwayHeader() const { return file_has_atac_kway_header_; }
    uint32_t AtacKwayReferenceIdFromFileHeader() const {
        return atac_kway_reference_id_from_file_header_;
    }

    // Check if reader is valid (file opened successfully)
    bool IsValid() const { return file_ != nullptr; }

    // errno from the failed open when !IsValid(); 0 otherwise.
    int OpenErrno() const { return open_errno_; }

    // Get current file path
    const std::string& GetPath() const { return path_; }

private:
    bool ConsumeAtacSpillFilePrefixIfPresent();

    std::string path_;
    FILE* file_;
    int open_errno_ = 0;
    bool prefix_checked_ = false;
    bool file_has_atac_spill_header_ = false;
    bool file_has_atac_kway_header_ = false;
    uint16_t atac_spill_schema_from_file_header_ = 0;
    uint32_t atac_kway_reference_id_from_file_header_ = 0;
    uint32_t atac_kway_block_records_remaining_ = 0;
    uint32_t atac_kway_block_bytes_remaining_ = 0;
};

#endif  // OVERFLOW_READER_H_
