#ifndef OVERFLOW_READER_H_
#define OVERFLOW_READER_H_

#include <string>
#include <cstdio>
#include <cstdint>

// Simple reader for overflow files
// Reads length-prefixed records sequentially
class OverflowReader {
public:
    explicit OverflowReader(const std::string& path);
    ~OverflowReader();

    // Read next record header and payload
    // Returns false on EOF, true on success
    bool ReadNext(uint32_t& out_rid, std::string& out_payload);

    // When the file begins with AtacSpillFileHeader, the prefix is consumed
    // on the first ReadNext and these reflect the header contents.
    bool FileHasAtacSpillHeader() const { return file_has_atac_spill_header_; }
    uint16_t AtacSpillSchemaFromFileHeader() const {
        return atac_spill_schema_from_file_header_;
    }

    // Check if reader is valid (file opened successfully)
    bool IsValid() const { return file_ != nullptr; }

    // Get current file path
    const std::string& GetPath() const { return path_; }

private:
    bool ConsumeAtacSpillFilePrefixIfPresent();

    std::string path_;
    FILE* file_;
    bool prefix_checked_ = false;
    bool file_has_atac_spill_header_ = false;
    uint16_t atac_spill_schema_from_file_header_ = 0;
};

#endif  // OVERFLOW_READER_H_
