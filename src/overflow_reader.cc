#include "overflow_reader.h"

#include <cassert>

#include "atac_spill_record.h"
#include "utils.h"

OverflowReader::OverflowReader(const std::string& path) 
    : path_(path), file_(nullptr) {
    file_ = fopen(path.c_str(), "rb");
}

OverflowReader::~OverflowReader() {
    if (file_) {
        fclose(file_);
        file_ = nullptr;
    }
}

bool OverflowReader::ConsumeAtacSpillFilePrefixIfPresent() {
    if (prefix_checked_) {
        return true;
    }
    prefix_checked_ = true;
    if (!file_) {
        return false;
    }
    uint32_t magic = 0;
    if (fread(&magic, sizeof(uint32_t), 1, file_) != 1) {
        return false;
    }
    if (magic != chromap::kAtacSpillFileMagic) {
        if (fseek(file_, 0, SEEK_SET) != 0) {
            return false;
        }
        file_has_atac_spill_header_ = false;
        return true;
    }
    chromap::AtacSpillFileHeader hdr;
    hdr.magic = magic;
    if (fread(reinterpret_cast<char*>(&hdr) + sizeof(uint32_t),
              sizeof(hdr) - sizeof(uint32_t), 1, file_) != 1) {
        return false;
    }
    if (hdr.format_version != 1 ||
        hdr.record_codec_version != chromap::kAtacSpillRecordCodecVersion) {
        chromap::ExitWithMessage(
            "Unsupported ATAC spill overflow file header version/codec");
    }
    file_has_atac_spill_header_ = true;
    atac_spill_schema_from_file_header_ = hdr.schema_mask;
    return true;
}

bool OverflowReader::ReadNext(uint32_t& out_rid, std::string& out_payload) {
    if (!file_) {
        return false;
    }
    if (!ConsumeAtacSpillFilePrefixIfPresent()) {
        return false;
    }
    
    // Read header: rid (4 bytes) + byte_len (4 bytes)
    uint32_t rid, byte_len;
    
    if (fread(&rid, sizeof(uint32_t), 1, file_) != 1) {
        // EOF or error
        return false;
    }
    
    if (fread(&byte_len, sizeof(uint32_t), 1, file_) != 1) {
        // Incomplete header - this is an error
        return false;
    }
    
    // Read payload
    out_payload.resize(byte_len);
    if (byte_len > 0) {
        if (fread(&out_payload[0], 1, byte_len, file_) != byte_len) {
            // Incomplete payload - this is an error
            return false;
        }
    }
    
    out_rid = rid;
    return true;
}
