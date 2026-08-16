#include "overflow_reader.h"

#include <cassert>
#include <cstring>

#include "atac_kway_spill.h"
#include "atac_spill_record.h"
#include "utils.h"

OverflowReader::OverflowReader(const std::string& path) 
    : path_(path), file_(nullptr) {
    file_ = fopen(path.c_str(), "rb");
    if (file_ != nullptr) {
        (void)setvbuf(file_, nullptr, _IOFBF, 8u * 1024u * 1024u);
        (void)ConsumeAtacSpillFilePrefixIfPresent();
    }
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
    char first_magic[8] = {};
    if (fread(first_magic, 1, sizeof(first_magic), file_) !=
        sizeof(first_magic)) {
        return false;
    }
    if (memcmp(first_magic, chromap::kAtacKwaySpillMagicV1,
               sizeof(first_magic)) == 0) {
        chromap::AtacKwaySpillFileHeaderV1 header = {};
        const uint16_t known_schema = static_cast<uint16_t>(
            chromap::kAtacSpillSchemaHasBamPair |
            chromap::kAtacSpillSchemaHasYHit |
            chromap::kAtacSpillSchemaIsBulk |
            chromap::kAtacSpillSchemaHasRawBarcodeEvidence);
        memcpy(header.magic, first_magic, sizeof(first_magic));
        if (fread(reinterpret_cast<char*>(&header) + sizeof(first_magic), 1,
                  sizeof(header) - sizeof(first_magic), file_) !=
                sizeof(header) - sizeof(first_magic) ||
            header.format_version != chromap::kAtacKwaySpillFormatVersion ||
            header.fixed_header_bytes != sizeof(header) ||
            header.record_codec_version !=
                chromap::kAtacKwaySpillRecordCodecVersion ||
            header.endian_marker != chromap::kAtacKwaySpillEndianMarker ||
            (header.schema_mask & ~known_schema) != 0 ||
            header.flags != 0 || header.reserved != 0) {
          chromap::ExitWithMessage(
              "Unsupported or invalid ATAC k-way spill file header");
        }
        file_has_atac_spill_header_ = true;
        file_has_atac_kway_header_ = true;
        atac_spill_schema_from_file_header_ = header.schema_mask;
        atac_kway_reference_id_from_file_header_ = header.reference_id;
        return true;
    }
    if (fseek(file_, 0, SEEK_SET) != 0) {
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
    if (hdr.format_version != chromap::kAtacSpillFileFormatVersion ||
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

    if (file_has_atac_kway_header_) {
        if (atac_kway_block_records_remaining_ == 0) {
            if (atac_kway_block_bytes_remaining_ != 0) {
                chromap::ExitWithMessage(
                    "ATAC k-way spill block has trailing bytes");
            }
            chromap::AtacKwaySpillBlockHeaderV1 block = {};
            const size_t got = fread(&block, 1, sizeof(block), file_);
            if (got == 0 && feof(file_)) {
                return false;
            }
            if (got != sizeof(block) ||
                block.magic != chromap::kAtacKwaySpillBlockMagic ||
                block.record_count == 0 || block.payload_bytes == 0 ||
                block.reserved != 0) {
                chromap::ExitWithMessage(
                    "Invalid or truncated ATAC k-way spill block header");
            }
            atac_kway_block_records_remaining_ = block.record_count;
            atac_kway_block_bytes_remaining_ = block.payload_bytes;
        }
        uint32_t byte_len = 0;
        if (atac_kway_block_bytes_remaining_ < sizeof(byte_len) ||
            fread(&byte_len, sizeof(byte_len), 1, file_) != 1 ||
            byte_len == 0 ||
            byte_len > atac_kway_block_bytes_remaining_ - sizeof(byte_len)) {
            chromap::ExitWithMessage(
                "Invalid ATAC k-way spill record length");
        }
        out_payload.resize(byte_len);
        if (fread(&out_payload[0], 1, byte_len, file_) != byte_len) {
            chromap::ExitWithMessage(
                "Truncated ATAC k-way spill record payload");
        }
        atac_kway_block_bytes_remaining_ -=
            static_cast<uint32_t>(sizeof(byte_len) + byte_len);
        --atac_kway_block_records_remaining_;
        if (atac_kway_block_records_remaining_ == 0 &&
            atac_kway_block_bytes_remaining_ != 0) {
            chromap::ExitWithMessage(
                "ATAC k-way spill block count/size mismatch");
        }
        out_rid = atac_kway_reference_id_from_file_header_;
        return true;
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
