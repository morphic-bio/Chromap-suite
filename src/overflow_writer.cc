#include "overflow_writer.h"

#include <sstream>
#include <unistd.h>
#include <cstdlib>
#include <cstring>
#include <limits>

#include "atac_kway_spill.h"
#include "atac_spill_record.h"
#include "utils.h"

namespace {
// Helper function to determine optimal temp directory for containers and regular systems
std::string GetOptimalTempDir(const std::string& user_specified = "") {
    if (!user_specified.empty()) {
        return user_specified;  // User override wins
    }
    
    // Check for Chromap-specific environment variable first
    const char* chromap_temp = getenv("CHROMAP_TEMP_DIR");
    if (chromap_temp && strlen(chromap_temp) > 0) {
        return chromap_temp;
    }
    
    // Check standard temp directory environment variables
    const char* tmpdir = getenv("TMPDIR");
    if (tmpdir && strlen(tmpdir) > 0) return tmpdir;
    
    const char* tmp = getenv("TMP");
    if (tmp && strlen(tmp) > 0) return tmp;
    
    const char* temp = getenv("TEMP");
    if (temp && strlen(temp) > 0) return temp;
    
    // Default fallback
    return "/tmp";
}
}

// Static member definitions
std::atomic<uint32_t> OverflowWriter::global_counter_{0};
thread_local std::unordered_map<uint32_t, FILE*> OverflowWriter::tls_files_;
thread_local std::unordered_map<uint32_t, std::string> OverflowWriter::tls_file_paths_;
thread_local std::unordered_map<uint32_t, std::vector<uint8_t>> OverflowWriter::tls_atac_block_buffers_;
thread_local std::unordered_map<uint32_t, uint32_t> OverflowWriter::tls_atac_block_record_counts_;
thread_local bool OverflowWriter::tls_initialized_{false};

OverflowWriter::OverflowWriter(const std::string& base_dir, const std::string& prefix)
    : base_dir_(base_dir), prefix_(prefix) {
    if (base_dir_.empty()) {
        base_dir_ = GetOptimalTempDir();
    }
}

void OverflowWriter::EnableAtacSpillFileHeader(uint16_t schema_mask) {
    atac_spill_header_enabled_ = true;
    atac_spill_schema_mask_ = schema_mask;
}

bool OverflowWriter::WriteAtacSpillFileHeaderIfNeeded(FILE* fp, uint32_t rid) {
    if (!atac_spill_header_enabled_ || !fp) {
        return true;
    }
    chromap::AtacKwaySpillFileHeaderV1 hdr = {};
    memcpy(hdr.magic, chromap::kAtacKwaySpillMagicV1, sizeof(hdr.magic));
    hdr.format_version = chromap::kAtacKwaySpillFormatVersion;
    hdr.fixed_header_bytes = sizeof(hdr);
    hdr.schema_mask = atac_spill_schema_mask_;
    hdr.record_codec_version = chromap::kAtacKwaySpillRecordCodecVersion;
    hdr.reference_id = rid;
    hdr.endian_marker = chromap::kAtacKwaySpillEndianMarker;
    return fwrite(&hdr, sizeof(hdr), 1, fp) == 1;
}

bool OverflowWriter::FlushAtacBlock(uint32_t rid) {
    auto file_it = tls_files_.find(rid);
    auto buffer_it = tls_atac_block_buffers_.find(rid);
    auto count_it = tls_atac_block_record_counts_.find(rid);
    if (file_it == tls_files_.end() || buffer_it == tls_atac_block_buffers_.end() ||
        count_it == tls_atac_block_record_counts_.end()) {
        return false;
    }
    if (count_it->second == 0) {
        return buffer_it->second.empty();
    }
    if (buffer_it->second.size() > std::numeric_limits<uint32_t>::max()) {
        return false;
    }
    chromap::AtacKwaySpillBlockHeaderV1 block_header = {};
    block_header.magic = chromap::kAtacKwaySpillBlockMagic;
    block_header.record_count = count_it->second;
    block_header.payload_bytes = static_cast<uint32_t>(buffer_it->second.size());
    if (fwrite(&block_header, sizeof(block_header), 1, file_it->second) != 1 ||
        fwrite(buffer_it->second.data(), 1, buffer_it->second.size(),
               file_it->second) != buffer_it->second.size()) {
        return false;
    }
    buffer_it->second.clear();
    count_it->second = 0;
    return true;
}

void OverflowWriter::WriteAtac(uint32_t rid,
                               const chromap::AtacSpillRecord& rec) {
    FILE* fp = GetFileForRid(rid);
    if (!fp) {
        chromap::ExitWithMessage("Cannot open ATAC k-way overflow file");
    }
    std::vector<uint8_t> encoded;
    std::string error;
    if (!chromap::EncodeAtacKwaySpillRecord(
            rec, atac_spill_schema_mask_, &encoded, &error) ||
        encoded.size() > std::numeric_limits<uint32_t>::max()) {
        chromap::ExitWithMessage(error.empty()
                                     ? "Cannot encode ATAC k-way record"
                                     : error);
    }
    std::vector<uint8_t>& block = tls_atac_block_buffers_[rid];
    const uint32_t encoded_bytes = static_cast<uint32_t>(encoded.size());
    const size_t old_size = block.size();
    block.resize(old_size + sizeof(encoded_bytes) + encoded.size());
    memcpy(block.data() + old_size, &encoded_bytes, sizeof(encoded_bytes));
    memcpy(block.data() + old_size + sizeof(encoded_bytes), encoded.data(),
           encoded.size());
    ++tls_atac_block_record_counts_[rid];
    if (block.size() >= chromap::kAtacKwaySpillTargetBlockBytes &&
        !FlushAtacBlock(rid)) {
        chromap::ExitWithMessage("Cannot flush ATAC k-way overflow block");
    }
}

OverflowWriter::~OverflowWriter() {
    // Close any remaining thread-local files
    for (auto it = tls_files_.begin(); it != tls_files_.end(); ++it) {
        if (it->second) {
            if (atac_spill_header_enabled_) {
                (void)FlushAtacBlock(it->first);
            }
            fclose(it->second);
        }
    }
    tls_files_.clear();
    tls_file_paths_.clear();
    tls_atac_block_buffers_.clear();
    tls_atac_block_record_counts_.clear();
}

FILE* OverflowWriter::GetFileForRid(uint32_t rid) {
    if (!tls_initialized_) {
        tls_initialized_ = true;
    }
    
    auto it = tls_files_.find(rid);
    if (it != tls_files_.end()) {
        return it->second;
    }
    
    // Create new file for this rid
    std::string filename = GenerateFilename(rid);
    FILE* file = fopen(filename.c_str(), "wb");
    if (!file) {
        return nullptr;
    }
    // ATAC already buffers complete multi-megabyte blocks. A second stdio
    // buffer per open reference file only multiplies memory use.
    if (atac_spill_header_enabled_) {
        (void)setvbuf(file, nullptr, _IONBF, 0);
    }
    if (!WriteAtacSpillFileHeaderIfNeeded(file, rid)) {
        fclose(file);
        return nullptr;
    }
    
    tls_files_[rid] = file;
    tls_file_paths_[rid] = filename;
    return file;
}

std::string OverflowWriter::GenerateFilename(uint32_t rid) {
    uint32_t counter = global_counter_.fetch_add(1);
    pid_t pid = getpid();
    
    // Use a simple thread identifier (address of thread-local variable)
    uintptr_t thread_id = reinterpret_cast<uintptr_t>(&tls_initialized_);
    
    std::ostringstream oss;
    oss << base_dir_;
    if (!base_dir_.empty() && base_dir_.back() != '/') {
        oss << "/";
    }
    oss << prefix_ << "_" << pid << "_" << thread_id << "_" << rid << "_" << counter << ".tmp";
    
    return oss.str();
}

std::vector<std::string> OverflowWriter::Close() {
    std::vector<std::string> file_paths;
    
    // Close all thread-local files and collect paths
    for (auto it = tls_files_.begin(); it != tls_files_.end(); ++it) {
        if (it->second) {
            if (atac_spill_header_enabled_ && !FlushAtacBlock(it->first)) {
                chromap::ExitWithMessage(
                    "Cannot finalize ATAC k-way overflow block");
            }
            fflush(it->second);
            fclose(it->second);
            
            auto path_it = tls_file_paths_.find(it->first);
            if (path_it != tls_file_paths_.end()) {
                file_paths.push_back(path_it->second);
            }
        }
    }
    
    tls_files_.clear();
    tls_file_paths_.clear();
    tls_atac_block_buffers_.clear();
    tls_atac_block_record_counts_.clear();
    
    return file_paths;
}
