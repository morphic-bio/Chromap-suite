#include "overflow_writer.h"

#include <sstream>
#include <unistd.h>
#include <cstdlib>
#include <cstring>

#include "atac_spill_record.h"

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

bool OverflowWriter::WriteAtacSpillFileHeaderIfNeeded(FILE* fp) {
    if (!atac_spill_header_enabled_ || !fp) {
        return true;
    }
    chromap::AtacSpillFileHeader hdr = {};
    hdr.magic = chromap::kAtacSpillFileMagic;
    hdr.format_version = 1;
    hdr.schema_mask = atac_spill_schema_mask_;
    hdr.record_codec_version = chromap::kAtacSpillRecordCodecVersion;
    hdr.reserved0 = 0;
    return fwrite(&hdr, sizeof(hdr), 1, fp) == 1;
}

OverflowWriter::~OverflowWriter() {
    // Close any remaining thread-local files
    for (auto it = tls_files_.begin(); it != tls_files_.end(); ++it) {
        if (it->second) {
            fclose(it->second);
        }
    }
    tls_files_.clear();
    tls_file_paths_.clear();
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
    if (!WriteAtacSpillFileHeaderIfNeeded(file)) {
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
    
    return file_paths;
}
