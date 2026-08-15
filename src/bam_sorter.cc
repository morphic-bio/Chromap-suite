#include "bam_sorter.h"
#include "utils.h"
#include <fstream>
#include <sstream>
#include <cstring>

namespace chromap {

// Static member initialization
std::atomic<uint64_t> BamSorter::monotonic_id_{0};

// SpillFileReader implementation for k-way merge
struct BamSorter::SpillFileReader {
    std::ifstream stream;
    BamRecord currentRecord;
    bool hasRecord;
    std::string filename;
    int sourceId;
    
    SpillFileReader(const std::string& fname, int id) 
        : hasRecord(false), filename(fname), sourceId(id) {
        stream.open(fname.c_str(), std::ios::binary);
    }
    
    bool readNext() {
        if (!stream.good()) {
            hasRecord = false;
            return false;
        }
        
        // Read totalSize first
        uint32_t totalSize;
        stream.read(reinterpret_cast<char*>(&totalSize), sizeof(uint32_t));
        if (stream.eof()) {
            hasRecord = false;
            return false;
        }
        
        const size_t METADATA_SIZE = sizeof(uint8_t) + sizeof(uint64_t);
        const size_t CORE_SIZE = sizeof(bam1_core_t);
        
        // Read the rest: [hasY][readId][core][data]
        uint32_t payloadSize = totalSize - sizeof(uint32_t);
        currentRecord.data.resize(totalSize);
        
        // Copy totalSize to start of data
        memcpy(currentRecord.data.data(), &totalSize, sizeof(uint32_t));
        
        // Read payload
        stream.read(currentRecord.data.data() + sizeof(uint32_t), payloadSize);
        // Check gcount() first - EOF can set good() false even when read succeeded
        if (stream.gcount() != static_cast<std::streamsize>(payloadSize)) {
            hasRecord = false;
            return false;
        }
        // Allow EOF after successful read (last record in file)
        if (stream.eof() && stream.gcount() == static_cast<std::streamsize>(payloadSize)) {
            // This is fine - we got all the data, EOF is expected on last record
        } else if (!stream.good() && !stream.eof()) {
            // Real error (not just EOF)
            hasRecord = false;
            return false;
        }
        
        // Extract fields for key and metadata
        uint8_t hasYFlag = currentRecord.data[sizeof(uint32_t)];
        uint64_t readId;
        memcpy(&readId, currentRecord.data.data() + sizeof(uint32_t) + sizeof(uint8_t), sizeof(uint64_t));
        
        bam1_core_t core;
        memcpy(&core, currentRecord.data.data() + sizeof(uint32_t) + METADATA_SIZE, CORE_SIZE);
        
        currentRecord.key = SortKey{core.tid, core.pos, core.flag, core.mtid, core.mpos, core.isize};
        currentRecord.hasY = (hasYFlag != 0);
        currentRecord.readId = readId;
        hasRecord = true;
        return true;
    }
    
    ~SpillFileReader() {
        if (stream.is_open()) stream.close();
    }
};

SortKey BamSorter::extractSortKey(const bam1_core_t* core) {
    return SortKey{core->tid, core->pos, core->flag, core->mtid, core->mpos, core->isize};
}

void BamSorter::serializeBam1(bam1_t* b, uint64_t readId, bool hasY, std::vector<char>& out) {
    const size_t METADATA_SIZE = sizeof(uint8_t) + sizeof(uint64_t);
    const size_t CORE_SIZE = sizeof(bam1_core_t);
    
    // Format: [totalSize:uint32][hasY:uint8][readId:uint64][core:bam1_core_t][data:l_data bytes]
    uint32_t totalSize = sizeof(uint32_t) + METADATA_SIZE + CORE_SIZE + b->l_data;
    out.resize(totalSize);
    
    uint8_t hasYFlag = hasY ? 1 : 0;
    char* ptr = out.data();
    
    // Write totalSize first
    memcpy(ptr, &totalSize, sizeof(uint32_t));
    ptr += sizeof(uint32_t);
    
    // Write metadata
    memcpy(ptr, &hasYFlag, sizeof(uint8_t));
    ptr += sizeof(uint8_t);
    memcpy(ptr, &readId, sizeof(uint64_t));
    ptr += sizeof(uint64_t);
    
    // Write core
    memcpy(ptr, &b->core, CORE_SIZE);
    ptr += CORE_SIZE;
    
    // Write data
    memcpy(ptr, b->data, b->l_data);
}

bam1_t* BamSorter::reconstructBam1(const char* storedData, uint32_t storedSize) {
    const size_t METADATA_SIZE = sizeof(uint8_t) + sizeof(uint64_t);
    const size_t CORE_SIZE = sizeof(bam1_core_t);
    
    // storedData format: [totalSize:uint32][hasY:uint8][readId:uint64][core:bam1_core_t][data:l_data]
    // Skip totalSize (already known), hasY, readId
    const char* corePtr = storedData + sizeof(uint32_t) + METADATA_SIZE;
    const char* dataPtr = corePtr + CORE_SIZE;
    uint32_t l_data = storedSize - sizeof(uint32_t) - METADATA_SIZE - CORE_SIZE;
    
    bam1_t* b = bam_init1();
    if (!b) {
        ExitWithMessage("Failed to initialize bam1_t");
    }
    
    // Copy core struct
    memcpy(&b->core, corePtr, CORE_SIZE);
    
    // Allocate/resize data buffer (same approach as ConvertToHtsBam)
    if (static_cast<int>(b->m_data) < static_cast<int>(l_data)) {
        b->data = (uint8_t*)realloc(b->data, l_data);
        if (!b->data) {
            bam_destroy1(b);
            ExitWithMessage("Failed to allocate BAM data buffer");
        }
        b->m_data = l_data;
    }
    b->l_data = l_data;
    memcpy(b->data, dataPtr, l_data);
    
    return b;
}

BamSorter::BamSorter(uint64_t maxRAM, const std::string& tmpDir)
    : maxRAM_(maxRAM), tmpDir_(tmpDir),
      currentRAM_(0), spillFileCounter_(0), finalized_(false),
      returnedBam_(nullptr), use_monotonic_id_(false) {
    records_.reserve(100000); // Pre-allocate space
}

BamSorter::~BamSorter() {
    cleanupSpillFiles();
    if (returnedBam_) {
        bam_destroy1(returnedBam_);
        returnedBam_ = nullptr;
    }
}

void BamSorter::addRecord(bam1_t* b, uint64_t readId, bool hasY) {
    if (!b || b->l_data == 0) return;
    
    // Use provided readId, or assign monotonic if needed
    uint64_t assigned_readId = readId;
    if (use_monotonic_id_) {
        assigned_readId = monotonic_id_.fetch_add(1, std::memory_order_relaxed);
    }
    
    BamRecord record;
    serializeBam1(b, assigned_readId, hasY, record.data);
    record.key = extractSortKey(&b->core);
    record.hasY = hasY;
    record.readId = assigned_readId;
    
    // Capture size BEFORE move (move empties record.data)
    uint64_t record_size = record.data.size();
    
    // Quick check without lock for common case
    bool needSpill = false;
    {
        std::lock_guard<std::mutex> lock(bufferMutex_);
        records_.push_back(std::move(record));
        // Estimate RAM usage: record data + overhead (use captured size, not moved record)
        currentRAM_ += record_size + sizeof(BamRecord) + 64; // 64 bytes overhead per record
        
        // Check if we need to spill to disk
        if (currentRAM_ > maxRAM_ && maxRAM_ > 0 && !records_.empty()) {
            needSpill = true;
        }
    }
    
    // Spill outside the lock to avoid blocking writers
    if (needSpill) {
        sortAndSpill();
    }
}

void BamSorter::sortAndSpill() {
    // Lock only to extract records, then release for sorting/writing
    std::vector<BamRecord> recordsToSpill;
    {
        std::lock_guard<std::mutex> lock(bufferMutex_);
        if (records_.empty()) {
            return;
        }
        // Move all records out (no copy due to move semantics)
        recordsToSpill = std::move(records_);
        records_.clear();
        currentRAM_ = 0;
    }
    
    // Sort records by SortKey + read_id (outside lock)
    std::sort(recordsToSpill.begin(), recordsToSpill.end(), BamRecordComparator());
    
    // Write sorted chunk to temp file
    int spillFileNum;
    {
        std::lock_guard<std::mutex> lock(spillFileCounterMutex_);
        spillFileNum = spillFileCounter_++;
    }
    std::string spillFile = tmpDir_ + "/chromap_sort_spill_" + std::to_string(spillFileNum) + ".dat";
    
    std::ofstream spillStream(spillFile.c_str(), std::ios::binary);
    if (!spillStream) {
        std::ostringstream errOut;
        errOut << "EXITING because of fatal ERROR: could not open spill file: " << spillFile << "\n";
        errOut << "SOLUTION: check that the disk is not full, increase the max number of open files with Linux command ulimit -n before running chromap";
        ExitWithMessage(errOut.str());
    }
    
    // Write records to spill file
    // Format: [totalSize:uint32][hasY:uint8][readId:uint64][core:bam1_core_t][data:l_data bytes]
    // Note: record.data already contains totalSize at the start (from serializeBam1)
    for (auto& record : recordsToSpill) {
        spillStream.write(record.data.data(), record.data.size());
        if (!spillStream.good()) {
            ExitWithMessage("Failed to write to spill file");
        }
    }
    
    spillStream.close();
    
    // Store spill file name (will be used for k-way merge)
    {
        std::lock_guard<std::mutex> lock(spillFileCounterMutex_);
        spillFiles_.push_back(spillFile);
    }
}

void BamSorter::finalize() {
    if (finalized_) return;
    
    std::lock_guard<std::mutex> lock(bufferMutex_);
    
    // Sort remaining in-memory records
    if (!records_.empty()) {
        std::sort(records_.begin(), records_.end(), BamRecordComparator());
    }
    
    // Initialize k-way merge with spill files
    initializeKWayMerge();
    
    finalized_ = true;
}

void BamSorter::initializeKWayMerge() {
    // Clear heap
    while (!mergeHeap_.empty()) {
        mergeHeap_.pop();
    }
    
    // Add in-memory records to heap (if any)
    if (!records_.empty()) {
        HeapEntry entry;
        entry.key = records_[0].key;
        entry.readId = records_[0].readId;
        entry.sourceId = -1;  // -1 indicates in-memory source
        entry.recordIdx = 0;   // First record in sorted vector
        entry.data = records_[0].data;
        entry.hasY = records_[0].hasY;
        mergeHeap_.push(std::move(entry));
    }
    
    // Open all spill files for reading and add first record from each to heap
    for (size_t i = 0; i < spillFiles_.size(); i++) {
        SpillFileReader* reader = new SpillFileReader(spillFiles_[i], static_cast<int>(i));
        if (reader->readNext()) {
            HeapEntry entry;
            entry.key = reader->currentRecord.key;
            entry.readId = reader->currentRecord.readId;
            entry.sourceId = reader->sourceId;
            entry.recordIdx = 0;  // Unused for spill files
            entry.data = reader->currentRecord.data;
            entry.hasY = reader->currentRecord.hasY;
            mergeHeap_.push(std::move(entry));
            spillReaders_.push_back(reader);
        } else {
            delete reader;
        }
    }
}

bool BamSorter::nextRecord(bam1_t** b, bool* hasY) {
    if (!finalized_) {
        finalize();
    }
    
    // Free previous returned buffer (contract: valid until next call)
    if (returnedBam_) {
        bam_destroy1(returnedBam_);
        returnedBam_ = nullptr;
    }
    
    // K-way merge using min-heap: SortKey then read_id
    while (!mergeHeap_.empty()) {
        // Get smallest entry from heap (key then read_id)
        HeapEntry top = mergeHeap_.top();
        mergeHeap_.pop();
        
        // Reconstruct bam1_t from stored data
        // top.data format: [totalSize:uint32][hasY:uint8][readId:uint64][core:bam1_core_t][data:l_data]
        uint32_t totalSize = top.data.size();
        returnedBam_ = reconstructBam1(top.data.data(), totalSize);
        *hasY = top.hasY;
        *b = returnedBam_;
        
        // Advance the source for the record we're returning
        if (top.sourceId == -1) {
            // Advance to next in-memory record
            size_t nextIdx = top.recordIdx + 1;
            
            // Push next in-memory record if available
            if (nextIdx < records_.size()) {
                HeapEntry nextEntry;
                nextEntry.key = records_[nextIdx].key;
                nextEntry.readId = records_[nextIdx].readId;
                nextEntry.sourceId = -1;
                nextEntry.recordIdx = nextIdx;
                nextEntry.data = records_[nextIdx].data;
                nextEntry.hasY = records_[nextIdx].hasY;
                mergeHeap_.push(std::move(nextEntry));
            }
        } else {
            // Read next record from this spill file
            SpillFileReader* reader = spillReaders_[top.sourceId];
            if (reader->readNext()) {
                // Push next record from this source
                HeapEntry nextEntry;
                nextEntry.key = reader->currentRecord.key;
                nextEntry.readId = reader->currentRecord.readId;
                nextEntry.sourceId = reader->sourceId;
                nextEntry.recordIdx = 0;  // Unused for spill files
                nextEntry.data = reader->currentRecord.data;
                nextEntry.hasY = reader->currentRecord.hasY;
                mergeHeap_.push(std::move(nextEntry));
            } else {
                // File exhausted, will be cleaned up later
            }
        }
        
        return true;
    }
    
    return false;
}

void BamSorter::cleanupSpillFiles() {
    // Close and delete all spill file readers
    for (auto* reader : spillReaders_) {
        delete reader;
    }
    spillReaders_.clear();
    
    // Delete spill files
    for (const std::string& spillFile : spillFiles_) {
        remove(spillFile.c_str());
    }
    spillFiles_.clear();
}

}  // namespace chromap
