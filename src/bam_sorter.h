#ifndef BAM_SORTER_H_
#define BAM_SORTER_H_

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <cstdint>
#include <climits>
#include <cstring>
#include <vector>
#include <string>
#include <mutex>
#include <queue>
#include <fstream>
#include <algorithm>
#include <atomic>

namespace chromap {

// Sort key fields extracted from bam1_core_t (STAR-Flex spec)
struct SortKey {
    int32_t tid;    // core.tid
    int32_t pos;    // core.pos
    uint16_t flag;  // core.flag
    int32_t mtid;   // core.mtid
    int32_t mpos;   // core.mpos
    int32_t isize;  // core.isize
    
    bool operator<(const SortKey& other) const {
        // Normalize unmapped to sort last
        int32_t t1 = (tid == -1) ? INT32_MAX : tid;
        int32_t t2 = (other.tid == -1) ? INT32_MAX : other.tid;
        int32_t p1 = (pos == -1) ? INT32_MAX : pos;
        int32_t p2 = (other.pos == -1) ? INT32_MAX : other.pos;
        return std::tie(t1, p1, flag, mtid, mpos, isize) <
               std::tie(t2, p2, other.flag, other.mtid, other.mpos, other.isize);
    }
};

struct BamRecord {
    std::vector<char> data;  // Serialized [core][data]
    SortKey key;             // Extracted sort key (for fast comparison)
    bool hasY;
    uint64_t readId;         // Final tie-break (deterministic within run)
};

// Entry for k-way merge heap
struct HeapEntry {
    SortKey key;
    uint64_t readId;
    int sourceId;           // -1 = in-memory, >=0 = spill file index
    size_t recordIdx;       // Index into records_ when sourceId == -1
    std::vector<char> data; // Owned copy for spill file records
    bool hasY;
};

// Comparator for std::sort: ascending by SortKey, then ascending by read_id
struct BamRecordComparator {
    bool operator()(const BamRecord& a, const BamRecord& b) const {
        if (!(a.key < b.key) && !(b.key < a.key)) {
            // Equal keys: tie-break by read_id (ascending)
            return a.readId < b.readId;
        }
        return a.key < b.key;
    }
};

// Comparator for std::priority_queue (min-heap): INVERTED for min-heap semantics
// Returns true if a should come AFTER b in sorted order (i.e., a has lower priority)
struct HeapLess {
    bool operator()(const HeapEntry& a, const HeapEntry& b) const {
        if (!(a.key < b.key) && !(b.key < a.key)) {
            // Equal keys: smaller read_id has HIGHER priority (comes first)
            return a.readId > b.readId;  // Inverted for min-heap
        }
        // Smaller key has HIGHER priority (comes first)
        return b.key < a.key;  // Inverted: b < a means a has lower priority
    }
};

class BamSorter {
public:
    BamSorter(uint64_t maxRAM, const std::string& tmpDir);
    
    // Add a BAM record to the sorter
    // b: fully populated bam1_t from ConvertToHtsBam()
    // readId: mapping.read_id_ for tie-breaking (or assigned monotonic if not deterministic)
    // hasY: Y-chromosome flag for routing
    // The sorter copies b internally; caller retains ownership
    void addRecord(bam1_t* b, uint64_t readId, bool hasY);
    
    // Finalize sorting (sort in-memory buffer, prepare k-way merge)
    void finalize();
    
    // Get next sorted record
    // Returns pointer to reconstructed bam1_t (valid until next call)
    // hasY: returns Y-chromosome flag for routing
    // Returns false when no more records
    bool nextRecord(bam1_t** b, bool* hasY);
    
    ~BamSorter();
    
private:
    // Extract SortKey from bam1_core_t
    static SortKey extractSortKey(const bam1_core_t* core);
    
    // Serialize bam1_t to storage format: [hasY][readId][core][data]
    static void serializeBam1(bam1_t* b, uint64_t readId, bool hasY, std::vector<char>& out);
    
    // Reconstruct bam1_t from storage format
    bam1_t* reconstructBam1(const char* data, uint32_t size);
    
    void sortAndSpill();
    void initializeKWayMerge();
    void cleanupSpillFiles();
    
    uint64_t maxRAM_;
    std::string tmpDir_;
    
    // Thread-safe buffer
    std::mutex bufferMutex_;
    std::vector<BamRecord> records_;
    uint64_t currentRAM_;
    
    // Spill files for when memory limit is exceeded
    std::vector<std::string> spillFiles_;
    std::mutex spillFileCounterMutex_;
    int spillFileCounter_;
    
    // Forward declaration - SpillFileReader defined in .cc file
    struct SpillFileReader;
    std::vector<SpillFileReader*> spillReaders_;
    bool finalized_;
    bam1_t* returnedBam_;  // Reconstructed bam1_t buffer returned to caller (valid until next call)
    
    // Min-heap for k-way merge: SortKey then read_id
    std::priority_queue<HeapEntry, std::vector<HeapEntry>, HeapLess> mergeHeap_;
    
    // Fallback: monotonic counter if read_id is not deterministic
    static std::atomic<uint64_t> monotonic_id_;
    bool use_monotonic_id_;
};

}  // namespace chromap

#endif  // BAM_SORTER_H_
