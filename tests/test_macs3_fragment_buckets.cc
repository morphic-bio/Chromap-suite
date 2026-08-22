#include <cassert>
#include <string>
#include <vector>

#include "macs3_fragment_buckets.h"

namespace {

macs3::FragmentRecord MakeRecord(int32_t chrom_id, int32_t start,
                                 int32_t end, uint32_t count) {
  macs3::FragmentRecord record;
  record.chrom_id = chrom_id;
  record.start = start;
  record.end = end;
  record.count = count;
  return record;
}

}  // namespace

int main() {
  using macs3::FragmentRecord;

  std::vector<std::string> names = {"empty0", "chr1", "empty2", "chr2"};
  std::vector<std::vector<FragmentRecord>> buckets(4);
  buckets[1].push_back(MakeRecord(1, 10, 20, 2));
  buckets[1].push_back(MakeRecord(1, 30, 40, 3));
  buckets[3].push_back(MakeRecord(3, 50, 60, 4));
  std::string error;
  assert(chromap::peaks::CompactNonEmptyMacs3FragmentBuckets(
      &buckets, &names, &error));
  assert(error.empty());
  assert(names.size() == 2);
  assert(names[0] == "chr1");
  assert(names[1] == "chr2");
  assert(buckets.size() == 2);
  assert(buckets[0].size() == 2);
  assert(buckets[1].size() == 1);
  assert(buckets[0][0].chrom_id == 0);
  assert(buckets[0][0].start == 10);
  assert(buckets[0][0].end == 20);
  assert(buckets[0][0].count == 2);
  assert(buckets[1][0].chrom_id == 1);
  assert(buckets[1][0].start == 50);
  assert(buckets[1][0].end == 60);
  assert(buckets[1][0].count == 4);

  names = {"empty0", "empty1"};
  buckets.assign(2, std::vector<FragmentRecord>());
  error.clear();
  assert(chromap::peaks::CompactNonEmptyMacs3FragmentBuckets(
      &buckets, &names, &error));
  assert(buckets.empty());
  assert(names.empty());

  names = {"chr1"};
  buckets.assign(2, std::vector<FragmentRecord>());
  error.clear();
  assert(!chromap::peaks::CompactNonEmptyMacs3FragmentBuckets(
      &buckets, &names, &error));
  assert(!error.empty());

  names = {"chr1"};
  buckets.assign(1, std::vector<FragmentRecord>());
  buckets[0].push_back(MakeRecord(1, 10, 20, 1));
  error.clear();
  assert(!chromap::peaks::CompactNonEmptyMacs3FragmentBuckets(
      &buckets, &names, &error));
  assert(!error.empty());

  return 0;
}
