// Unit tests for peak_caller/frag_compact_store (packing, policy, materialize).
#include "libmacs3/frag_compact_store.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

using chromap::peaks::CeilLog2U32;
using chromap::peaks::CompactFrag64;
using chromap::peaks::FragPackingPolicy;
using chromap::peaks::FragPeakMemoryAccumulator;
using chromap::peaks::FragPeakStorageClass;
using chromap::peaks::WideFrag;

void Assert(bool cond, const char* what) {
  if (!cond) {
    std::cerr << "FAIL: " << what << std::endl;
    std::exit(1);
  }
}

void TestCeilLog2() {
  Assert(CeilLog2U32(0) == 0u, "ceil_log2(0)");
  Assert(CeilLog2U32(1) == 0u, "ceil_log2(1)");
  Assert(CeilLog2U32(2) == 1u, "ceil_log2(2)");
  Assert(CeilLog2U32(2001) == 11u, "ceil_log2(2001)");
}

void TestPolicy2000Compact() {
  FragPackingPolicy p;
  std::string err;
  Assert(FragPackingPolicy::FromMaxInsertSize(2000, 16, &p, &err), "policy 2000");
  Assert(p.storage == FragPeakStorageClass::kCompact64, "2000->compact");
  Assert(p.length_bits == 11, "len_bits@2000");
  Assert(p.count_bits == 21, "count_bits@2000");
  Assert(p.max_length == 2047u, "max_len@2000");
  Assert(p.max_count == (1u << 21) - 1u, "max_count@2000");
}

void TestPolicyHugeWide() {
  FragPackingPolicy p;
  std::string err;
  // Force count_bits < 16: need length_bits > 16, e.g. max_insert ~2^17
  int huge = 400000;  // ceil_log2(400001) = 19, count_bits=13 < 16
  Assert(FragPackingPolicy::FromMaxInsertSize(huge, 16, &p, &err), "policy huge");
  Assert(p.storage == FragPeakStorageClass::kWide, "huge->wide");
}

void TestCompactPackUnpackRoundtrip() {
  FragPackingPolicy p;
  std::string err;
  Assert(FragPackingPolicy::FromMaxInsertSize(2000, 16, &p, &err), "policy for pack");
  for (uint32_t len = 1u; len < 200u; len += 17u) {
    for (uint32_t cnt = 1u; cnt < 100u; cnt += 9u) {
      uint32_t packed = (cnt << p.length_bits) | (len & ((1u << p.length_bits) - 1u));
      uint32_t low = packed & ((1u << p.length_bits) - 1u);
      uint32_t c = packed >> p.length_bits;
      Assert(len == low, "len roundtrip");
      Assert(cnt == c, "cnt roundtrip");
      (void)sizeof(CompactFrag64);
    }
  }
}

void TestOverflow() {
  std::string e;
  FragPeakMemoryAccumulator a(2000, 20, &e);
  Assert(a.IsValid(), "acc valid");
  std::vector<std::string> n = {"chr1"};
  a.ResizeForReferenceNames(1, n);
  e.clear();
  Assert(!a.Add(0, 0, 10000, 1, &e), "length overflow");
  Assert(!e.empty(), "length overflow msg");
  e.clear();
  Assert(a.Add(0, 0, 10, 1, &e), "valid row");
  e.clear();
  uint32_t too_many = a.policy().max_count + 1u;
  Assert(!a.Add(0, 100, 200, too_many, &e), "count overflow");
  Assert(!e.empty(), "count err msg");
}

void TestMaterializeOrder() {
  std::string e;
  FragPeakMemoryAccumulator a(2000, 16, &e);
  Assert(a.IsValid() && e.empty(), "mat acc");
  // chrB rid=0, chrA rid=1 -> materialize lex order A then B
  std::vector<std::string> names = {"chrB", "chrA"};
  a.ResizeForReferenceNames(2, names);
  Assert(a.Add(0, 10, 20, 2, &e), "r0");
  e.clear();
  Assert(a.Add(1, 5, 15, 1, &e), "r1a");
  e.clear();
  Assert(a.Add(1, 3, 8, 3, &e), "r1b");
  e.clear();
  std::vector<chromap::peaks::ChromFragments> chs;
  Assert(a.MaterializeToChromFragments(&chs, &e), "mat");
  Assert(chs.size() == 2u, "two chroms");
  Assert(chs[0].name == "chrA" && chs[1].name == "chrB", "lex chrom order");
  Assert(chs[0].frags.size() == 2u, "chrA nfrag");
  Assert(chs[0].frags[0].start == 3 && chs[0].frags[0].end == 8, "sort start");
  Assert(chs[0].frags[1].start == 5 && chs[0].frags[1].end == 15, "sort end");
  Assert(chs[1].frags[0].start == 10 && chs[1].frags[0].count == 2, "chrB");
}

}  // namespace

int main() {
  TestCeilLog2();
  TestPolicy2000Compact();
  TestPolicyHugeWide();
  TestCompactPackUnpackRoundtrip();
  TestOverflow();
  TestMaterializeOrder();
  static_assert(sizeof(CompactFrag64) == 8, "size CompactFrag64");
  static_assert(sizeof(WideFrag) == 12, "size WideFrag");
  std::cout << "OK test_frag_compact_store" << std::endl;
  return 0;
}
