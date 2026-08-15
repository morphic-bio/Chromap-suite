#ifndef BARCODE_CORRECTION_H_
#define BARCODE_CORRECTION_H_

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

namespace chromap {

enum class BarcodeCorrectionStatus : uint8_t {
  kRejected = 0,
  kInWhitelist = 1,
  kCorrected = 2,
};

struct BarcodeCorrectionCandidate {
  uint32_t corrected_base_index1 = 0;
  uint8_t correct_base1 = 0;
  uint32_t corrected_base_index2 = 0;
  uint8_t correct_base2 = 0;
  double score = 0.0;
  uint64_t corrected_key = 0;

  BarcodeCorrectionCandidate(uint32_t index1, uint8_t base1,
                             uint32_t index2, uint8_t base2,
                             double candidate_score, uint64_t key)
      : corrected_base_index1(index1),
        correct_base1(base1),
        corrected_base_index2(index2),
        correct_base2(base2),
        score(candidate_score),
        corrected_key(key) {}

  bool operator>(const BarcodeCorrectionCandidate &other) const {
    return std::tie(score, corrected_base_index1, correct_base1,
                    corrected_base_index2, correct_base2) >
           std::tie(other.score, other.corrected_base_index1,
                    other.correct_base1, other.corrected_base_index2,
                    other.correct_base2);
  }
};

// Correct one packed barcode without mutating a SequenceBatch. `lookup` must
// return the exact-whitelist abundance for a key and set `found` when the key
// is present (including present keys whose abundance is zero). N-mask bits use
// the packed key's little-endian base positions.
template <typename BarcodeQuality, typename AbundanceLookup>
BarcodeCorrectionStatus CorrectPackedBarcode(
    uint64_t barcode_key, uint32_t barcode_length, uint32_t barcode_n_mask,
    const BarcodeQuality &barcode_qual, int error_threshold,
    double confidence_threshold, uint64_t num_sample_barcodes,
    const AbundanceLookup &lookup, uint64_t *corrected_key) {
  if (corrected_key == nullptr || barcode_length == 0 || barcode_length > 32 ||
      barcode_qual.size() != barcode_length || error_threshold < 0 ||
      error_threshold > 2) {
    return BarcodeCorrectionStatus::kRejected;
  }
  *corrected_key = barcode_key;

  std::vector<uint32_t> n_positions;
  for (uint32_t i = 0; i < barcode_length; ++i) {
    if ((barcode_n_mask & (uint32_t{1} << i)) != 0) {
      n_positions.push_back(i);
    }
  }
  if (n_positions.size() > static_cast<size_t>(error_threshold)) {
    return BarcodeCorrectionStatus::kRejected;
  }

  bool exact_found = false;
  lookup(barcode_key, &exact_found);
  if (n_positions.empty() && exact_found) {
    return BarcodeCorrectionStatus::kInWhitelist;
  }
  if (error_threshold == 0 || num_sample_barcodes == 0) {
    return BarcodeCorrectionStatus::kRejected;
  }

  std::vector<BarcodeCorrectionCandidate> candidates;
  const uint64_t base_mask = 3;
  uint32_t i_start = 0;
  uint32_t i_end = barcode_length;
  uint32_t i_alternatives = 3;
  if (!n_positions.empty()) {
    i_start = n_positions[0];
    i_end = n_positions[0] + 1;
    i_alternatives = 4;
  }

  for (uint32_t i = i_start; i < i_end; ++i) {
    const uint64_t key_without_i =
        barcode_key & ~(base_mask << (2 * i));
    uint64_t base1 = (barcode_key >> (2 * i)) & base_mask;
    for (uint32_t ti = 0; ti < i_alternatives; ++ti) {
      base1 = (base1 + 1) & base_mask;
      const uint64_t key1 = key_without_i | (base1 << (2 * i));
      bool key1_found = false;
      const uint64_t abundance1 = lookup(key1, &key1_found);
      if (key1_found) {
        int adjusted_qual = barcode_qual[barcode_length - 1 - i] - 33;
        adjusted_qual = std::max(3, std::min(40, adjusted_qual));
        const double score =
            std::pow(10.0, -adjusted_qual / 10.0) *
            (abundance1 / static_cast<double>(num_sample_barcodes));
        candidates.push_back(BarcodeCorrectionCandidate{
            barcode_length - 1 - i, static_cast<uint8_t>(base1), 0, 0,
            score, key1});
      }

      if (error_threshold != 2) {
        continue;
      }
      uint32_t j_start = i + 1;
      uint32_t j_end = barcode_length;
      uint32_t j_alternatives = 3;
      if (n_positions.size() == 2) {
        j_start = n_positions[1];
        j_end = n_positions[1] + 1;
        j_alternatives = 4;
      }
      for (uint32_t j = j_start; j < j_end; ++j) {
        const uint64_t key_without_j = key1 & ~(base_mask << (2 * j));
        uint64_t base2 = (key1 >> (2 * j)) & base_mask;
        for (uint32_t tj = 0; tj < j_alternatives; ++tj) {
          base2 = (base2 + 1) & base_mask;
          const uint64_t key2 = key_without_j | (base2 << (2 * j));
          bool key2_found = false;
          const uint64_t abundance2 = lookup(key2, &key2_found);
          if (!key2_found) {
            continue;
          }
          int q1 = barcode_qual[barcode_length - 1 - i] - 33;
          int q2 = barcode_qual[barcode_length - 1 - j] - 33;
          q1 = std::max(3, std::min(40, q1));
          q2 = std::max(3, std::min(40, q2));
          const double score =
              std::pow(10.0, -(q1 + q2) / 10.0) *
              (abundance2 / static_cast<double>(num_sample_barcodes));
          candidates.push_back(BarcodeCorrectionCandidate{
              barcode_length - 1 - i, static_cast<uint8_t>(base1),
              barcode_length - 1 - j, static_cast<uint8_t>(base2), score,
              key2});
        }
      }
    }
  }

  if (candidates.empty()) {
    return BarcodeCorrectionStatus::kRejected;
  }
  if (candidates.size() == 1) {
    *corrected_key = candidates[0].corrected_key;
    return BarcodeCorrectionStatus::kCorrected;
  }

  std::sort(candidates.begin(), candidates.end(),
            std::greater<BarcodeCorrectionCandidate>());
  double sum_score = 0.0;
  for (const auto &candidate : candidates) {
    sum_score += candidate.score;
  }
  if (candidates[0].score / sum_score > confidence_threshold) {
    *corrected_key = candidates[0].corrected_key;
    return BarcodeCorrectionStatus::kCorrected;
  }
  return BarcodeCorrectionStatus::kRejected;
}

}  // namespace chromap

#endif  // BARCODE_CORRECTION_H_
