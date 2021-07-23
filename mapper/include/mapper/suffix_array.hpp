#ifndef MAPPER_SUFFIX_ARRAY_HPP_
#define MAPPER_SUFFIX_ARRAY_HPP_

#include <unordered_set>
#include <iostream>
#include <limits>

#include "lcp.hpp"
#include "nucleic_acid.hpp"
#include "suffix_arr_algo.hpp"

namespace mapper {

struct Match {
  uint32_t query_id, target_id;
  uint32_t query_position, target_position;
  uint32_t match_size;
  bool reversed;
};

namespace detail {

template <typename T>
struct Cmp {
  T& sequence;  // fine since the object's lifetime is equalt or longer than cmp
  std::uint32_t offset;
  Cmp(T& sequence, uint32_t offset) : sequence(sequence), offset(offset) {}

  bool operator()(const uint32_t i, const std::uint64_t pattern) const noexcept {
    if (i + offset >= sequence.size()) return true;
    return sequence.Code(i + offset) < pattern;
  }
  bool operator()(const std::uint64_t pattern, const uint32_t i) const noexcept {
    if (i + offset >= sequence.size()) return false;
    return pattern < sequence.Code(i + offset);
  }
};
}  // namespace detail

template <typename T>
void LookSa(T& query, const std::uint32_t position, const std::uint32_t size, std::vector<Match>& result,
            const std::vector<std::uint32_t>& sa, const T& sequence, const int quality, const double min_match,
            const std::uint32_t range_size, const bool discard, const bool extended_search) {
  std::uint32_t lo = 0, hi = static_cast<std::uint32_t>(sa.size());
  std::uint32_t matched_size = 0;

  for (; matched_size < size || (extended_search && (hi - lo > range_size)); ++matched_size) {
    if (position + matched_size >= query.size()) {
      break;
    }

    if (quality != 0 && !query.Score(position + matched_size)) {
      continue;
    }

    const auto pattern = query.Code(position + matched_size);
    const auto er = std::equal_range(sa.begin() + lo, sa.begin() + hi, pattern, detail::Cmp(sequence, matched_size));

    if (er.first == er.second) {
      break;
    }

    lo = er.first - sa.begin();
    hi = er.second - sa.begin();
  }

  if ((discard && (hi - lo >= range_size)) || (matched_size < size * min_match)) {
    return;
  }

  for (; lo < hi; lo++) {
    Match m = {query.Id(), sequence.Id(), position, sa[lo], matched_size, query.is_reverse_complement};
    result.emplace_back(m);
  }
}

// TODO: consider return types
void FindCharacteristic(std::vector<std::unordered_set<std::uint32_t>>& positions, std::vector<uint32_t>& lcp,
                        std::vector<uint32_t>& sa, const std::uint32_t max_offset, const std::uint32_t num_positions,
                        const std::uint32_t max_size, const std::uint32_t min_size);

template <typename T>
class SequencesCollectionSuffixArray {
 public:
  SequencesCollectionSuffixArray(const SequencesCollectionSuffixArray&) = default;
  SequencesCollectionSuffixArray& operator=(const SequencesCollectionSuffixArray&) = delete;

  SequencesCollectionSuffixArray(SequencesCollectionSuffixArray&&) = default;
  SequencesCollectionSuffixArray& operator=(SequencesCollectionSuffixArray&&) = delete;

  virtual ~SequencesCollectionSuffixArray() {}

  SequencesCollectionSuffixArray(std::vector<T>& sequences, bool create_lcp = false, int sample_length = 75,
                                 int lcp_search_size = 100)
      : sequences(sequences) {
    number_of_sequences = sequences.size();
    suffix_arrays = std::vector<std::vector<std::uint32_t>>(number_of_sequences);

#pragma omp parallel for
    for (std::uint32_t i = 0; i < number_of_sequences; ++i) {
      suffix_arrays[i] = ConstructSuffixArray(sequences[i]);
    }

    if (create_lcp) {
      lcps = std::vector<std::vector<std::uint32_t>>(number_of_sequences);
      positions = std::vector<std::vector<std::unordered_set<std::uint32_t>>>(number_of_sequences);

#pragma omp parallel for
      for (std::uint32_t i = 0; i < number_of_sequences; ++i) {
        lcps[i] = CreateLCParray(sequences[i].InflateData(), suffix_arrays[i]);
        positions[i] = std::vector<std::unordered_set<std::uint32_t>>(suffix_arrays[i].size(),
                                                                      std::unordered_set<std::uint32_t>());
        FindCharacteristic(positions[i], lcps[i], suffix_arrays[i], lcp_search_size, 5, 10000, sample_length);
      }
    }
  }

  void Look(T& query, std::uint32_t position, std::uint32_t size, std::vector<std::vector<Match>>& result,
            double min_match, int quality, int RANGE_SIZE, bool discard, bool extended_search) {
    if (result.size() == 0) result = std::vector<std::vector<Match>>(number_of_sequences, std::vector<Match>());

    for (std::uint32_t sequence_index = 0; sequence_index < number_of_sequences; sequence_index++)
      LookSa<T>(query, position, size, result[sequence_index], suffix_arrays[sequence_index], sequences[sequence_index],
                quality, min_match, RANGE_SIZE, discard, extended_search);
  }

  std::vector<std::vector<std::unordered_set<std::uint32_t>>> positions;

 private:
  std::vector<std::vector<std::uint32_t>> lcps;
  std::vector<std::vector<std::uint32_t>> suffix_arrays;
  std::vector<T>& sequences;
  std::uint32_t number_of_sequences;
};

}  // namespace mapper

#endif  // MAPPER_SUFFIX_ARRAY_HPP_
