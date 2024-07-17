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
  Cmp(T& sequence_, uint32_t offset_) : sequence(sequence_), offset(offset_) {}

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

std::vector<std::uint32_t> compute_lcp_segment_tree(const std::vector<std::uint32_t>& lcp);

std::uint32_t segment_tree_query_min(const std::vector<std::uint32_t>& segment_tree, const std::uint32_t qlow,
                                     const std::uint32_t qhigh, const std::uint32_t low, const std::uint32_t high,
                                     const std::uint32_t pos);

// https://github.com/tatiana/suffix_array

std::tuple<std::vector<std::uint32_t>, std::vector<std::uint32_t>, std::vector<std::uint32_t>> ComputeBinarySearchPath(
    const std::uint32_t n);

template <typename T>
std::uint32_t longest_common_prefix_by_index(const T& text1, const std::uint32_t pos1, const T& text2,
                                             const std::uint32_t pos2) {
  std::uint32_t i = 0;
  while (true) {
    if (pos1 + i >= text1.size() || pos2 + i >= text2.size()) {
      break;
    }
    if (text1.Code(pos1 + i) != text2.Code(pos2 + i)) {
      break;
    }
    i += 1;
  }
  return i;
}

template <typename T>
std::vector<std::uint32_t> ComputeXlcp(const std::vector<std::uint32_t>& xm, const std::vector<std::uint32_t>& m,
                                       const std::vector<std::uint32_t>& sa, const T& sequence,
                                       const std::vector<std::uint32_t>& segment_tree) {
  // TODO DELETUS DELETUS TODO TODO DELETUS DELETUS TODO TODO DELETUS DELETUS TODO TODO DELETUS DELETUS TODO TODO
  // DELETUS DELETUS TODO
  // std::vector<std::uint32_t> xlcp_array(sequence.size());
  // for (std::uint32_t i = 0; i < xm.size(); ++i) {
  //   std::uint32_t pos1 = xm[i];
  //   std::uint32_t pos2 = m[i];
  //   std::uint32_t lcp = longest_common_prefix_by_index(sequence, sa[pos1], sequence, sa[pos2]);
  //   // stderr lcp
  //   std::cerr << "lcp: " << lcp << std::endl;
  //   // xlcp_array[pos2] = lcp;
  // }
  // return xlcp_array;
  // TODO DELETUS DELETUS TODO TODO DELETUS DELETUS TODO TODO DELETUS DELETUS TODO TODO DELETUS DELETUS TODO TODO
  // DELETUS DELETUS TODO

  uint32_t n = sequence.size();
  std::uint32_t n2 = 1;
  while (n2 < n) {
    n2 *= 2;
  }

  std::vector<std::uint32_t> xlcp_array(sequence.size());
  for (std::uint32_t i = 0; i < xm.size(); ++i) {
    std::uint32_t pos1 = xm[i];
    std::uint32_t pos2 = m[i];
    std::uint32_t pos_min = std::min(pos1, pos2);
    std::uint32_t pos_max = std::max(pos1, pos2);
    pos_max = pos_max - 1;
    std::uint32_t xlcp = segment_tree_query_min(segment_tree, pos_min, pos_max, 0, n2 - 1, 0);
    xlcp_array[pos2] = xlcp;
  }
  return xlcp_array;
}

template <typename T>
std::uint32_t lcp_with_substring(const T& sequence, const std::vector<std::uint32_t>& sa, const T& substring,
                                 const std::uint32_t position, const std::uint32_t size, const std::uint32_t l,
                                 const std::uint32_t r, const std::vector<std::uint32_t>& llcp,
                                 const std::vector<std::uint32_t>& rlcp, std::uint32_t M) {
  std::uint32_t lcp = 0;
  if (l >= r) {
    if (llcp[M] == l) {
      lcp = l + longest_common_prefix_by_index(sequence, sa[M] + l, substring, l + position);
    } else if (llcp[M] > l) {
      lcp = l;
    } else {
      lcp = llcp[M];
    }
  } else {
    if (rlcp[M] == r) {
      lcp = r + longest_common_prefix_by_index(sequence, sa[M] + r, substring, r + position);
    } else if (rlcp[M] > r) {
      lcp = r;
    } else {
      lcp = rlcp[M];
    }
  }
  return std::min(lcp, size);
}

template <typename T>
std::uint32_t lcp_search_lo(const T& sequence, const std::vector<std::uint32_t>& sa, const T& substring,
                            const std::uint32_t position, const std::uint32_t size,
                            const std::vector<std::uint32_t>& llcp, const std::vector<std::uint32_t>& rlcp) {
  std::uint32_t substring_size = size;
  std::uint32_t n = static_cast<std::uint32_t>(sequence.size());
  std::uint32_t L = 0;
  std::uint32_t R = n - 1;
  std::uint32_t l = longest_common_prefix_by_index(sequence, sa[L], substring, position);
  std::uint32_t r = longest_common_prefix_by_index(sequence, sa[R], substring, position);
  while ((R - L) > 0) {
    std::uint32_t M = (L + R) / 2;
    std::uint32_t m = lcp_with_substring(sequence, sa, substring, position, size, l, r, llcp, rlcp, M);
    if (m == substring_size) {
      R = M;
      r = m;
    } else if (sequence.Code(sa[M] + m) > substring.Code(m + position)) {
      R = M;
      r = m;
    } else {
      if (L == M) {
        ++L;
        break;
      }
      L = M;
      l = m;
    }
  }
  return L;
}

template <typename T>
std::uint32_t lcp_search_hi(const T& sequence, const std::vector<std::uint32_t>& sa, const T& substring,
                            const std::uint32_t position, const std::uint32_t size,
                            const std::vector<std::uint32_t>& llcp, const std::vector<std::uint32_t>& rlcp) {
  std::uint32_t substring_size = substring.size();
  std::uint32_t n = static_cast<std::uint32_t>(sequence.size());
  std::uint32_t L = 0;
  std::uint32_t R = n - 1;
  std::uint32_t l = longest_common_prefix_by_index(sequence, sa[L], substring, position);
  std::uint32_t r = longest_common_prefix_by_index(sequence, sa[R], substring, position);
  while ((R - L) > 0) {
    std::uint32_t M = (L + R) / 2;
    std::uint32_t m = lcp_with_substring(sequence, sa, substring, position, size, l, r, llcp, rlcp, M);
    if (m == substring_size) {
      if (L == M) {
        if (std::min(longest_common_prefix_by_index(sequence, sa[M + 1], substring, 0), substring_size) ==
            substring_size)
          ++L;
        break;
      }
      L = M;
      l = m;
    } else if (sequence.Code(sa[M] + m) > substring.Code(m + position)) {
      R = M;
      r = m;
    } else {
      if (L == M) {
        ++L;
        break;
      }
      L = M;
      l = m;
    }
  }
  return R;
}

template <typename T>
void LookSa(T& query, const std::uint32_t position, const std::uint32_t size, std::vector<Match>& result,
            const std::vector<std::uint32_t>& sa, const T& sequence, const int quality, const double min_match,
            const std::uint32_t range_size, const bool discard, const bool extended_search,
            const std::vector<std::uint32_t>& llcp, const std::vector<std::uint32_t>& rlcp) {
  std::uint32_t lo = 0, hi = static_cast<std::uint32_t>(sequence.size() - 1);
  std::uint32_t matched_size = 0;

  // TODO use quality, extended_search, min_match etc

  // for (; matched_size < size || (extended_search && (hi - lo > range_size)); ++matched_size) {
  //   if (position + matched_size >= query.size()) {
  //     break;
  //   }

  //   if (quality != 0 && !query.Score(position + matched_size)) {
  //     continue;
  //   }

  //   const auto pattern = query[position + matched_size];
  //   // const auto er = std::equal_range(sa.begin() + lo, sa.begin() + hi, pattern, detail::Cmp(sequence,
  //   matched_size));

  //   // const auto first = std::lower_bound(sa.begin() + lo, sa.begin() + hi, pattern, detail::Cmp(sequence,
  //   matched_size));
  //   // const auto last = std::upper_bound(sa.begin() + lo, sa.begin() + hi, pattern, detail::Cmp(sequence,
  //   matched_size));

  //   if (first == last) {
  //     break;
  //   }

  //   lo = first;
  //   hi = last;
  // }

  lo = lcp_search_lo(sequence, sa, query, position, size, llcp, rlcp);
  hi = lcp_search_hi(sequence, sa, query, position, size, llcp, rlcp);

  if (discard && (hi - lo >= range_size)) {
    return;
  }

  matched_size = size;

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

  SequencesCollectionSuffixArray(std::vector<T>& sequences_, bool create_lcp = false, int sample_length = 75,
                                 int lcp_search_size = 100)
      : sequences(sequences_) {
    number_of_sequences = static_cast<std::uint32_t>(sequences_.size());
    suffix_arrays = std::vector<std::vector<std::uint32_t>>(number_of_sequences);

#pragma omp parallel for
    for (std::uint32_t i = 0; i < number_of_sequences; ++i) {
      suffix_arrays[i] = ConstructSuffixArray(sequences_[i]);
    }

    if (create_lcp) {
      lcps = std::vector<std::vector<std::uint32_t>>(number_of_sequences);
      llcp = std::vector<std::vector<std::uint32_t>>(number_of_sequences);
      rlcp = std::vector<std::vector<std::uint32_t>>(number_of_sequences);
      positions = std::vector<std::vector<std::unordered_set<std::uint32_t>>>(number_of_sequences);

#pragma omp parallel for
      for (std::uint32_t i = 0; i < number_of_sequences; ++i) {
        lcps[i] = CreateLCParray(sequences_[i].InflateData(), suffix_arrays[i]);

        std::uint32_t n = sequences_[i].size();
        auto [lm, m, rm] = ComputeBinarySearchPath(n);
        std::vector<std::uint32_t> segment_tree = compute_lcp_segment_tree(lcps[i]);
        rlcp[i] = ComputeXlcp(rm, m, suffix_arrays[i], sequences_[i], segment_tree);
        llcp[i] = ComputeXlcp(lm, m, suffix_arrays[i], sequences_[i], segment_tree);

        positions[i] = std::vector<std::unordered_set<std::uint32_t>>(suffix_arrays[i].size(),
                                                                      std::unordered_set<std::uint32_t>());
        FindCharacteristic(positions[i], lcps[i], suffix_arrays[i], lcp_search_size, 5, 10000, sample_length);
      }
    }
  }

  void Look(T& query, std::uint32_t position, std::uint32_t size, std::vector<std::vector<Match>>& result,
            double min_match, int quality, int RANGE_SIZE, bool discard, bool extended_search) {
    for (std::uint32_t sequence_index = 0; sequence_index < number_of_sequences; sequence_index++)
      LookSa<T>(query, position, size, result[sequence_index], suffix_arrays[sequence_index], sequences[sequence_index],
                quality, min_match, RANGE_SIZE, discard, extended_search, llcp[sequence_index], rlcp[sequence_index]);
  }

  std::vector<std::vector<std::unordered_set<std::uint32_t>>> positions;

 private:
  std::vector<std::vector<std::uint32_t>> lcps;
  std::vector<std::vector<std::uint32_t>> suffix_arrays;
  std::vector<std::vector<std::uint32_t>> llcp;
  std::vector<std::vector<std::uint32_t>> rlcp;
  std::vector<T>& sequences;
  std::uint32_t number_of_sequences;
};

}  // namespace mapper

#endif  // MAPPER_SUFFIX_ARRAY_HPP_
