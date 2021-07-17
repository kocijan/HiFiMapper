#ifndef SUFFIX_ARRAY_HPP_
#define SUFFIX_ARRAY_HPP_

#include "lcp.hpp"
#include "nucleic_acid.hpp"
#include "suffix_array_linear.hpp"

namespace {

struct match {
  uint32_t query_id, target_id;
  uint32_t query_position, target_position;
  uint32_t match_size;
  bool reversed;
  friend std::ostream &operator<<(std::ostream &out, const match &m) {
    out << "(" << m.query_id << ", " << m.query_position << ") -> ";
    out << "(" << m.target_id << ", " << m.target_position << ") ";
    out << "size: " << m.match_size;
    if (!m.reversed)
      out << " strain: +";
    else
      out << " strain: -";
    return out;
  }
};

template <typename T> struct Comp {
  T &sequence;
  std::uint32_t offset;
  Comp(T &_sequence, uint32_t _offset) : sequence(_sequence), offset(_offset) {}

  bool operator()(const uint32_t i, const std::uint64_t pattern) const {
    if (i + offset >= sequence.size())
      return true;
    return sequence[i + offset] < pattern;
  }
  bool operator()(const std::uint64_t pattern, const uint32_t i) const {
    if (i + offset >= sequence.size())
      return false;
    return pattern < sequence[i + offset];
  }
};

template <typename T>
void lookSa(T &query, std::uint32_t position, std::uint32_t size,
            std::vector<match> &result, std::vector<std::uint32_t> &sa,
            T &sequence, int quality, double min_match, int RANGE_SIZE,
            bool discard, bool extended_search) {

  std::uint32_t lo = 0, hi = sa.size();
  std::uint32_t matched_size = 0;

  for (; matched_size < size || (extended_search && (hi - lo > RANGE_SIZE));
       ++matched_size) {
    if (position + matched_size >= query.size())
      break;

    if (quality != 0 && !query.Score(position + matched_size)) {
      continue;
    }

    auto pattern = query[position + matched_size];
    auto er = std::equal_range(sa.begin() + lo, sa.begin() + hi, pattern,
                               Comp<T>(sequence, matched_size));

    if (er.first == er.second)
      break;

    lo = er.first - sa.begin();
    hi = er.second - sa.begin();
  }

  if ((discard && (hi - lo >= RANGE_SIZE)) || (matched_size < size * min_match))
    return;

  for (; lo < hi; lo++) {
    match m = {query.Id(), sequence.Id(), position,
               sa[lo],     matched_size,          query.is_reverse_complement};
    result.emplace_back(m);
  }
}

void findCharacteristic(
    std::vector<std::unordered_set<std::uint32_t>> &positions,
    std::vector<uint32_t> &lcp, std::vector<uint32_t> &sa, int max_offset,
    int num_positions, int max_size, int min_size) {
#pragma omp parallel for
  for (int i = 0; i < sa.size(); i++) {
    auto &poss = positions[sa[i]];
    std::uint32_t current_prefix = UINT32_MAX;
    int counter = 0;
    for (int j = i; j < i + max_offset && j < sa.size(); ++j) {
      if (counter >= num_positions)
        break;
      current_prefix = std::min(current_prefix, lcp[j]);
      if (current_prefix < min_size)
        break;
      if (current_prefix < max_size) {
        bool include = true;
        for (auto pref : poss) {
          if (abs(static_cast<int>(pref) - static_cast<int>(current_prefix)) <
              min_size) {
            include = false;
            break;
          }
        }
        if (include) {
          poss.insert(current_prefix);
          counter++;
        }
      }
    }

    counter = 0;
    current_prefix = UINT32_MAX;
    for (int j = i - 1; j >= 0 && j >= i - max_offset; --j) {
      if (counter >= num_positions)
        break;
      current_prefix = std::min(current_prefix, lcp[j]);
      if (current_prefix < min_size)
        break;
      if (current_prefix < max_size) {
        bool include = true;
        for (auto pref : poss) {
          if (abs(static_cast<int>(pref) - static_cast<int>(current_prefix)) <
              min_size) {
            include = false;
            break;
          }
        }
        if (include) {
          poss.insert(current_prefix);
          counter++;
        }
      }
    }

    if (poss.size() > num_positions) {
      vector<std::uint32_t> v(poss.begin(), poss.end());
      sort(v.begin(), v.end());
      poss.clear();
      std::copy(v.begin(), v.begin() + num_positions,
                std::inserter(poss, poss.end()));
    }
  }
}

}; // namespace

namespace suffix_array {

template <typename T> class SequencesCollectionSuffixArray {
public:
  SequencesCollectionSuffixArray(const SequencesCollectionSuffixArray &) =
      default;
  SequencesCollectionSuffixArray &
  operator=(const SequencesCollectionSuffixArray &) = delete;

  SequencesCollectionSuffixArray(SequencesCollectionSuffixArray &&) = default;
  SequencesCollectionSuffixArray &
  operator=(SequencesCollectionSuffixArray &&) = delete;

  virtual ~SequencesCollectionSuffixArray() {}

  SequencesCollectionSuffixArray(std::vector<T> &sequences,
                                 bool create_lcp = false,
                                 int sample_length = 75,
                                 int lcp_search_size = 100)
      : sequences(sequences) {
    number_of_sequences = sequences.size();
    suffix_arrays =
        std::vector<std::vector<std::uint32_t>>(number_of_sequences);

#pragma omp parallel for
    for (std::uint32_t i = 0; i < number_of_sequences; ++i)
      suffix_arrays[i] = suffix_array_linear(sequences[i]);

    if (create_lcp) {
      lcps = std::vector<std::vector<std::uint32_t>>(number_of_sequences);
      positions = std::vector<std::vector<std::unordered_set<std::uint32_t>>>(
          number_of_sequences);

#pragma omp parallel for
      for (std::uint32_t i = 0; i < number_of_sequences; ++i) {
        lcps[i] = repetitive_analysis::createLCParray(
            sequences[i].InflateData(), suffix_arrays[i]);
        positions[i] = std::vector<std::unordered_set<std::uint32_t>>(
            suffix_arrays[i].size(), std::unordered_set<std::uint32_t>());
        findCharacteristic(positions[i], lcps[i], suffix_arrays[i],
                           lcp_search_size, 5, 10000, sample_length);
      }
    }
  }

  void look(T &query, std::uint32_t position, std::uint32_t size,
            std::vector<std::vector<match>> &result, double min_match,
            int quality, int RANGE_SIZE, bool discard, bool extended_search) {

    if (result.size() == 0)
      result = std::vector<std::vector<match>>(number_of_sequences,
                                               std::vector<match>());

    for (std::uint32_t sequence_index = 0; sequence_index < number_of_sequences;
         sequence_index++)
      lookSa<T>(query, position, size, result[sequence_index],
                suffix_arrays[sequence_index], sequences[sequence_index],
                quality, min_match, RANGE_SIZE, discard, extended_search);
  }

  std::vector<std::vector<std::unordered_set<std::uint32_t>>> positions;

private:
  std::vector<std::vector<std::uint32_t>> lcps;
  std::vector<std::vector<std::uint32_t>> suffix_arrays;
  std::vector<T> &sequences;
  std::uint32_t number_of_sequences;
};

} // namespace suffix_array

#endif // SUFFIX_ARRAY_HPP_