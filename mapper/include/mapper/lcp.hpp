#ifndef MAPPER_LCP_HPP_
#define MAPPER_LCP_HPP_

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <vector>
#include <cmath>

namespace mapper {

template <typename T>
std::vector<std::uint32_t> CreateLCParray(const T& target, const std::vector<std::uint32_t>& suffix_array) {
  const auto target_length = target.size();

  std::vector<std::uint32_t> rank(target_length, 0U);
  for (std::uint32_t i = 0; i < target_length; i++) {
    rank[suffix_array[i]] = i;
  }

  std::vector<std::uint32_t> LCP(target_length, 0U);
  std::uint32_t prefix_length = 0;
  for (std::uint32_t suffix_index = 0U; suffix_index < target_length; suffix_index++) {
    if (rank[suffix_index] == target_length - 1U) {
      prefix_length = 0;
      continue;
    }

    std::uint32_t next_rank = rank[suffix_index] + 1;
    std::uint32_t next_suffix_index = suffix_array[next_rank];
    while (suffix_index + prefix_length < target_length && next_suffix_index + prefix_length < target_length &&
           target[suffix_index + prefix_length] == target[next_suffix_index + prefix_length]) {
      prefix_length++;
    }

    LCP[rank[suffix_index]] = prefix_length;
    prefix_length = prefix_length > 0 ? prefix_length - 1 : 0;
  }
  return LCP;
}

std::vector<std::uint32_t> CreateForwardArray(const std::vector<std::uint32_t>& suffix_array,
                                              const std::vector<std::uint32_t>& lcp);

}  // namespace mapper

#endif /* MAPPER_LCP_HPP_ */
