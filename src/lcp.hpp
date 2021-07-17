#pragma once

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

using namespace std;

namespace repetitive_analysis {

template <typename T>
std::vector<uint32_t>
createLCParray(const T &target, const std::vector<uint32_t> &suffix_array) {
  auto target_length = target.size();

  std::vector<uint32_t> rank(target_length, 0);
  for (uint32_t i = 0; i < target_length; i++)
    rank[suffix_array[i]] = i;

  std::vector<uint32_t> LCP(target_length, 0);
  int prefix_length = 0;
  for (uint32_t suffix_index = 0; suffix_index < target_length;
       suffix_index++) {
    if (rank[suffix_index] == target_length - 1) {
      prefix_length = 0;
      continue;
    }

    uint32_t next_rank = rank[suffix_index] + 1;
    uint32_t next_suffix_index = suffix_array[next_rank];
    while (suffix_index + prefix_length < target_length &&
           next_suffix_index + prefix_length < target_length &&
           target[suffix_index + prefix_length] ==
               target[next_suffix_index + prefix_length])
      prefix_length++;

    LCP[rank[suffix_index]] = prefix_length;
    prefix_length = prefix_length > 0 ? prefix_length - 1 : 0;
  }
  return LCP;
}

std::vector<uint32_t>
createForwardArray(const std::vector<uint32_t> &suffix_array,
                   const std::vector<uint32_t> &lcp) {
  // koliko unaprijed ima istih slova
  auto n = suffix_array.size(); // velicina SA
  std::vector<uint32_t> fa(n);
  for (int i = 0; i < n; i++) {
    // na indeksu i u fa je na mjestu i vrijednost koja odgovara stvarnom slovu
    // i
    auto &val = fa[suffix_array[i]];
    val = std::max(val, lcp[i]); // stavimo na mjesto i duljinu prefiksa
    if (i < n - 1) {
      auto &val2 = fa[suffix_array[i + 1]];
      val2 = std::max(val2, lcp[i]);
    }
  }
  return fa;
}

} // namespace repetitive_analysis