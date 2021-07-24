#include "mapper/suffix_array.hpp"

namespace mapper {
std::ostream& operator<<(std::ostream& out, const Match& m) {
  out << "(" << m.query_id << ", " << m.query_position << ") -> ";
  out << "(" << m.target_id << ", " << m.target_position << ") ";
  out << "size: " << m.match_size;
  if (!m.reversed)
    out << " strain: +";
  else
    out << " strain: -";
  return out;
}
void FindCharacteristic(std::vector<std::unordered_set<std::uint32_t>>& positions, std::vector<uint32_t>& lcp,
                        std::vector<uint32_t>& sa, int max_offset, int num_positions, int max_size, int min_size) {
#pragma omp parallel for
  for (int i = 0; i < sa.size(); i++) {
    auto& poss = positions[sa[i]];
    std::uint32_t current_prefix = UINT32_MAX;
    int counter = 0;
    for (int j = i; j < i + max_offset && j < sa.size(); ++j) {
      if (counter >= num_positions) break;
      current_prefix = std::min(current_prefix, lcp[j]);
      if (current_prefix < min_size) break;
      if (current_prefix < max_size) {
        bool include = true;
        for (auto pref : poss) {
          if (abs(static_cast<int>(pref) - static_cast<int>(current_prefix)) < min_size) {
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
      if (counter >= num_positions) break;
      current_prefix = std::min(current_prefix, lcp[j]);
      if (current_prefix < min_size) break;
      if (current_prefix < max_size) {
        bool include = true;
        for (auto pref : poss) {
          if (abs(static_cast<int>(pref) - static_cast<int>(current_prefix)) < min_size) {
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
      std::vector<std::uint32_t> v(poss.begin(), poss.end());
      std::sort(v.begin(), v.end());
      poss.clear();

      std::copy(v.cbegin(), std::next(v.cbegin(), num_positions), std::inserter(poss, poss.end()));
    }
  }
}
}  // namespace mapper
