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
                        std::vector<uint32_t>& sa, const std::uint32_t max_offset, const std::uint32_t num_positions,
                        const std::uint32_t max_size, const std::uint32_t min_size) {
#pragma omp parallel for
  for (std::size_t i = 0U; i < sa.size(); i++) {
    auto& poss = positions[sa[i]];
    std::uint32_t current_prefix = std::numeric_limits<std::uint32_t>::max();
    std::uint32_t counter = 0;
    for (std::size_t j = i; j < i + max_offset && j < sa.size(); ++j) {
      if (counter >= num_positions) break;
      current_prefix = std::min(current_prefix, lcp[j]);
      if (current_prefix < min_size) break;
      if (current_prefix < max_size) {
        bool include = true;
        for (auto pref : poss) {
          if (static_cast<std::uint32_t>(abs(static_cast<int>(pref) - static_cast<int>(current_prefix))) < min_size) {
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

    counter = 0U;
    current_prefix = std::numeric_limits<std::uint32_t>::max();
    if (i > 0UL) {
      for (std::size_t j = i - 1UL; j >= 0UL && j >= i - max_offset; --j) {
        if (counter >= num_positions) break;
        current_prefix = std::min(current_prefix, lcp[j]);
        if (current_prefix < min_size) break;
        if (current_prefix < max_size) {
          bool include = true;
          for (auto pref : poss) {
            // TODO: figure out how to without type casts?
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
    }

    if (poss.size() > num_positions) {
      std::vector<std::uint32_t> v(poss.begin(), poss.end());
      std::sort(v.begin(), v.end());
      poss.clear();

      std::copy(v.begin(), v.begin() + num_positions, std::inserter(poss, poss.end()));
    }
  }
}

}  // namespace mapper
