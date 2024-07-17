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

// https://github.com/tatiana/suffix_array

std::vector<std::uint32_t> compute_lcp_segment_tree(const std::vector<std::uint32_t>& lcp) {
  std::uint32_t n = lcp.size();
  // into n2 calculate the next power of 2
  std::uint32_t n2 = 1;
  while (n2 < n) {
    n2 *= 2;
  }
  std::uint32_t m = 2 * n2 - 1;
  std::vector<std::uint32_t> segment_tree(m);
  for (std::uint32_t i = 0; i < n; ++i) {
    segment_tree[n2 - 1 + i] = lcp[i];
  }
  for (std::uint32_t i = n; i < n2; ++i) {
    segment_tree[n2 - 1 + i] = std::numeric_limits<std::uint32_t>::max();
  }
  for (std::int32_t i = n2 - 2; i >= 0; --i) {
    segment_tree[i] = std::min(segment_tree[2 * i + 1], segment_tree[2 * i + 2]);
  }
  return segment_tree;
}

std::uint32_t segment_tree_query_min(const std::vector<std::uint32_t>& segment_tree, const std::uint32_t qlow,
                                     const std::uint32_t qhigh, const std::uint32_t low, const std::uint32_t high,
                                     const std::uint32_t pos) {
  if (qlow <= low && qhigh >= high) {
    return segment_tree[pos];
  }
  if (qlow > high || qhigh < low) {
    return std::numeric_limits<std::uint32_t>::max();
  }
  std::uint32_t mid = (low + high) / 2;
  return std::min(segment_tree_query_min(segment_tree, qlow, qhigh, low, mid, 2 * pos + 1),
                  segment_tree_query_min(segment_tree, qlow, qhigh, mid + 1, high, 2 * pos + 2));
}

std::tuple<std::vector<std::uint32_t>, std::vector<std::uint32_t>, std::vector<std::uint32_t>> ComputeBinarySearchPath(
    const std::uint32_t n) {
  std::vector<std::uint32_t> lm(n);
  std::vector<std::uint32_t> rm(n);
  std::vector<std::uint32_t> m(n);

  std::uint32_t l_index = 0;
  std::uint32_t r_index = n - 1;

  std::vector<std::pair<std::uint32_t, std::uint32_t>> previous_layer;
  previous_layer.push_back(std::make_pair(l_index, r_index));
  std::uint32_t tree_level = static_cast<std::uint32_t>(std::log2(n));
  while (tree_level > 0) {
    std::vector<std::pair<std::uint32_t, std::uint32_t>> current_layer;
    for (auto it = previous_layer.begin(); it != previous_layer.end(); ++it) {
      l_index = it->first;
      r_index = it->second;
      std::uint32_t m_index = (l_index + r_index) / 2;
      if (m[m_index] == 0) {
        m[m_index] = m_index;
        rm[m_index] = r_index;
        lm[m_index] = l_index;
      }
      current_layer.push_back(std::make_pair(l_index, m_index));
      current_layer.push_back(std::make_pair(m_index, r_index));
    }
    previous_layer = current_layer;
    tree_level -= 1;
  }

  for (auto it = previous_layer.begin(); it != previous_layer.end(); ++it) {
    l_index = it->first;
    r_index = it->second;
    std::uint32_t m_index = (l_index + r_index) / 2;
    if (m[m_index] == 0) {
      m[m_index] = m_index;
      rm[m_index] = r_index;
      lm[m_index] = l_index;
    }
  }

  lm.erase(lm.begin());
  lm.erase(lm.end() - 1);
  m.erase(m.begin());
  m.erase(m.end() - 1);
  rm.erase(rm.begin());
  rm.erase(rm.end() - 1);

  return std::make_tuple(lm, m, rm);
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
