#pragma once

#include "suffix_array.hpp"
#include <math.h>
#include <vector>

struct overlap {
  std::uint32_t query_id, target_id;
  std::uint32_t query_start_position, target_start_position;
  std::uint32_t query_end_position, target_end_position;
  double score;
  bool reversed;
  friend std::ostream &operator<<(std::ostream &out, const overlap &o) {
    out << "q_id : " << o.query_id << " t_id : " << o.target_id << " ";
    out << "q_pos (" << o.query_start_position << ", " << o.query_end_position
        << ") -> ";
    out << "t_pos (" << o.target_start_position << ", " << o.target_end_position
        << ") -> ";
    out << "score: " << o.score;
    if (!o.reversed)
      out << " strain: +";
    else
      out << " strain: -";
    out << std::endl;
    return out;
  }
};

std::vector<overlap> chain(std::vector<match> &matches, int GAP, int BANDWIDTH,
                           int minimal_anchors,
                           double secondary_to_primary_ratio,
                           int secondary_alignements) {
  // ova funkcija napisana je pod pretpostavkom da su matchevi izmedu jednog
  // upita i jedne reference

  std::vector<overlap> result;

  if (matches.size() == 0) {
    overlap o = {0, 0, 0, 0, 0, 0, -1, false};
    result.emplace_back(o);
    return result;
  }

  // sortirati po zavrsnoj poziciji matcha
  sort(matches.begin(), matches.end(), [](auto m1, auto m2) -> bool {
    if (m1.target_position + m1.match_size - 1 !=
        m2.target_position + m2.match_size - 1)
      return m1.target_position + m1.match_size - 1 <
             m2.target_position + m2.match_size - 1;
    return m1.query_position + m1.match_size - 1 <
           m2.query_position + m2.match_size - 1;
  });

  double avg_match_size;
  for (auto m : matches) {
    avg_match_size += m.match_size;
  }
  avg_match_size /= matches.size();

  std::vector<double> scores(matches.size());
  std::vector<uint32_t> predecessors(matches.size());
  std::vector<bool> used(matches.size(), false);

  int best_end_index = 0;

  for (int i = 0; i < matches.size(); i++) {
    const auto &m1 = matches[i];
    scores[i] = m1.match_size;
    predecessors[i] = i;

    for (int j = i - 1; j >= 0 && i - j <= BANDWIDTH; j--) {
      const auto &m2 = matches[j];

      std::uint32_t query_dist = (m1.query_position + m1.match_size - 1) -
                                 (m2.query_position + m2.match_size - 1);
      std::uint32_t target_dist = (m1.target_position + m1.match_size - 1) -
                                  (m2.target_position + m2.match_size - 1);
      std::uint32_t matching_bases =
          std::min(std::min(query_dist, target_dist), m1.match_size);

      double gap_cost = 0;
      if ((m2.query_position + m2.match_size - 1) >
              (m1.query_position + m1.match_size - 1) ||
          std::max(query_dist, target_dist) > GAP) {
        gap_cost = HUGE_VAL;
      } else {
        std::uint32_t l = query_dist >= target_dist ? query_dist - target_dist
                                                    : target_dist - query_dist;
        if (l != 0)
          gap_cost = avg_match_size * 0.01 * l + 0.5 * log2(l);
      }

      double score = scores[j] + matching_bases - gap_cost;
      if (score > scores[i] && score > m1.match_size) {
        scores[i] = score;
        predecessors[i] = j;
      }
    }
  }

  std::vector<uint32_t> anchors_sorted(matches.size());
  std::iota(anchors_sorted.begin(), anchors_sorted.end(), 0);
  std::sort(anchors_sorted.begin(), anchors_sorted.end(),
            [&](uint32_t i, uint32_t j) { return scores[i] > scores[j]; });

  best_end_index = anchors_sorted[0];
  int included = 0;

  while (best_end_index != -1 || included < secondary_alignements) {
    used[best_end_index] = true;
    int anchor_counter = 0;
    int best_start_index = best_end_index;
    while (predecessors[best_start_index] != best_start_index) {
      best_start_index = predecessors[best_start_index];
      ++anchor_counter;
      used[best_start_index] = true;
    }
    used[best_start_index] = true;

    overlap res = {matches[best_start_index].query_id,
                   matches[best_start_index].target_id,
                   matches[best_start_index].query_position,
                   matches[best_start_index].target_position,
                   matches[best_end_index].query_position +
                       matches[best_end_index].match_size - 1,
                   matches[best_end_index].target_position +
                       matches[best_end_index].match_size - 1,
                   scores[best_end_index],
                   matches[best_start_index].reversed};

    bool include = true;
    for (auto o : result) {
      if (abs(static_cast<int>(o.target_start_position) -
              static_cast<int>(res.target_start_position)) <= 1000 &&
          abs(static_cast<int>(o.target_end_position) -
              static_cast<int>(res.target_end_position)) <= 1000)
        include = false;
    }
    if (include && anchor_counter >= minimal_anchors) {
      ++included;
      result.emplace_back(res);
    }

    best_end_index = -1;
    for (auto i = 0; i < anchors_sorted.size(); ++i) {
      if (!used[anchors_sorted[i]]) {
        best_end_index = anchors_sorted[i];
        break;
      }
    }

    if (scores[best_end_index] <
        secondary_to_primary_ratio * scores[anchors_sorted[0]])
      break;
  }

  return result;
}
