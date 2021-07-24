#include "mapper/chaining.hpp"

#include <cassert>
#include <limits>

namespace mapper {
std::vector<Overlap> Chain(std::vector<Match>& matches, int GAP, int BANDWIDTH, int minimal_anchors,
                           double secondary_to_primary_ratio, int secondary_alignements) {
  // ova funkcija napisana je pod pretpostavkom da su matchevi izmedu jednog
  // upita i jedne reference

  std::vector<Overlap> result;

  if (matches.size() == 0) {
    Overlap o = {0, 0, 0, 0, 0, 0, -1, false};
    result.emplace_back(o);
    return result;
  }

  // sortirati po zavrsnoj poziciji matcha
  std::sort(matches.begin(), matches.end(), [](const Match& m1, const Match& m2) -> bool {
    if (m1.target_position + m1.match_size - 1 != m2.target_position + m2.match_size - 1) {
      return m1.target_position + m1.match_size - 1 < m2.target_position + m2.match_size - 1;
    } else {
      return m1.query_position + m1.match_size - 1 < m2.query_position + m2.match_size - 1;
    }
  });

  double avg_match_size;
  for (const auto& m : matches) {
    avg_match_size += m.match_size;
  }
  avg_match_size /= matches.size();

  std::vector<double> scores(matches.size());
  std::vector<std::uint32_t> predecessors(matches.size());
  std::vector<std::uint8_t> used(matches.size(), false);

  int best_end_index = 0;

  for (int i = 0; i < matches.size(); i++) {
    const auto& m1 = matches[i];
    scores[i] = m1.match_size;
    predecessors[i] = i;

    for (int j = i - 1; j >= 0 && i - j <= BANDWIDTH; j--) {
      const auto& m2 = matches[j];

      std::uint32_t query_dist = (m1.query_position + m1.match_size - 1) - (m2.query_position + m2.match_size - 1);
      std::uint32_t target_dist = (m1.target_position + m1.match_size - 1) - (m2.target_position + m2.match_size - 1);
      std::uint32_t matching_bases = std::min(std::min(query_dist, target_dist), m1.match_size);

      double gap_cost = 0;
      if ((m2.query_position + m2.match_size - 1) > (m1.query_position + m1.match_size - 1) ||
          std::max(query_dist, target_dist) > GAP) {
        gap_cost = std::numeric_limits<double>::max() / 4.0;
      } else {
        std::uint32_t l = query_dist >= target_dist ? query_dist - target_dist : target_dist - query_dist;
        if (l != 0) {
          gap_cost = avg_match_size * 0.01 * l + 0.5 * log2(l);
        }
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
            [&](const std::uint32_t i, const std::uint32_t j) -> bool { return scores.at(i) > scores.at(j); });

  best_end_index = anchors_sorted[0];
  int included = 0;

  while (best_end_index != -1 || included < secondary_alignements) {
    used[best_end_index] = true;
    int anchor_counter = 0;
    int best_start_index = best_end_index;
    while (predecessors[best_start_index] != best_start_index) {
      best_start_index = predecessors[best_start_index];

      assert(best_end_index < predecessors.size());

      ++anchor_counter;
      used[best_start_index] = true;
    }
    used[best_start_index] = true;

    assert(best_end_index < matches.size());
    assert(best_end_index < scores.size());

    Overlap res = {matches[best_start_index].query_id,
                   matches[best_start_index].target_id,
                   matches[best_start_index].query_position,
                   matches[best_start_index].target_position,
                   matches[best_end_index].query_position + matches[best_end_index].match_size - 1,
                   matches[best_end_index].target_position + matches[best_end_index].match_size - 1,
                   scores[best_end_index],
                   matches[best_start_index].reversed};

    bool include = true;
    for (auto o : result) {
      if (abs(static_cast<int>(o.target_start_position) - static_cast<int>(res.target_start_position)) <= 1000 &&
          abs(static_cast<int>(o.target_end_position) - static_cast<int>(res.target_end_position)) <= 1000)
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

    assert(best_end_index != -1);
    if (scores[best_end_index] < secondary_to_primary_ratio * scores[anchors_sorted[0]]) {
      break;
    }
  }

  return result;
}
}  // namespace mapper
