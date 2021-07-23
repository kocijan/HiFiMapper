#include "mapper/chaining.hpp"

namespace mapper {

std::ostream& operator<<(std::ostream& out, const Overlap& o) {
  out << "q_id : " << o.query_id << " t_id : " << o.target_id << " ";
  out << "q_pos (" << o.query_start_position << ", " << o.query_end_position << ") -> ";
  out << "t_pos (" << o.target_start_position << ", " << o.target_end_position << ") -> ";
  out << "score: " << o.score;
  if (!o.reversed)
    out << " strain: +";
  else
    out << " strain: -";
  out << std::endl;
  return out;
}

std::vector<Overlap> Chain(std::vector<Match>& matches, const std::uint32_t gap, const std::uint32_t bandwidth,
                           std::uint32_t minimal_anchors, const double secondary_to_primary_ratio,
                           const std::uint32_t secondary_alignements) {
  // assumes matches are exclusive pairs of {target_i, query_j}

  std::vector<Overlap> result;

  if (matches.size() == 0) {
    Overlap o = {0, 0, 0, 0, 0, 0, -1, false};
    result.emplace_back(o);
    return result;
  }

  // sortirati po zavrsnoj poziciji matcha
  sort(matches.begin(), matches.end(), [](auto m1, auto m2) -> bool {
    if (m1.target_position + m1.match_size - 1 != m2.target_position + m2.match_size - 1)
      return m1.target_position + m1.match_size - 1 < m2.target_position + m2.match_size - 1;
    return m1.query_position + m1.match_size - 1 < m2.query_position + m2.match_size - 1;
  });

  double avg_match_size = 0.0;
  for (auto m : matches) {
    avg_match_size += m.match_size;
  }
  avg_match_size /= static_cast<double>(matches.size());

  std::vector<double> scores(matches.size());
  std::vector<uint32_t> predecessors(matches.size());
  std::vector<std::uint8_t> used(matches.size(), 0U);

  std::uint32_t best_end_index = 0;
  for (std::uint32_t i = 0; i < matches.size(); i++) {
    const auto& m1 = matches[i];
    scores[i] = m1.match_size;
    predecessors[i] = i;

    if (i > 0U) {
      for (std::uint32_t band = 1U; band <= i && band <= bandwidth; ++band) {
        auto const j = i - band;
        const auto& m2 = matches[j];
        std::uint32_t query_dist = (m1.query_position + m1.match_size - 1U) - (m2.query_position + m2.match_size - 1U);
        std::uint32_t target_dist =
            (m1.target_position + m1.match_size - 1U) - (m2.target_position + m2.match_size - 1U);
        std::uint32_t matching_bases = std::min(std::min(query_dist, target_dist), m1.match_size);

        double gap_cost = 0.0;
        if ((m2.query_position + m2.match_size - 1U) > (m1.query_position + m1.match_size - 1U) ||
            std::max(query_dist, target_dist) > gap) {
          gap_cost = std::numeric_limits<double>::max() / 4.0;
        } else {
          std::uint32_t l = query_dist >= target_dist ? query_dist - target_dist : target_dist - query_dist;
          if (l != 0) gap_cost = avg_match_size * 0.01 * l + 0.5 * log2(l);
        }

        const double score = scores[j] + matching_bases - gap_cost;
        if (score > scores[i] && score > m1.match_size) {
          scores[i] = score;
          predecessors[i] = j;
        }
      }
    }
  }

  std::vector<std::uint32_t> anchors_sorted(matches.size());
  std::iota(anchors_sorted.begin(), anchors_sorted.end(), 0U);
  std::sort(anchors_sorted.begin(), anchors_sorted.end(),
            [&scores](const std::uint32_t i, const std::uint32_t j) -> bool { return scores[i] > scores[j]; });

  best_end_index = anchors_sorted[0];
  std::uint32_t included = 0;

  bool terminate = false;
  while (terminate || included < secondary_alignements) {
    used[best_end_index] = true;
    std::uint32_t anchor_counter = 0;
    std::uint32_t best_start_index = best_end_index;
    while (predecessors[best_start_index] != best_start_index) {
      best_start_index = predecessors[best_start_index];
      ++anchor_counter;
      used[best_start_index] = true;
    }
    used[best_start_index] = true;

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

    terminate = true;
    for (std::uint32_t i = 0U; i < anchors_sorted.size(); ++i) {
      if (!used[anchors_sorted[i]]) {
        best_end_index = anchors_sorted[i];
        terminate = false;
        break;
      }
    }

    if (scores[best_end_index] < secondary_to_primary_ratio * scores[anchors_sorted[0]]) {
      break;
    }
  }

  return result;
}
}  // namespace mapper
