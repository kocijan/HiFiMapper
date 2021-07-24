#ifndef MAPPER_CHAINING_HPP_
#define MAPPER_CHAINING_HPP_

#include "mapper/suffix_array.hpp"

#include <iostream>
#include <vector>
#include <cmath>

namespace mapper {

struct Overlap {
  std::uint32_t query_id, target_id;
  std::uint32_t query_start_position, target_start_position;
  std::uint32_t query_end_position, target_end_position;
  double score;
  bool reversed;
  friend std::ostream& operator<<(std::ostream& out, const Overlap& o) {
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
};

std::vector<Overlap> Chain(std::vector<Match>& matches, int GAP, int BANDWIDTH, int minimal_anchors,
                           double secondary_to_primary_ratio, int secondary_alignements);

}  // namespace mapper

#endif /* MAPPER_CHAINING_HPP_ */
