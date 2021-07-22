#ifndef MAPPER_CHAINING_HPP_
#define MAPPER_CHAINING_HPP_

#include <vector>
#include <cmath>

#include "mapper/suffix_array.hpp"

namespace mapper {

struct Overlap {
  std::uint32_t query_id, target_id;
  std::uint32_t query_start_position, target_start_position;
  std::uint32_t query_end_position, target_end_position;
  double score;
  bool reversed;
};

std::ostream& operator<<(std::ostream& out, const Overlap& o);

std::vector<Overlap> Chain(std::vector<Match>& matches, const std::uint32_t gap, const std::uint32_t bandwidth,
                           std::uint32_t minimal_anchors, const double secondary_to_primary_ratio,
                           const std::uint32_t secondary_alignements);  // namespace mapper
}  // namespace mapper

#endif /* MAPPER_CHAINING_HPP_ */
