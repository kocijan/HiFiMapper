#include "mapper/lcp.hpp"

namespace mapper {

std::vector<uint32_t> CreateForwardArray(const std::vector<uint32_t>& suffix_array, const std::vector<uint32_t>& lcp) {
  // koliko unaprijed ima istih slova
  auto n = suffix_array.size();  // velicina SA
  std::vector<uint32_t> fa(n);
  for (int i = 0; i < n; i++) {
    // na indeksu i u fa je na mjestu i vrijednost koja odgovara stvarnom slovu
    // i
    auto& val = fa[suffix_array[i]];
    val = std::max(val, lcp[i]);  // stavimo na mjesto i duljinu prefiksa
    if (i < n - 1) {
      auto& val2 = fa[suffix_array[i + 1]];
      val2 = std::max(val2, lcp[i]);
    }
  }
  return fa;
}

}  // namespace mapper
