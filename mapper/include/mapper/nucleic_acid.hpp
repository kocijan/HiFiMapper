// Copyright (c) 2020 Robert Vaser

// taken from Biosoup

#ifndef MAPPER_NUCLEIC_ACID_HPP_
#define MAPPER_NUCLEIC_ACID_HPP_

#include <memory_resource>
#include <algorithm>
#include <stdexcept>
#include <cstdint> #include <numeric>
#include <string>
#include <memory>
#include <vector>
#include <atomic>
#include <limits>
#include <cmath>

namespace biosoup {
constexpr static std::uint8_t kNucleotideCoder[] = {
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
    255, 0,   255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0,
    1,   1,   0,   255, 255, 2,   3,   255, 255, 2,   255, 1,   0,   255, 255, 255, 0,   1,   3,   3,   2,   0,
    255, 3,   255, 255, 255, 255, 255, 255, 255, 0,   1,   1,   0,   255, 255, 2,   3,   255, 255, 2,   255, 1,
    0,   255, 255, 255, 0,   1,   3,   3,   2,   0,   255, 3,   255, 255, 255, 255, 255, 255};

constexpr static char kNucleotideDecoder[] = {'A', 'C', 'G', 'T'};

class NucleicAcid {
 public:
  NucleicAcid() = default;

  NucleicAcid(const std::string& name, const std::string& data)
      : NucleicAcid(name.c_str(), static_cast<std::uint32_t>(name.size()), data.c_str(),
                    static_cast<std::uint32_t>(data.size())) {}

  NucleicAcid(const char* name, std::uint32_t name_len, const char* data, std::uint32_t data_len)
      : id(num_objects++),
        name(name, name_len),
        deflated_data(),
        deflated_quality(),
        inflated_len(data_len),
        is_reverse_complement(0) {
    deflated_data.reserve(static_cast<std::size_t>(data_len / 32. + .999));
    std::uint64_t block = 0;
    for (std::uint32_t i = 0; i < data_len; ++i) {
      std::uint64_t c = kNucleotideCoder[static_cast<std::uint8_t>(data[i])];
      if (c == 255ULL) {
        throw std::invalid_argument("[biosoup::NucleicAcid::NucleicAcid] error: not a nucleotide");
      }
      block |= c << ((i << 1) & 63);
      if (((i + 1) & 31) == 0 || i == data_len - 1) {
        deflated_data.emplace_back(block);
        block = 0;
      }
    }
  }

  NucleicAcid(const std::string& name, const std::string& data, const std::string& quality)
      : NucleicAcid(name.c_str(), static_cast<std::uint32_t>(name.size()), data.c_str(),
                    static_cast<std::uint32_t>(data.size()), quality.c_str(),
                    static_cast<std::uint32_t>(quality.size())) {}

  NucleicAcid(const char* name, std::uint32_t name_len, const char* data, std::uint32_t data_len, const char* quality,
              std::uint32_t quality_len)
      : NucleicAcid(name, name_len, data, data_len) {
    // this is changed
    deflated_quality.reserve(static_cast<std::size_t>(1. * quality_len / 64 + .999));
    std::uint64_t block = 0;
    for (std::uint32_t i = 0; i < quality_len; ++i) {
      std::uint32_t c = (quality[i] - '!') < kQualityThreshold ? 0 : 1;
      block = block << 1UL | c;
      if ((i + 1) % 64 == 0 || i == quality_len - 1) {
        deflated_quality.emplace_back(block);
        block = 0;
      }
    }
  }

  NucleicAcid(const NucleicAcid&) = default;
  NucleicAcid& operator=(const NucleicAcid&) = default;

  NucleicAcid(NucleicAcid&&) = default;
  NucleicAcid& operator=(NucleicAcid&&) = default;

  ~NucleicAcid() = default;

  std::uint32_t size() const { return inflated_len; }

  std::uint64_t Code(std::uint32_t i) const {
    std::uint64_t x = 0;
    if (is_reverse_complement) {
      i = inflated_len - i - 1;
      x = 3;
    }
    // i >> 5 the block in which the value is stored, one block contains 32
    // nucleotides
    return ((deflated_data[i >> 5] >> ((i << 1) & 63)) & 3) ^ x;
  }

  bool Score(std::uint32_t i) const {
    if (is_reverse_complement) {
      i = inflated_len - i - 1;
    }
    return ((deflated_quality[i >> 6] << (i % 64)) >> 63 & 1) ? true : false;
  }

  std::string InflateData(std::uint32_t i = 0, std::uint32_t len = std::numeric_limits<std::uint32_t>::max()) const {
    if (i >= inflated_len) {
      return std::string{};
    }
    len = std::min(len, inflated_len - i);

    std::string dst{};
    dst.reserve(len);
    for (; len; ++i, --len) dst += kNucleotideDecoder[Code(i)];
    return dst;
  }

  std::string InflateQuality(std::uint32_t i = 0,
                             std::uint32_t len = -1) const {  // NOLINT
    if (deflated_quality.empty() || i >= inflated_len) {
      return std::string{};
    }

    len = std::min(len, inflated_len - i);

    std::string dst{};
    dst.reserve(len);
    for (; len; ++i, --len) {
      dst += Score(i) ? '~' : '!';
    }
    return dst;
  }

  void ReverseAndComplement() {  // Watson-Crick base pairing
    is_reverse_complement ^= 1;
  }

  std::uint32_t Id() const noexcept { return id; }

  std::string const& Name() const noexcept { return name; }

  static std::atomic<std::uint32_t> num_objects;
  static int kQualityThreshold;  // Phred quality treshold

  std::uint32_t id;
  std::string name;

  std::pmr::vector<std::uint64_t> deflated_data;
  std::pmr::vector<std::uint64_t> deflated_quality;

  std::uint32_t inflated_len;
  bool is_reverse_complement;
};
}  // namespace biosoup

#endif
