// Copyright (c) 2021 Suzana Pratljacic

#ifndef SIMULATOR_UTILITY_HPP_
#define SIMULATOR_UTILITY_HPP_

#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/sequence.hpp"


namespace simulator_utility {

inline bool endsWith(const std::string& fullString, const std::string& suffix){
  if (fullString.size() < suffix.size()) return false;
  return fullString.compare(fullString.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline bool endsWithAny(const std::string& path, std::vector<std::string> &validSuffixes){
  for (auto const &suffix : validSuffixes){
    if (endsWith(path, suffix)) return true;
  }
  return false;
}

std::unique_ptr<bioparser::Parser<biosoup::Sequence>>
    CreateParser(const std::string& path) {

  //the file can be parsed with fasta parser
  std::vector<std::string> fastaValidSuffixes (
    {".fasta", ".fasta.gz", ".fna", ".fna.gz", ".fa", ".fa.gz"});
  if (endsWithAny(path, fastaValidSuffixes)) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(path);
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  //the file can be parsed with fastq parser
  std::vector<std::string> fastqValidSuffixes (
    {".fastq", ".fastq.gz", ".fq", ".fq.gz"});
  if (endsWithAny(path, fastqValidSuffixes)) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastqParser>(path);
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  std::cerr << "[HiFiMapper::CreateParser] error: file " << path
            << " has unsupported format extension (valid extensions: .fasta, "
            << ".fasta.gz, .fna, .fna.gz, .fa, .fa.gz, .fastq, .fastq.gz, "
            << ".fq, .fq.gz)"
            << std::endl;
  return nullptr;
}
}  // namespace simulator_utility

#endif  // SIMULATOR_UTILITY_HPP_