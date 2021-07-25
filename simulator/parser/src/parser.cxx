#include "parser/parser.hpp"

namespace sim::parser {

namespace detail {

inline bool EndsWith(const std::string& fullString, const std::string& suffix) {
  if (fullString.size() < suffix.size()) return false;
  return fullString.compare(fullString.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline bool EndsWithAny(const std::string& path, std::vector<std::string>& validSuffixes) {
  for (auto const& suffix : validSuffixes) {
    if (EndsWith(path, suffix)) return true;
  }
  return false;
}
}  // namespace detail

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(const std::string& path) {
  // the file can be parsed with fasta parser
  std::vector<std::string> fastaValidSuffixes({".fasta", ".fasta.gz", ".fna", ".fna.gz", ".fa", ".fa.gz"});
  if (detail::EndsWithAny(path, fastaValidSuffixes)) {
    try {
      return bioparser::Parser<biosoup::Sequence>::Create<bioparser::FastaParser>(path);
    } catch (const std::invalid_argument& exception) {
      std::cerr << exception.what() << std::endl;
      return nullptr;
    }
  }

  // the file can be parsed with fastq parser
  std::vector<std::string> fastqValidSuffixes({".fastq", ".fastq.gz", ".fq", ".fq.gz"});
  if (detail::EndsWithAny(path, fastqValidSuffixes)) {
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
            << ".fq, .fq.gz)" << std::endl;
  return nullptr;
}
}  // namespace sim::parser
