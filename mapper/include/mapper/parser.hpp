#ifndef MAPPER_PARSER_HPP_
#define MAPPER_PARSER_HPP_

#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"

#include "nucleic_acid.hpp"

namespace mapper {

std::unique_ptr<bioparser::Parser<biosoup::NucleicAcid>> CreateParser(const std::string &path);

}  // namespace mapper

#endif  // MAPPER_PARSER_HPP_
