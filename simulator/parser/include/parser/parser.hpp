// Copyright (c) 2021 Suzana Pratljacic

// TODO: move to separate target, shared between simulator and the mapper
#ifndef SIMULATOR_PARSER_HPP_
#define SIMULATOR_PARSER_HPP_

#include <iostream>

#include "bioparser/fasta_parser.hpp"
#include "bioparser/fastq_parser.hpp"
#include "biosoup/sequence.hpp"

namespace sim::parser {

std::unique_ptr<bioparser::Parser<biosoup::Sequence>> CreateParser(const std::string& path);  // namespace sim::parser

}

#endif  // SIMULATOR_PARSER_HPP_
