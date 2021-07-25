// Copyright (c) 2021 Suzana Pratljacic

#include <unordered_set>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdint>
#include <random>
#include <thread>
#include <mutex>

#include <getopt.h>
#include <unistd.h>
#include <omp.h>

#include "biosoup/timer.hpp"

#include "mapper/nucleic_acid.hpp"
#include "mapper/suffix_array.hpp"
#include "mapper/chaining.hpp"
#include "mapper/parser.hpp"

namespace {

void Help() {
  std::cout << "usage: ./HiFimapper [options ...] <target> [<sequences>]\n"
               "  # default output is stdout\n"
               "  <target>\n"
               "    path to the targets in FASTA/FASTQ format\n" "  <sequences>\n" "    path to the queries in FASTA/FASTQ format\n"
               "\n"
               "  options:\n"
               "    -t, --threads <int>l\n"
               "      defaul: std::thread::hardware_concurency() / 2"
               "      number of threads\n "
               "    -l, --sample-_length <int>\n"
               "      default: 50\n"
               "      the length of the samples\n"
               "    -c, --sample_count <int>\n"
               "      default: 20\n"
               "      the number of samples extracted from each query\n"
               "    -m, --min_match <double>\n"
               "      default: 0.8\n"
               "      percentage of the sample that must be mapped for the "
               "match to be valid\n"
               "    -q, --quality <int>\n"
               "      default: 90\n"
               "      phred quality\n"
               "    -f, --frequency <int>\n"
               "      default: 10\n"
               "      maximum number of matches\n"
               "    -b, --bandwidth <int>\n"
               "      default: 10\n"
               "      size of bandwidth in which sample hits can be chained\n"
               "    -g, --gap <int>\n"
               "      default: 10000\n"
               "      maximal gap between sample hits in a chain\n"
               "    -d, --discard <bool>\n"
               "      default: false\n"
               "      discarding matches that occur more times than the "
               "default frequency\n"
               "    -e, --extended_search <bool>\n"
               "      default: false\n"
               "      allows the extended search heuristics\n"
               "    -a, --sequential <bool>\n"
               "      default: false\n"
               "      blast like sampling algorithm\n"
               "    -i, --lcp_information <bool>\n"
               "      deafult: false\n"
               "      use lcp information for sampling\n"
               "    -N, --secondary_alignements <int>\n"
               "      default: 5\n"
               "      number of secondary alignements\n"
               "    -p, --ratio <double>\n"
               "      deafult: 0.8\n"
               "      secondary to primary alignements ratio\n"
               "    -n, --minimal_anchors <int>\n"
               "      default: 3\n"
               "      minimal number of anchors on chain\n"
               "    -x, --lcp_search_size <int>\n"
               "      defult: 100\n"
               "      number of suffixes / 2 to calculate lcp\n"
               "    -h, --help\n"
               "      prints the usage\n";
}

bool FindChains(std::uint32_t query_id, std::vector<std::vector<mapper::Match>>& matches,
                std::vector<std::vector<mapper::Match>>& matches_complement, std::vector<mapper::Overlap>& result,
                int secondary_alignements, double secondary_to_primary_ratio, int gap, int bandwidth,
                int minimal_anchors) {
  std::vector<mapper::Overlap> overlaps;

  // for (auto v : matches) for (auto m : v) std::cerr << m << std::endl;
  // for(auto v : matches_complement) for(auto m : v) std::cerr << m <<
  // std::endl;

  for (auto& m : matches) {
    auto chains = mapper::Chain(m, gap, bandwidth, minimal_anchors, secondary_to_primary_ratio, secondary_alignements);
    overlaps.insert(overlaps.end(), chains.begin(), chains.end());
  }

  for (auto& m : matches_complement) {
    auto chains = mapper::Chain(m, gap, bandwidth, minimal_anchors, secondary_to_primary_ratio, secondary_alignements);
    overlaps.insert(overlaps.end(), chains.begin(), chains.end());
  }

  if (overlaps.size() == 0) {
    return false;
  }

  std::sort(overlaps.begin(), overlaps.end(),
            [](const mapper::Overlap& a, const mapper::Overlap& b) -> bool { return a.score > b.score; });
  result.insert(result.end(), overlaps.begin(),
                overlaps.size() >= secondary_alignements ? overlaps.begin() + secondary_alignements : overlaps.end());
  return true;
}

}  // namespace

std::atomic<std::uint32_t> biosoup::NucleicAcid::num_objects{0};
int biosoup::NucleicAcid::QUALITY_TRESHOLD = 90;

int main(int argc, char** argv) {
  static struct option options[] = {{"threads", required_argument, nullptr, 't'},
                                    {"sample-length", required_argument, nullptr, 'l'},
                                    {"sample-count", required_argument, nullptr, 'c'},
                                    {"min-match", required_argument, nullptr, 'm'},
                                    {"quality", required_argument, nullptr, 'q'},
                                    {"seed", required_argument, nullptr, 's'},
                                    {"frequency", required_argument, nullptr, 'f'},
                                    {"discard", no_argument, nullptr, 'd'},
                                    {"extended-search", no_argument, nullptr, 'e'},
                                    {"sequential", no_argument, nullptr, 'a'},
                                    {"lcp-information", no_argument, nullptr, 'i'},
                                    {"secondary-alignements", required_argument, nullptr, 'N'},
                                    {"ratio", required_argument, nullptr, 'p'},
                                    {"bandwidth", required_argument, nullptr, 'b'},
                                    {"gap", required_argument, nullptr, 'g'},
                                    {"minimal-anchors", required_argument, nullptr, 'n'},
                                    {"lcp-search-size", required_argument, nullptr, 'x'},
                                    {"help", no_argument, nullptr, 'h'},
                                    {nullptr, 0, nullptr, 0}};

  srand(time(NULL));

  int threads = std::thread::hardware_concurrency() / 2;
  int sample_length = 75;
  int sample_count = 20;
  double min_match = 0;
  int quality = 0;
  int seed = rand();
  int frequency = 10;
  bool discard = false;
  bool extended_search = false;
  bool sequential = false;
  bool lcp_information = false;
  int minimal_anchors = 3;
  int bandwidth = 10;
  int gap = 10000;
  int secondary_alignements = 5;
  double secondary_to_primary_ratio = 0.8;
  int lcp_search_size = 100;

  std::string tartget_path = "";
  std::string queries_path = "";

  const char* optstr = "t:l:c:m:q:s:f:deaiuN:p:b:g:n:x:";
  char arg;
  while ((arg = getopt_long(argc, argv, optstr, options, 0)) != -1) {
    switch (arg) {
      case 'd':
        discard = true;
        break;
      case 'a':
        sequential = true;
        break;
      case 'i':
        lcp_information = true;
        break;
      case 'e':
        extended_search = true;
        break;
      case 'b':
        bandwidth = std::atoi(optarg);
        break;
      case 'N':
        secondary_alignements = std::atoi(optarg);
        break;
      case 'n':
        minimal_anchors = std::atoi(optarg);
        break;
      case 'p':
        secondary_to_primary_ratio = std::stod(optarg);
        break;
      case 'g':
        gap = std::atoi(optarg);
        break;
      case 't':
        threads = std::atoi(optarg);
        break;
      case 'l':
        sample_length = std::atoi(optarg);
        break;
      case 'c':
        sample_count = std::atoi(optarg);
        break;
      case 's':
        seed = std::atoi(optarg);
        break;
      case 'm':
        min_match = std::stod(optarg);
        break;
      case 'q':
        quality = std::atoi(optarg);
        break;
      case 'f':
        frequency = std::atoi(optarg);
        break;
      case 'x':
        lcp_search_size = std::atoi(optarg);
        break;
      case 'h':
        Help();
        return 0;
      default:
        return 1;
    }
  }
    
  std::pmr::pool_options pool_options;
  pool_options.largest_required_pool_block = 25001 / 32. + .999;
  pool_options.max_blocks_per_chunk = 100000;

  std::pmr::unsynchronized_pool_resource pool(pool_options);
  biosoup::set_sequence_mem_src(std::pmr::get_default_resource());

  // TODO: check for number of threads
  const auto omp_threads = threads - 1;

  omp_set_dynamic(0);
  omp_set_num_threads(threads);
  biosoup::NucleicAcid::QUALITY_TRESHOLD = quality;

  if (argc == 1 || optind == argc) {
    std::cerr << "[HiFimapper::] error: missing target file" << std::endl;
    return 0;
  }
  tartget_path = argv[optind++];

  if (optind == argc) {
    std::cerr << "[HiFimapper::] error: missing queries file" << std::endl;
    return 0;
  }
  queries_path = argv[optind];

  auto tparser = mapper::CreateParser(tartget_path);
  if (tparser == nullptr) {
    return 1;
  }

  auto qparser = mapper::CreateParser(queries_path);
  if (qparser == nullptr) {
    return 1;
  }

  biosoup::Timer timer{};
  timer.Start();

  biosoup::NucleicAcid::num_objects = 0;
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> targets;
  try {
    targets = tparser->Parse(-1);
  } catch (std::invalid_argument& exception) {
    std::cerr << exception.what() << std::endl;
    return 1;
  }

  std::vector<biosoup::NucleicAcid> references(targets.size());
  for (auto id = 0; id < targets.size(); id++) {
    references[id] = *(targets[id].get());
  }

  std::cerr << "[HiFimapper::] Loaded reference in " << timer.Stop() << std::endl;
  timer.Start();
  std::cerr << "[HiFimapper::] Number of references: " << references.size() << std::endl;

  std::cerr << "[HiFimapper::] Creating suffix array...\n";

  auto suffixArray = mapper::SequencesCollectionSuffixArray<biosoup::NucleicAcid>(references, lcp_information,
                                                                                  sample_length, lcp_search_size);

  std::cerr << "[HiFimapper::] Suffix array constructed in " << timer.Stop() << std::endl;
  timer.Start();

  std::default_random_engine random_engine;
  random_engine.seed(seed);

  std::mutex mtx;
  std::atomic<bool> READ_ALL{false};
  std::vector<std::unique_ptr<biosoup::NucleicAcid>> queries;

  auto func = [&]() {
    while (true) {
      biosoup::NucleicAcid::num_objects = 0;
      auto records = qparser->Parse(1ULL << 28);
      if (records.size() == 0) {
        READ_ALL = true;
        return;
      }
      while (true) {
        if (queries.empty()) break;
        sleep(0.001); // doesn't really sleep
      }
      mtx.lock();
      queries = std::move(records);
      mtx.unlock();
    }
  };

  std::thread reader(func);

  int unmapped = 0;
  double total_time = 0;

  while (true) {
    while (true) {
      mtx.lock();
      if (!queries.empty()) break;
      mtx.unlock();
      sleep(0.001);
    }
    mtx.unlock();

    std::vector<biosoup::NucleicAcid> reads(queries.size());
#pragma omp parallel for num_threads(omp_threads) 
    for (auto id = 0; id < queries.size(); id++) {
      reads[id] = *(queries[id].get());
    }

    std::cerr << "[HiFimapper::] Loaded queries in " << timer.Stop() << std::endl;
    std::cerr << "[HiFimapper::] Number of queries: " << queries.size() << std::endl;
    timer.Start();

    const std::vector<std::vector<mapper::Match>> helper_vec(references.size(), std::vector<mapper::Match>());

    std::vector<std::vector<std::vector<mapper::Match>>> matches(queries.size(), helper_vec);
    std::vector<std::vector<std::vector<mapper::Match>>> matches_complement(queries.size(), helper_vec);

#pragma omp parallel for num_threads(omp_threads)
    for (int q_id = 0; q_id < queries.size(); q_id++) {
      auto& n = reads[q_id];

      if (!sequential) {
        std::uniform_int_distribution<int> position_distribution;
        if (lcp_information)
          position_distribution = std::uniform_int_distribution<int>(0, n.size() - sample_length - 10000);
        else
          position_distribution = std::uniform_int_distribution<int>(0, n.size() - sample_length - 1);

        for (uint32_t i = 0; i < sample_count; i++) {
          std::uint32_t position = position_distribution(random_engine);
          suffixArray.look(n, position, sample_length, matches[q_id], min_match, quality, frequency, discard,
                           extended_search);

          n.ReverseAndComplement();
          suffixArray.look(n, position, sample_length, matches_complement[q_id], min_match, quality, frequency, discard,
                           extended_search);
          n.ReverseAndComplement();
        }

        if (lcp_information) {
          std::vector<std::uint32_t> new_positions;
          for (auto m : matches[q_id][0]) {
            for (auto p : suffixArray.positions[m.target_id][m.target_position])
              if (p + m.query_position < n.size() - sample_length - 1) new_positions.emplace_back(p + m.query_position);
          }

          sort(new_positions.begin(), new_positions.end());
          for (int i = 0; i < new_positions.size(); i++) {
            suffixArray.look(n, new_positions[i] - 2, sample_length, matches[q_id], min_match, quality, frequency,
                             discard, extended_search);
          }

          std::vector<std::uint32_t> new_positions_c;
          for (auto m : matches_complement[q_id][0]) {
            for (auto p : suffixArray.positions[m.target_id][m.target_position])
              if (p + m.query_position < n.size() - sample_length - 1)
                new_positions_c.emplace_back(p + m.query_position);
          }

          sort(new_positions_c.begin(), new_positions_c.end());
          n.ReverseAndComplement();
          for (int i = 0; i < new_positions_c.size(); i++) {
            suffixArray.look(n, new_positions_c[i] - 2, sample_length, matches_complement[q_id], min_match, quality,
                             frequency, discard, extended_search);
          }
          n.ReverseAndComplement();
        }

      } else {
        for (uint32_t position = 0; position < n.size() - sample_length; position += sample_length) {
          suffixArray.look(n, position, sample_length, matches[q_id], min_match, quality, frequency, discard,
                           extended_search);
          n.ReverseAndComplement();

          suffixArray.look(n, position, sample_length, matches_complement[q_id], min_match, quality, frequency, discard,
                           extended_search);
          n.ReverseAndComplement();
        }
      }
    }

    auto searching_time = timer.Stop();
    std::cerr << "[HiFimapper::] Batch search time: " << searching_time << std::endl;
    timer.Start();
    total_time += searching_time;

    std::vector<std::vector<mapper::Overlap>> overlaps(queries.size());
    std::cerr << "[HiFimapper::] Start chaining...." << std::endl;

    std::vector<std::uint8_t> unmapped_ids(queries.size());

#pragma omp parallel for num_threads(omp_threads)
    for (int q_id = 0; q_id < queries.size(); q_id++) {
      if (!FindChains(q_id, matches[q_id], matches_complement[q_id], overlaps[q_id], secondary_alignements,
                      secondary_to_primary_ratio, gap, bandwidth, minimal_anchors)) {
        unmapped_ids[q_id] = true;
        ++unmapped;
      }
    }

    const auto chaining_time = timer.Stop();

    std::cerr << "[HiFimapper::] Reads chained in " << chaining_time << std::endl;
    total_time += chaining_time;
    timer.Start();

    for (const auto& a : overlaps) {
      for (const auto& o : a) {
        std::string strain = o.reversed ? "-" : "+";
        std::cout << reads[o.query_id].Name() << "\t" << reads[o.query_id].size() << "\t" << 0 << "\t"
                  << reads[o.query_id].size() - 1 << "\t" << strain << "\t" << references[o.target_id].Name() << "\t"
                  << references[o.target_id].size() << "\t" << o.target_start_position << "\t" << o.target_end_position
                  << "\t"
                  << "0"
                  << "\t"
                  << std::max(o.target_end_position - o.target_start_position + 1,
                              o.query_end_position - o.query_start_position + 1)
                  << "\t" << 255 << "\n";
      }
    }

    mtx.lock();
    queries.erase(queries.begin(), queries.end());
    mtx.unlock();
    if (READ_ALL) break;
  }

  reader.join();

  timer.Stop();
  std::cerr << std::fixed << "[HiFiMapper::] Mapping time: " << total_time << std::endl;
  std::cerr << "[HiFiMapper::] Unmapped: " << unmapped << std::endl;
  std::cerr << "[HiFiMapper::] Real time: " << timer.elapsed_time() << "s" << std::endl;
  return 0;
}
