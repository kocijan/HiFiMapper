#include <bits/stdc++.h>
#include <getopt.h>
#include <regex>
#include <sys/time.h>

#include "biosoup/overlap.hpp"
#include "bioparser/paf_parser.hpp"
#include "bioparser/parser.hpp"
#include "parser.hpp"

struct PafOverlap {
 public:
  PafOverlap(
      const char* q_name, std::uint32_t q_name_len,
      std::uint32_t q_len,
      std::uint32_t q_begin,
      std::uint32_t q_end,
      char orientation,
      const char* t_name, std::uint32_t t_name_len,
      std::uint32_t t_len,
      std::uint32_t t_begin,
      std::uint32_t t_end,
      std::uint32_t score,
      std::uint32_t overlap_len,
      std::uint32_t quality)
      : q_name(q_name, q_name_len),
        q_len(q_len),
        q_begin(q_begin),
        q_end(q_end),
        orientation(orientation),
        t_name(t_name, t_name_len),
        t_len(t_len),
        t_begin(t_begin),
        t_end(t_end),
        score(score),
        overlap_len(overlap_len),
        quality(quality) {}

  std::uint32_t q_begin;
  std::uint32_t q_end;
  char orientation;
  std::uint32_t t_begin;
  std::uint32_t t_end;
  std::uint32_t score;
  std::string q_name;
  std::uint32_t q_len;
  std::string t_name;
  std::uint32_t t_len;
  std::uint32_t overlap_len;
  std::uint32_t quality;
};

// taken from Petar Velickovic github
inline int needleman_wunsch(const std::string &sequence1,
                            const std::string &sequence2,
                            std::vector<std::vector<int>> &dp, int match_score,
                            int mismatch_score, int gap_score) {
  for (auto i = 0; i <= sequence1.size(); i++)
    dp[i][0] = dp[0][i] = -i * gap_score;

  for (int i = 1; i <= sequence1.size(); i++) {
    for (int j = 1; j <= sequence2.size(); j++) {
      int score =
          sequence1[i - 1] == sequence2[j - 1] ? match_score : -mismatch_score;
      dp[i][j] = std::max(
          dp[i - 1][j - 1] + score,
          std::max(dp[i - 1][j] - gap_score, dp[i][j - 1] - gap_score));
    }
  }
  return dp[sequence1.size()][sequence2.size()];
}

inline std::pair<std::string, std::string>
get_optimal_alignment(const std::string &sequence1,
                      const std::string &sequence2,
                      std::vector<std::vector<int>> &dp, int match_score,
                      int mismatch_score, int gap_score) {

  std::string retA, retB;
  std::stack<char> SA, SB;
  int ii = sequence1.size(), jj = sequence2.size();
  while (ii != 0 || jj != 0) {
    if (ii == 0) {
      SA.push('-');
      SB.push(sequence2[jj - 1]);
      jj--;
    } else if (jj == 0) {
      SA.push(sequence1[ii - 1]);
      SB.push('-');
      ii--;
    } else {
      int S = (sequence1[ii - 1] == sequence2[jj - 1]) ? match_score
                                                       : -mismatch_score;
      if (dp[ii][jj] == dp[ii - 1][jj - 1] + S) {
        SA.push(sequence1[ii - 1]);
        SB.push(sequence2[jj - 1]);
        ii--;
        jj--;
      } else if (dp[ii - 1][jj] > dp[ii][jj - 1]) {
        SA.push(sequence1[ii - 1]);
        SB.push('-');
        ii--;
      } else {
        SA.push('-');
        SB.push(sequence2[jj - 1]);
        jj--;
      }
    }
  }
  while (!SA.empty()) {
    retA += SA.top();
    retB += SB.top();
    SA.pop();
    SB.pop();
  }
  return make_pair(retA, retB);
}

void train(biosoup::Sequence& reference, std::string& subreads_path, std::string& paf_path){
  
  std::unordered_map<std::string, biosoup::Sequence> reads;
  biosoup::Sequence::num_objects = 0;
  auto mparser = simulator_utility::CreateParser(subreads_path);
  auto reads_pointers = mparser->Parse(-1);
  for (auto id=0; id < reads_pointers.size(); id++){
    biosoup::Sequence s = *(reads_pointers[id].get());
    reads[s.name] = s;
  }


  std::unordered_map<std::string, std::vector<PafOverlap>> paf_lines;
  auto pafparser = bioparser::Parser<PafOverlap>::Create<bioparser::PafParser>(paf_path);
  auto paf_pointers = pafparser->Parse(-1);
  for (auto id=0; id < paf_pointers.size(); id++){
      PafOverlap o = *(paf_pointers[id].get());
      auto it = paf_lines.find(o.q_name);
      if (it == paf_lines.end()) paf_lines[o.q_name] = std::vector<PafOverlap>();
      paf_lines[o.q_name].push_back(o);
  }

  std::vector<int> quality(reference.data.size(), 0);
  std::vector<int> counter(reference.data.size(), 0);

  for (const auto & [ key, value ] : reads) {
    double treshold = 0.95;
    
    for (const auto & overlap : paf_lines[key]){
      std::string sequence1 = reference.data.substr(overlap.t_begin, overlap.t_end - overlap.t_begin);
      std::string sequence2 = value.data.substr(overlap.q_begin, overlap.q_end - overlap.q_begin);
      int vel = std::max(sequence1.size(), sequence2.size()) + 1;
      std::vector<std::vector<int>> dp (vel, std::vector<int>(vel, 0));
      double score = needleman_wunsch(sequence1, sequence2, dp, 1, 0, 0) / (double) (vel-1);
    
      if (score > treshold){
        std::pair<std::string, std::string> align = get_optimal_alignment(sequence1, sequence2, dp, 1, 0, 0);
        std::uint32_t ref_pos = overlap.t_begin;
        std::uint32_t query_pos = overlap.q_begin;
        for(auto i = 0; i < align.first.size() && i < align.second.size(); i++){
          if (align.first[i] == '-') {
            //INSERTION
            quality[ref_pos] += 33;
            counter[ref_pos] ++;
            query_pos++;
          } else if (align.second[i] == '-') {
            //DELETION
            quality[ref_pos] += 33;
            counter[ref_pos] ++;
            ref_pos++;
          } else {
            if (align.first[i] == align.second[i]){
              //MATCH
              quality[ref_pos] += value.quality[query_pos];
              counter[ref_pos]++;
              ref_pos++;
              query_pos++;
            } else {
              //MISMATCH
              quality[ref_pos] += 33;
              counter[ref_pos]++;
              ref_pos++;
              query_pos++;
            }
          }
        }

      }
    }
  }
  std::string retQuality;
  for (auto i = 0; i < quality.size(); i++){
    int q = counter[i] > 0 ? static_cast<int>(round((double) quality[i] / counter[i])) : 126;
    retQuality += (char)q; 
  }
  reference.quality = retQuality;
}

void complement(std::string &read, std::string &quality) {
  for (auto &it : read) {
    switch (static_cast<char>(std::toupper(static_cast<unsigned char>(it)))) {
    case 'A':
      it = 'T';
      break;
    case 'C':
      it = 'G';
      break;
    case 'G':
      it = 'C';
      break;
    case 'T':
    case 'U':
      it = 'A';
      break;
    case 'R':
      it = 'Y';
      break; // A || G
    case 'Y':
      it = 'R';
      break; // C || T (U)
    case 'K':
      it = 'M';
      break; // G || T (U)
    case 'M':
      it = 'K';
      break; // A || C
    case 'S':
      break; // C || G
    case 'W':
      break; // A || T (U)
    case 'B':
      it = 'V';
      break; // !A
    case 'D':
      it = 'H';
      break; // !C
    case 'H':
      it = 'D';
      break; // !G
    case 'V':
      it = 'B';
      break; // !T (!U)
    default:
      break; // N || -
    }
  }
  std::reverse(read.begin(), read.end());
  std::reverse(quality.begin(), quality.end());
}

void extractErrorModel(biosoup::Sequence &reference, std::uint32_t position,
                       std::uint32_t read_length,
                       std::default_random_engine &gene, std::string &read,
                       std::string &quality) {

  std::uniform_real_distribution<double> error_distribution(0.0, 1.0);
  std::uniform_int_distribution<int> bases_distribution(0, 3);
  std::vector<char> bases = {'A', 'C', 'G', 'T'};

  std::uint32_t sequenced = 0;
  std::uint32_t sequence_index = position;
  while (sequenced < read_length && sequence_index < reference.data.size()) {

    // DELETION
    // deletion probability is defined by the average quality score of the
    // 5’neighbors
    int NEIGH_SIZE = 2;
    int quality_sum = reference.quality[sequence_index];

    if (sequence_index > 0)
      for (auto i = sequence_index - 1;
           i > 0 && i >= sequence_index - NEIGH_SIZE; i--)
        quality_sum += reference.quality[i];

    if (sequence_index < reference.data.size())
      for (auto i = sequence_index + 1;
           i < reference.data.size() && i <= sequence_index + NEIGH_SIZE; i++)
        quality_sum += reference.quality[i];

    double deletion_probability =
        pow(10, -1 * (quality_sum / (NEIGH_SIZE * 2 + 1) - 33) / 10);

    if (error_distribution(gene) <= deletion_probability) {
      sequence_index++;
      continue;
    }

    // substitution and insertion error probability is defined by the quality
    // score at the position
    double error_probability =
        pow(10, -1 * (reference.quality[sequence_index] - 33) / 10);
    // NO ERROR
    if (error_distribution(gene) > error_probability) {
      read += reference.data[sequence_index];
      quality += reference.quality[sequence_index];
      sequenced++;
      sequence_index++;
      continue;
    }

    // ERROR
    double error_type = error_distribution(gene);
    if (error_type < 0.5) {
      // SUBSTITUTION
      while (true) {
        char base = bases[bases_distribution(gene)];
        if (base != reference.data[sequence_index]) {
          read += base;
          quality += reference.quality[sequence_index];
          sequenced++;
          sequence_index++;
          break;
        }
      }
    } else {
      // INSERTION
      if (sequenced != 0 && error_distribution(gene) < 0.5)
        read += read[sequenced - 1]; // INSERT PREVIOUS LETTER
      else
        read += bases[bases_distribution(gene)]; // INSERT RANDOM LETTER
      quality += reference.quality[sequence_index];
      sequenced++;
    }
  }
}

void extractErrorUniform(biosoup::Sequence &reference, std::uint32_t position,
                         std::uint32_t read_length,
                         std::default_random_engine &gene, std::string &read,
                         std::string &quality, double error_probability) {

  std::uniform_real_distribution<double> error_distribution(0.0, 1.0);
  std::uniform_int_distribution<int> bases_distribution(0, 3);
  std::vector<char> bases = {'A', 'C', 'G', 'T'};

  std::uint32_t sequenced = 0;
  std::uint32_t sequence_index = position;
  while (sequenced < read_length && sequence_index < reference.data.size()) {
    // NO ERROR
    if (error_distribution(gene) > error_probability) {
      read += reference.data[sequence_index];
      quality += '~';
      sequenced++;
      sequence_index++;
      continue;
    }

    // ERROR
    double error_type = error_distribution(gene);
    if (error_type < 0.33) {
      // SUBSTITUTION
      while (true) {
        char base = bases[bases_distribution(gene)];
        if (base != reference.data[sequence_index]) {
          read += base;
          quality += '!';
          sequenced++;
          sequence_index++;
          break;
        }
      }
    } else if (error_type < 0.66) {
      // INSERTION
      if (sequenced != 0 && error_distribution(gene) < 0.5)
        read += read[sequenced - 1]; // INSERT PREVIOUS LETTER
      else
        read += bases[bases_distribution(gene)]; // INSERT RANDOM LETTER
      quality += '!';
      sequenced++;
    } else {
      // DELETION
      sequence_index++;
    }
  }
}

void generateReads(biosoup::Sequence &reference, std::uint32_t read_length,
                   std::uint32_t num_reads, float complement_probability,
                   std::uint32_t SEED, double error_probability = 0) {

  std::default_random_engine gene;
  gene.seed(SEED);

  bool error_model = reference.quality.size() != 0;

  std::uniform_int_distribution<int> position_distribution(
      0, reference.data.size() - read_length - 1);
  std::uniform_real_distribution<float> complement_distribution(0.0, 1.0);

  bool reversed = false;
  for (int i = 0; i < num_reads; i++) {
    reversed = false;
    int position = position_distribution(gene);

    std::string read;
    std::string quality;

    if (error_model)
      extractErrorModel(reference, position, read_length, gene, read, quality);
    else
      extractErrorUniform(reference, position, read_length, gene, read, quality,
                          error_probability);

    if (complement_distribution(gene) < complement_probability){
      complement(read, quality);
      reversed = true;
    }

    std::cout << "@" << reference.name;
    if (reversed) printf("_%d_TARGET_SIZE_%d_STRAIN_-_POSITION_%09d\n", i, read_length, position);
    else printf("_%d_TARGET_SIZE_%d_STRAIN_+_POSITION_%09d\n", i, read_length, position);
    std::cout << read << std::endl;
    std::cout << "+\n";
    std::cout << quality << std::endl;
  }
}

void Help() {
  std::cout << "usage: ./Readsgen [options ...] <reference_path>"
               "[<subreads_path>, <paf_path>]\n"
               "\n"
               "  # default output is stdout\n"
               "  <reference_path>\n"
               "    path to the reference\n"
               "  <subreads_path>\n"
               "    path to the subreads of the reference, requiren when model is set\n"
               "  <alignment_path>\n"
               "    path to the alignements paf file\n"
               "\n"
               "  options:\n"
               "    -m, --model\n"
               "      error model will be trained from the given reads and alignements\n "
               "    -l, --read_length <int>\n"
               "      default: 25000\n"
               "      the length of the reads\n"
               "    -n, --num_reads <int>\n"
               "      default: 10000\n"
               "      the number of readings to be generated\n"
               "    -c, --complement\n"
               "      default: true\n"
               "      reverse complement\n"
               "    -p, --error_probability <double>\n"
               "      default: 0.01\n"
               "      probability of base sequencing error\n"
               "    -s, --seed <int>\n"
               "      random seed\n"
               "    -h, --help\n"
               "      prints the usage\n";
}

std::atomic<std::uint32_t> biosoup::Sequence::num_objects{0};

int main(int argc, char **argv) {

  std::uint32_t read_length = 25000;
  std::uint32_t num_reads = 10000;
  float complement_probability = 0.5;
  float error_probability = 0.01;
  std::string reference_path = "";
  std::uint32_t seed = random();
  bool model = false;

  static struct option options[] = {
      {"read_length", required_argument, nullptr, 'l'},
      {"num_reads", required_argument, nullptr, 'n'},
      {"complement_probability", required_argument, nullptr, 'c'},
      {"error_probability", required_argument, nullptr, 'p'},
      {"seed", required_argument, nullptr, 's'},
      {"model", no_argument, nullptr, 'm'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0}};

  const char *optstr = "l:n:c:p:s:m";
  char arg;
  while ((arg = getopt_long(argc, argv, optstr, options, 0)) != -1) {
    switch (arg) {
    case 'l':
      read_length = std::atoi(optarg);
      break;
    case 'n':
      num_reads = std::atoi(optarg);
      break;
    case 'c':
      complement_probability = std::stod(optarg);
      break;
    case 'p':
      error_probability = std::stod(optarg);
      break;
    case 'm':
      model = true;
      break;
    case 's':
      seed = std::atoi(optarg);
      break;
    case 'h':
      Help();
      return 0;
    default:
      return 1;
    }
  }

  if (argc == 1 || optind == argc) {
    Help();
    return 0;
  }

  reference_path = std::string(argv[optind++]);

  std::string subreads_path = "";
  std::string paf_path = "";

  if (model) {
    if (argc - optind < 2) {
      Help();
      return 0;
    }
    subreads_path = std::string(argv[optind++]);
    paf_path = std::string(argv[optind]);
  }

  biosoup::Sequence::num_objects = 0;
  auto rparser = simulator_utility::CreateParser(reference_path);
  auto ref_pointers = rparser->Parse(-1);
  biosoup::Sequence reference = *(ref_pointers[0].get());

  if (read_length > reference.data.size()) {
    std::cerr << "[Readsgen::] error: wrong read length" << std::endl;
    return 1;
  }

  if (model) train(reference, subreads_path, paf_path);

  generateReads(reference, read_length, num_reads, complement_probability, seed,
                error_probability);
}
