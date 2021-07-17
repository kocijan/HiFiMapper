#include <bits/stdc++.h>
#include <getopt.h>
#include <sys/time.h>

namespace reference_utils {
void generate_reference(std::string name, int REFERENCE_SIZE, int SEED,
                        const std::map<int, int> &additional_seeds,
                        double probability) {
  std::default_random_engine generator;
  generator.seed(SEED);
  std::default_random_engine generator_err;
  generator_err.seed(SEED * SEED + 1);

  std::uniform_real_distribution<float> fl;

  std::uniform_int_distribution<int> bases_distribution(0, 3);
  std::vector<char> bases = {'A', 'C', 'G', 'T'};
  std::cout << ">" << name << std::endl;
  int print_new_line = 0;
  for (int i = 0; i < REFERENCE_SIZE; i++) {
    print_new_line++;
    if (additional_seeds.count(i)) {
      generator.seed(additional_seeds.at(i));
    }
    int ind = bases_distribution(generator);
    if (fl(generator_err) < probability)
      ind = (ind + bases_distribution(generator_err)) % 4;
    std::cout << bases[ind];
    if (print_new_line == 80) {
      std::cout << std::endl;
      print_new_line = 0;
    }
  }
  printf("\n");
}

void Help() {
  std::cout
      << "usage: ./Refgen [options ...] <reference_size> [instructions ..]\n"
         "\n"
         "  # default output is stdout\n"
         "  <reference_size>\n"
         "    the size of the reference you want to generate\n"
         "  instruction\n"
         "   list of [index seed] values, allows the creation of repetitive segments\n "
         "\n"
         "  options:\n"
         "    -r, --repetitive-counter <int>\n"
         "      default: 0\n"
         "      the number of positions in for which you want to specify a "
         "seed to simulate repetitive segments\n"
         "    -s, --seed <int>\n"
         "      random seed that is set at the beginning of the reference "
         "generation\n"
         "    -p, --probability <double>\n"
         "      the probability that the base will be mutated, to make "
         "repetitive regions distinguishable\n"
         "    -n, --name <string>\n"
         "      reference name\n"
         "    -h, --help\n"
         "      prints the usage\n";
}

} // namespace reference_utils

int main(int argc, char **argv) {
  std::uint32_t r = 0;
  std::uint32_t s = rand();
  double p = 0.01;
  std::string name = "ARTIFICIAL_REF";
  std::uint64_t reference_length = 0;

  static struct option options[] = {
      {"repetitive-counter", required_argument, nullptr, 'r'},
      {"seed", required_argument, nullptr, 's'},
      {"probability", required_argument, nullptr, 'p'},
      {"name", required_argument, nullptr, 'n'},
      {"help", no_argument, nullptr, 'h'},
      {nullptr, 0, nullptr, 0}};

  const char *optstr = "r:s:p:n:";
  char arg;
  while ((arg = getopt_long(argc, argv, optstr, options, 0)) != -1) {
    switch (arg) {
    case 'r':
      r = std::atoi(optarg);
      break;
    case 's':
      s = std::atoi(optarg);
      break;
    case 'p':
      s = std::stod(optarg);
      break;
    case 'n':
      name = std::string(optarg);
      break;
    case 'h':
      reference_utils::Help();
      return 0;
    default:
      return 1;
    }
  }

  if (argc == 1 || optind == argc) {
    reference_utils::Help();
    return 0;
  }

  reference_length = atoi(argv[optind++]);

  std::map<int, int> m;

  if (r > 0) {
    if (argc - optind != r * 2) {
      std::cerr << "[Refgen::] error: wrong number of position seed pairs"
                << std::endl;
      return 1;
    }

    for (auto i = optind; i < argc; i += 2) {
      int index = atoi(argv[i]);
      int seed = atoi(argv[i + 1]);
      m[index] = seed;
    }
  }

  reference_utils::generate_reference(name, reference_length, s, m, p);
}