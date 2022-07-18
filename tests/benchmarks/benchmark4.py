from benchmark_utils import Reference, Reads, ResultsHiFiMapper, ResultsMinimap2, ResultsHiFiMapper
from benchmark_utils import benchmark_parameters
from benchmark_utils import PafEvaluator
from benchmark_utils import get_instructions

REF_SIZE = 25000000

instructions = get_instructions(REF_SIZE, 30000, 0.1, seed = 49)
ref = Reference.generate(reference_size=REF_SIZE, seed = 7, instructions = instructions, probability=0.0001)

reads = Reads.generate(ref, read_length = 25000, error_probability = 0.001, num_reads = 100000, seed = 7)


def evalf(sample_length, sample_count):
    resh,comm = ResultsHiFiMapper.create(ref, 
    reads, sample_count=sample_count, sample_length=sample_length, threads = 256,
    min_match=0, frequency=1, lcp_information=True, lcp_search_size = 5000)
    evalutor = PafEvaluator(str(resh), str(reads))
    print(comm)
    print([evalutor.correct / evalutor.total, int(resh.get_unmapped()), round(float(resh.get_mapping_time()), 2)])
    return [evalutor.correct / evalutor.total, int(resh.get_unmapped()), round(float(resh.get_mapping_time()), 2)]


print( 
    benchmark_parameters(
    [150, 500, 1000],
    [20, 40],
    "l",
    'c',
    evalf,
    tablefmt='latex',
    num_experiments=1)
)
