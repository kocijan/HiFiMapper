from benchmark_utils import Reference, Reads, ResultsHiFiMapper, ResultsMinimap2, ResultsWinnowmap
from benchmark_utils import benchmark_parameters
from benchmark_utils import PafEvaluator
from benchmark_utils import get_instructions

REF_SIZE = 25000000

instructions = get_instructions(REF_SIZE, 30000, 0.1, seed = 49)
ref = Reference.generate(reference_size=REF_SIZE, seed = 7, instructions = instructions, probability=0.0001)

reads = Reads.generate(ref, read_length = 25000, error_probability = 0.001, num_reads = 100000, seed = 7)


def evalf(sample_length, sample_count):
    resh, comm = ResultsHiFiMapper.create(
        ref, reads, sample_count=sample_count, sample_length=sample_length,
        threads = 256, min_match=0, secondary_alignements = 1, minimal_anchors=1, discard=True, frequency=3)
    evalutor = PafEvaluator(str(resh), str(reads))
    return [evalutor.correct / (evalutor.total + 0.001), int(resh.get_unmapped()), round(float(resh.get_mapping_time()), 2)]


print( 
    benchmark_parameters(
    [25, 50, 75, 100, 150],
    [10, 20, 50, 70, 100],
    "l",
    'c',
    evalf,
    tablefmt='latex',
    num_experiments=1)
)
