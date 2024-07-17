from benchmark_utils import Reference, Reads, ResultsHiFiMapper, ResultsMinimap2, ResultsWinnowmap
from benchmark_utils import benchmark_parameters
from benchmark_utils import PafEvaluator
from benchmark_utils import get_instructions

REF_SIZE = 2500000

instructions = get_instructions(REF_SIZE, 3000, 0.1, seed = 49)
ref = Reference.generate(reference_size=REF_SIZE, seed = 7, instructions = instructions, probability=0.001)

reads = Reads.generate(ref, read_length = 2500, error_probability = 0.01, num_reads = 10000, seed = 7)


def evalfDFalse(frequency, min_match):
    resh, comm = ResultsHiFiMapper.create(
        ref, reads, sample_count=10, sample_length=100,
        threads = 256, secondary_alignements = 1, minimal_anchors=1, min_match=min_match, frequency=frequency, discard=False)
    evalutor = PafEvaluator(str(resh), str(reads))
    return [evalutor.correct / (evalutor.total + 0.001), int(resh.get_unmapped()), round(float(resh.get_mapping_time()), 2)]

def evalfDTrue(frequency, min_match):
    resh, comm = ResultsHiFiMapper.create(
        ref, reads, sample_count=10, sample_length=100,
        threads = 256, secondary_alignements = 1, minimal_anchors=1, min_match=min_match, frequency=frequency, discard=True)
    evalutor = PafEvaluator(str(resh), str(reads))
    return [evalutor.correct / (evalutor.total + 0.001), int(resh.get_unmapped()), round(float(resh.get_mapping_time()), 2)]

print("d=True")
print( 
    benchmark_parameters(
    [0, 3, 10, 30, 100],
    [0, 0.1, 0.5, 0.9, 1],
    "f",
    'm',
    evalfDTrue,
    tablefmt='latex',
    num_experiments=1)
)

print("d=False")
print( 
    benchmark_parameters(
    [0, 3, 10, 30, 100],
    [0, 0.1, 0.5, 0.9, 1],
    "f",
    'm',
    evalfDFalse,
    tablefmt='latex',
    num_experiments=1)
)
