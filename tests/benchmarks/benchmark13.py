from benchmark_utils import Reference, Reads, ResultsHiFiMapper, ResultsMinimap2, ResultsWinnowmap
from benchmark_utils import compare_methods
from benchmark_utils import PafEvaluator
from benchmark_utils import CHROMOSOME_19_REFERENCE

ref = Reference(CHROMOSOME_19_REFERENCE)
reads = Reads.generate(ref, read_length = 25000, error_probability = 0.001, num_reads = 600000, seed = 7)


data = []

for sample_length in [50, 100, 250, 500]:
    sample_count = 5 * (10000 // sample_length)
    res, comm = ResultsHiFiMapper.create(
        ref, reads, sample_length=sample_length, sample_count = int(sample_count),
        threads = 256, min_match = 0, quality = 0, frequency = 15,
        discard = True, extended_search=True, sequential = False, secondary_alignements = 1, seed = 175, minimal_anchors = 1)
    evalutor = PafEvaluator(str(res), str(reads))
    data.append([comm, res.get_time(), evalutor.correct, evalutor.wrong, evalutor.unmapped])
    print(comm)
    print(data[-1])

for sample_length in [50, 100, 250, 500]:
    sample_count = 5 * (10000 // sample_length)
    res, comm = ResultsHiFiMapper.create(
        ref, reads, sample_length=sample_length, sample_count = int(sample_count),
        threads = 256, min_match = 0, quality = 90, frequency = 15,
        discard = True, extended_search=True, sequential = False, secondary_alignements = 1, seed = 175, minimal_anchors = 1)
    evalutor = PafEvaluator(str(res), str(reads))
    data.append([comm, res.get_time(), evalutor.correct, evalutor.wrong, evalutor.unmapped])
    print(comm)
    print(data[-1])

print(compare_methods(data, tablefmt='latex'))
