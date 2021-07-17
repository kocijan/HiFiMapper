
from benchmark_utils import Reference, Reads, ResultsHiFiMapper, ResultsMinimap2, ResultsWinnowmap
from benchmark_utils import compare_methods
from benchmark_utils import PafEvaluator
from benchmark_utils import CHROMOSOME_19_REFERENCE

ref = Reference(CHROMOSOME_19_REFERENCE)
reads = Reads.generate(ref, read_length = 25000, error_probability = 0.001, num_reads = 600000, seed = 7)


data = []

for sample_length in [100, 500, 1000]:
    for sample_count in [20, 50]:
        for frequency in [15, 100]:
            for lcp_search_size in [1000, 5000]:
                res,comm = ResultsHiFiMapper.create(
                    ref, reads, sample_length=sample_length, sample_count=sample_count,
                    threads = 256, min_match = 0.9, quality = 90, frequency = frequency,
                    discard = True, extended_search=False, sequential = False, lcp_information = True,
                    seed = 115, secondary_alignements = 1, secondary_to_primary_ratio=0.95, minimal_anchors=1, lcp_search_size = lcp_search_size)
                evalutor = PafEvaluator(str(res), str(reads))
                data.append([comm, res.get_time(), evalutor.correct, evalutor.wrong, evalutor.unmapped])
                print(data[-1])


for sequential in [True, False]:
    for extended_search in [True, False]:
        for discard in [True, False]:
            res,comm = ResultsHiFiMapper.create(
                ref, reads, sample_length=1000, sample_count=20,
                threads = 256, min_match = 0.9, quality = 90, frequency = 15,
                discard = discard, extended_search=extended_search, sequential = sequential, lcp_information = False,
                seed = 115, secondary_alignements = 1, secondary_to_primary_ratio=0.95, minimal_anchors=1)
            evalutor = PafEvaluator(str(res), str(reads))
            data.append([comm, res.get_time(), evalutor.correct, evalutor.wrong, evalutor.unmapped])
            print(data[-1])


for sequential in [True, False]:
    for extended_search in [True, False]:
        for discard in [True, False]:
            res,comm = ResultsHiFiMapper.create(
                ref, reads, sample_length=500, sample_count=20,
                threads = 256, min_match = 0.9, quality = 90, frequency = 15,
                discard = discard, extended_search=extended_search, sequential = sequential, lcp_information = False,
                seed = 115, secondary_alignements = 1, secondary_to_primary_ratio=0.95, minimal_anchors=1)
            evalutor = PafEvaluator(str(res), str(reads))
            data.append([comm, res.get_time(), evalutor.correct, evalutor.wrong, evalutor.unmapped])
            print(data[-1])



resm, comm = ResultsMinimap2.create(ref, reads, threads = 256, ssecondary_alignments = 0, min_ratio=1)
evalutorm = PafEvaluator(str(resm), str(reads), filter_mappings = True)
data.append([comm, resm.get_time(), evalutorm.correct, evalutorm.wrong, evalutorm.unmapped])


resw, comm = ResultsWinnowmap.create(ref, reads, threads = 256, min_ratio=1)
evalutorw = PafEvaluator(str(resw), str(reads),  filter_mappings = True)
data.append([comm, resw.get_time(), evalutorw.correct, evalutorw.wrong, evalutorw.unmapped])


print(compare_methods(data, tablefmt='latex'))
