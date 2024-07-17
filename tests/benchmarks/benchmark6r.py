from benchmark_utils import Reference, Reads, ResultsHiFiMapper, ResultsMinimap2, ResultsWinnowmap
from benchmark_utils import compare_methods
from benchmark_utils import RealEvaluator
from benchmark_utils import CHROMOSOME_13_REFERENCE
from benchmark_utils.table import compare_real_methods
from benchmark_utils import benchmark_parameters

ref = Reference(CHROMOSOME_13_REFERENCE)
# reads = Reads.generate(ref, read_length = 25000, error_probability = 0.001, num_reads = 600000, seed = 7)
reads = Reads.generate(ref, read_length = 25000, error_probability = 0.001, num_reads = 1, seed = 7)


data = []

# for sample_length in [500, 1000]:
#     for sample_count in [20, 50]:
#         for frequency in [75, 150]:
#             for lcp_search_size in [5000]:
#                 res,comm = ResultsHiFiMapper.create(
#                     ref, reads, sample_length=sample_length, sample_count=sample_count,
#                     threads = 256, min_match = 0.9, quality = 90, frequency = frequency,
#                     discard = True, extended_search=False, sequential = False, lcp_information = True,
#                     seed = 115, secondary_alignements = 1, minimal_anchors=1, lcp_search_size = lcp_search_size)
#                 evalutor = PafEvaluator(str(res), str(reads))
#                 data.append([comm, res.get_time(), evalutor.correct, evalutor.wrong, evalutor.unmapped])
#                 print(data[-1])


# for sequential in [True, False]:
#     for extended_search in [True, False]:
#         for discard in [True, False]:
#             res,comm = ResultsHiFiMapper.create(
#                 ref, reads, sample_length=1000, sample_count=20,
#                 threads = 256, min_match = 0.9, quality = 90, frequency = 75,
#                 discard = discard, extended_search=extended_search, sequential = sequential, lcp_information = False,
#                 seed = 115, secondary_alignements = 1, minimal_anchors=1)
#             evalutor = PafEvaluator(str(res), str(reads))
#             data.append([comm, res.get_time(), evalutor.correct, evalutor.wrong, evalutor.unmapped])
#             print(data[-1])

# def evalf(sample_length, sample_count):
#     resh,comm = ResultsHiFiMapper.create(
#                 ref, reads, sample_length=sample_length, sample_count=sample_count,
#                 threads = 256, min_match = 0.9, quality = 90, frequency = 75,
#                 discard = True, extended_search=True, sequential = False, lcp_information = False,
#                 seed = 115, secondary_alignements = 1, minimal_anchors=1)
#     evalutor = RealEvaluator(str(resh), str(reads), filter_mappings=True)
#     return [evalutor.unmapped, round(float(resh.get_mapping_time()), 2)]

# print( 
#     benchmark_parameters(
#     [20, 30, 40, 50, 60, 80, 100, 200, 400],
#     [25, 50, 100, 200, 400],
#     "l",
#     'c',
#     evalf,
#     tablefmt='latex',
#     num_experiments=1)
# )

# exit

# for sequential in [True, False]:
for sequential in [False]:
# for sample_count in [10]
    # for extended_search in [True, False]:
    for extended_search in [True]:
        # for discard in [True, False]:
        for discard in [True]:
            res,comm = ResultsHiFiMapper.create(
                ref, reads, sample_length=50, sample_count=200,
                threads = 256, min_match = 0.9, quality = 90, frequency = 75,
                discard = discard, extended_search=extended_search, sequential = sequential, lcp_information = False,
                seed = 115, secondary_alignements = 1, minimal_anchors=1)
            evalutor = RealEvaluator(str(res), str(reads), filter_mappings=True)
            data.append([comm, res.get_time(), evalutor.total, evalutor.unmapped])
            evalutor.draw_coverage(113566686, "images/chromosome13_real_hifi6_1.png", ylim = 80)
            print(data[-1])



resm, comm = ResultsMinimap2.create(ref, reads, threads = 256, ssecondary_alignments = 0, min_ratio=1)
evalutorm = RealEvaluator(str(resm), str(reads), filter_mappings=True)
data.append([comm, resm.get_time(), evalutorm.total, evalutorm.unmapped])
evalutorm.draw_coverage(113566686, "images/chromosome13_real_minimap6_1.png", ylim = 140)


# resw, comm = ResultsWinnowmap.create(ref, reads, threads = 256, min_ratio=1)
# evalutorw = PafEvaluator(str(resw), str(reads), filter_mappings=True)
# data.append([comm, resw.get_time(), evalutorw.correct, evalutorw.wrong, evalutorw.unmapped])

print(compare_real_methods(data, tablefmt='latex'))

