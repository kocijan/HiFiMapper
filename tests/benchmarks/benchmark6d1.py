from benchmark_utils import (
    Reference,
    Reads,
    ResultsHiFiMapper,
    ResultsMinimap2,
    ResultsWinnowmap,
)
from benchmark_utils import compare_methods
from benchmark_utils import PafEvaluator
from benchmark_utils import CHROMOSOME_13_REFERENCE, CHROMOSOME_13_READS

ref = Reference(CHROMOSOME_13_REFERENCE)
# reads = Reads.generate(
#     ref, read_length=25000, error_probability=0.001, num_reads=600000, seed=7
# )
# reads = Reads(CHROMOSOME_13_READS)
# reads = "/home/mjkocijan/hifimapper/HiFiMapper/tests/benchmarks/data/sd_0001.fixed.limited.fastq"
reads = "/home/mjkocijan/hifimapper/HiFiMapper/tests/benchmarks/data/sd_0001.fixed.fastq"

# THIS BENCHMARK IS THE BEST ONE SO FAR !!!

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


# for sequential in [True, False]:
for sequential in [False]:
    # for extended_search in [True, False]:
    for extended_search in [False]:
        # for discard in [True, False]:
        # for discard in [True]:
        # for sample_length in [150, 200, 250]:
        # for bandwidth in [10, 100, 1000]:
        for bandwidth in [100]:
            res, comm = ResultsHiFiMapper.create(
                ref,
                reads,
                # sample_length=sample_length,
                # sample_length=200,
                sample_length=20,
                # sample_count=200,
                sample_count=50,
                threads=1,
                # threads=256,
                # min_match=0.8,
                min_match=0,
                quality=0,
                frequency=400,
                # discard=discard,
                discard=True,
                extended_search=extended_search,
                sequential=sequential,
                lcp_information=False,
                seed=115,
                secondary_alignements=1000,
                # bandwidth=1000,
                bandwidth=bandwidth,
                minimal_anchors=1,
            )
            evalutor = PafEvaluator(str(res), str(reads))
            # draw ne radi dok je secondary alignments visok
            # evalutor.draw_coverage(113566686, "/home/mjkocijan/images/chromosome13_benchmark6d1_real_hifi_1.png", ylim = 80)
            data.append(
                [
                    comm,
                    res.get_time(),
                    evalutor.correct,
                    evalutor.wrong,
                    evalutor.unmapped,
                ]
            )
            print(data[-1])

# for f in [1.0, 0.1, 0.0]:
# for f in [1.0, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0]:
# for f in [0, 1, 2, 4, 8]:
# for f in [0.5, 0.1, 0.02, 0.004, 0.0008]:
# for f in [0.1, 0.0001, 0.0000001]:
# resm, comm = ResultsMinimap2.create(ref, reads, threads = 256, ssecondary_alignments = 0, min_ratio=1, filter=f)
# resm, comm = ResultsMinimap2.create(ref, reads, threads = 256, ssecondary_alignments = 0, min_ratio=1)
# evalutorm = PafEvaluator(str(resm), str(reads), filter_mappings=True)
# data.append([comm, resm.get_time(), evalutorm.correct, evalutorm.wrong, evalutorm.unmapped])

# for f in [0.1, 0.001, 0.00001]:
# resw, comm = ResultsWinnowmap.create(ref, reads, threads = 256, min_ratio=1, filter=f)
# resw, comm = ResultsWinnowmap.create(ref, reads, threads = 256, min_ratio=1)
# evalutorw = PafEvaluator(str(resw), str(reads), filter_mappings=True)
# data.append([comm, resw.get_time(), evalutorw.correct, evalutorw.wrong, evalutorw.unmapped])

resm, comm = ResultsMinimap2.create(ref, reads, threads = 256, ssecondary_alignments = 1000, min_ratio=1)
#filter mapping je false da bismo dobili sve sekundarne?
evalutorm = PafEvaluator(str(resm), str(reads), filter_mappings=False)
data.append([comm, resm.get_time(), evalutorm.correct, evalutorm.wrong, evalutorm.unmapped])
# draw ne radi dok je secondary alignments visok
# evalutorm.draw_coverage(113566686, "/home/mjkocijan/images/chromosome13_benchmark6d1_real_minimap_1.png", ylim = 140)

print(compare_methods(data, tablefmt="latex"))
