from benchmark_utils import Reference, Reads, ResultsHiFiMapper, ResultsMinimap2, ResultsWinnowmap
from benchmark_utils import compare_methods
from benchmark_utils import PafEvaluator
from benchmark_utils import CHROMOSOME_13_REFERENCE

ref = Reference(CHROMOSOME_13_REFERENCE)
reads = Reads.generate(ref, read_length = 25000, error_probability = 0.001, num_reads = 600000, seed = 7)


data = []

res,comm = ResultsHiFiMapper.create(
    ref, reads, sample_length=500, sample_count=20,
    threads = 256, min_match = 0.9, quality = 90, frequency = 75,
    discard = True, extended_search=True, sequential = True, lcp_information = False,
    # seed = 115,
    secondary_alignements = 1, minimal_anchors=1)
evalutor = PafEvaluator(str(res), str(reads))
evalutor.draw_corrcet_coverage(113566686, "images/chromosome13_correct_hifi.png")
evalutor.draw_wrong_coverage(113566686, "images/chromosome13_wrong_hifi.png", ylim = 250)
data.append([comm, res.get_time(), evalutor.correct, evalutor.wrong, evalutor.unmapped])
print(data[-1])


resm, comm = ResultsMinimap2.create(ref, reads, threads = 256, ssecondary_alignments = 0, min_ratio=1)
evalutorm = PafEvaluator(str(resm), str(reads), filter_mappings=True)
evalutorm.draw_corrcet_coverage(113566686, "images/chromosome13_correct_minimap.png")
evalutorm.draw_wrong_coverage(113566686, "images/chromosome13_wrong_minimap.png", ylim = 250)
data.append([comm, resm.get_time(), evalutorm.correct, evalutorm.wrong, evalutorm.unmapped])


# resw, comm = ResultsWinnowmap.create(ref, reads, threads = 256, min_ratio=1)
# evalutorw = PafEvaluator(str(resw), str(reads), filter_mappings=True)
# evalutorw.draw_corrcet_coverage(113566686, "images/chromosome13_correct_winnowmap.png")
# evalutorw.draw_wrong_coverage(113566686, "images/chromosome13_wrong_winnowmap.png", ylim = 250)
# data.append([comm, resw.get_time(), evalutorw.correct, evalutorw.wrong, evalutorw.unmapped])


print(compare_methods(data, tablefmt='latex'))

