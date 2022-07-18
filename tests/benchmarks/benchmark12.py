from benchmark_utils import Reference, Reads, ResultsHiFiMapper, ResultsMinimap2, ResultsWinnowmap
from benchmark_utils import compare_real_methods
from benchmark_utils import RealEvaluator
from benchmark_utils import CHROMOSOME_13_REFERENCE, CHROMOSOME_13_READS

ref = Reference(CHROMOSOME_13_REFERENCE)
reads = Reads(CHROMOSOME_13_READS)


data = []


res, comm = ResultsHiFiMapper.create(
    ref, reads, sample_length=250, sample_count = 200,
    threads = 256, min_match = 0, quality = 0, frequency = 75,
    discard = True, extended_search=True, sequential = False, secondary_alignements = 1, seed = 175, minimal_anchors = 1)
evalutor = RealEvaluator(str(res), str(reads))
data.append([comm, res.get_time(), evalutor.unmapped])
evalutor.draw_coverage(113566686, "images/chromosome13_real_hifi.png", ylim = 800)

res, comm = ResultsMinimap2.create(ref, reads, threads = 256, ssecondary_alignments=1, min_ratio=1)
evalutor = RealEvaluator(str(res), str(reads), filter_mappings = True)
evalutor.draw_coverage(113566686, "images/chromosome13_real_minimap.png", ylim = 800)
data.append([comm, res.get_time(), evalutor.unmapped])


res, comm = ResultsWinnowmap.create(ref, reads, threads = 256, min_ratio = 1)
evalutor = RealEvaluator(str(res), str(reads), filter_mappings = True)
evalutor.draw_coverage(113566686, "images/chromosome13_real_winnowmap.png", ylim = 800)
data.append([comm, res.get_time(), evalutor.unmapped])

print(compare_real_methods(data, tablefmt='latex'))
