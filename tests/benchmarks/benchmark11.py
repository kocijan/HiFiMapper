from benchmark_utils import Reference, Reads, ResultsHiFiMapper, ResultsMinimap2, ResultsWinnowmap
from benchmark_utils import compare_real_methods
from benchmark_utils import RealEvaluator
from benchmark_utils import CHROMOSOME_19_REFERENCE, CHROMOSOME_19_READS

ref = Reference(CHROMOSOME_19_REFERENCE)
reads = Reads(CHROMOSOME_19_READS)


data = []

res, comm = ResultsHiFiMapper.create(
    ref, reads, sample_length=250, sample_count = 200,
    threads = 256, min_match = 0, quality = 0, frequency = 15,
    discard = True, extended_search=True, sequential = False, secondary_alignements = 1, seed = 175, minimal_anchors = 1)
evalutor = RealEvaluator(str(res), str(reads))
evalutor.draw_coverage(61707364, "images/chromosome19_real_hifi.png", ylim = 100)
data.append([comm, res.get_time(), evalutor.unmapped])
print(comm)
print(data[-1])

res, comm = ResultsMinimap2.create(ref, reads, threads = 256, ssecondary_alignments=1, min_ratio=1)
evalutor = RealEvaluator(str(res), str(reads),  filter_mappings = True)
evalutor.draw_coverage(61707364, "images/chromosome19_real_minimap.png", ylim = 100)
data.append([comm, res.get_time(), evalutor.unmapped])
print(data[-1])


res, comm = ResultsWinnowmap.create(ref, reads, threads = 256, min_ratio = 1)
evalutor = RealEvaluator(str(res), str(reads),  filter_mappings = True)
evalutor.draw_coverage(61707364, "images/chromosome19_real_winnowmap.png", ylim = 100)
data.append([comm, res.get_time(), evalutor.unmapped])
print(data[-1])

print(compare_real_methods(data, tablefmt='latex'))