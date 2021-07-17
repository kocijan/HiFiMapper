from benchmark_utils.proxies import (
    Reference,
    Reads,
    ResultsHiFiMapper,
    ResultsMinimap2,
    ResultsWinnowmap,
)

from benchmark_utils.table import benchmark_parameters, compare_methods, compare_real_methods

from benchmark_utils.evaluator import PafEvaluator, RealEvaluator

from benchmark_utils.utils import get_instructions

from benchmark_utils.data import CHROMOSOME_13_READS, CHROMOSOME_13_REFERENCE, CHROMOSOME_19_READS, CHROMOSOME_19_REFERENCE
