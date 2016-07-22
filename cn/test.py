from cn.cn_methods import TestVoronoiCoordFinder
from cn.core import Benchmark

tvcf_params = {"cutoff": 5.0}

methods = [ TestVoronoiCoordFinder(params=tvcf_params)]
x = Benchmark(methods)
x.benchmark()
x.report()