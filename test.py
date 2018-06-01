import pandas as pd
pd.set_option("display.max_columns",999) # display all columns

from pymatgen.analysis.local_env import BrunnerNN_reciprocal, BrunnerNN_relative, BrunnerNN_real, EconNN, JMolNN, \
                                        MinimumDistanceNN, MinimumOKeeffeNN, MinimumVIRENN, \
                                        VoronoiNN, VoronoiNN_modified, CrystalNN
from materialscoord.core import Benchmark, HumanInterpreter

methods = [BrunnerNN_reciprocal(), HumanInterpreter()]

structure_groups = ["common_binaries"]

unique_sites = 4

algo = ["BrunnerNN_reciprocal"]

uw_bm = Benchmark(methods=methods, structure_groups=structure_groups, unique_sites=True,
                  use_weights=False, cations=True)
uw_bm.benchmark()

uw = uw_bm.report(totals=False, separate_columns=True, max_sites=unique_sites)
print(uw)