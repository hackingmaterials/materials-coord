import pandas as pd
pd.set_option("display.max_columns",50) # display all columns

from pymatgen.analysis.local_env import VoronoiNN, VoronoiNN_modified, JMolNN, MinimumDistanceNN, \
                                        MinimumOKeeffeNN, MinimumVIRENN, \
                                        BrunnerNN, EconNN
from materialscoord.cn_methods import HumanInterpreter, TestVoronoiCoordFinder
from materialscoord.core import Benchmark

methods = [MinimumDistanceNN(), HumanInterpreter()]
structure_groups = ["zeolites"]

bm = Benchmark(methods=methods, structure_groups=structure_groups)
bm.benchmark()

print "done"

num_sites = 6
p = bm.report(totals=False, separate_columns=True, max_sites=num_sites)

import json
from collections import Counter

sub_p = None # diff between element-wise cn and human interpreter
sub_p_avg = pd.DataFrame()

algo = ["MinimumDistanceNN"]

cn_dict = {}
for i in range(len(algo)):
    for j in range(num_sites):
        site = p[algo[i]+str(j)]
        hi_site = p["HumanInterpreter"+str(j)]
        for k in range(len(site)):
            temp = Counter(site[k])
            temp.subtract(hi_site[k])
            cn_dict[site.keys()[k]] = json.dumps(dict(temp))
        df = pd.DataFrame.from_dict(cn_dict, orient='index')
        df = df.rename(columns={0:algo[i]+str(j)})
        sub_p = pd.concat([sub_p, df], axis=1)

def deserialize(dictionary):
    for key, val in dictionary.iteritems():
        for i, j in val.iteritems():
            j = json.loads(j)
            j = dict(sorted(j.items()))
            val[i] = j
        dictionary[key] = val
    return dictionary

de_sub_p = deserialize(sub_p)

new_df = num_us(de_sub_p, num_sites)

print new_df