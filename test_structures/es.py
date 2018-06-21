import os
import glob
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

find_structure = glob.glob(os.path.join(module_dir, "misc", "Rb3Cr_mp-975045_computed"+"*"))
print(find_structure)
s = Structure.from_file(find_structure[0])
es = SpacegroupAnalyzer(s).get_symmetrized_structure().equivalent_sites

all = []
for i in es:
    l = []
    for j in i:
        l.append(j.species_string)
    all.append(l)
print(all)
print([x[0] for x in es])


def report(self, totals=False, separate_columns=True):
    Reports
    the
    benchmark as a
    pandas
    dataframe.This is the
    recommended
    method
    for pulling the
        CNs
        obtained
        by
        each
        method.Default
        settings
        are
        used in NBFuncs

        class to calculate
    benchmark
    score.Other
    combinations(ie
    totals = True and separate_columns = True) are
    for
        visual
        purposes
        only.

    :param
    totals(bool): option
    to
    report
    only
    total
    CNs
    of
    a
    site.
    Defaults
    to
    False, meaning
    element - wise
    CN is listed.

:param
separate_columns(bool): option
to
format
the
dataframe
such
that
each
total
CN is in a
separate
column.
Defaults
to
True.
:returns
cn
benchmarks as pandas
dataframe.

TODO: there
's definitely some repeat code... not sure if its worth to go through and create separate functions?
data = {}
for m in self.methods:
    if totals:
        s_dict = {}
        for k in m._cns:
            rev_cns = []
            for i, j in m._cns[k]:
                rev_cns.append((i, sum(j.values())))
            s_dict[k] = rev_cns

        if separate_columns:
            for i in range(self.nsites):
                data[m.__class__.__name__ + str(i)] = {}

            for kc, kv in s_dict.items():
                l = len(kv)
                if l < self.nsites:
                    for i in range(self.nsites - l):
                        kv.append(("null", 0))
                for i in range(self.nsites):
                    data[m.__class__.__name__ + str(i)][kc] = kv[i][1]
        else:
            data[m.__class__.__name__] = s_dict

    else:
        if separate_columns:
            s_dict = {}
            for i in range(self.nsites):
                data[m.__class__.__name__ + str(i)] = {}
            for j in m._cns:
                temp = []
                for s, t in m._cns[j]:
                    if isinstance(t, dict):
                        temp.append((s, t))
                s_dict[j] = temp

                for a, b in s_dict.items():
                    l = len(b)
                    if l < self.nsites:
                        for i in range(self.nsites - l):
                            b.append(("null", {}))
                    for z in range(self.nsites):
                        data[m.__class__.__name__ + str(z)][a] = b[z][1]
        else:
            data[m.__class__.__name__] = m._cns

return pd.DataFrame(data=data)

