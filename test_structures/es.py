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

    def benchmark(self):
        """
        Calculates CN for each structure site using NN method(s).
        Dictionary of calculated CNs are stored in m._cns

        nsites (int) is used to determine the number of sites each structure has
        and uses the max number of sites as the number of columns in the framework.
        """
        nsites = []
        for m in self.methods:
            for name, structure in self.test_structures.items():
                cns = []
                if self.unique_sites:
                    es = SpacegroupAnalyzer(structure).get_symmetrized_structure().equivalent_sites
                    sites = [structure.index(x[0]) for x in es]
                else:
                    sites = range(len(structure))

                for key, val in self.hi.items():
                    if name == key:
                        for j in sites:
                            if isinstance(m, NearNeighbors):
                                tmpcn = m.get_cn_dict(structure, j, self.use_weights)
                            else:
                                tmpcn = m.compute(structure, j)
                                if tmpcn == "null":
                                    continue
                            if self.nround:
                                self._roundcns(tmpcn, self.nround)
                            cns.append((structure[j].species_string, tmpcn))
                        if self.cation_anion:
                            for mat, cat in self.cations.items():
                                if (name == mat) and cat:
                                    cns = self._popel(cns, cat)
                        elif self.anion_cation:
                            for mat, an in self.anions.items():
                                if name == mat:
                                    cns = self._popel(cns, an)
                    m._cns[name] = cns
                nsites.append(len(cns))
        self.nsites = max(nsites)

