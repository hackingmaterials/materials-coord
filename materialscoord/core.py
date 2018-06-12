import yaml
import abc
import os
import glob
import pandas as pd
from collections import Counter
from collections import OrderedDict
import seaborn as sns
import matplotlib.pyplot as plt
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.local_env import NearNeighbors

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

class Benchmark(object):
    """
    Class for performing CN benchmarks on a set of structures using the selected set of methods.
    :param methods: (list) CN methods. All methods must be subclassed from CNBase.
    :param structure_groups: (str) or (list) groups of test structures. Defaults to "elemental"
            Current options include "elemental", "common_binaries", "laves", but will be
            significantly expanded in future.
    :param custom_set (str): Full path to custom set of external structures to be loaded. Can be used to
            apply CN methods on user specified structures.
    :param unique_sites: (bool) Only calculate CNs of symmetrically unique sites in structures.
            This is essential to get a cleaner output. Defaults to True.
    :param nround: (int) Rounds CNs to given number of decimals. Defaults to 3. nround=0 means
            no rounding.
    """

    def __init__(self, methods, structure_groups="elemental", custom_set=None,
                 unique_sites=True, nround=2, use_weights=False,
                 cation_anion=False, anion_cation=False):

        self.methods = methods
        self.structure_groups = structure_groups if isinstance(structure_groups, list) else [structure_groups]

        self.test_structures = OrderedDict()
        self.cations = OrderedDict()
        self.anions = OrderedDict()
        if custom_set:
            self.structure_groups = None
            self.custom = custom_set
            self._load_test_structures(None)
        else:
            for g in self.structure_groups:
                self._load_test_structures(g)

        self.unique_sites = unique_sites
        self.nround = nround
        self.use_weights = use_weights
        self.cation_anion = cation_anion
        self.anion_cation = anion_cation

        for m in self.methods:
            assert isinstance(m, (NearNeighbors, CNBase))
            m._cns = {}
        print("Initialization successful.")

        p = os.path.join(module_dir, "..", "test_structures", "human_interpreter.yaml")

        with open(p) as f:
            hi = yaml.load(f)
        self.hi = hi

    def _load_test_structures(self, group):
        """
        Loads the structure group from test_structures
        :param group: (str) group name, options: "elemental". Defaults to "elemental"
        """
        if self.structure_groups:
            p = os.path.join(module_dir, "..", "test_structures", group, "*")
        else:
            if os.path.isdir(self.custom):
                p = os.path.join(self.custom, "*")
            else:
                p = self.custom
        str_files = glob.glob(p)
        for s in str_files:
            name = os.path.basename(s).split(".")[0]
            structure = Structure.from_file(s)

            cations = []
            anions = []
            for i in structure:
                i = str(i).split(']', 1)[1]
                if i.endswith('+'):
                    if '0' in i: # metals
                        pass
                    else:
                        el = ''.join(x for x in str(i) if x.isalpha())
                        if el not in cations:
                            cations.append(el)
                if str(i).endswith('-'):
                    el = ''.join(x for x in str(i) if x.isalpha())
                    if el not in anions:
                        anions.append(el)
            self.cations[name] = cations
            self.anions[name] = anions

            structure.remove_oxidation_states()
            self.test_structures[name] = structure

    def benchmark(self):
        """
        Performs the benchmark calculations.
        """
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
                                if name == mat:
                                    if not cat: # metals
                                        pass
                                    else:
                                        cns = [i for i in cns if i[0] in cat]
                                        for tup in cns:
                                            if tup[0] in tup[1].keys():
                                                tup[1].pop(tup[0])
                        if self.anion_cation:
                            for mat, an in self.anions.items():
                                if name == mat:
                                    cns = [i for i in cns if i[0] in an]
                                    for tup in cns:
                                        if tup[0] in tup[1].keys():
                                            tup[1].pop(tup[0])
                    m._cns[name] = cns

    def report(self, totals=False, separate_columns=True, max_sites=5):
        """
        Reports the benchmark as a pandas DataFrame. This is the recommended method for pulling the
        CNs obtained by each method.
        :param totals: (bool) option to report only total CNs of a site. Defaults to False, meaning element-wise CN
            is listed.
        :param separate_columns: (bool) option to format the data-frame such that each total CN is in a separate column.
            totals must be set True.
        :param max_sites: (int) maximum number of unique sites to have in report for each method. Only used if totals
            and separate_columns are set True. Defaults to 5

        """
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
                    for i in range(max_sites):
                        data[m.__class__.__name__ + str(i)] = {}

                    for kc, kv in s_dict.items():
                        l = len(kv)
                        if l < max_sites:
                            for i in range(max_sites - l):
                                kv.append(("null", 0))
                        for i in range(max_sites):
                            data[m.__class__.__name__ + str(i)][kc] = kv[i][1]
                else:
                    data[m.__class__.__name__] = s_dict

            else:
                if separate_columns:
                    s_dict = {}
                    for i in range(max_sites):
                        data[m.__class__.__name__ + str(i)] = {}
                    for j in m._cns:
                        temp = []
                        for s, t in m._cns[j]:
                            if isinstance(t, dict):
                                temp.append((s, t))
                        s_dict[j] = temp

                        for a, b in s_dict.items():
                            l = len(b)
                            if l < max_sites:
                                for i in range(max_sites - l):
                                    b.append(("null", {}))
                            for z in range(max_sites):
                                data[m.__class__.__name__ + str(z)][a] = b[z][1]
                else:
                    data[m.__class__.__name__] = m._cns

        return pd.DataFrame(data=data)

        #index = self.test_structures.keys()
        #return pd.DataFrame(data=data, index=list(index)) <-- doesn't work for some reason???

    @staticmethod
    def _roundcns(d, ndigits):
        """
        rounds all values in a dict to ndigits
        """
        for k,v in d.items():
            if isinstance(v, list):
                pass
            else:
                d[k]=round(v,ndigits)

    #@staticmethod
    #def _ionlist(el):


class NbFuncs(Benchmark):

    def __init__(self, Benchmark, unique_sites):

        self.df = Benchmark.report()
        self.methods = Benchmark.methods
        self.test_structures = Benchmark.test_structures
        self.cation_anion = Benchmark.cation_anion
        self.anion_cation = Benchmark.anion_cation
        self.nround = Benchmark.nround
        self.unique_sites = unique_sites

        self.nohi = [i for i in list(self.df.columns) if 'HumanInterpreter' not in i]
        self.hi = [i for i in list(self.df.columns) if 'HumanInterpreter' in i]

    def sub_hi(self):
        """
        element-wise subtraction of human interpreted cn
        from nn algo calculated nn = error in nn algo-calculated cn
        """

        df = {}
        for i in range(len(self.nohi)):
            cndict = {}
            for j in range(self.unique_sites):
                if str(j) in self.nohi[i]:
                    site = self.df[self.nohi[i]]
                    hisite = self.df[self.hi[j]]
                    for k in range(len(site)):
                        coord = [z for z in hisite[k].values()]
                        if all(isinstance(z, float) for z in coord) or not coord:
                            temp = Counter(site[k])
                            temp.subtract(hisite[k])
                            cndict[site.keys()[k]] = dict(temp)
                        else:
                            t = list(site[k].values())[0]
                            lsub = []
                            dsub = {}
                            for cn in list(hisite[k].values())[0]:
                                tsub = t - cn
                                lsub.append(tsub)
                            minsub = min(map(abs, lsub))
                            dsub[list(hisite[k].keys())[0]] = minsub
                            cndict[site.keys()[k]] = dict(dsub)
                    df[self.nohi[i]] = dict(cndict)

        return pd.DataFrame(df)

    def abs_df(self):
        """
        abs value of error of nn-algo calculated cn
        """
        df = self.sub_hi()

        for key, val in df.items():
            for i, j in val.items():
                if j != {}:
                    absval = map(abs, j.values())
                else:
                    absval = {}
                zip_absval = dict(zip(j.keys(), absval))
                val[i] = zip_absval
            df[key] = val

        return df

    def cif_stats(self):

        df = self.abs_df()

        ts = os.path.join(module_dir, "..", "test_structures")

        structures = []
        for i in df.index:
            find_structure = glob.glob(os.path.join(ts, "*", i + "*"))
            s = Structure.from_file(find_structure[0])
            structures.append(s)

        us = []
        for i in structures:
            es = SpacegroupAnalyzer(i).get_symmetrized_structure().equivalent_sites
            if self.cation_anion:
                ces = [x for x in es if '+' in x[0].species_string]
                sites = [len(x) for x in ces]
            else:
                aes = [x for x in es if '-' in x[0].species_string]
                sites = [len(x) for x in aes]
            if len(sites) < self.unique_sites:
                sites.extend([0] * (self.unique_sites - len(sites)))
            us.append(sites)

        summed = []
        for i in us:
            summed.append(sum(i))

        if self.cation_anion:
            df['unique site cations'] = us
            df['total cations'] = summed
        else:
            df['unique site anions'] = us
            df['total anions'] = summed

        return df

    def mult_equiv(self):

        df = self.cif_stats()

        counter = 0
        for nn in df.keys()[:-2]:
            if self.cation_anion:
                num_equiv = [[num for num in equiv] for equiv in df['unique site cations']]
            else:
                num_equiv = [[num for num in equiv] for equiv in df['unique site anions']]
            other_counter = 0
            for j in df[nn]:
                j = j.update((x, [y]*num_equiv[other_counter][counter]) for x, y in j.items())
                other_counter += 1
                if other_counter == int(len(df.index)):
                    other_counter = 0
            counter += 1
            if counter == self.unique_sites:
                counter = 0

        return df

    def merge(self):

        df = self.mult_equiv()

        merged = {}
        for m in self.methods:
            if m.__class__.__name__ == 'HumanInterpreter':
                pass
            else:
                algo = df[[i for i in list(df.columns) if m.__class__.__name__ in i]]
                extended = {}
                count = 0
                for a in algo[:len(algo.index)].values:
                    extension = {}
                    for i in a:
                        for key, val in i.items():
                            if key in extension.keys():
                                extension[key].extend(val)
                            else:
                                extension[key] = val
                    extended[algo.index[count]] = extension
                    count += 1
                merged[m.__class__.__name__] = dict(extended)

        return pd.DataFrame(merged)

    def total(self):

        df = self.merge()

        totsum = {}
        for algo, mats in df.items():
            matsum = {}
            for mat, coord in mats.items():
                for coords, li in coord.items():
                    sumli = sum(li)
                    matsum[mat] = sumli
                totsum[algo] = dict(matsum)

        return pd.DataFrame(totsum)

    def div(self):

        df = pd.concat([self.total(), self.cif_stats()[self.cif_stats().columns[-1:]]], axis=1)

        for m in self.methods:
            if m.__class__.__name__ == "HumanInterpreter":
                pass
            else:
                algo = df[m.__class__.__name__]
                if self.cation_anion:
                    div = algo.divide(df['total cations'])
                else:
                    div = algo.divide(df['total anions'])
                df[m.__class__.__name__] = div

        if self.cation_anion:
            df = df.drop('total cations', axis=1)
        else:
            df = df.drop('total anions', axis=1)

        df = df.round(self.nround)

        return df

    def final(self):

        df = self.div()

        fig, ax = plt.subplots(figsize=(20, 10))

        sns.set(font='Times New Roman')
        sns.set(font_scale=1)

        hm = sns.heatmap(df, annot=True, cmap="BuPu", vmax=10)

        ax.set_xticklabels(df.columns.tolist(), rotation=60)
        ax.set_yticklabels(df.index.tolist())

        plt.show()

class CNBase:
    __metaclass__ = abc.ABCMeta
    """
    This is an abstract base class for implementation of CN algorithms. All CN methods
    must be subclassed from this class, and have a compute method that returns CNs as
    a dict.
    """

    def __init__(self, params=None):
        """
        :param params: (dict) of parameters to pass to compute method.
        """
        self._params = params if params else {}
        self._cns = {}

    @abc.abstractmethod
    def compute(self, structure, n):
        """
        :param structure: (Structure) a pymatgen Structure
        :param n: (int) index of the atom in structure that the CN will be calculated
            for.
        :return: Dict of CN's for the site n. (e.g. {'O': 4.4, 'F': 2.1})
        """
        pass

class HumanInterpreter(CNBase):
    """
    This is a special CN method that reads a yaml file where "human interpretations" of coordination
    numbers are given.
    """
    def __init__(self, custom_interpreter=None, custom_test_structures=None):

        p = os.path.join(module_dir, "..", "test_structures", "human_interpreter.yaml")

        t = os.path.join(module_dir, "..", "test_structures")
        interpreter = custom_interpreter if custom_interpreter else p
        test_structures = custom_test_structures if custom_test_structures else t

        with open(interpreter) as f:
            hi = yaml.load(f)

        for k in hi.keys():

            hi[k].append(len(hi[k]))

            if custom_test_structures:
                find_structure = glob.glob(os.path.join(test_structures, k + "*"))
            else:
                find_structure = glob.glob(os.path.join(test_structures, "*", k+"*"))
            if len(find_structure)==0:
                continue

            s = Structure.from_file(find_structure[0])
            s.remove_oxidation_states()
            hi[k].append(s)

        super(HumanInterpreter, self).__init__(params=hi)

    def compute(self, structure, n):
        for v in self._params.values():
            #print(v)
            if structure == v[-1]:
                if len(v[:-1]) != len(v[-1]):
                    # means possibly reduced structure is used by human interpreter
                    # therefore get the equivalent sites and replace n with its
                    # index in the list of unique sites
                    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
                    es = SpacegroupAnalyzer(structure).get_symmetrized_structure().equivalent_sites
                    #print [x[0] for x in es]
                    sites = [structure.index(x[0]) for x in es]
                    n = sites.index(n)
                if n < v[-2]:
                    cn = list(v[n].values())[0]
                    return cn
        return "null"