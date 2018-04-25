import pandas as pd
from collections import Counter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import yaml
import glob
import os
from pymatgen.core.structure import Structure
import seaborn as sns
import matplotlib.pyplot as plt

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

class nb_funcs(object):

    def __init__(self, df, algo_names, unique_sites=24, elemental=False, cations=True, anions=False):

        self.df = df
        self.algo_names = sorted(algo_names)
        self.unique_sites = unique_sites
        self.elemental = elemental
        self.cations = cations

        c = os.path.join(module_dir, "..", "test_structures", "cat_an.yaml")

        with open(c) as t:
            cats = yaml.load(t)

        self.cats = cats



    def order_cols(self, hi=True):

        if not hi:
            try:
                self.algo_names.remove("HumanInterpreter")
            except:
                pass

        cols = self.df.columns.tolist()

        all_cols = []
        for a in self.algo_names:
            ones_cols = []
            tens_cols = []
            for i in range(len(cols)):
                if cols[i][:-1] == a:
                    ones_cols.append(cols[i])
                if cols[i][:-2] == a:
                    tens_cols.append(cols[i])
            ones_cols.extend(tens_cols)
            all_cols.extend(ones_cols)
        return all_cols

    def mv_df(self):

        mv_df = {}
        for key, val in self.df.items():
            for i in range(self.unique_sites):
                if key == "MinimumVIRENN" + str(i):
                    site = self.df["MinimumVIRENN" + str(i)]
                    for mv_key, mv_val in site.items():
                        ndict = {}
                        for x, y in mv_val.items():
                            x = ''.join(a for a in x if not a.isdigit())
                            x = ''.join(b for b in x if b != "-")
                            x = ''.join(c for c in x if c != "+")
                            ndict[x] = y
                        val[mv_key] = dict(ndict)
            mv_df[key] = val

        mv_df = pd.DataFrame(mv_df)
        return mv_df

    def sub_hi(self):
        """
        element-wise subtraction of human interpreted cn
        from nn algo calculated nn = error in nn algo-calculated cn
        """

        mv_df = self.mv_df()

        sub_hi = {} # initializing dict for diff between element-wise cn and human interpreter
        for i in range(len(self.algo_names)):
            cn_dict = {}
            for j in range(self.unique_sites):
                site = mv_df[self.algo_names[i]+str(j)]
                hi_site = mv_df["HumanInterpreter"+str(j)]
                for k in range(len(site)):
                    ifli = [z for z in hi_site[k].values()]
                    if all(isinstance(z, float) for z in ifli) or len(ifli) == 0:
                        temp = Counter(site[k])
                        temp.subtract(hi_site[k])
                        cn_dict[site.keys()[k]] = dict(temp)
                    else:
                        t = list(site[k].values())[0]
                        lsub = []
                        dsub = {}
                        for coord in list(hi_site[k].values())[0]:
                            tsub = t - coord
                            lsub.append(tsub)
                        min_sub = min(map(abs, lsub))
                        dsub[list(hi_site[k].keys())[0]] = min_sub
                        cn_dict[site.keys()[k]] = dict(dsub)
                sub_hi[self.algo_names[i]+str(j)] = dict(cn_dict)

        sub_hi = pd.DataFrame(sub_hi)
        ordered_cols = self.order_cols(hi=False)
        sub_hi = sub_hi[ordered_cols]
        sub_hi = sub_hi.reindex(self.df.index.tolist())

        return sub_hi

    def abs_df(self):
        """
        abs value of error of nn-algo calculated cn
        """
        sub_hi = self.sub_hi()

        for key, val in sub_hi.items():
            for i, j in val.items():
                if j != {}:
                    abs_val = map(abs, j.values())
                zip_abs_val = dict(zip(j.keys(), abs_val))
                val[i] = zip_abs_val
            sub_hi[key] = val

        abs_df = sub_hi
        abs_df = abs_df.reindex(self.df.index.tolist())
        return abs_df

    def cif_stats(self):

        abs_df = self.abs_df()

        t = os.path.join(module_dir, "..", "test_structures")

        struc = []
        for i in abs_df.index:
            find_structure = glob.glob(os.path.join(t, "*", i + "*"))
            s = Structure.from_file(find_structure[0])
            struc.append(s)

        num_us_list = []
        for j in struc:
            es = SpacegroupAnalyzer(j).get_symmetrized_structure().equivalent_sites
            sites = [len(x) for x in es]
            if len(sites) < self.unique_sites:
                sites.extend([0] * (self.unique_sites - len(sites)))
            num_us_list.append(sites)

        abs_df['num equiv site atoms'] = num_us_list

        mat_dict = {}
        mat_list = []
        counter = 0
        for mat in struc:
            cn_dict = {}
            for line in mat:
                element = line.species_string
                if element not in cn_dict:
                    cn_dict[element] = 1
                else:
                    cn_dict[element] += 1
            mat_dict[abs_df.index[counter]] = dict(cn_dict)
            counter += 1
            mat_list.append(dict(cn_dict))

        if self.cations:
            test = {}
            for key, val in mat_dict.items():
                if self.elemental:
                    for i, j in val.items():
                        test[key] = j
                else:
                    for mat, cat in self.cats.items():
                        if key == mat:
                            for w in list(val):
                                if self.elemental:
                                    test[mat] = val[w]
                                else:
                                    if w in cat:
                                        test[mat] = val[w]
                                        del val[w]

        test = pd.Series(test)
        abs_df['num unit cell atoms'] = test

        cs_df = abs_df
        cs_df = cs_df.reindex(self.df.index.tolist())
        return cs_df

    def mult_equiv(self):

        cs_df = self.cif_stats()

        counter = 0
        for nn in cs_df.keys()[:-2]:
            num_equiv = [[num for num in equiv] for equiv in cs_df['num equiv site atoms']]
            other_counter = 0
            for j in cs_df[nn]:
                j = j.update((x, y*num_equiv[other_counter][counter]) for x, y in j.items())
                other_counter += 1
                if other_counter == int(len(cs_df.index)):
                    other_counter = 0
            counter += 1
            if counter == self.unique_sites:
                counter = 0

        me_df = cs_df
        me_df = me_df.reindex(self.df.index.tolist())
        return me_df

    def tot_algo_sum(self):

        me_df = self.mult_equiv()

        algo_sum = {}
        for a in range(len(self.algo_names)):
            each_algo = me_df.loc[:, self.algo_names[a] + str(0): self.algo_names[a] + str(self.unique_sites - 1)]
            sum_each_algo = {}
            count = 0
            for d in each_algo[:len(each_algo.index)].values:
                if self.cations:
                    l = []
                    for i in range(len(each_algo.keys())):
                        if d[i] != {}:
                            l.append(d[i])
                        #print(l[0])
                        #for k in l[0]:
                        #    print(k)
                    summed = {k: sum(di[k] for di in l) for k in l[0]}
                    #summed = {k: sum(di for di in l[1]) for k in l[0]}
                else:
                    for i in range(len(each_algo.keys())):
                        if d[i] != {}:
                            test = Counter(d[i])
                    summed = {k: sum(di[k] for di in l) for k in l}
                sum_each_algo[each_algo.index[count]] = summed
                count += 1
            algo_sum[self.algo_names[a]] = dict(sum_each_algo)

        totsum = {}
        for mat, su in algo_sum.items():
            new_sum = {}
            for k, v in su.items():
                for i, j in v.items():
                    new_sum[k] = j
            totsum[mat] = dict(new_sum)

        tas_df = pd.DataFrame(totsum)
        tas_df['num unit cell atoms'] = me_df['num unit cell atoms']

        tas_df = tas_df.reindex(self.df.index.tolist())

        return tas_df

    def div_df(self):

        tas_df = self.tot_algo_sum()

        for i in range(len(tas_df.columns)-1):
            mat = tas_df[self.algo_names[i]]
            hi_mat = tas_df["num unit cell atoms"]
            c = mat.divide(hi_mat)

            tas_df[self.algo_names[i]] = c

        div_df = tas_df
        div_df = div_df.reindex(self.df.index.tolist())
        return div_df

    def final(self):

        div_df = self.div_df()

        final_df = div_df.drop('num unit cell atoms', axis=1)
        final_df = final_df.reindex(self.df.index.tolist())

        return final_df