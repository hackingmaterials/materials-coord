from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import glob
import os
from pymatgen.core.structure import Structure
from pymatgen.core.sites import Site

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

def num_us(df, num_sites):

    t = os.path.join(module_dir, "..", "test_structures")

    struc = []
    for i in df.index:
        find_structure = glob.glob(os.path.join(t, "*", i+"*"))
        s = Structure.from_file(find_structure[0])
        struc.append(s)

    num_us_list = []
    for j in struc:
        es = SpacegroupAnalyzer(j).get_symmetrized_structure().equivalent_sites
        #print [x.species_string for x in es]
        sites = [len(x) for x in es]
        if len(sites) < num_sites:
            sites.extend([0] * (num_sites - len(sites)))
        num_us_list.append(sites)

    mat_list = []
    for mat in struc:
        cn_dict = {}
        for line in mat:
            element = line.species_string
            if element not in cn_dict:
                cn_dict[element] = 1
            else:
                cn_dict[element] += 1
        mat_list.append(dict(cn_dict))

    df['mat_list'] = mat_list
    df['test'] = num_us_list
    return df


