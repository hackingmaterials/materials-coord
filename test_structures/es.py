import os
import glob
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

find_structure = glob.glob(os.path.join(module_dir, "clus", "Al1H10O5"+"*"))
#print(find_structure)
s = Structure.from_file(find_structure[0])
print(s)
es = SpacegroupAnalyzer(s).get_symmetrized_structure().equivalent_sites

all = []
for i in es:
    l = []
    for j in i:
        l.append(j.species_string)
    all.append(l)
#print(all)
#print(len(all))
#print([x[0] for x in es])

