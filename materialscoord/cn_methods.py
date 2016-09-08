import yaml
import os
import glob
from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from materialscoord.core import CNBase
from materialscoord.external_src.imported import EffectiveCoordFinder_modified, \
    VoronoiCoordFinder_modified
from pymatgen.core.structure import Structure

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))

# Sample implementations of pymatgens VoronoiCoordFinder (O'Keeffe's method), a modification
# of that method, and Hoppe's ECoN method.

class TestVoronoiCoordFinder(CNBase):
    """
    O'Keefe's CN's as implemented in pymatgen.

    Note: We only need to define the compute method using the CNBase class
        that returns a dict object.
        params can be accessed via ._params and passed to actual function.
    """
    def compute(self, structure, n):
        params = self._params
        vor = VoronoiCoordFinder(structure, **params)
        vorp = vor.get_voronoi_polyhedra(n)
        cdict = {}
        for i in vorp:
            if i.species_string not in cdict:
                cdict[i.species_string] = vorp[i]
            else:
                cdict[i.species_string] += vorp[i]
        return cdict


class TestECoN(CNBase):
    """
    ECoN
    """
    def compute(self, structure, n):
        params = self._params
        x = EffectiveCoordFinder_modified(structure, n)
        return x.get_cns(**params)


class TestVoronoiCoordFinder_mod(CNBase):
    """
    Modified VoronoiCoordFinder that considers only neighbors
    with at least 50% weight of max(weight).
    """
    def compute(self, structure, n):
        params = self._params
        x = VoronoiCoordFinder_modified(structure, n)
        return x.get_cns(**params)


class TestBrunnerReciprocal(CNBase):
    """
    Brunner's CN described as counting the atoms that are within the largest reciprocal interatomic distance.

    From the paper: G.O. Brunner, A definition of coordination and its relevance in structure types AlB2 and NiAs.
        Acta Crys. A33 (1977) 226.

    """
    def compute(self, structure, n):
        params = self._params
        tol = params.get("tol", 1.0e-4)
        r = params.get("radius", 8.0)
        nl = structure.get_neighbors(structure.sites[n], r)
        ds = [i[-1] for i in nl]
        ds.sort()
        ns = [1.0/ds[i] - 1.0/ds[i+1] for i in range(len(ds) - 1)]
        d_max = ds[ ns.index(max(ns)) ]
        cn = {}
        for i in nl:
            if i[-1] < d_max + tol:
                el = i[0].species_string
                if el in cn:
                    cn[el] += 1.0
                else:
                    cn[el] = 1.0
        return cn


class HumanInterpreter(CNBase):
    """
    This is a special CN method that reads a yaml file where "human interpretations" of coordination
    numbers are given.
    """
    def __init__(self, custom_interpreter=None, custom_test_structures=None, average_ranges=True):

        p = os.path.join(module_dir, "..", "test_structures", "human_interpreter.yaml")
        t = os.path.join(module_dir, "..", "test_structures")
        interpreter = custom_interpreter if custom_interpreter else p
        test_structures = custom_test_structures if custom_test_structures else t

        with open(interpreter) as f:
            hi = yaml.load(f)

        for k in hi.keys():
            if custom_test_structures:
                find_structure = glob.glob(os.path.join(test_structures, k + "*"))
            else:
                find_structure = glob.glob(os.path.join(test_structures, "*", k+"*"))
            if len(find_structure)==0:
                continue

            s = Structure.from_file(find_structure[0])
            hi[k].append(s)

        self._average_ranges = average_ranges
        super(HumanInterpreter, self).__init__(params=hi)

    def compute(self, structure, n):

        for v in self._params.values():
            if structure == v[-1]:
                if len(v[:-1]) != len(v[-1]):
                    # means possibly reduced structure is used by human interpreter
                    # therefore get the equivalent sites and replace n with its
                    # index in the list of unique sites
                    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
                    es = SpacegroupAnalyzer(structure).get_symmetrized_structure().equivalent_sites
                    sites = [structure.index(x[0]) for x in es]
                    n = sites.index(n)
                cn = v[n].values()[0]
                for i,j in cn.items():
                    if isinstance(j,list) and self._average_ranges:
                        cn[i] = sum(j) / float(len(j))
                return cn
        return "null"
