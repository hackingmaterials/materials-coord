from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from materialscoord.core import CNBase
from materialscoord.external_src.imported import EffectiveCoordFinder_modified, \
    VoronoiCoordFinder_modified


# Sample implementation of pymatgens VoronoiCoordFinder

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