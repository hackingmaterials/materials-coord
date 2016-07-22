from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from cn.core import CNBase


# Sample implementation of pymatgens VoronoiCoordFinder

class TestVoronoiCoordFinder(CNBase):
    """
    We only need to define a compute method using the CNBase class
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