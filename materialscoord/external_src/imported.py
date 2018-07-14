from pymatgen.analysis.structure_analyzer import VoronoiCoordFinder
from collections import defaultdict
import math

"""
This module contains CN methods imported from external sources.
See each class for the original authors.
"""


class VoronoiCoordFinder_modified(object):
    """
    Author: S. Bajaj (LBL)
    Modified: M. Aykol (LBL)
    """
    def __init__(self, structure, n):
        self._structure = structure
        self.n = n

    def get_cns(self):
        siteno = self.n
        try:
            vor = VoronoiCoordFinder(self._structure).get_voronoi_polyhedra(siteno)
            weights = VoronoiCoordFinder(self._structure).get_voronoi_polyhedra(siteno).values()
        except RuntimeError as e:
            print e

        coordination = {}
        max_weight = max(weights)
        for v in vor:
            el = v.species_string
            if vor[v] > 0.50 * max_weight:
                if el in coordination:
                    coordination[el]+=1
                else:
                    coordination[el] = 1
        return coordination



class EffectiveCoordFinder_modified(object):

    """
    Author: S. Bajaj (LBL)
    Modified: M. Aykol (LBL)
    Finds the average effective coordination number for each cation in a given structure. It
    finds all cation-centered polyhedral in the structure, calculates the bond weight for each peripheral ion in the
    polyhedral, and sums up the bond weights to obtain the effective coordination number for each polyhedral. It then
    averages the effective coordination of all polyhedral with the same cation at the central site.
    We use the definition from Hoppe (1979) to calculate the effective coordination number of the polyhedrals:
    Hoppe, R. (1979). Effective coordination numbers (ECoN) and mean Active fictive ionic radii (MEFIR).
    Z. Kristallogr. , 150, 23-52.
    ECoN = sum(exp(1-(l_i/l_av)^6)), where l_av = sum(l_i*exp(1-(1_i/l_min)))/sum(exp(1-(1_i/l_min)))
    """

    def __init__(self, structure, n):
        self._structure = structure
        self.n = n

    def get_cns(self, radius=10.0):
        """
        Get a specie-centered polyhedra for a structure
        :param radius: (float) distance in Angstroms for bond cutoff
        :return: (dict) A dictionary with keys corresponding to different ECoN coordination numbers for site n.
        """
        site = self._structure.sites[self.n]

        all_bond_lengths = []
        neighbor_list = []
        bond_weights = []
        for neighbor in self._structure.get_neighbors(site, radius):  # entry = (site, distance)
            if neighbor[1] < radius:
                all_bond_lengths.append(neighbor[1])
                neighbor_list.append(neighbor[0].species_string)

        weighted_avg = calculate_weighted_avg(all_bond_lengths)
        cns = {}
        for i in neighbor_list:
            cns[i]=0.0

        for bond in range(len(all_bond_lengths)):
            bond_weight = math.exp(1-(all_bond_lengths[bond]/weighted_avg)**6)
            cns[ neighbor_list[bond]]+=bond_weight
        return cns


def calculate_weighted_avg(bonds):
    """
    Author: S. Bajaj (LBL)
    Get the weighted average bond length given by the effective coordination number formula in Hoppe (1979)
    :param bonds: (list) list of floats that are the bond distances between a cation and its peripheral ions
    :return: (float) exponential weighted average
    """
    minimum_bond = min(bonds)
    weighted_sum = 0.0
    total_sum = 0.0
    for entry in bonds:
        weighted_sum += entry*math.exp(1 - (entry/minimum_bond)**6)
        total_sum += math.exp(1-(entry/minimum_bond)**6)
    return weighted_sum/total_sum