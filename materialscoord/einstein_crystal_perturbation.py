# coding: utf-8

from __future__ import division, unicode_literals

import numpy as np

from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import Site, PeriodicSite


def perturb_einstein_crystal_style(sites, sqrt_kBT_over_kspring=1.0):
    """
    Perturb the position of each site so that the distribution around
    the equilibrium position yields a normal distribution for each
    Cartesian component.  The perturbation complies thus with the
    expectation for an Einstein crystal, in which the potential is
    given by V(dr) = 1/2 * kspring * (dr)^2.  kspring denotes
    the spring constant with which the sites are tethered to their
    equilibrium position, and dr is the distance of the site under
    consideration from its equilibrium position.  The displacements
    are obtained with numpy's random.randn() method, the values of
    which are scaled by sqrt_kBT_over_kspring.
    TODO: move to Structure class!

    Args:
        sites ([Site]): list of Site objects to be perturbed.
        sqrt_kBT_over_kspring (float): width of the underlying normal
            distribution, which is sigma = (kB*T/kspring)^0.5.
    Returns:
        [Site]: copy of input site list that is perturbed.
    """
    sites_perturbed = []
    for s in sites:
        displ = np.random.randn(3)
        new_coords = sqrt_kBT_over_kspring * displ + s.coords
        if isinstance(s, Site):
            sites_perturbed.append(Site(s.species_string, new_coords))
        elif isinstance(s, PeriodicSite):
            sites_perturbed.append(PeriodicSite(
                s.species_string, new_coords, s.lattice,
                coords_are_cartesian=True))
        else:
            raise RuntimeError('expected Site or PeriodicSite.')
    return sites_perturbed


