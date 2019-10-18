"""Functions for perturbing structures according to the Einstein crystal test rig."""

from copy import deepcopy

import numpy as np

from pymatgen import Structure


def perturb_einstein_crystal(structure: Structure, sigma: float = 1.0):
    """
    Perturb structure according to the Einstein crystal model.

    Each site is perturbed so that the distribution around the equilibrium position
    yields a normal distribution for each Cartesian component.

    The perturbation complies thus with the expectation for an Einstein crystal,
    in which the potential is given by V(dr) = 1/2 * kspring * (dr)^2.
    kspring denotes the spring constant with which the sites are tethered to
    their equilibrium position, and dr is the distance of the site under
    consideration from its equilibrium position.

    Args:
        structure: A structure.
        sigma: Width of the underlying normal distribution, equivalent to
            (kB*T/kspring)^0.5.

    Returns:
        A crystal structure with the atoms displaced.
    """
    structure = deepcopy(structure)

    for i in range(len(structure)):
        displacement = np.random.randn(3) * sigma
        structure.translate_sites([i], displacement, frac_coords=False)

    return structure
