import unittest

import numpy as np
from pymatgen.core.structure import Structure

from materialscoord.einstein_crystal_perturbation import perturb_einstein_crystal


class EinsteinTest(unittest.TestCase):
    """Test einstein perturbation functionality."""

    def test_perturb(self):
        # basic test to check structures are not equal after perturbation
        structure = Structure.from_spacegroup(
            225,
            [[5.7, 0, 0], [0, 5.7, 0], [0, 0, 5.7]],
            ["Na1+", "Cl1-"],
            [[0, 0, 0], [0.5, 0, 0]],
        )
        perturb_structure = perturb_einstein_crystal(structure)

        orig_sites = [s.frac_coords for s in structure]
        perturb_sites = [s.frac_coords for s in perturb_structure]

        # check coord arrays are not almost equal
        self.assertRaises(
            AssertionError,
            np.testing.assert_array_almost_equal,
            orig_sites,
            perturb_sites,
        )
