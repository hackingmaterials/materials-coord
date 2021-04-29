import unittest
from copy import deepcopy

import numpy as np
from pymatgen.analysis.local_env import EconNN, MinimumVIRENN, VoronoiNN
from pymatgen.core.structure import Structure

from materialscoord.core import Benchmark


class BenchmarkTest(unittest.TestCase):
    """Test the Benchmark class."""

    def setUp(self):
        # set up a test structure, the coordinations are not correct and are only
        # for test purposes.
        structure = Structure.from_spacegroup(
            225,
            [[5.7, 0, 0], [0, 5.7, 0], [0, 0, 5.7]],
            ["Na1+", "Cl1-"],
            [[0, 0, 0], [0.5, 0, 0]],
        )
        all_sites_coordination = [{"Cl": 6}] * 4 + [{"Na": 8, "Cl": [6, 8]}] * 4
        structure.add_site_property("coordination", all_sites_coordination)
        self.structures = {"test_structure": structure}
        self.nn_methods = [MinimumVIRENN(), EconNN()]

    def test_from_structure_group(self):
        # test catching bad structure groups
        self.assertRaises(ValueError, Benchmark.from_structure_group, "invalid_group")

        # test a single structure group
        bm = Benchmark.from_structure_group("elemental")
        self.assertEqual(len(bm.structures), 16)
        self.assertTrue("Ga_12174" in bm.structures)

        # test multiple structure groups
        bm = Benchmark.from_structure_group(["elemental", "common_binaries"])
        self.assertEqual(len(bm.structures), 27)
        self.assertTrue("Ga_12174" in bm.structures)
        self.assertTrue("ZnS_sphalerite_651455" in bm.structures)

    def test_initialization(self):
        # test basic initialization
        bm = Benchmark(self.structures)
        self.assertEqual(len(bm.structures), 1)
        self.assertTrue("test_structure" in bm.structures)

        # check other parameters initialized correctly
        self.assertEqual(bm.max_nsites, 2)
        self.assertTrue(bm.all_structures_have_oxi)

        # check site information correct
        self.assertEqual(len(bm.site_information), 1)
        site_info = bm.site_information["test_structure"]
        self.assertEqual(tuple(site_info["unique_idxs"]), (0, 4))
        self.assertEqual(tuple(site_info["all_idxs"]), (0, 1))
        self.assertEqual(tuple(site_info["cation_idxs"]), (0,))
        self.assertEqual(tuple(site_info["anion_idxs"]), (1,))
        self.assertEqual(tuple(site_info["all_degens"]), (4, 4))
        self.assertEqual(tuple(site_info["cation_degens"]), (4,))
        self.assertEqual(tuple(site_info["anion_degens"]), (4,))
        self.assertEqual(site_info["all_total"], 8)
        self.assertEqual(site_info["cation_total"], 4)
        self.assertEqual(site_info["anion_total"], 4)
        self.assertEqual(site_info["cations"], {"Na"})
        self.assertEqual(site_info["anions"], {"Cl"})

        # Now try testing without symmetry
        bm = Benchmark(self.structures, symprec=None)
        site_info = bm.site_information["test_structure"]

        self.assertEqual(bm.max_nsites, 8)
        self.assertEqual(len(bm.site_information), 1)
        self.assertEqual(tuple(site_info["unique_idxs"]), (0, 1, 2, 3, 4, 5, 6, 7))
        self.assertEqual(tuple(site_info["all_idxs"]), (0, 1, 2, 3, 4, 5, 6, 7))
        self.assertEqual(tuple(site_info["cation_idxs"]), (0, 1, 2, 3))
        self.assertEqual(tuple(site_info["anion_idxs"]), (4, 5, 6, 7))
        self.assertEqual(tuple(site_info["all_degens"]), (1, 1, 1, 1, 1, 1, 1, 1))
        self.assertEqual(tuple(site_info["cation_degens"]), (1, 1, 1, 1))
        self.assertEqual(tuple(site_info["anion_degens"]), (1, 1, 1, 1))
        self.assertEqual(site_info["all_total"], 8)
        self.assertEqual(site_info["cation_total"], 4)
        self.assertEqual(site_info["anion_total"], 4)
        self.assertEqual(site_info["cations"], {"Na"})
        self.assertEqual(site_info["anions"], {"Cl"})

        # test structures without site property
        structures = deepcopy(self.structures)
        structures["test_structure"].remove_site_property("coordination")
        self.assertRaises(AttributeError, Benchmark, structures)

        # test perturb structure initialization, i.e. check sites are not equal
        structures = deepcopy(self.structures)
        bm = Benchmark(structures, perturb_sigma=0.1)
        orig_sites = structures["test_structure"][0].frac_coords
        perturb_sites = bm.structures["test_structure"][0].frac_coords
        self.assertRaises(
            AssertionError,
            np.testing.assert_array_almost_equal,
            orig_sites,
            perturb_sites,
        )

    def test_benchmark(self):
        bm = Benchmark(self.structures)

        # test dataframe output
        results = bm.benchmark(self.nn_methods)
        expected_results = {
            "EconNN0": {"test_structure": {"Cl": 6}},
            "EconNN1": {"test_structure": {"Na": 6}},
            "MinimumVIRENN0": {"test_structure": {"Cl": 6}},
            "MinimumVIRENN1": {"test_structure": {"Na": 6}},
        }
        self.assertEqual(results.to_dict(), expected_results)

        # test dict output
        results = bm.benchmark(self.nn_methods, return_dataframe=False)
        vire_results = {"test_structure": [{"Cl": 6}, {"Na": 6}]}
        econ_results = {"test_structure": [{"Cl": 6}, {"Na": 6}]}
        self.assertEqual(results[self.nn_methods[0]], vire_results)
        self.assertEqual(results[self.nn_methods[1]], econ_results)

        # test no symmetry
        bm = Benchmark(self.structures, symprec=None)
        no_sym_results = bm.benchmark(self.nn_methods)
        expected_no_sym_results = {
            "EconNN0": {"test_structure": {"Cl": 6}},
            "EconNN1": {"test_structure": {"Cl": 6}},
            "EconNN2": {"test_structure": {"Cl": 6}},
            "EconNN3": {"test_structure": {"Cl": 6}},
            "EconNN4": {"test_structure": {"Na": 6}},
            "EconNN5": {"test_structure": {"Na": 6}},
            "EconNN6": {"test_structure": {"Na": 6}},
            "EconNN7": {"test_structure": {"Na": 6}},
            "MinimumVIRENN0": {"test_structure": {"Cl": 6}},
            "MinimumVIRENN1": {"test_structure": {"Cl": 6}},
            "MinimumVIRENN2": {"test_structure": {"Cl": 6}},
            "MinimumVIRENN3": {"test_structure": {"Cl": 6}},
            "MinimumVIRENN4": {"test_structure": {"Na": 6}},
            "MinimumVIRENN5": {"test_structure": {"Na": 6}},
            "MinimumVIRENN6": {"test_structure": {"Na": 6}},
            "MinimumVIRENN7": {"test_structure": {"Na": 6}},
        }

        self.assertEqual(no_sym_results.to_dict(), expected_no_sym_results)

    def test_score(self):
        # test all scores
        bm = Benchmark(self.structures)
        scores = bm.score(self.nn_methods)
        expected_scores = {
            "EconNN": {"test_structure": 4.0, "Total": 4.0},
            "MinimumVIRENN": {"test_structure": 4.0, "Total": 4.0},
        }
        self.assertEqual(scores.to_dict(), expected_scores)

        # test cation scores
        scores = bm.score(self.nn_methods, site_type="cation")
        expected_scores = {
            "EconNN": {"Total": 0.0, "test_structure": 0.0},
            "MinimumVIRENN": {"Total": 0.0, "test_structure": 0.0},
        }
        self.assertEqual(scores.to_dict(), expected_scores)

        # test anion scores
        scores = bm.score(self.nn_methods, site_type="anion")
        expected_scores = {
            "EconNN": {"Total": 8.0, "test_structure": 8.0},
            "MinimumVIRENN": {"Total": 8.0, "test_structure": 8.0},
        }
        self.assertEqual(scores.to_dict(), expected_scores)

        # test cation-anion filtering for all sites
        scores = bm.score(self.nn_methods, cation_anion=True)
        expected_scores = {
            "EconNN": {"test_structure": 1.0, "Total": 1.0},
            "MinimumVIRENN": {"test_structure": 1.0, "Total": 1.0},
        }
        self.assertEqual(scores.to_dict(), expected_scores)

        # test cation-anion filtering for cation sites
        scores = bm.score(self.nn_methods, cation_anion=True, site_type="cation")
        expected_scores = {
            "EconNN": {"test_structure": 0.0, "Total": 0.0},
            "MinimumVIRENN": {"test_structure": 0.0, "Total": 0.0},
        }
        self.assertEqual(scores.to_dict(), expected_scores)

        # test cation-anion filtering for anion sites
        scores = bm.score(self.nn_methods, cation_anion=True, site_type="anion")
        expected_scores = {
            "EconNN": {"test_structure": 2.0, "Total": 2.0},
            "MinimumVIRENN": {"test_structure": 2.0, "Total": 2.0},
        }
        self.assertEqual(scores.to_dict(), expected_scores)

        # test returning raw sites
        scores = bm.score(self.nn_methods, return_raw_site_scores=True)
        expected_scores = {
            "MinimumVIRENN": {"test_structure": [0, -8]},
            "EconNN": {"test_structure": [0, -8]},
        }
        self.assertEqual(scores.to_dict(), expected_scores)

    def test_multiple_same_methods(self):
        # test that if running the benchmark on multiple NN methods of the same class,
        # that each NN method is named differently and appears in the benchmark and
        # score results
        bm = Benchmark(self.structures)
        nn_methods = [VoronoiNN(), VoronoiNN(tol=0.5)]

        results = bm.benchmark(nn_methods)
        expected_results = {
            "VoronoiNN(0)0": {"test_structure": {"Cl": 6}},
            "VoronoiNN(0)1": {"test_structure": {"Na": 6}},
            "VoronoiNN(1)0": {"test_structure": {"Cl": 6}},
            "VoronoiNN(1)1": {"test_structure": {"Na": 6}},
        }
        self.assertEqual(results.to_dict(), expected_results)

        scores = bm.score(nn_methods)
        expected_scores = {
            "VoronoiNN(0)": {"Total": 4.0, "test_structure": 4.0},
            "VoronoiNN(1)": {"Total": 4.0, "test_structure": 4.0},
        }
        self.assertEqual(scores.to_dict(), expected_scores)
