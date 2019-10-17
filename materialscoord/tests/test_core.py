import unittest

from materialscoord.core import Benchmark
from pymatgen import Structure
from pymatgen.analysis.local_env import EconNN, VoronoiNN


class BenchmarkTest(unittest.TestCase):

    def setUp(self):
        # set up a test structure
        structure = Structure.from_spacegroup(
            225,
            [[5.7, 0, 0], [0, 5.7, 0], [0, 0, 5.7]],
            ["Na1+", "Cl1-"],
            [[0, 0, 0], [0.5, 0, 0]]
        )
        all_sites_coordination = [{"Cl": 6}] * 4 + [{"Na": 6}] * 4
        structure.add_site_property("coordination", all_sites_coordination)
        self.structures = {"test_structure": structure}
        self.nn_methods = [EconNN(), VoronoiNN()]

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
        self.assertEqual(tuple(site_info["cation_idxs"]), (0, ))
        self.assertEqual(tuple(site_info["anion_idxs"]), (1, ))
        self.assertEqual(tuple(site_info["all_degens"]), (4, 4))
        self.assertEqual(tuple(site_info["cation_degens"]), (4, ))
        self.assertEqual(tuple(site_info["anion_degens"]), (4, ))
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

    def test_benchmark(self):
        bm = Benchmark(self.structures)

        # test dataframe output
        results = bm.benchmark(self.nn_methods)
        expected_results = {'EconNN0': {'test_structure': {'Na': 12, 'Cl': 6}},
                            'EconNN1': {'test_structure': {'Na': 6, 'Cl': 12}},
                            'VoronoiNN0': {'test_structure': {'Cl': 6}},
                            'VoronoiNN1': {'test_structure': {'Na': 6}}}
        self.assertEqual(results.to_dict(), expected_results)

        # test dict output
        results = bm.benchmark(self.nn_methods, return_dataframe=False)
        econ_results = {'test_structure': [{'Na': 12, 'Cl': 6}, {'Na': 6, 'Cl': 12}]}
        voronoi_results = {'test_structure': [{'Cl': 6}, {'Na': 6}]}
        self.assertEqual(results[self.nn_methods[0]], econ_results)
        self.assertEqual(results[self.nn_methods[1]], voronoi_results)

        # test no symmetry
        bm = Benchmark(self.structures, symprec=None)
        no_sym_results = bm.benchmark(self.nn_methods)
        expected_no_sym_results = {
            'EconNN0': {'test_structure': {'Na': 12, 'Cl': 6}},
            'EconNN1': {'test_structure': {'Cl': 6, 'Na': 12}},
            'EconNN2': {'test_structure': {'Cl': 6, 'Na': 12}},
            'EconNN3': {'test_structure': {'Cl': 6, 'Na': 12}},
            'EconNN4': {'test_structure': {'Na': 6, 'Cl': 12}},
            'EconNN5': {'test_structure': {'Na': 6, 'Cl': 12}},
            'EconNN6': {'test_structure': {'Na': 6, 'Cl': 12}},
            'EconNN7': {'test_structure': {'Na': 6, 'Cl': 12}},
            'VoronoiNN0': {'test_structure': {'Cl': 6}},
            'VoronoiNN1': {'test_structure': {'Cl': 6}},
            'VoronoiNN2': {'test_structure': {'Cl': 6}},
            'VoronoiNN3': {'test_structure': {'Cl': 6}},
            'VoronoiNN4': {'test_structure': {'Na': 6}},
            'VoronoiNN5': {'test_structure': {'Na': 6}},
            'VoronoiNN6': {'test_structure': {'Na': 6}},
            'VoronoiNN7': {'test_structure': {'Na': 6}}}

        self.assertEqual(no_sym_results.to_dict(), expected_no_sym_results)

    def test_score(self):
        # test all scores
        bm = Benchmark(self.structures)
        scores = bm.score(self.nn_methods)
        expected_scores = {'EconNN': {'test_structure': 12.0, 'Total': 12.0},
                           'VoronoiNN': {'test_structure': 0.0, 'Total': 0.0}}
        self.assertEqual(scores.to_dict(), expected_scores)

        # test cation scores
        bm = Benchmark(self.structures)
        scores = bm.score(self.nn_methods, site_type="cation")
        self.assertEqual(scores.to_dict(), expected_scores)

        # test anion scores
        bm = Benchmark(self.structures)
        scores = bm.score(self.nn_methods, site_type="anion")
        self.assertEqual(scores.to_dict(), expected_scores)

        # test cation-anion filtering
        bm = Benchmark(self.structures)
        scores = bm.score(self.nn_methods, cation_anion=True)
        expected_scores = {'EconNN': {'test_structure': 0.0, 'Total': 0.0},
                           'VoronoiNN': {'test_structure': 0.0, 'Total': 0.0}}
        self.assertEqual(scores.to_dict(), expected_scores)

    def test_multiple_same_methods(self):
        # test that if running the benchmark on multiple NN methods of the same class,
        # that each NN method is named differently and appears in the benchmark and
        # score results
        bm = Benchmark(self.structures)
        nn_methods = [VoronoiNN(), VoronoiNN(tol=0.5)]

        results = bm.benchmark(nn_methods)
        expected_results = {'VoronoiNN(0)0': {'test_structure': {'Cl': 6}},
                            'VoronoiNN(0)1': {'test_structure': {'Na': 6}},
                            'VoronoiNN(1)0': {'test_structure': {'Cl': 6}},
                            'VoronoiNN(1)1': {'test_structure': {'Na': 6}}}
        self.assertEqual(results.to_dict(), expected_results)

        scores = bm.score(nn_methods)
        expected_scores = {'VoronoiNN(0)': {'test_structure': 0.0, 'Total': 0.0},
                           'VoronoiNN(1)': {'test_structure': 0.0, 'Total': 0.0}}
        self.assertEqual(scores.to_dict(), expected_scores)
