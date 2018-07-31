from __future__ import division, unicode_literals
import unittest

from materialscoord.core import Benchmark, HumanInterpreter
from pymatgen.analysis.local_env import MinimumDistanceNN, VoronoiNN, CrystalNN

class CoreTest(unittest.TestCase):

    def test_benchmark(self):
        #TODO: test custom_set

        bm = Benchmark()

        rd = {1: 1.0005, 2: 1.0006}
        rd = bm._roundcns(rd, 3)
        self.assertEqual(rd[1], 1)
        self.assertEqual(rd[2], 1.001)

        cns = [('Ca', {'O': 8.0}), ('W', {'O': 4.0}), ('O', {'W': 1.0})]
        ions = ['Ca', 'W']
        exp_res = [('Ca', {u'O': 8.0}), ('W', {u'O': 4.0})]
        res = bm._popel(cns, ions)
        self.assertEqual(len(res), 2)
        self.assertTrue(res[0][0] in [er[0] for er in exp_res])
        self.assertTrue(res[1][0] in [er[0] for er in exp_res])
        for r in res:
            for er in exp_res:
                if r[0] == er[0]:
                    self.assertTrue(len(list(r[1].keys())), 1)
                    k = list(r[1].keys())[0]
                    ek = list(er[1].keys())[0]
                    self.assertTrue(k == ek)
                    self.assertTrue(r[1][k] == er[1][k])

        mdnn = MinimumDistanceNN()
        vnn = VoronoiNN()
        cnn = CrystalNN()
        df = bm.benchmark((mdnn, vnn, cnn))
        self.assertTrue(type(df.loc['As_alpha_16518']['CrystalNN0']), dict)
        self.assertEqual(list(df.loc['As_alpha_16518']['CrystalNN0'].keys())[0], 'As')
        self.assertEqual(df.loc['As_alpha_16518']['CrystalNN0']['As'], 3)
        self.assertEqual(df.loc['C_diamond_52054']['CrystalNN0']['C'], 4)
        self.assertEqual(df.loc['C_graphite_76767']['CrystalNN0']['C'], 3)
        self.assertEqual(df.loc['C_graphite_76767']['CrystalNN1']['C'], 3)
        self.assertEqual(df.loc['Cu_52256']['CrystalNN0']['Cu'], 12)
        self.assertEqual(df.loc['Ga_12174']['CrystalNN0']['Ga'], 12)
        self.assertEqual(df.loc['Hg_alpha_104296']['CrystalNN0']['Hg'], 6)
        self.assertEqual(df.loc['La_43573']['CrystalNN0']['La'], 12)
        self.assertEqual(df.loc['La_43573']['CrystalNN1']['La'], 12)
        self.assertEqual(df.loc['Mg_52260']['CrystalNN0']['Mg'], 12)
        self.assertEqual(df.loc['Mn_alpha_42743']['CrystalNN1']['Mn'], 16)
        self.assertEqual(df.loc['Mn_alpha_42743']['CrystalNN2']['Mn'], 13)
        self.assertEqual(df.loc['Mn_alpha_42743']['CrystalNN3']['Mn'], 12)
        self.assertEqual(df.loc['Mn_beta_41775']['CrystalNN1']['Mn'], 12)
        self.assertEqual(df.loc['P_black_23836']['CrystalNN0']['P'], 3)
        self.assertEqual(df.loc['Se_trigonal_23068']['CrystalNN0']['Se'], 2)
        self.assertEqual(df.loc['Sm_76031']['CrystalNN1']['Sm'], 12)
        self.assertEqual(df.loc['Sn_beta_106072']['CrystalNN0']['Sn'], 6)
        self.assertEqual(df.loc['U_alpha_16056']['CrystalNN0']['U'], 12)
        self.assertEqual(df.loc['W_alpha_43667']['CrystalNN0']['W'], 8)
        self.assertEqual(df.loc['As_alpha_16518']['MinimumDistanceNN0']['As'], 3)
        self.assertEqual(df.loc['Mn_alpha_42743']['MinimumDistanceNN1']['Mn'], 10)
        self.assertEqual(df.loc['Mn_alpha_42743']['MinimumDistanceNN2']['Mn'], 4)
        self.assertEqual(df.loc['Mn_alpha_42743']['MinimumDistanceNN3']['Mn'], 4)
        self.assertEqual(df.loc['Mn_beta_41775']['MinimumDistanceNN0']['Mn'], 6)
        self.assertEqual(df.loc['Mn_beta_41775']['MinimumDistanceNN1']['Mn'], 12)
        self.assertEqual(df.loc['P_black_23836']['MinimumDistanceNN0']['P'], 3)
        self.assertEqual(df.loc['Se_trigonal_23068']['MinimumDistanceNN0']['Se'], 2)
        self.assertEqual(df.loc['Sm_76031']['MinimumDistanceNN1']['Sm'], 12)
        self.assertEqual(df.loc['U_alpha_16056']['MinimumDistanceNN0']['U'], 4)
        self.assertEqual(df.loc['As_alpha_16518']['VoronoiNN0']['As'], 19)
        self.assertEqual(df.loc['C_diamond_52054']['VoronoiNN0']['C'], 16)
        self.assertEqual(df.loc['C_graphite_76767']['VoronoiNN0']['C'], 13)
        self.assertEqual(df.loc['C_graphite_76767']['VoronoiNN1']['C'], 17)
        self.assertEqual(df.loc['Hg_alpha_104296']['VoronoiNN0']['Hg'], 12)
        self.assertEqual(df.loc['Mg_52260']['VoronoiNN0']['Mg'], 14)
        self.assertEqual(df.loc['Mn_beta_41775']['VoronoiNN1']['Mn'], 14)
        self.assertEqual(df.loc['P_black_23836']['VoronoiNN0']['P'], 15)
        self.assertEqual(df.loc['Se_trigonal_23068']['VoronoiNN0']['Se'], 16)
        self.assertEqual(df.loc['Sn_beta_106072']['VoronoiNN0']['Sn'], 18)
        self.assertEqual(df.loc['W_alpha_43667']['VoronoiNN0']['W'], 14)

    def test_human_interpreter(self):
        hi = HumanInterpreter()
        self.assertEqual(hi.compute(None, 0), "null")

if __name__ == "__main__":
    unittest.main()
