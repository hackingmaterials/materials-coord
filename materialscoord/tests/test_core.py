from __future__ import division, unicode_literals
import unittest

from materialscoord.core import Benchmark, HumanInterpreter
from pymatgen.analysis.local_env import MinimumDistanceNN, VoronoiNN, CrystalNN

class CoreTest(unittest.TestCase):

    def test_benchmark(self):
        #TODO: test custom_set

        mdnn = MinimumDistanceNN()
        vnn = VoronoiNN()
        cnn = CrystalNN()
        bm = Benchmark((mdnn, vnn, cnn))

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

    def test_human_interpreter(self):
        hi = HumanInterpreter()
        self.assertEqual(hi.compute(None, 0), "null")

if __name__ == "__main__":
    unittest.main()
