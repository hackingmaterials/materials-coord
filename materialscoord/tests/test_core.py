from __future__ import division, unicode_literals
import unittest

from materialscoord.core import HumanInterpreter

class CoreTest(unittest.TestCase):

    def test_human_interpreter(self):
        hi = HumanInterpreter()
        self.assertEqual(hi.compute(None, 0), "null")

if __name__ == "__main__":
    unittest.main()
