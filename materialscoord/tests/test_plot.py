import unittest
from pandas import DataFrame

from materialscoord.plot import plot_benchmark_scores


class PlotTest(unittest.TestCase):
    """Test plotting functions."""

    def test_plot(self):
        """Simple test to check the plot function doesn't error."""
        data = {
            "EconNN": {"test_structure": 2.0, "Total": 2.0},
            "MinimumVIRENN": {"test_structure": 2.0, "Total": 2.0},
        }
        scores = DataFrame(data=data)
        plt = plot_benchmark_scores(scores)
        self.assertNotEqual(plt, None)
