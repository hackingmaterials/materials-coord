# MaterialsCoord Benchmark

MaterialsCoord provides an infrastructure for comparing the performance of
different crystal structure bonding algorithms against a benchmark set of
materials.

Atomic coordination numbers (CNs) are an important descriptor for the local
environment of condensed materials. A range of algorithms for calculating
coordination numbers have been proposed, the performance of which are often
unknown outside the chemical domain they were designed for. MaterialsCoord
allows for the benchmarking of different CN algorithms against a set of crystal
structures for which the bonding environments have been experimentally
determined.

## Usage

The basic API of MaterialsCoord can be summarized as:

- `materialscoord.Benchmark` provides the primary interface.
  It takes a set of structures that have been decorated with the
  `"coordination"` site property (see the [Preparing crystal structures for benchmarking](#preparing-crystal-structures-for-benchmarking)
  section).
- The `Benchmark.from_structure_group()` method provides a convenience function
  to initialize the benchmark from a list of pre-prepared structure groups.
  More information on the available groups is provided in the
  [Benchmarking structure sets](#benchmarking-structure-sets)
  section.
- `Benchmark.benchmark()` accepts a list of coordination algorithms and runs them
  on the set of test structures.
- `Benchmark.score()` accepts a list of coordination algorithms and compares
  their performance against a human interpretation of crystal structure bonding.

A full introduction to using MaterialsCoord is provided in the
[introduction notebook](https://github.com/hackingmaterials/MaterialsCoord/tree/master/examples/introduction-to-MaterialsCoord.ipynb).

### Coordination number algorithms

MaterialsCoord has been designed to interface with the `NearNeighbors` methods
implemented in [`pymatgen.analysis.local_env`](https://pymatgen.org/pymatgen.analysis.local_env.html#module-pymatgen.analysis.local_env).
A number of methods exist, including:

- [`CrystalNN`](https://pymatgen.org/pymatgen.analysis.local_env.html#pymatgen.analysis.local_env.CrystalNN)
- [`VoronoiNN`](https://pymatgen.org/pymatgen.analysis.local_env.html#pymatgen.analysis.local_env.VoronoiNN)
- [`BrunnerNN_reciprocal`](https://pymatgen.org/pymatgen.analysis.local_env.html#pymatgen.analysis.local_env.BrunnerNN_reciprocal)
- [`MinimumDistanceNN`](https://pymatgen.org/pymatgen.analysis.local_env.html#pymatgen.analysis.local_env.MinimumDistanceNN)

### Benchmarking structure sets

The benchmark includes around 95 pre-prepared test structures for which the
correct coordination has been determined by a human interpreter (taken from the
scientific literature). The structures have been grouped into a number of
material classes, including:

- `"elemental"`: Simple elemental materials, including diamond, graphite, Ga,
  and α-As.
- `"common_binaries"`: Simple and more complex binary structures, including
  rocksalt NaCl, rutile TiO<sub>2</sub>, and γ-brass.
- `"ABX3"`: ABX<sub>3</sub> structured ternary materials, including perovskite
  SrTiO<sub>3</sub> and argonite CaCO<sub>3</sub>.
- `"ABX4"`: ABX<sub>4</sub> structured ternary materials, including zircon,
  (ZrSiO<sub>4</sub>) and wolframite (FeWO<sub>4</sub>).
- `"A2BX4"`: A<sub>2</sub>BX<sub>4</sub> structured ternary materials, including
  olivine Fe<sub>2</sub>SiO<sub>4</sub>.

The full set of materials classes and crystal structures can be found
in the [structures directory](https://github.com/hackingmaterials/MaterialsCoord/tree/master/materialscoord/structures).

MaterialsCoord supports benchmarking on custom structures. See the
[Preparing crystal structures for benchmarking](#preparing-crystal-structures-for-benchmarking)
section for more details.

### Running the Benchmark

The MaterialsCoord benchmark should be run using the python API. For example,
we can compare the Effective Coordination Number (ECoN) and O'Keeffe's Voronoi
CN method on the `"elemental"` and `"common_binaries"` structure groups as follows:

```python
from pymatgen.analysis.local_env import EconNN, VoronoiNN
from materialscoord.core import Benchmark

nn_methods = [EconNN(), VoronoiNN()]

bm = Benchmark.from_structure_group(["elemental", "common_binaries"])
bm.score(nn_methods)
```

The `score` function will return the results as a [Pandas](https://pandas.pydata.org)
`DataFrame` object. Further details can be found in the
[introduction notebook](https://github.com/hackingmaterials/MaterialsCoord/tree/master/examples).

### Preparing crystal structures for benchmarking

MaterialsCoord uses pymatgen `Structure` objects to internally store the crystal
structure data. The structures provided in the default structure groups contain
an additional site property called `"coordination"`. This property contains
information on the human interpreted coordination environment for each site.

For example, if a site is bonded to two oxygen and three chlorine atoms, the
coordination property should be a dictionary of ``{"O": 2, "Cl": 3}``. If
multiple coordination numbers are deemed correct these should be provided as a
list. For example, if the coordination to gallium could be 4 or 12 coordinate,
the coordination property should a dictionary of ``{"Ga": [4, 12]}``.

A full guide to preparing crystal structures for use in MaterialsCoord is in
provided in the
[benchmarking using custom structures notebook](https://github.com/hackingmaterials/MaterialsCoord/tree/master/examples/benchmarking-custom-structures.ipynb).

## How to cite MaterialsCoord

*A research paper has been submitted. This section will be updated when the
paper has been published online.*

## Installation

MaterialsCoord can be installed from source using:

```bash
git clone https://github.com/hackingmaterials/MaterialsCoord.git
cd MaterialsCoord
pip install .
```

MaterialsCoord requires Python 3.6+.

## What’s new?

Track changes to MaterialsCoord through the
[Changelog](https://github.com/hackingmaterials/MaterialsCoord/blob/master/CHANGELOG.rst).

## Contributing

MaterialsCoord is in early development but we still welcome your
contributions. Please read our [contribution guidelines](https://github.com/hackingmaterials/MaterialsCoord/blob/master/CONTRIBUTING.rst)
for more information. We maintain a list of all
contributors [here](https://github.com/hackingmaterials/MaterialsCoord/blob/master/CONTRIBUTORS.rst).

## License

MaterialsCoord is released under a modified BSD license;
the full text can be found
[here](https://github.com/hackingmaterials/MaterialsCoord/blob/master/LICENSE).
