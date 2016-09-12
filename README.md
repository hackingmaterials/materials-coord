# MaterialsCoord

MaterialsCoord provides an infrastructure for comparing coordination numbers produced by different
approaches against a benchmark set of material structures.

Atomic coordination numbers (CNs) are one of the most important descriptors for the local environments in condensed materials, but they
have no universal definition. Therefore, researchers have come up with a wide range of algorithms, performance of which are often
unknown outside the chemical domain they were particularly designed for. MaterialsCoord aims at providing the necessary infrastructure
to host, interface, compare and benchmark different CN algorithms against a selected set of crystal structures. It further allows
human interpretations of CN environments to be incorporated into benchmarking (in progress).

[1. How do I benchmark Coordination Number algorithms with MaterialsCoord?](#how_do_i_benchmark_cn_algos)  
[2. Which coordination number algorithms are currently available?](#cn_algos)    
[3. How do I implement a new Coordination Number algorithm?](#how_do_i_implement_cn_algos)    
[4. How can I use MaterialsCoord on my own structures?](#how_can_i_use_my_own_structures)  
[5. What is the HumanInterpreter?](#humaninterpreter)  
[6. Installation and requirements](#install)  

<a name="how_do_i_benchmark_cn_algos"/>
## How do I benchmark Coordination Number algorithms with MaterialsCoord?

Simply put:
* `Benchmark` class provides the necessary infrastructure to perform the comparison of CN algorithms. When being initialized, 
it takes the group of test structures (or your own structures) and a list of `CNBase` methods (CN calculation methods) as arguments.
* `Benchmark.benchmark()` performs the CN calculations on the selected test structures
* `Benchmark.report()` provides different types of reports that summarizes the results of the benchmarking.

For example, we can compare the Effective Coordination Number (ECoN) and O'Keeffe's Voronoi CN method using a set of unique elemental
crystal structures:

```python
from materialscoord.cn_methods import TestECoN, TestVoronoiCoordFinder
from materialscoord.core import Benchmark
bm = Benchmark([TestECoN(), TestVoronoiCoordFinder()], "elemental")
bm.benchmark()
bm.report()
```

Further details can be found in [examples](https://github.com/aykol/MaterialsCoord/tree/master/examples) provided in [benchmark_examples.ipynb](https://github.com/aykol/MaterialsCoord/blob/master/examples/benchmark_examples.ipynb).

<a name="cn_algos"/>
## Which coordination number algorithms are currently available?
Currently, MaterialsCoord has the following coordination number algorithms implemented in `materialscoord.cn_methods`:
- `TestVoronoiCoordFinder`: Weighted Voronoi CNs (O'Keeffe's method).
- `TestECoN`: Effective Coordination Numbers, ECoN, (Hoppe's method.
- `TestVoronoiCoordFinder_mod`: A modified version of `TestVoronoiCoordFinder`.
- `TestVoronoiLegacy`: Basic Voronoi facet counting.
- `TestBrunnerReciprocal`: Brunner's method of largest reciprocal gap in interactomic distances.
- `TestBrunnerRelative`: Brunner's method of largest relative gap in interactomic distances.
- `TestBrunnerReal`: Brunner's method of largest gap in interactomic distances.
- `TestDelaunay`: David Mrdjenovich et al.'s Delaunay triangulation based algorithm (under development at LBNL).

You can see the details of the available algorithms [here](https://github.com/aykol/MaterialsCoord/blob/master/materialscoord/cn_methods.py). New algorithms are welcome; simply submit a pull request on github.

<a name="how_do_i_implement_cn_algos"/>
## How do I implement a new Coordination Number algorithm?

This is fairly simple:

1. Define a new class that is subclassed from `materialscoord.core.CNBase`.
2. The class must have a method named `compute` which takes a pymatgen Structure and site-index as input,
and returns a dictionary of CNs for that site; e.g. `{'O': 4.4, 'F': 2.1}`.

For example:
```python
class MyCoordinationNumberAlgorithm(CNBase):
    """
    My new algorithm
    """
    def compute(self, structure, n):
        params = self._params

        # ... here your algorithm finds CNs of site n.
        # e.g. cns = my_algorithm(structure, n, **params)

        return cns
```

Any parameters the algorithm needs can be passed to the class when initializing using the params keyword as a dictionary. These can later
be accessed from `compute` method as the `_params` attribute of the class.

In method `compute` you can do whatever is necessary to interface the algorithm with MaterialsCoord. Options include:
* Simplest: Implementing the entire algorithm within the `compute` method.
* Recommended: Add the necessary "bulky" code to a relevant module in external_src package and import as you define
  `compute`.
* Not recommended: call or import a library/program outside of MaterialsCoord within the `compute` method.
This is not recommended as external dependencies will restrict portability,
but maybe unavoidable if the external algorithm is part of another python package,
or is written in some other language such as Java. In that case `compute` can simply serve as a wrapper that calls
and post processes the output of the external code.

<a name="how_can_i_use_my_own_structures"/>
## How can I use MaterialsCoord on my own structures?
There are different ways of how one can do this.
* Add a new "group" folder that includes your structures into the test_structures folder. Then you can initialize Benchmark with
your `structure_groups = group_name` (i.e. the name of your folder). You can then call your new structure group any time you want.
* If you don't want to permanently add the structures to MaterialsCoord but rather run CN methods on some external structures,
you can use the `custom_set` argument of `Benchmark` to provide the path to a set of structure files. If `custom_set` is given,
`Benchmark` will ignore any `structure_group` provided. 

An example can be found in the [custom_tests](https://github.com/aykol/MaterialsCoord/blob/master/examples/custom_tests.ipynb) notebook.

Note that in any case, the structures provided can be of any type that pymatgen can automatically interpret (cif, POSCAR, etc.)

<a name="humaninterpreter"/>
## What is the HumanInterpreter?
`HumanInterpreter` is a special CN method that provides CNs from the "human" interpreted coordination environments stored in the
human_interpreter.yaml file. It can be added as a CN method along with other CN methods. Currently only the `common_binaries` structure group
has human interpreted CNs, but more will be added soon.

If you interpret the CN environment in a structure and you want to add it to MaterialsCoord to compare against availble CN methods, you can basically add the CN numbers to the human_interpreter.yaml as a dictionary that matches
the name of the structure file in `test_structures` or `custom_set` you provided (without the file extension, if it has any). For example:
```yaml
Fe3O4_spinel:
    - Fe:
        Fe: 0.0
        O: 6.0
    - Fe:
        Fe: 0.0
        O: 4.0
    - O:
        Fe: 4.0
        O: 0.0
```
describes the CN environment in the spinel Fe3O4 spinel structure, where there are two different unique Fe sites, and one O site. Sites are given as a list. And for each site a dictionary
of neighboring elements are provided. In this sub-dictionary of surroundings of a given site, we do not differentiate between sites and only list the total CN number for each chemical element (i.e. there is one Fe in each) opposite to the unique sites list (where we had two Fe sites).

Since it is not always easy to interpret the exact coordination number, MaterialsCoord will accept a range. Let's assume hypothetically we aren't sure how many Fe neighbors the O site has,
but we guess it's between 3.0 and 5.0, we can write:
```yaml
Fe3O4_spinel:
    - Fe:
        Fe: 0.0
        O: 6.0
    - Fe:
        Fe: 0.0
        O: 4.0
    - O:
        Fe:
            - 3.0
            - 5.0
        O: 0.0
```



<a name="install"/>
## Installation and requirements
For now, cloning the repo and executing `pyhton setup.py develop` should work fine. MaterialsCoord requires a number of other pyhton packages installed, such as pymatgen. Successfully
installing [pymatgen](http://pymatgen.org) with all of its dependencies, should also satisfy the dependencies of MaterialsCoord.
