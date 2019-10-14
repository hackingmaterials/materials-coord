import re
import os
import glob
import warnings

from pathlib import Path
from typing import List, Optional, Dict, Union

import numpy as np
import pandas as pd

from collections import defaultdict

from pkg_resources import resource_filename

from pymatgen import Specie
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.local_env import NearNeighbors

from materialscoord.einstein_crystal_perturbation import perturb_einstein_crystal

module_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)))


class Benchmark(object):
    """
    Class for performing coordination number benchmarks on a set of structures
    using different nearest neighbor methods.

    Args:
        structures: A set of structures. Should be given as a dictionary of
            ``{"name": Structure}``. The structures should be decorated
            with a "coordination" property on each site specifying the correct
            coordination environment. For example, if a site is bonded to
            two oxygen and three chlorine atoms, the coordination property
            should be a dictionary of ``{"O": 2, "Cl": 3}``. If multiple
            coordination numbers are deemed correct these should be provided
            as a list. For example, if the coordination to gallium could be
            4 or 12 coordinate, the coordination property should a dictionary
            of ``{"Ga": [4, 12]}``.
        symprec: If not None, the benchmark will use symmetry to reduce
            the number of sites for which the coordinate number is calculated.
            symprec is the symmetry precision in Angstrom.
        use_weights: Whether or not to use the coordination number
            method weighting scheme. Defaults to False.
        perturb_sigma: If not None, this will enable the
            Einstein crystal test rig mode. Each site will be displaced
            according a normal distribution with the width equal to
            perturb_sigma. The perturbation complies thus with the expectation
            for an Einstein crystal, in which the potential is given by
            V(dr) = 1/2 * kspring * (dr)^2. kspring denotes the spring constant
            with which the sites are tethered to their equilibrium position, and
            dr is the distance of the site under consideration from its
            equilibrium position.

    TODO:
        - Use predict_oxidation_states() for structures without oxidation
            states present? Otherwise, score doesn't work with MP cif files.
    """

    def __init__(
        self,
        structures: Dict[str, Structure],
        symprec: float = 0.01,
        use_weights: bool = False,
        perturb_sigma: Optional[float] = None,
    ):
        self.structures = structures
        self.symprec = symprec
        self.use_weights = use_weights

        # use this to cache benchmark results
        self._benchmark: Dict[NearNeighbors, Dict[str, List]] = defaultdict(dict)

        for name, structure in structures.items():
            if "coordination" not in structure.site_properties:
                raise AttributeError(
                    "{} structure does not have a 'coordination' site property"
                )

        if perturb_sigma:
            for name, structure in self.structures.items():
                self.structures[name] = perturb_einstein_crystal(structure, perturb_sigma)

        # precompute the symmetrized structures to save time during the
        # benchmark. Also, determine the total number of unique
        # cations/anions each structure.
        self.site_information = {}
        self.max_nsites = 0

        n_structures_with_oxi = 0
        for name, structure in self.structures.items():
            if self.symprec:
                sga = SpacegroupAnalyzer(structure, symprec=self.symprec)

                equiv_sites = sga.get_symmetrized_structure().equivalent_sites
                unique_site_idxs = np.unique(
                    sga.get_symmetry_dataset()["equivalent_atoms"]
                )
            else:
                equiv_sites = [[s] for s in structure]
                unique_site_idxs = list(range(len(structure)))

            # cation_idxs and anion_idxs are the indexes of the cations/anions
            # in the list of unique sites (not the list of ALL sites).
            cation_degens = []
            cation_idxs = []
            anion_degens = []
            anion_idxs = []
            cations = set()
            anions = set()
            for i, sites in enumerate(equiv_sites):
                # Specie check needed to see if site has an oxidation state at all.
                if isinstance(sites[0].specie, Specie) and sites[0].specie.oxi_state >= 0:
                    # based on previous implementation, neutral ions will be scored as
                    # cations, however, neutral to neutral bonding is allowed, hence
                    # don't add the element to the cations set as this will be used to
                    # prevent cation-cation bonding later on.
                    cation_degens.append(len(sites))
                    cation_idxs.append(i)
                    if sites[0].specie.oxi_state > 0:
                        cations.add(sites[0].specie.element.name)
                elif isinstance(sites[0].specie, Specie) and sites[0].specie.oxi_state < 0:
                    anion_degens.append(len(sites))
                    anion_idxs.append(i)
                    anions.add(sites[0].specie.element.name)

            all_degens = [len(x) for x in equiv_sites]
            total_all = sum(all_degens)

            total_cations = sum(cation_degens)
            total_anions = sum(anion_degens)

            self.site_information[name] = {
                "unique_idxs": unique_site_idxs,
                "all_idxs": list(range(len(unique_site_idxs))),
                "cation_idxs": cation_idxs,
                "anion_idxs": anion_idxs,
                "all_degens": all_degens,
                "cation_degens": cation_degens,
                "anion_degens": anion_degens,
                "all_total": total_all,
                "cation_total": total_cations,
                "anion_total": total_anions,
                "cations": cations,
                "anions": anions,
            }

            self.max_nsites = max(self.max_nsites, len(unique_site_idxs))

            if cation_idxs or anion_idxs:
                n_structures_with_oxi += 1

        # useful to know if all structures have oxidation states when calculating scores
        self.all_structures_have_oxi = n_structures_with_oxi == len(structures)

        # once we have calculated anion and cation type, remove oxidation
        # states from the structures. Not sure if this is a good idea but I'm
        # just copying from the previous implementation
        for structure in self.structures.values():
            structure.remove_oxidation_states()

    @classmethod
    def from_structure_group(cls, structure_groups: Union[str, List[str]], **kwargs):
        """
        Initialises the benchmark from a list of test structure classes.

        Args:

            structure_groups: One or more test structure groups. Options
                include: "elemental", "common_binaries", "ABX3", "ABX4",
                "A2BX4", "laves". See the "test_structures" folder for
                the full list of options. Defaults to "elemental".
            **kwargs: Additional keyword arguments that will be passed
                to the Benchmark constructor.
        """
        if isinstance(structure_groups, str):
            structure_groups = [structure_groups]

        str_path = resource_filename("materialscoord", "structures")

        filenames = []
        for materials_class in structure_groups:
            path = os.path.join(str_path, materials_class, "*.json")
            filenames.extend(glob.glob(path))

        structures = {}
        for filename in filenames:
            name = Path(filename).stem
            structures[name] = Structure.from_file(filename)

        return cls(structures, **kwargs)

    def benchmark(self, methods: List[NearNeighbors], return_dataframe: bool = True):
        """
        Calculates the coordination numbers for all sites in all structures
        using each nearest neighbor method.

        Args:
            methods: A list of NearNeighbors methods. E.g., from
                ``pymatgen.analysis.local_env``.
            return_dataframe: Whether to return the results as a pandas
                dataframe.

        Returns:
            If ``return_dataframe``. The benchmark results as a pandas
            DataFrame, else a dictionary formatted as::

                {method: {structure_name: List[cn_dicts]}}

            Where the cn_dicts are given for each unique site in the structure.
            See the docstring for `NearNeighbors.get_cn_dict` for the format of
            cn_dict.
        """
        for method in methods:
            for name in self.structures:
                if method not in self._benchmark or name not in self._benchmark[method]:
                    self._benchmark[method][name] = self._benchmark_structure(
                        name, method)

        if not return_dataframe:
            return self._benchmark

        df_data = defaultdict(dict)

        for method in methods:
            for name in self.structures:
                for site_idx in range(self.max_nsites):
                    column = method.__class__.__name__ + str(site_idx)

                    if site_idx < len(self.site_information[name]["unique_idxs"]):
                        val = self._benchmark[method][name][site_idx]
                    else:
                        val = None

                    df_data[column][name] = val

        return pd.DataFrame(data=df_data)

    def score(self, methods: List[NearNeighbors], site_type: str = "all") -> pd.DataFrame:
        r"""
        Assigns a score for each near neighbor method for each structure.

        The score is calculated by taking the summation of the absolute value
        of the error in CN prediction (CN_calc - CN_expected) multiplied by the
        degeneracy over the total number of sites. I.e.,

        score = \sum_i^N_sites^unique | CN_i^calc - CN_i^expected | * N_i^degen)
                / N_sites

        Args:
            methods: A list of NearNeighbors methods. E.g., from
                ``pymatgen.analysis.local_env``.
            site_type: Which sites to include when calculating the score.
                Options are "all", "cation", "anion". For example, if "cation"
                is chosen, the scores will only reflect the coordination numbers
                of the cation sites. In addition, only cation to anion
                bonds will be considered and any cation to cation bonds will
                be ignored. Note that cation and anion can only be used on
                input structures that have oxidation states.

        Returns:
            The scores as a Pandas DataFrame.
        """
        if site_type != "all" and not self.all_structures_have_oxi:
            warnings.warn(
                "Not all structures have oxidation states, "
                "{} based scoring is not recommended".format(site_type)
            )

        results = self.benchmark(methods, return_dataframe=False)

        scores: Dict[str, Dict[str, float]] = defaultdict(dict)
        for method in methods:
            for name in self.structures:
                scores[method.__class__.__name__][name] = self._score_structure(
                    name, results[method][name], site_type=site_type
                )

        df = pd.DataFrame(data=scores)
        df.loc["Total"] = df.sum(axis=0)

        return df

    def _benchmark_structure(self, name: str, nn: NearNeighbors):
        """
        Run a near neighbor algorithm on a structure.

        Args:
            name: A structure name.
            nn: A near neighbor method.

        Returns:
            A list of dictionaries for each unique site in the structure. See
            the docstring for `NearNeighbors.get_cn_dict` for the format of the
            dictionaries.
        """
        results = []
        for i in self.site_information[name]["unique_idxs"]:
            cn_dict = nn.get_cn_dict(self.structures[name], i, self.use_weights)

            if nn.__class__.__name__ == "MinimumVIRENN":
                tmp_cn_dict = {}
                for k, v in cn_dict.items():
                    k = re.sub("[^a-zA-Z]+", "", k)
                    tmp_cn_dict[k] = v
                cn_dict = tmp_cn_dict

            results.append(cn_dict)

        return results

    def _score_structure(
        self, name: str, predictions: List[Dict[str, float]], site_type: str = "all"
    ) -> float:
        r"""
        Calculate the score based on a set of predictions.

        The score is calculated by taking the summation of the absolute value
        of the error in CN prediction (CN_calc - CN_expected) multiplied by the
        degeneracy over the total number of sites. I.e.,

        score = \sum_i^N_sites^unique | CN_i^calc - CN_i^expected | * N_i^degen)
                / N_sites

        Args:
            name: A structure name.
            predictions: The coordination predictions for each unique site in
                the structure. Formatted as a list of coordination dictionaries.
                See the docstring for `NearNeighbors.get_cn_dict` for the format
                of the dictionaries.
            site_type: Which sites to include when calculating the score.
                Options are "all", "cation", "anion". For example, if "cation"
                is chosen, the scores will only reflect the coordination numbers
                of the cation sites. In addition, only cation to anion
                bonds will be considered and any cation to cation bonds will
                be ignored. Note that cation and anion can only be used on
                input structures that have oxidation states.

        Returns:
            The score for the structure as a float.
        """
        structure = self.structures[name]
        idxs = self.site_information[name]["{}_idxs".format(site_type)]
        degens = self.site_information[name]["{}_degens".format(site_type)]
        total = self.site_information[name]["{}_total".format(site_type)]

        if site_type != "all":
            # only consider bonds to oppositely charged ions
            # TODO: This seems wrong to me as it could make an algorithm appear
            #  to be doing better than it actually is, but I'm this behaviour
            #  to be consistent with the previous implementation -- AG
            filter_elements = self.site_information[name]["{}s".format(site_type)]
        else:
            filter_elements = set()

        # as site_idx below is the index of the site in the unique sites NOT the
        # index in the overall structure, we first get actual coordinations
        # of the unique sites.
        coordinations = [
            structure[i].properties["coordination"]
            for i in self.site_information[name]["unique_idxs"]
        ]

        score = 0
        for site_idx, site_degen in zip(idxs, degens):
            prediction = predictions[site_idx]
            coordination = coordinations[site_idx]

            # create a list of possible bonding elements (these are the species
            # in both the known coordination dict and the predicted coordination
            # dict)
            elements = set(list(coordination.keys()) + list(prediction.keys()))

            # exclude species that are of the opposite charge. See note above.
            elements = elements.difference(filter_elements)

            site_score = 0
            for element in elements:
                pred_val = prediction.get(element, 0)
                actual_val = coordination.get(element, 0)

                if isinstance(actual_val, list):
                    # there are multiple options for coordination, choose the
                    # value that results in the smallest score
                    site_score += min([abs(pred_val - x) for x in actual_val])
                else:
                    site_score += abs(pred_val - actual_val)

            score += site_score * site_degen

        if total == 0:
            return np.nan

        return score / total
