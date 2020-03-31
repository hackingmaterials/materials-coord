"""This module contains the primary benchmarking class in MaterialsCoord."""

import re
import warnings
from copy import deepcopy

from pathlib import Path
from typing import List, Optional, Dict, Union, Any

import numpy as np
import pandas as pd

from collections import defaultdict, Counter

from pkg_resources import resource_filename

from pymatgen import Specie
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.local_env import NearNeighbors

from materialscoord.einstein_crystal_perturbation import perturb_einstein_crystal

_resource_dir = resource_filename("materialscoord", "structures")
_vire_re = re.compile("[^a-zA-Z]+")
_el_re = re.compile(r"[\d+-.]*")

CN_dict = Dict[str, float]  # define coordination dictionary type


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
        perturb_sigma: If not None, this will enable the
            Einstein crystal test rig mode. Each site will be displaced
            according a normal distribution with the width equal to
            perturb_sigma. The perturbation complies thus with the expectation
            for an Einstein crystal, in which the potential is given by
            V(dr) = 1/2 * kspring * (dr)^2. kspring denotes the spring constant
            with which the sites are tethered to their equilibrium position, and
            dr is the distance of the site under consideration from its
            equilibrium position.
        remove_oxidation_states: Remove oxidation states from structures before
            performing the benchmark. As oxidation states will not always be known,
            any near neighbor methods that require oxidation states (i.e.,
            MinimumVIRENN) should include a method to assign oxidation states
            automatically. Setting this option to False allows testing whether the
            automatic assignment of oxidation states will itself impact performance.
            Furthermore, some methods, such as CrystalNN, have slightly different
            behaviour if oxidation states are present. Again this option enables
            testing this behaviour.
        reciprocal_coordination: Several near neighbor methods are not reciprocal. I.e.,
            if site A is bonded to site B, it is not guaranteed that site B is bonded
            to site A. Enabling this option ensures that coordination is reciprocal by
            evaluating the coordination of all sites and including all bonds. This
            behaviour is the same as that provided by NearNeighbor.get_bonded_structure.
    """

    all_structure_groups: List[str] = [
        filename.stem for filename in Path(_resource_dir).iterdir() if filename.is_dir()
    ]

    def __init__(
        self,
        structures: Dict[str, Structure],
        symprec: Optional[float] = 0.01,
        perturb_sigma: Optional[float] = None,
        remove_oxidation_states: bool = True,
        reciprocal_coordination: bool = True,
    ):
        # make a deep copy to avoid modifying structures in place
        self.structures = deepcopy(structures)
        self.symprec = symprec
        self.reciprocal_coordination = reciprocal_coordination

        # use this to cache benchmark results
        self._benchmark: Dict[NearNeighbors, Dict[str, List]] = defaultdict(dict)

        for name, structure in structures.items():
            if "coordination" not in structure.site_properties:
                raise AttributeError(
                    "{} structure does not have the 'coordination' site "
                    "property".format(name)
                )

        if perturb_sigma:
            for name, structure in self.structures.items():
                self.structures[name] = perturb_einstein_crystal(
                    structure, perturb_sigma
                )

        # precompute the symmetrized structures to save time during the benchmark. Also,
        # determine the total number of unique cations/anions each structure.
        self.site_information: Dict[str, Dict[str, Any]] = {}
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
                if (
                    isinstance(sites[0].specie, Specie)
                    and sites[0].specie.oxi_state >= 0
                ):
                    # based on previous implementation, neutral ions will be scored as
                    # cations, however, neutral to neutral bonding is allowed, hence
                    # don't add the element to the cations set as this will be used to
                    # prevent cation-cation bonding later on.
                    cation_degens.append(len(sites))
                    cation_idxs.append(i)
                    if sites[0].specie.oxi_state > 0:
                        cations.add(sites[0].specie.name)
                elif (
                    isinstance(sites[0].specie, Specie)
                    and sites[0].specie.oxi_state < 0
                ):
                    anion_degens.append(len(sites))
                    anion_idxs.append(i)
                    anions.add(sites[0].specie.name)

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

        if remove_oxidation_states:
            for structure in self.structures.values():
                structure.remove_oxidation_states()

    @classmethod
    def from_structure_group(
        cls, structure_groups: Union[str, List[str]], **kwargs
    ) -> "Benchmark":
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

        filenames: List[Path] = []
        for structure_group in structure_groups:
            if structure_group not in Benchmark.all_structure_groups:
                raise ValueError(
                    '"{}" is not a valid structure group'.format(structure_group)
                )

            filenames.extend(Path(_resource_dir, structure_group).glob("*.json"))

        structures = {}
        for filename in filenames:
            name = Path(filename).stem
            structures[name] = Structure.from_file(str(filename))

        return cls(structures, **kwargs)

    def benchmark(
        self, methods: List[NearNeighbors], return_dataframe: bool = True
    ) -> Union[pd.DataFrame, Dict[NearNeighbors, Dict[str, List[CN_dict]]]]:
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
                        name, method
                    )

        if not return_dataframe:
            return self._benchmark

        method_names = _get_method_names(methods)
        df_data: Dict[str, Dict[str, float]] = defaultdict(dict)
        for method_name, method in zip(method_names, methods):
            for name in self.structures:
                for site_idx in range(self.max_nsites):
                    column = method_name + str(site_idx)

                    if site_idx < len(self.site_information[name]["unique_idxs"]):
                        val = self._benchmark[method][name][site_idx]
                    else:
                        val = None

                    df_data[column][name] = val

        return pd.DataFrame(data=df_data)

    def score(
        self,
        methods: List[NearNeighbors],
        site_type: str = "all",
        cation_anion: bool = False,
        return_raw_site_scores: bool = False,
    ) -> pd.DataFrame:
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
            cation_anion: If True, the score will only include bonding to ions of
                opposing charge. I.e., cation-cation or anion-anion bonding is ignored.
                This option will only affect the scores for input structures that have
                oxidation states.
            return_raw_site_scores: If true, the "raw site scores" are returned. The
                raw scores are given as a list of CN^calc - CN^expected for each
                inequivalent site. Note that the raw_site_score can sometimes be zero
                even if the correct bonding is not not determined. This can arise when
                the coordination numbers themselves are correct but the bonds themselves
                aren't. For example, if the human interpreted bonding is
                {"Cl": 2, "Br": 3} and the predicted bonding is {"Cl": 3, "Br": 2},
                the raw score for that site will be 0.

        Returns:
            The scores as a Pandas DataFrame.
        """
        if site_type != "all" and not self.all_structures_have_oxi:
            warnings.warn(
                "Not all structures have oxidation states, "
                "{} based scoring is not recommended".format(site_type)
            )

        results = self.benchmark(methods, return_dataframe=False)

        method_names = _get_method_names(methods)
        scores: Dict[str, Dict[str, Union[float, List[float]]]] = defaultdict(dict)
        for method_name, method in zip(method_names, methods):
            for name in self.structures:
                scores[method_name][name] = self._score_structure(
                    name,
                    results[method][name],
                    site_type=site_type,
                    cation_anion=cation_anion,
                    return_raw_site_scores=return_raw_site_scores,
                )

        df = pd.DataFrame(data=scores)

        if not return_raw_site_scores:
            df.loc["Total"] = df.sum(axis=0)

        return df

    def _benchmark_structure(self, name: str, nn: NearNeighbors) -> List[CN_dict]:
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
        if self.reciprocal_coordination:
            bonded_structure = nn.get_bonded_structure(self.structures[name])

        for i in self.site_information[name]["unique_idxs"]:
            if self.reciprocal_coordination:
                connected_sites = bonded_structure.get_connected_sites(i)
                cn_dict = _connected_sites_to_cn_dict(connected_sites)
            else:
                cn_dict = nn.get_cn_dict(self.structures[name], i)

            if nn.__class__.__name__ == "MinimumVIRENN":
                cn_dict = {_vire_re.sub("", k): v for k, v in cn_dict.items()}

            # remove oxidation states from the element label if they are present
            # sum the coordinations for the same elements but different oxi states
            tmp_cn_dict: Dict[str, float] = defaultdict(int)
            for k, v in cn_dict.items():
                tmp_cn_dict[_el_re.sub("", k)] += v
            cn_dict = dict(tmp_cn_dict)

            results.append(cn_dict)
        return results

    def _score_structure(
        self,
        name: str,
        all_predictions: List[Dict[str, float]],
        site_type: str = "all",
        cation_anion: bool = False,
        return_raw_site_scores: bool = False,
    ) -> Union[float, List[float]]:
        r"""
        Calculate the score based on a set of predictions.

        The score is calculated by taking the summation of the absolute value
        of the error in CN prediction (CN_calc - CN_expected) multiplied by the
        degeneracy over the total number of sites. I.e.,

        score = \sum_i^N_sites^unique | CN_i^calc - CN_i^expected | * N_i^degen)
                / N_sites

        Args:
            name: A structure name.
            all_predictions: The coordination predictions for each unique site in
                the structure. Formatted as a list of coordination dictionaries.
                See the docstring for `NearNeighbors.get_cn_dict` for the format
                of the dictionaries.
            site_type: Which sites to include when calculating the score.
                Options are "all", "cation", "anion". For example, if "cation"
                is chosen, the scores will only reflect the coordination numbers
                of the cation sites. Note that cation and anion can only be used on
                input structures that have oxidation states.
            cation_anion: If True, the score will only include bonding to ions of
                opposing charge. I.e., cation-cation or anion-anion bonding is ignored.
                This option will only affect the scores for input structures that have
                oxidation states.
            return_raw_site_scores: If true, the "raw site scores" are returned. The
                raw scores are given as a list of CN^calc - CN^expected for each
                inequivalent site. Note that the raw_site_score can sometimes be zero
                even if the correct bonding is not not determined. This can arise when
                the coordination numbers themselves are correct but the bonds themselves
                aren't. For example, if the human interpreted bonding is
                {"Cl": 2, "Br": 3} and the predicted bonding is {"Cl": 3, "Br": 2},
                the raw score for that site will be 0.

        Returns:
            If ``return_raw_site_scores`` is ``False`` (the default), the score will be
            the overall score for the structure as a float.
            If ``return_raw_site_scores`` is ``True``, the score will be a list of the
            raw site scores. See the docstring for ``return_raw_site_scores`` for more
            details.
        """
        structure = self.structures[name]

        # idxs are the indexes of the sites we are interested in IN the list of unique
        # sites. I.e., not in the list of ALL structural sites
        idxs = self.site_information[name]["{}_idxs".format(site_type)]
        degens = self.site_information[name]["{}_degens".format(site_type)]
        total = self.site_information[name]["{}_total".format(site_type)]
        cations = self.site_information[name]["cations"]
        anions = self.site_information[name]["anions"]

        # as site_idx below is the index of the site in the unique sites NOT the
        # index in the overall structure, we first get actual coordinations
        # of the unique sites and then get the index of the sites were are interested in
        coordinations = [
            structure[i].properties["coordination"]
            for i in np.array(self.site_information[name]["unique_idxs"])[idxs]
        ]

        # similarly we want to know the site species types
        elements = [
            structure[i].specie.name
            for i in np.array(self.site_information[name]["unique_idxs"])[idxs]
        ]

        # finally, as the predictions are already provided only for each unique site
        # select the predictions we are interested in
        predictions = [all_predictions[i] for i in idxs]

        if return_raw_site_scores:
            score: Union[float, List[float]] = []
        else:
            score = 0.0

        for site_degen, site_element, prediction, coordination in zip(
            degens, elements, predictions, coordinations
        ):
            # create a list of possible bonding elements (these are the species
            # in both the known coordination dict and the predicted coordination
            # dict)
            bond_elements = set(list(coordination.keys()) + list(prediction.keys()))

            # exclude species that are of the opposite charge
            if cation_anion and site_element in cations:
                bond_elements = bond_elements.difference(cations)
            elif cation_anion and site_element in anions:
                bond_elements = bond_elements.difference(anions)

            site_score = 0
            for bond_element in bond_elements:
                pred_val = prediction.get(bond_element, 0)
                actual_val = coordination.get(bond_element, 0)

                if isinstance(actual_val, list):
                    # there are multiple options for coordination, choose the
                    # value that results in the smallest score
                    diffs = [pred_val - x for x in actual_val]
                    min_idx = np.argmin(list(map(abs, diffs)))
                    raw_site_score = diffs[min_idx]
                else:
                    raw_site_score = pred_val - actual_val

                if return_raw_site_scores:
                    site_score += raw_site_score
                else:
                    site_score += abs(raw_site_score)

            if not return_raw_site_scores:
                score += site_score * site_degen
            else:
                score.append(site_score)  # type: ignore

        if return_raw_site_scores:
            return score

        if total == 0:
            return np.nan

        return score / total


def _get_method_names(methods: List[NearNeighbors]) -> List[str]:
    """
    Get the names of near neighbor methods from a list of methods.

    If multiple near neighbor methods of the same class are present in the list,
    the names will be numbered to differentiate them.

    Args:
        methods: A list of near neighbor methods.

    Returns:
        A list of method names. If each near neighbor method is only present once,
        the method names will just be the class names. However, if a method is present
        more than once, the name will the class name + (the index)
    """
    str_names = [method.__class__.__name__ for method in methods]

    if len(set(str_names)) == len(str_names):
        return str_names

    method_names_counter: Dict[str, int] = defaultdict(int)
    method_names = []
    for name in str_names:
        method_names.append("{}({})".format(name, method_names_counter[name]))
        method_names_counter[name] += 1

    return method_names


def _connected_sites_to_cn_dict(connected_sites):
    counts = Counter([x.site.specie.name for x in connected_sites])
    return dict(counts)
