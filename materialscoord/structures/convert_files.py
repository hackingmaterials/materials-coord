from glob import glob
from pathlib import Path

import numpy as np
from monty.serialization import loadfn, dumpfn

from pymatgen import Structure
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.util.provenance import StructureNL

answers = loadfn("human_interpreter.yaml")
folders = glob('*/')

for folder in folders:
    structure_files = glob('{}/*.cif'.format(folder))
    print(folder)

    if folder == "problematic_structures/":
        continue

    for structure_file in structure_files:
        cif = CifParser(structure_file)
        # structure = cif.get_structures()[0]
        reference = cif.get_bibtex_string()
        structure = Structure.from_file(structure_file)

        name = Path(structure_file).stem
        print(name)

        # es = SpacegroupAnalyzer(structure).get_symmetrized_structure().equivalent_sites
        # reordered_sites = [structure.index(y) for x in es for y in x]
        #
        # structure = Structure(
        #     structure.lattice,
        #     [structure.species[x] for x in reordered_sites],
        #     [structure.frac_coords[x] for x in reordered_sites],
        #     )
        # # es = SpacegroupAnalyzer(structure).get

        # es = SpacegroupAnalyzer(structure).get_symmetrized_structure().equivalent_sites
        # structure = SpacegroupAnalyzer(structure).get_symmetrized_structure()
        # print(structure)
        sga = SpacegroupAnalyzer(structure, symprec=1)
        equivalent_atoms = sga.get_symmetry_dataset()["equivalent_atoms"]
        _, inverse_idx = np.unique(equivalent_atoms, return_inverse=True)

        cn_dicts = [answers[name][idx] for idx in inverse_idx]

        cns = []
        for cn_dict in cn_dicts:
            values = list(cn_dict.values())
            cns.append(values[0])

        structure.add_site_property("coordination", cns)
        snl = StructureNL(structure, ["Alex Ganose <aganose@lbl.gov>"],
                          references=reference,
                          remarks="Structure is part of MaterialsCoord test set")
        dumpfn(snl, "{}/{}.json".format(folder, name), indent=4)

