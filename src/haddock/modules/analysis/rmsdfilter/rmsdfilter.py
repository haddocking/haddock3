"""RMSD calculation for the rmsdfilter module."""

import copy
from pathlib import Path
import numpy as np
from haddock import log
from haddock.libs.libalign import (
    ALIGNError,
    calc_rmsd,
    centroid,
    get_atoms,
    kabsch,
    load_coords,
)
from haddock.libs.libontology import PDBFile


class RMSDFilter:
    """Compute RMSD between model and reference."""

    def __init__(
        self,
        identificator: int,
        model: PDBFile,
        reference: Path,
        path: Path,
        params: dict,
        align_func,
        ref_id: int = 1,
    ) -> None:
        self.identificator = identificator
        self.model = model
        self.reference = reference
        self.path = path
        self.params = params
        self.align_func = align_func
        self.ref_id = ref_id
        self.rmsd = float("nan")

    def run(self) -> "RMSDFilter":
        """Compute RMSD and return self."""
        allatoms = self.params["allatoms"]
        keep_hetatm = self.params["keep_hetatm"]

        # map model residue numbers onto reference residue numbers
        try:
            model2ref_numbering, model2ref_chain_dict = self.align_func(
                self.reference, self.model, self.path
            )
        except ALIGNError:
            log.warning(
                f"Alignment failed between {self.reference} "
                f"and {self.model}, skipping..."
            )
            # deepcopy so the scheduler gets an independent object with rmsd=nan
            return copy.deepcopy(self)

        atoms = get_atoms(self.model, full=allatoms)
        atoms.update(get_atoms(self.reference, full=allatoms))

        ref_coord_dic, _ = load_coords(
            self.reference,
            atoms,
            keep_hetatm=keep_hetatm,
        )
        try:
            # numbering_dic remaps model residues to reference numbering so that
            # coordinate keys are directly comparable between the two dicts
            mod_coord_dic, _ = load_coords(
                self.model,
                atoms,
                numbering_dic=model2ref_numbering,
                model2ref_chain_dict=model2ref_chain_dict,
                keep_hetatm=keep_hetatm,
            )
        except ALIGNError as e:
            log.warning(e)
            return copy.deepcopy(self)

        # Only superpose on shared atoms between model and reference
        common_keys = ref_coord_dic.keys() & mod_coord_dic.keys()
        if not common_keys:
            log.warning(
                f"No common atoms found between {self.reference} and {self.model}"
            )
            return copy.deepcopy(self)

        Q = np.asarray([ref_coord_dic[k] for k in common_keys])
        P = np.asarray([mod_coord_dic[k] for k in common_keys])

        # Kabsch superposition: centre both structures Q (ref) and P (model), 
        # find optimal rotation U, apply it to P, then compute RMSD on the superposed coordinates.
        Q = Q - centroid(Q)
        P = P - centroid(P)
        U = kabsch(P, Q)
        P = np.dot(P, U)

        self.rmsd = calc_rmsd(P, Q)
        return copy.deepcopy(self)
