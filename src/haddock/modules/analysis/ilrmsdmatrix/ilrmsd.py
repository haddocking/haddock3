"""ilRMSD calculations."""
import os
from pathlib import Path

import numpy as np

from haddock import log
from haddock.core.typing import NDFloat, Any
from haddock.libs.libalign import (
    get_atoms,
    )
from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.caprieval.capri import load_contacts


class ContactJob:
    """A Job dedicated to the fast retrieval of receptor contacts."""

    def __init__(
            self,
            output,
            params,
            contact_obj):

        log.info(f"core {contact_obj.core}, initialising Contact...")
        self.params = params
        self.output = output
        self.contact_obj = contact_obj
        log.info(f"core {contact_obj.core}, Contact initialised")

    def run(self):
        """Run this ContactJob."""
        log.info(f"core {self.contact_obj.core}, running Contact...")
        self.contact_obj.run()
        self.contact_obj.output()
        return


class Contact:
    """Contact class."""

    def __init__(
            self,
            model_list: list[PDBFile],
            output_name: str,
            core: int,
            path: Path,
            contact_distance_cutoff: float = 5.0,
            **params: dict[str, Any],
            ):
        """Initialise Contact class."""
        self.model_list = model_list
        self.output_name = output_name
        self.core = core
        self.path = path
        self.contact_distance_cutoff = contact_distance_cutoff
        self.atoms = {}
        self.receptor_chain = params["params"]["receptor_chain"]
        self.ligand_chains = params["params"]["ligand_chains"]
        self.unique_rec_res: NDFloat = []
        self.unique_lig_res: NDFloat = []
        for m in model_list:
            self.atoms.update(get_atoms(m))
        
    def run(self):
        """Run contact calculations."""
        nmodels = len(self.model_list)
        # create empty set of receptor interface residues
        receptor_interface_residues = []
        ligand_interface_residues = []
        for n in range(nmodels):
            contacts = load_contacts(
                self.model_list[n],
                cutoff=self.contact_distance_cutoff,
                )
            rec_resids = [
                con[1] if con[0] == self.receptor_chain
                else con[3] for con in contacts
                ]
            lig_resids = [
                con[1] if con[0] in self.ligand_chains
                else con[3] for con in contacts
                ]
            if rec_resids != []:
                receptor_interface_residues.append(np.unique(rec_resids))
            if lig_resids != []:
                ligand_interface_residues.append(np.unique(lig_resids))
        # concatenate all the receptor residues
        if receptor_interface_residues != []:
            rec_np_arr = np.concatenate(receptor_interface_residues)
            self.unique_rec_res = np.unique(rec_np_arr)
        # now the ligand residues
        if ligand_interface_residues != []:
            lig_np_arr = np.concatenate(ligand_interface_residues)
            self.unique_lig_res = np.unique(lig_np_arr)
       
    def output(self):
        """Write down unique contacts to file."""
        output_fname = Path(self.path, self.output_name)
        with open(output_fname, "w") as out_fh:
            # first the receptor
            out_fh.write(f"{self.receptor_chain} ")
            out_fh.write(" ".join([str(res) for res in self.unique_rec_res]))
            out_fh.write(f"{os.linesep}")
            # now the ligand
            for chain in self.ligand_chains:
                out_fh.write(f"{chain} ")
                res_str = " ".join([str(res) for res in self.unique_lig_res])
                out_fh.write(res_str)
                out_fh.write(f"{os.linesep}")
