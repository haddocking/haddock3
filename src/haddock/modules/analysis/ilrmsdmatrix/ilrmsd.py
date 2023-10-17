"""ilRMSD calculations."""
import os
from pathlib import Path

import numpy as np

from haddock import log
from haddock.libs.libalign import (
    calc_rmsd,
    centroid,
    get_atoms,
    kabsch,
    load_coords,
    make_range,
    )
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
        """Run this RMSDJob."""
        log.info(f"core {self.contact_obj.core}, running ilRMSD...")
        self.contact_obj.run()
        self.contact_obj.output()
        return


class Contact:
    """Contact class."""

    def __init__(
            self,
            model_list,
            output_name,
            core,
            path,
            **params,
            ):
        """Initialise Contact class."""
        self.model_list = model_list
        self.output_name = output_name
        self.core = core
        self.path = path
        self.atoms = {}
        self.receptor_chain = params["params"]["receptor_chain"]
        self.ligand_chains = params["params"]["ligand_chains"]
        self.unique_rec_res = []
        self.unique_lig_res = []
        for m in model_list:
            self.atoms.update(get_atoms(m))
        
    def run(self):
        """Run contact calculations."""
        nmodels = len(self.model_list)
        # create empty set of receptor interface residues
        receptor_interface_residues = []
        ligand_interface_residues = []
        for n in range(nmodels):
            contacts = load_contacts(self.model_list[n], cutoff=5.0)
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
        # concatenate all the receptor residues
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
        

class ilRMSDJob:
    """A Job dedicated to the fast rmsd calculation."""

    def __init__(
            self,
            output,
            params,
            ilrmsd_obj):

        log.info(f"core {ilrmsd_obj.core}, initialising ilRMSD...")
        log.info(f"core {ilrmsd_obj.core}, # of pairs : {ilrmsd_obj.npairs}")
        self.output = output
        self.params = params
        self.ilrmsd_obj = ilrmsd_obj
        log.info(f"core {ilrmsd_obj.core}, ilRMSD initialised")

    def run(self):
        """Run this RMSDJob."""
        log.info(f"core {self.ilrmsd_obj.core}, running ilRMSD...")
        self.ilrmsd_obj.run()
        self.ilrmsd_obj.output()
        return


class ilRMSD:
    """RMSD class."""

    def __init__(
            self,
            model_list,
            core,
            npairs,
            start_ref,
            start_mod,
            filter_resdic,
            output_name,
            path,
            **params,
            ):
        """
        Initialise RMSD class.

        Parameters
        ----------

        model_list : list
            List of models

        core : int
            index of the current core

        npairs : int
            the number of pairs of structures

        start_ref : int
            the index of the first reference structure

        start_mod : int
            the index of the first mobile structure. The class performs npairs
            RMSD calculations starting from the pair (start_ref, start_mod)

        output_name : str
            name of the core-specific output file

        path : pathlib.Path
            path to the current directory

        **params : dict
            additional parameters
        """
        self.model_list = model_list
        self.core = core
        self.npairs = npairs
        self.start_ref = start_ref
        self.start_mod = start_mod
        self.receptor_chain = params["params"]["receptor_chain"]
        self.filter_resdic = filter_resdic
        self.output_name = output_name
        self.path = path
        self.atoms = {}
        for m in model_list:
            self.atoms.update(get_atoms(m))
        # data array
        self.data = np.zeros((self.npairs, 3))

    def calc_ilrmsd(self, ref, mod):
        """
        Calculate ilRMSD.

        Parameters
        ----------
        ref : int
            index of the reference structure
        
        mod : int
            index of the mobile structure
        
        Returns
        -------
        ilrmsd : float
            ilRMSD value
        """
        ref_coord_dic, _ = load_coords(
            self.model_list[ref], self.atoms, self.filter_resdic
        )
        
        mod_coord_dic, _ = load_coords(
            self.model_list[mod], self.atoms, self.filter_resdic
        )
        # calculating ilRMSD
        # find atoms present in both interfaces
        Q_int = []
        P_int = []
        common_keys = ref_coord_dic.keys() & mod_coord_dic.keys()
        for k in sorted(common_keys):
            ref_xyz = ref_coord_dic[k]
            mod_xyz = mod_coord_dic[k]
            Q_int.append(ref_xyz)
            P_int.append(mod_xyz)
        Q_int = np.asarray(Q_int)
        P_int = np.asarray(P_int)

        chain_ranges = {}
        for i, segment in enumerate(sorted(common_keys)):
            chain, _, _ = segment
            if chain not in chain_ranges:
                chain_ranges[chain] = []
            chain_ranges[chain].append(i)
        chain_ranges = make_range(chain_ranges)
        obs_chains = list(chain_ranges.keys())  # observed chains
        if len(obs_chains) < 2:
            log.warning("Not enough chains for calculating ilrmsd")
            ilrmsd = None
        else:
            r_chain = self.receptor_chain
            l_chains = [c for c in obs_chains if c != self.receptor_chain]
            r_start, r_end = chain_ranges[r_chain]
            l_starts = [chain_ranges[l_chain][0] for l_chain in l_chains]
            l_ends = [chain_ranges[l_chain][1] for l_chain in l_chains]
            # put system at origin of the receptor interface
            Q_r_int = Q_int[r_start: r_end + 1]
            P_r_int = P_int[r_start: r_end + 1]
            Q_int = Q_int - centroid(Q_r_int)
            P_int = P_int - centroid(P_r_int)
            # put interfaces at the origin
            # find the rotation that minimizes the receptor interface rmsd
            Q_r_int = Q_int[r_start: r_end + 1]
            P_r_int = P_int[r_start: r_end + 1]
            U_int = kabsch(P_r_int, Q_r_int)
            P_int = np.dot(P_int, U_int)
            
            # Identify ligand coords concatenating all the ligand chains
            Q_l_int = np.empty((0, 3))
            P_l_int = np.empty((0, 3))
            for l_st, l_end in zip(l_starts, l_ends):
                Q_l_int = np.concatenate((Q_l_int, Q_int[l_st: l_end + 1]))
                P_l_int = np.concatenate((P_l_int, P_int[l_st: l_end + 1]))
            # # this will be the interface-ligand-rmsd
            ilrmsd = calc_rmsd(P_l_int, Q_l_int)
            return ilrmsd


    def run(self):
        """Run calculations."""
        # initialising the number of pairs
        ref = self.start_ref
        mod = self.start_mod
        nmodels = len(self.model_list)
        for n in range(self.npairs):
            # calculating RMSD
            ilrmsd = self.calc_ilrmsd(ref, mod)
            if ilrmsd is None:
                break
            # saving output (adding one for consistency with clusterfcc)
            self.data[n, 0] = ref + 1
            self.data[n, 1] = mod + 1
            self.data[n, 2] = ilrmsd
            # updating indices
            if mod == (nmodels - 1):
                ref += 1
                mod = ref + 1
            else:
                mod += 1

    def output(
            self,
            ):
        """Write down the RMSD matrix."""
        output_fname = Path(self.path, self.output_name)
        # check if there are very low values in the RMSD vector
        check_low_values = np.isclose(
            self.data[:, 2],
            np.zeros(self.npairs),
            atol=0.1
            ).any()
        if check_low_values:
            log.warning(f"core {self.core}: low values of RMSD detected.")
        with open(output_fname, "w") as out_fh:
            for data in list(self.data):
                data_str = f"{data[0]:.0f} {data[1]:.0f} {data[2]:.3f}"
                data_str += os.linesep
                out_fh.write(data_str)


def get_pair(nmodels, idx):
    """Get the pair of structures given the 1D matrix index."""
    if (nmodels < 0 or idx < 0):
        err = "get_pair cannot accept negative numbers"
        err += f"Input is {nmodels} , {idx}"
        raise ValueError(err)
    # solve the second degree equation
    b = 1 - (2 * nmodels)
    i = (-b - np.sqrt(b ** 2 - 8 * idx)) // 2
    j = idx + i * (b + i + 2) // 2 + 1
    return (int(i), int(j))
