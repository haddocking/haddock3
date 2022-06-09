"""RMSD calculations."""
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
    )


class RMSDJob:
    """A Job dedicated to the fast rmsd calculation."""

    def __init__(
            self,
            output,
            params,
            rmsd_obj):

        log.info(f"core {rmsd_obj.core}, initialising RMSD...")
        log.info(f"core {rmsd_obj.core}, # of pairs : {rmsd_obj.npairs}")
        self.output = output
        self.params = params
        self.rmsd_obj = rmsd_obj
        log.info(f"core {rmsd_obj.core}, RMSD initialised")

    def run(self):
        """Run this RMSDJob."""
        log.info(f"core {self.rmsd_obj.core}, running RMSD...")
        self.rmsd_obj.run()
        self.rmsd_obj.output()
        return


class RMSD:
    """RMSD class."""

    def __init__(
            self,
            model_list,
            core,
            npairs,
            start_ref,
            start_mod,
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
        # choice of atoms
        if "params" in params.keys():
            self.filter_resdic = {
                key[-1]: value for key, value
                in params["params"].items()
                if key.startswith("resdic")
                }
        else:
            self.filter_resdic = {}
        if self.filter_resdic:
            log.info(f"Using filtering dictionary {self.filter_resdic}")
        else:
            log.info("No filtering dictionary, using all residues")
        self.output_name = output_name
        self.path = path
        self.atoms = {}
        for m in model_list:
            self.atoms.update(get_atoms(m))
        # data array
        self.data = np.zeros((self.npairs, 3))

    def run(self):
        """Run calculations."""
        # initialising the number of pairs
        ref = self.start_ref
        mod = self.start_mod
        nmodels = len(self.model_list)
        for n in range(self.npairs):

            ref_coord_dic, _ = load_coords(
                self.model_list[ref], self.atoms, self.filter_resdic
                )

            mod_coord_dic, _ = load_coords(
                self.model_list[mod], self.atoms, self.filter_resdic
                )
            P = []
            Q = []
            for k in ref_coord_dic.keys() & mod_coord_dic.keys():
                ref_xyz = ref_coord_dic[k]
                mod_xyz = mod_coord_dic[k]
                Q.append(ref_xyz)
                P.append(mod_xyz)
            Q = np.asarray(Q)
            P = np.asarray(P)
            Q = Q - centroid(Q)
            P = P - centroid(P)
            U = kabsch(P, Q)
            P = np.dot(P, U)
            rmsd = calc_rmsd(P, Q)
            # saving output (adding one for consistency with clusterfcc)
            self.data[n, 0] = ref + 1
            self.data[n, 1] = mod + 1
            self.data[n, 2] = rmsd
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
        raise Exception(err)
    # solve the second degree equation
    b = 1 - (2 * nmodels)
    i = (-b - np.sqrt(b ** 2 - 8 * idx)) // 2
    j = idx + i * (b + i + 2) // 2 + 1
    return (int(i), int(j))


def rmsd_dispatcher(nmodels, tot_npairs, ncores):
    """Optimal dispatching of rmsd jobs."""
    base_pairs = tot_npairs // ncores
    modulo = tot_npairs % ncores
    npairs = []
    for core in range(ncores):
        if core < modulo:
            npairs.append(base_pairs + 1)
        else:
            npairs.append(base_pairs)
    # each core must know how many pairs and where to start
    index = 0
    start_structures = [0]
    end_structures = [1]
    for el in npairs[:-1]:
        index += el
        pair = get_pair(nmodels, index)
        start_structures.append(pair[0])
        end_structures.append(pair[1])
    return npairs, start_structures, end_structures
