"""RMSD calculations."""
import os
import numpy as np
from pathlib import Path

from haddock import log
from haddock.core.typing import AtomsDict
from haddock.libs.libalign import get_atoms, load_coords
from haddock.libs.libsubprocess import BaseJob


class RMSDJob(BaseJob):
    """
    Instantiate a subprocess job with inverted args and input.

    Runs with the following scheme, INPUT comes first:

        $ cmd INPUT ARGS
    """

    def make_cmd(self) -> None:
        """Execute job in subprocess."""
        self.cmd = " ".join([
            os.fspath(self.executable),
            os.fspath(self.input),
            ' '.join(map(str, self.args)),  # empty string if no args
            ])
        return


def get_pair(nmodels: int, idx: int) -> tuple[int, int]:
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


def rmsd_dispatcher(nmodels: int, tot_npairs: int,
                    ncores: int) -> tuple[list[int], list[int], list[int]]:
    """Optimal dispatching of rmsd jobs."""
    base_pairs = tot_npairs // ncores
    modulo = tot_npairs % ncores
    npairs: list[int] = []
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


class XYZWriterJob:
    """A Job dedicated to the parallel writing of xyz files."""

    def __init__(
            self,
            xyzwriter_obj):
        """Initialise XYZWriterJob."""
        self.xyzwriter_obj = xyzwriter_obj
        self.output = xyzwriter_obj.output_name

    def run(self):
        """Run this XYZWriterJob."""
        log.info(f"core {self.xyzwriter_obj.core}, running XYZWriter...")
        self.xyzwriter_obj.run()
        return


class XYZWriter:
    """XYZWriter class."""

    def __init__(
            self,
            model_list,
            output_name,
            core,
            n_atoms,
            common_keys,
            filter_resdic,
            allatoms=False,
            ):
        """Initialise Contact class."""
        self.model_list = model_list
        self.output_name = output_name
        self.core = core
        self.n_atoms = n_atoms
        self.common_keys = common_keys
        self.filter_resdic = filter_resdic
        self.allatoms = allatoms
        
    def run(self) -> None:
        """write xyz coordinates."""
        with open(self.output_name, "w") as traj_xyz:
            for mod in self.model_list:
                atoms: AtomsDict = get_atoms(mod, self.allatoms)

                ref_coord_dic, _ = load_coords(
                mod, atoms, self.filter_resdic
                )
                # now we filter the dictionary with the common keys
                common_coord_dic = {k: v for k, v in ref_coord_dic.items() if k in self.common_keys}  
                # write header
                traj_xyz.write(f"{self.n_atoms}{os.linesep}{os.linesep}")
                # write the coordinates
                for k in self.common_keys:
                    v = common_coord_dic[k]
                    at_string = ''.join([str(el) for el in k])
                    traj_xyz.write(f"{at_string} {v[0]} {v[1]} {v[2]}{os.linesep}")
        return
