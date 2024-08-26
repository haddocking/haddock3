"""RMSD calculations."""

import os
import numpy as np
from pathlib import Path
import copy
import subprocess
import shlex

from haddock.core.typing import AtomsDict, Optional
from haddock.libs.libontology import PDBFile
from haddock.libs.libalign import get_atoms, load_coords


class RMSDJob:
    """RMSDJob class.

    This is a very similar class to `BaseJob` in `libsubprocess.py`.
    It is re-defined here because we need a special implementation for the `run` method
    that is not compatible with the `BaseJob` class.
    """

    def __init__(self, input, output, executable, *args) -> None:
        self.input = input
        self.output = output
        self.executable = executable
        self.args = args
        self.output_data = ""

    def make_cmd(self) -> None:
        """Execute job in subprocess."""

        # =======================================================================
        # IMPORTANT: `fast-rmsdmatrix` will always write a file!
        #
        # here it will be `self.args[0]`, for example if `self.args[0]` is "42",
        #  it will produce a file named `rmsd_42.matrix` and write the output there.
        #
        # =======================================================================

        self.cmd = " ".join(
            [
                os.fspath(self.executable),
                os.fspath(self.input),
                " ".join(map(str, self.args)),  # empty string if no args
            ]
        )
        return

    def run(self) -> "RMSDJob":
        """Execute job in subprocess."""

        self.make_cmd()

        p = subprocess.Popen(
            shlex.split(self.cmd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            close_fds=True,
        )

        # FIXME: Add error handling
        _, _ = p.communicate()

        # =============================================================================
        # IMPORTANT: `fast-rmsdmatrix` will always write a file!
        #
        # This part here can only be removed if the `fast-rmsdmatrix` code is changed to
        #  write to stdout instead of a file!
        #
        assert Path(self.output).exists(), f"Output file {self.output} does not exist"
        with open(self.output, "r") as f:
            stdout = f.read()

        Path(self.output).unlink()
        # =============================================================================

        self.output_data = stdout

        return copy.deepcopy(self)


def get_pair(nmodels: int, idx: int) -> tuple[int, int]:
    """Get the pair of structures given the 1D matrix index."""
    if nmodels < 0 or idx < 0:
        err = "get_pair cannot accept negative numbers"
        err += f"Input is {nmodels} , {idx}"
        raise ValueError(err)
    # solve the second degree equation
    b = 1 - (2 * nmodels)
    i = (-b - np.sqrt(b**2 - 8 * idx)) // 2
    j = idx + i * (b + i + 2) // 2 + 1
    return (int(i), int(j))


def rmsd_dispatcher(
    nmodels: int, tot_npairs: int, ncores: int
) -> tuple[list[int], list[int], list[int]]:
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


# TODO: Move this class into a library, since it's used by other modules
class XYZWriter:
    """XYZWriter class."""

    def __init__(
        self,
        model: PDBFile,
        n_atoms: int,
        common_keys: list[str],
        allatoms: bool = False,
        filter_resdic: Optional[dict] = None,
        output_name: Optional[Path] = None,
    ):
        """Initialise Contact class."""
        self.model = model
        self.output_name = output_name
        self.n_atoms = n_atoms
        self.common_keys = common_keys
        self.filter_resdic = filter_resdic
        self.allatoms = allatoms
        self.output_data = ""

    def run(self) -> "XYZWriter":
        """write xyz coordinates."""
        atoms: AtomsDict = get_atoms(self.model, self.allatoms)

        ref_coord_dic, _ = load_coords(self.model, atoms, self.filter_resdic)
        # now we filter the dictionary with the common keys
        common_coord_dic = {
            k: v for k, v in ref_coord_dic.items() if k in self.common_keys
        }
        # write header
        self.output_data += f"{self.n_atoms}{os.linesep}{os.linesep}"
        # write the coordinates
        for k in self.common_keys:
            v = common_coord_dic[k]
            at_string = "".join([str(el) for el in k])
            self.output_data += f"{at_string} {v[0]} {v[1]} {v[2]}{os.linesep}"

        return copy.deepcopy(self)

    def write_output(self):
        """Write the output data to the output file."""
        if self.output_name:
            with open(self.output_name, "w") as f:
                f.write(self.output_data)
