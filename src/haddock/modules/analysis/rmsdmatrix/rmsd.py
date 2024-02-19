"""RMSD calculations."""
import os
import numpy as np

from haddock import log
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
        print(f"executing command {self.cmd}")
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
