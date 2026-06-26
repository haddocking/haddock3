"""Module in charge of MPI execution of tasks."""
import pickle
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional
from shutil import which

from haddock import log
from haddock.core.exceptions import HaddockTermination


class MPIScheduler:
    """Schedules tasks to be executed via MPI."""

    def __init__(self, tasks: list[Any], ncores: Optional[int] = None) -> None:
        self.tasks = tasks
        self.cwd = Path.cwd()
        self.ncores = ncores

    def run(self) -> None:
        """Send it to the haddock3-mpitask runner."""
        pkl_tasks = self._pickle_tasks()
        log.info(
            f"Executing tasks with the haddock3-mpitask runner using "
            f"{self.ncores} processors..."
            )
        if which("mpirun") is not None:
            cmd = f"mpirun -np {self.ncores} haddock3-mpitask {pkl_tasks}"
        elif which("srun") is not None:
            cmd = f"srun haddock3-mpitask {pkl_tasks}"
        else:
            log.error("mpirun or srun are not available on the system")
            log.error("Terminating run!")
            raise HaddockTermination

        log.debug(f"MPI cmd is {cmd}")

        p = subprocess.run(
            shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )

        # out = p.stdout.decode("utf-8")
        err = p.stderr.decode("utf-8")

        if err:
            log.error(err)
            log.error("Terminating run!")
            raise HaddockTermination

    def _pickle_tasks(self) -> Path:
        """Pickle the tasks."""
        fpath = Path(self.cwd, "mpi.pkl")
        log.debug(f"Pickling the tasks at {fpath}")
        with open(fpath, "wb") as output_handler:
            pickle.dump(self.tasks, output_handler)
        return fpath
