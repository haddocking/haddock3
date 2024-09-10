import tempfile
from pathlib import Path
from typing import Generator

import pytest

from haddock.clis.cli_mpi import main as cli_mpi_main
from haddock.libs.libmpi import MPIScheduler
from haddock.libs.libsubprocess import CNSJob
import subprocess
import os

from . import CNS_EXEC


class Task:
    """A task to be executed by the MPIScheduler."""

    def __init__(self, output_fname: Path):
        self.output_fname = output_fname

    def run(self):
        self.output_fname.touch()


@pytest.fixture
def cnsjob() -> Generator[CNSJob, None, None]:
    """A simple CNSJob that writes a message to a file."""
    with tempfile.NamedTemporaryFile(suffix=".out") as out_f:
        cns_inp = f"""
set display={out_f.name} end
display hello_from_test_mpi
close {out_f.name} end
stop"""

        yield CNSJob(
            input_file=cns_inp,
            output_file=Path(out_f.name),
            cns_exec=CNS_EXEC,
        )


@pytest.fixture
def task() -> Generator[Task, None, None]:
    """A simple task that creates a file."""
    with tempfile.NamedTemporaryFile(suffix=".out") as out_f:
        yield Task(output_fname=Path(out_f.name))


def test_run_cli_mpi_main_cnsjob(cnsjob):
    """Test the cli_mpi_main function with a CNSJob."""
    with tempfile.TemporaryDirectory() as tempdir:
        scheduler = MPIScheduler(tasks=[cnsjob], ncores=1)
        scheduler.cwd = Path(tempdir)
        scheduler._pickle_tasks()

        pickled_tasks = Path(scheduler.cwd, "mpi.pkl")
        assert pickled_tasks.exists()
        assert pickled_tasks.stat().st_size > 0

        cli_mpi_main(pickled_tasks)

        scheduler.tasks[0].output_file.exists()


def test_run_cli_mpi_main_task(task):
    """Test the cli_mpi_main function with a Task."""
    with tempfile.TemporaryDirectory() as tempdir:
        scheduler = MPIScheduler(tasks=[task], ncores=1)
        scheduler.cwd = Path(tempdir)
        scheduler._pickle_tasks()

        pickled_tasks = Path(scheduler.cwd, "mpi.pkl")
        assert pickled_tasks.exists()
        assert pickled_tasks.stat().st_size > 0

        cli_mpi_main(pickled_tasks)

        assert scheduler.tasks[0].output_fname.exists()


@pytest.fixture(autouse=True)
def cleanup():
    """Cleanup the MPI environment after each test."""
    # IMPORTANT! #####################################################################
    # `mpi4py` when imported, sets up the MPI environment. This environment is
    # shared across all tests, so we need to clean it up after each test. to
    # avoid interference between tests.
    # IMPORTANT! #####################################################################
    yield
    import haddock.clis.cli_mpi

    haddock.clis.cli_mpi.MPI = None
    haddock.clis.cli_mpi.COMM = None

    result = subprocess.run(["env"], capture_output=True, text=True)
    mpi_vars = [
        line.split("=")[0] for line in result.stdout.split("\n") if "MPI" in line
    ]
    for var in mpi_vars:
        os.environ[var] = ""
