import pickle
import random
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from haddock.libs.libmpi import MPIScheduler


@pytest.fixture
def mpischeduler():
    with tempfile.TemporaryDirectory() as tempdir:

        scheduler = MPIScheduler(tasks=[1, 2, 3], ncores=4)
        scheduler.cwd = Path(tempdir)

        yield scheduler


def test_mpischeduler_init():

    cores = random.randint(1, 99)
    mpi = MPIScheduler(tasks=[1, 2, 3], ncores=cores)

    assert mpi.tasks == [1, 2, 3]
    assert mpi.cwd == Path.cwd()
    assert mpi.ncores == cores


def test_mpischduler_run(mocker, mpischeduler):
    # Mock the necessary methods and objects
    mock_pickle_tasks = mocker.patch.object(
        mpischeduler, "_pickle_tasks", return_value="mocked_pkl_tasks"
    )
    mock_subprocess_run = mocker.patch("subprocess.run")
    mock_sys_exit = mocker.patch("sys.exit")

    # Set up the mock subprocess.run return value
    mock_process = MagicMock()
    mock_process.stderr.decode.return_value = ""
    mock_subprocess_run.return_value = mock_process

    # Call the run method
    mpischeduler.run()

    # Assertions
    mock_pickle_tasks.assert_called_once()

    mock_subprocess_run.assert_called_once_with(
        [
            "mpirun",
            "-np",
            str(mpischeduler.ncores),
            "haddock3-mpitask",
            "mocked_pkl_tasks",
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    mock_sys_exit.assert_not_called()

    # Test error case
    mock_process.stderr.decode.return_value = "Error occurred"
    mpischeduler.run()
    mock_sys_exit.assert_called_once()


def test__pickle_tasks(mpischeduler):

    result = mpischeduler._pickle_tasks()

    expected_path = Path(mpischeduler.cwd, "mpi.pkl")
    assert result == expected_path
    assert expected_path.exists()

    # Verify the content of the pickled file
    with open(expected_path, "rb") as f:
        unpickled_tasks = pickle.load(f)
    assert unpickled_tasks == mpischeduler.tasks
