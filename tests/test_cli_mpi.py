from unittest import mock

import pytest

from haddock.clis import cli_mpi
from haddock.clis.cli_mpi import get_mpi, split_tasks


def test_cli_has_maincli():
    """
    Test maincli func in CLI.

    maincli is used in setup.py.
    """
    assert cli_mpi.maincli


def test_cli_has_comm():
    """
    Test COMM object in CLI.

    COMM is used in main function.
    """
    assert cli_mpi.COMM is None


def test_split_tasks():

    t = split_tasks([1, 2, 3, 4, 5, 6, 7, 8, 9], 3)

    assert t == [[1, 4, 7], [2, 5, 8], [3, 6, 9]]


def test_get_mpi_success():
    # Mock the import of mpi4py.MPI
    with mock.patch.dict(
        "sys.modules", {"mpi4py": mock.Mock(), "mpi4py.MPI": mock.Mock()}
    ):
        from mpi4py import MPI as _MPI

        MPI, COMM = get_mpi()
        assert MPI == _MPI
        assert COMM == _MPI.COMM_WORLD


def test_get_mpi_import_error():
    # Mock the absence of mpi4py.MPI to simulate ImportError
    with mock.patch.dict("sys.modules", {"mpi4py": None}):
        with mock.patch("sys.exit") as mock_exit:
            get_mpi()
            mock_exit.assert_called_once()


# Cleanup fixture to reset global state after each test
@pytest.fixture(autouse=True)
def cleanup():
    yield
    import haddock.clis.cli_mpi

    haddock.clis.cli_mpi.MPI = None
    haddock.clis.cli_mpi.COMM = None
