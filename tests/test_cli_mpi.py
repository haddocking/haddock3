from mpi4py import MPI

from haddock.clis import cli_mpi
from haddock.clis.cli_mpi import split_tasks


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
    assert cli_mpi.COMM

    assert isinstance(cli_mpi.COMM, MPI.Intracomm)


def test_split_tasks():

    t = split_tasks([1, 2, 3, 4, 5, 6, 7, 8, 9], 3)

    assert t == [[1, 4, 7], [2, 5, 8], [3, 6, 9]]
