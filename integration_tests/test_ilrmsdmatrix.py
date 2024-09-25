"""Intergration tests for the IL-RMSD matrix module."""

import os
import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.ilrmsdmatrix import \
    DEFAULT_CONFIG as DEFAULT_ILRMSD_CONFIG
from haddock.modules.analysis.ilrmsdmatrix import \
    HaddockModule as IlrmsdmatrixModule

from integration_tests import GOLDEN_DATA


@pytest.fixture
def ilrmsdmatrix_module():
    """Provide a parametrized IL-RMSD matrix module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield IlrmsdmatrixModule(
            order=0, path=tmpdir, initial_params=DEFAULT_ILRMSD_CONFIG
        )


class MockPreviousIO:
    """Mock of the ModuleIO class."""

    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        """Mock of the io.retrive_models from ModuleIO."""
        shutil.copy(
            Path(GOLDEN_DATA, "protglyc_complex_1.pdb"),
            Path(self.path, "protglyc_complex_1.pdb"),
        )
        shutil.copy(
            Path(GOLDEN_DATA, "protglyc_complex_2.pdb"),
            Path(self.path, "protglyc_complex_2.pdb"),
        )

        model_list = [
            PDBFile(file_name="protglyc_complex_1.pdb", path=self.path),
            PDBFile(file_name="protglyc_complex_2.pdb", path=self.path),
        ]

        return model_list


class MockPreviousIO_protprot:
    """Mock previous IO class."""

    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        """Retrieve models."""
        shutil.copy(
            Path(GOLDEN_DATA, "protprot_complex_1.pdb"),
            Path(self.path, "protprot_complex_1.pdb"),
        )
        shutil.copy(
            Path(GOLDEN_DATA, "protprot_complex_2.pdb"),
            Path(self.path, "protprot_complex_2.pdb"),
        )
        model_list = [
            PDBFile(file_name="protprot_complex_1.pdb", path=self.path),
            PDBFile(file_name="protprot_complex_2.pdb", path=self.path),
        ]
        return model_list


def test_ilrmsdmatrix_default(ilrmsdmatrix_module, mocker):
    """Test the topoaa module."""
    ilrmsdmatrix_module.previous_io = MockPreviousIO(path=ilrmsdmatrix_module.path)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )
    ilrmsdmatrix_module.run()

    exp_ilrmsd_matrix = Path(ilrmsdmatrix_module.path, "ilrmsd.matrix")
    exp_contacts_file = Path(ilrmsdmatrix_module.path, "receptor_contacts.con")

    assert exp_ilrmsd_matrix.exists(), "ilrmsd.matrix does not exist"
    assert exp_contacts_file.exists(), "receptor_contacts.con does not exist"

    with open(exp_ilrmsd_matrix) as f:
        assert f.readline() == f"1 2 11.877{os.linesep}"

    with open(exp_contacts_file) as f:
        lines = f.readlines()
        assert (
            lines[0]
            == f"A 35 44 46 48 50 52 56 57 58 59 61 62 63 73 75 98 101 102 103 104 106 107 108 109 110 112{os.linesep}"
        )  # noqa : E501
        assert lines[1] == f"B 1 2 3 4{os.linesep}"


def test_ilrmsdmatrix_run(ilrmsdmatrix_module, mocker):
    """Test ilrmsdmatrix run method."""
    ilrmsdmatrix_module.previous_io = MockPreviousIO_protprot(
        path=ilrmsdmatrix_module.path
    )
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )

    ilrmsdmatrix_module.run()

    ilrmsd_matrix_f = Path(ilrmsdmatrix_module.path, "ilrmsd.matrix")
    assert ilrmsd_matrix_f.exists()

    receptor_con_f = Path(ilrmsdmatrix_module.path, "receptor_contacts.con")
    assert receptor_con_f.exists()

    with open(ilrmsd_matrix_f, encoding="utf-8", mode="r") as f:
        assert f.readline() == f"1 2 16.715{os.linesep}"

    with open(
        receptor_con_f,
        encoding="utf-8",
        mode="r",
    ) as f:
        lines = f.readlines()
        assert (
            lines[0]
            == f"A 37 38 39 40 43 44 45 69 71 72 75 90 93 94 96 132{os.linesep}"
        )
        assert lines[1] == f"B 10 11 12 16 17 48 51 52 53 54 56 57{os.linesep}"


def test_ilrmsdmatrix_run_swappedchains(ilrmsdmatrix_module, mocker):
    """Test ilrmsdmatrix run method."""
    ilrmsdmatrix_module.previous_io = MockPreviousIO_protprot(
        path=ilrmsdmatrix_module.path
    )
    ilrmsdmatrix_module.params["receptor_chain"] = "B"
    ilrmsdmatrix_module.params["ligand_chains"] = ["A"]

    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )

    ilrmsdmatrix_module.run()

    ilrmsd_matrix_f = Path(ilrmsdmatrix_module.path, "ilrmsd.matrix")
    assert ilrmsd_matrix_f.exists()

    receptor_con_f = Path(ilrmsdmatrix_module.path, "receptor_contacts.con")
    assert receptor_con_f.exists()

    with open(ilrmsd_matrix_f, encoding="utf-8", mode="r") as f:
        assert f.readline() == f"1 2 15.166{os.linesep}"

    with open(
        receptor_con_f,
        encoding="utf-8",
        mode="r",
    ) as f:
        lines = f.readlines()
        assert lines[0] == f"B 10 11 12 16 17 48 51 52 53 54 56 57{os.linesep}"
        assert (
            lines[1]
            == f"A 37 38 39 40 43 44 45 69 71 72 75 90 93 94 96 132{os.linesep}"
        )  # noqa : E501
