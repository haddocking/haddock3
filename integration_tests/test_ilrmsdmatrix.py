"""Intergration tests for the IL-RMSD matrix module."""

import os
import tempfile
from pathlib import Path
import shutil
import pytest
import pytest_mock  # noqa : F401

from haddock.libs.libontology import PDBFile

from haddock.modules.analysis.ilrmsdmatrix import (
    DEFAULT_CONFIG as DEFAULT_ILRMSD_CONFIG,
    HaddockModule as IlrmsdmatrixModule,
    )

from . import golden_data
from tests import golden_data as tests_golden_data


@pytest.fixture
def ilrmsdmatrix_module():
    """Provide a parametrized IL-RMSD matrix module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ilrmsdmatrix = IlrmsdmatrixModule(
            order=0, path=tmpdir, initial_params=DEFAULT_ILRMSD_CONFIG
            )
        yield ilrmsdmatrix


class MockPreviousIO():
    """Mock of the ModuleIO class."""

    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        """Mock of the io.retrive_models from ModuleIO."""
        shutil.copy(
            Path(golden_data, "protglyc_complex_1.pdb"),
            Path(".", "protglyc_complex_1.pdb"),
            )
        shutil.copy(
            Path(golden_data, "protglyc_complex_2.pdb"),
            Path(".", "protglyc_complex_2.pdb"),
            )
        
        model_list = [
            PDBFile(file_name="protglyc_complex_1.pdb", path="."),
            PDBFile(file_name="protglyc_complex_2.pdb", path="."),
            ]

        return model_list
    

class MockPreviousIO_protprot:
    """Mock previous IO class."""

    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        """Retrieve models."""
        shutil.copy(
            Path(tests_golden_data, "protprot_complex_1.pdb"),
            Path(".", "protprot_complex_1.pdb"),
            )
        shutil.copy(
            Path(tests_golden_data, "protprot_complex_2.pdb"),
            Path(".", "protprot_complex_2.pdb"),
            )
        model_list = [
            PDBFile(file_name="protprot_complex_1.pdb", path="."),
            PDBFile(file_name="protprot_complex_2.pdb", path="."),
            ]
        return model_list


def test_ilrmsdmatrix_default(ilrmsdmatrix_module, mocker):
    """Test the topoaa module."""
    ilrmsdmatrix_module.previous_io = MockPreviousIO(
        path=ilrmsdmatrix_module.path
        )
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
        )
    ilrmsdmatrix_module.run()
    # expected paths
    exp_ilrmsd_matrix = Path(ilrmsdmatrix_module.path, "ilrmsd.matrix")
    exp_contacts_file = Path(ilrmsdmatrix_module.path, "receptor_contacts.con")
    assert exp_ilrmsd_matrix.exists(), "ilrmsd.matrix does not exist"
    assert exp_contacts_file.exists(), "receptor_contacts.con does not exist"
    # open files and check content
    with open(exp_ilrmsd_matrix) as f:
        assert f.readline() == "1 2 11.877\n"
    with open(exp_contacts_file) as f:
        lines = f.readlines()
        assert lines[0] == f"A 35 44 46 48 50 52 56 57 58 59 61 62 63 73 75 98 101 102 103 104 106 107 108 109 110 112{os.linesep}"  # noqa : E501
        assert lines[1] == f"B 1 2 3 4{os.linesep}"


def test_ilrmsdmatrix_run(ilrmsdmatrix_module, mocker):
    """Test ilrmsdmatrix run method."""
    ilrmsdmatrix_module.previous_io = MockPreviousIO_protprot(path=ilrmsdmatrix_module.path)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
        )
    ilrmsdmatrix_module.run()
    assert Path(ilrmsdmatrix_module.path, "ilrmsd.matrix").exists()
    assert Path(ilrmsdmatrix_module.path, "receptor_contacts.con").exists()
    with open(Path(ilrmsdmatrix_module.path, "ilrmsd.matrix")) as f:
        assert f.readline() == f"1 2 16.715{os.linesep}"
    with open(Path(ilrmsdmatrix_module.path, "receptor_contacts.con")) as f:
        lines = f.readlines()
        assert lines[0] == f"A 37 38 39 40 43 44 45 69 71 72 75 90 93 94 96 132{os.linesep}"  # noqa : E501
        assert lines[1] == f"B 10 11 12 16 17 48 51 52 53 54 56 57{os.linesep}"


def test_ilrmsdmatrix_run_swappedchains(ilrmsdmatrix_module, mocker):
    """Test ilrmsdmatrix run method."""
    ilrmsdmatrix_module.previous_io = MockPreviousIO_protprot(path=ilrmsdmatrix_module.path)
    ilrmsdmatrix_module.params["receptor_chain"] = "B"
    ilrmsdmatrix_module.params["ligand_chains"] = ["A"]
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
        )
    ilrmsdmatrix_module.run()
    assert Path(ilrmsdmatrix_module.path, "ilrmsd.matrix").exists()
    assert Path(ilrmsdmatrix_module.path, "receptor_contacts.con").exists()
    with open(Path(ilrmsdmatrix_module.path, "ilrmsd.matrix")) as f:
        assert f.readline() == f"1 2 15.166{os.linesep}"
    with open(Path(ilrmsdmatrix_module.path, "receptor_contacts.con")) as f:
        lines = f.readlines()
        assert lines[0] == f"B 10 11 12 16 17 48 51 52 53 54 56 57{os.linesep}"
        assert lines[1] == f"A 37 38 39 40 43 44 45 69 71 72 75 90 93 94 96 132{os.linesep}"  # noqa : E501
