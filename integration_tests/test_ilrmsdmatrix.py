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


@pytest.fixture
def ilrmsdmatrix_module():
    """Provide a parametrized IL-RMSD matrix module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        ilrmsdmatrix = IlrmsdmatrixModule(
            order=1, path=Path(tmpdir), initial_params=DEFAULT_ILRMSD_CONFIG
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
