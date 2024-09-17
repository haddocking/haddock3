import os
import tempfile
from pathlib import Path
import shutil
import pytest
import pytest_mock

from haddock.libs.libontology import PDBFile

from haddock.modules.analysis.rmsdmatrix import DEFAULT_CONFIG as DEFAULT_RMSD_CONFIG
from haddock.modules.analysis.rmsdmatrix import HaddockModule as rmsdmatrixModule
from integration_tests import GOLDEN_DATA


@pytest.fixture
def rmsdmatrix_module():
    with tempfile.TemporaryDirectory() as tmpdir:
        ilrmsdmatrix = rmsdmatrixModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_RMSD_CONFIG
        )
        yield ilrmsdmatrix


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
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


def test_rmsdmatrix_default(rmsdmatrix_module, mocker):
    """Test the rmsdmatrix module."""
    rmsdmatrix_module.previous_io = MockPreviousIO(path=rmsdmatrix_module.path)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models", return_value=None
    )
    rmsdmatrix_module.run()
    # expected paths
    exp_rmsd_matrix = Path(rmsdmatrix_module.path, "rmsd.matrix")
    assert exp_rmsd_matrix.exists(), "rmsd.matrix does not exist"
    # open files and check content
    with open(exp_rmsd_matrix) as f:
        assert f.readline() == "1 2 3.326\n"
