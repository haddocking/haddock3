import tempfile
from pathlib import Path

import pytest
import shutil

from haddock.modules.analysis.alascan import DEFAULT_CONFIG as DEFAULT_ALASCAN_CONFIG
from haddock.modules.analysis.alascan import HaddockModule as AlascanModule
from haddock.libs.libontology import PDBFile
from . import CNS_EXEC, DATA_DIR, has_cns

@pytest.fixture
def alascan_module():
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        alascan = AlascanModule(
            order=0, path=tmpdir, initial_params=DEFAULT_ALASCAN_CONFIG
        )
        alascan.params["int_cutoff"] = 3.5
        alascan.params["output"] = True
        yield alascan
    
class MockPreviousIO:
    def retrieve_models(self, individualize: bool = False):
        shutil.copy(Path("..", "tests", "golden_data", "protprot_complex_1.pdb"), Path(".", "protprot_complex_1.pdb"))
        shutil.copy(Path("..", "tests", "golden_data", "protprot_complex_2.pdb"), Path(".", "protprot_complex_2.pdb"))
        model_list = [
            PDBFile(file_name="protprot_complex_1.pdb", path="."),
            PDBFile(file_name="protprot_complex_2.pdb", path="."),
        ]
        return model_list

@has_cns
def test_alascan_default(alascan_module, mocker):
    """Test the topoaa module."""
    alascan_module.previous_io = MockPreviousIO()
    mocker.patch("haddock.modules.BaseHaddockModule.export_io_models", return_value = None)
    alascan_module.run()

    expected_csv1 = Path(alascan_module.path, "scan_protprot_complex_1.csv")
    expected_csv2 = Path(alascan_module.path, "scan_protprot_complex_2.csv")
    expected_clt_csv = Path(alascan_module.path, "scan_clt_-.csv")

    assert expected_csv1.exists(), f"{expected_csv1} does not exist"
    assert expected_csv2.exists(), f"{expected_csv2} does not exist"
    assert expected_clt_csv.exists(), f"{expected_clt_csv} does not exist"
    