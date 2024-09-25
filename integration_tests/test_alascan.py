import tempfile
from pathlib import Path

import pytest
import shutil
import pandas as pd

from haddock.modules.analysis.alascan import DEFAULT_CONFIG as DEFAULT_ALASCAN_CONFIG
from haddock.modules.analysis.alascan import HaddockModule as AlascanModule
from haddock.libs.libontology import PDBFile
from tests import golden_data


@pytest.fixture
def alascan_module():
    """Return a default alascan module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        alascan = AlascanModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_ALASCAN_CONFIG
        )
        alascan.params["int_cutoff"] = 3.5
        yield alascan


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(
            Path(golden_data, "protprot_complex_1.pdb"),
            Path(self.path, "protprot_complex_1.pdb"),
        )
        shutil.copy(
            Path(golden_data, "protprot_complex_2.pdb"),
            Path(self.path, "protprot_complex_2.pdb"),
        )
        model_list = [
            PDBFile(file_name="protprot_complex_1.pdb", path=self.path),
            PDBFile(file_name="protprot_complex_2.pdb", path=self.path),
        ]

        return model_list

    def output(self):
        return None


def test_alascan_default(alascan_module, mocker):
    """Test the alascan module."""
    alascan_module.previous_io = MockPreviousIO(path=alascan_module.path)
    alascan_module.run()

    expected_csv1 = Path(alascan_module.path, "scan_protprot_complex_1.csv")
    expected_csv2 = Path(alascan_module.path, "scan_protprot_complex_2.csv")
    expected_clt_csv = Path(alascan_module.path, "scan_clt_-.csv")

    assert expected_csv1.exists(), f"{expected_csv1} does not exist"
    assert expected_csv2.exists(), f"{expected_csv2} does not exist"
    assert expected_clt_csv.exists(), f"{expected_clt_csv} does not exist"

    # check single complex csv
    df = pd.read_csv(expected_csv1, sep="\t", comment="#")
    assert df.shape == (10, 15), f"{expected_csv1} has wrong shape"
    # ARG 17 B should have a negative delta_score
    assert df.loc[df["ori_resname"] == "ARG"].iloc[0, :]["delta_score"] < 0.0

    # check clt csv
    df_clt = pd.read_csv(expected_clt_csv, sep="\t", comment="#")
    assert df_clt.shape == (18, 11), f"{expected_clt_csv} has wrong shape"
    # average delta score of A-38-ASP should be negative
    assert (
        df_clt.loc[df_clt["full_resname"] == "A-38-ASP"].iloc[0, :]["delta_score"] < 0.0
    )
