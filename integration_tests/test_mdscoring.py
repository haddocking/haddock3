"""integration test for mdscoring module."""

import shutil
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from haddock.libs.libontology import PDBFile, TopologyFile
from haddock.modules.scoring.mdscoring import \
    DEFAULT_CONFIG as DEFAULT_MDSCORING_CONFIG
from haddock.modules.scoring.mdscoring import HaddockModule as mdscoringModule

from integration_tests import GOLDEN_DATA


@pytest.fixture
def mdscoring_module():
    """Return a default mdscoring module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mdscoring_module = mdscoringModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_MDSCORING_CONFIG
        )
        # lower number of steps for faster testing
        mdscoring_module.params["watersteps"] = 200
        mdscoring_module.params["watercoolsteps"] = 200
        # enable per interface scoring
        mdscoring_module.params["per_interface_scoring"] = True
        yield mdscoring_module


class MockPreviousIO:
    """Mock the previous IO module."""

    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(
            Path(GOLDEN_DATA, "protglyc_complex_1.pdb"),
            Path(self.path, "protglyc_complex_1.pdb"),
        )

        shutil.copy(
            Path(GOLDEN_DATA, "protglyc_complex_1.psf"),
            Path(self.path, "protglyc_complex_1.psf"),
        )

        # add the topology to the models
        model_list = [
            PDBFile(
                file_name="protglyc_complex_1.pdb",
                path=self.path,
                topology=TopologyFile(
                    path=self.path,
                    file_name="protglyc_complex_1.psf",
                ),
            ),
        ]
        return model_list

    def output(self) -> None:
        return None


def test_mdscoring_default(mdscoring_module, calc_fnat):
    """Test the mdscoring module."""
    mdscoring_module.previous_io = MockPreviousIO(path=mdscoring_module.path)
    mdscoring_module.run()

    expected_pdb1 = Path(mdscoring_module.path, "mdscoring_1.pdb")
    expected_csv = Path(mdscoring_module.path, "mdscoring.tsv")

    assert expected_pdb1.exists(), f"{expected_pdb1} does not exist"

    assert expected_csv.exists(), f"{expected_csv} does not exist"
    df = pd.read_csv(expected_csv, sep="\t", comment="#")
    assert df.columns.tolist() == ["structure", "original_name", "md5", "score"]
    assert df.shape == (1, 4)
    assert df["score"].dtype == float
    # the model should have highly negative score
    assert all(df["score"] < -10)
    # the model should have a Fnat close to the starting structure
    fnat = calc_fnat(
        model=Path(mdscoring_module.path, "mdscoring_1.pdb"),
        native=Path(GOLDEN_DATA, "protglyc_complex_1.pdb"),
    )
    assert fnat == pytest.approx(0.90, abs=0.1)

    # check the interface scoring
    expected_interface_csv = Path(mdscoring_module.path, "mdscoring_A_B.tsv")
    assert expected_interface_csv.exists(), f"{expected_interface_csv} does not exist"
    df_perint = pd.read_csv(expected_interface_csv, sep="\t", comment="#")
    # check that the score is equal to the global score (it's a dimer!)
    for perint, sc in zip(df_perint["score"].tolist(), df["score"].tolist()):
        assert perint == pytest.approx(sc, abs=0.01)
