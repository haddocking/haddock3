"""integration test for emscoring module"""

import shutil
import tempfile
from pathlib import Path

import pandas as pd
import pytest

from haddock.libs.libontology import PDBFile, TopologyFile
from haddock.modules.scoring.emscoring import \
    DEFAULT_CONFIG as DEFAULT_EMSCORING_CONFIG
from haddock.modules.scoring.emscoring import HaddockModule as EmscoringModule

from integration_tests import GOLDEN_DATA


@pytest.fixture
def emscoring_module():
    """Return a default emscoring module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        emscoring_module = EmscoringModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_EMSCORING_CONFIG
        )
        # lower number of steps for faster testing
        emscoring_module.params["nemsteps"] = 5
        emscoring_module.params["per_interface_scoring"] = True
        yield emscoring_module


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(
            Path(GOLDEN_DATA, "protglyc_complex_1.pdb"),
            Path(self.path, "protglyc_complex_1.pdb"),
        )

        # add the topology to the models
        psf_file = Path(GOLDEN_DATA, "protglyc_complex_1.psf")
        model_list = [
            PDBFile(
                file_name="protglyc_complex_1.pdb",
                path=self.path,
                topology=TopologyFile(psf_file),
            ),
        ]
        return model_list

    def output(self) -> None:
        return None


class MockPreviousIO_3chains:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        file_bn = "ab_ag_BHL"
        shutil.copy(
            Path(GOLDEN_DATA, f"{file_bn}.pdb"),
            Path(self.path, f"{file_bn}.pdb"),
        )

        # add the topology to the models
        psf_file = Path(GOLDEN_DATA, f"{file_bn}.psf")
        model_list = [
            PDBFile(
                file_name=f"{file_bn}.pdb",
                path=self.path,
                topology=TopologyFile(psf_file),
            ),
        ]
        return model_list

    def output(self) -> None:
        return None


def test_emscoring_default(emscoring_module, calc_fnat):
    """Test the emscoring module."""
    emscoring_module.previous_io = MockPreviousIO(path=emscoring_module.path)
    emscoring_module.run()

    expected_pdb1 = Path(emscoring_module.path, "emscoring_1.pdb")
    expected_csv = Path(emscoring_module.path, "emscoring.tsv")

    assert expected_pdb1.exists(), f"{expected_pdb1} does not exist"

    assert expected_csv.exists(), f"{expected_csv} does not exist"
    df = pd.read_csv(expected_csv, sep="\t", comment="#")
    assert df.columns.tolist() == ["structure", "original_name", "md5", "score"]
    assert df.shape == (1, 4)
    assert df["score"].dtype == float
    # the model should have a negative score (lower than 10)
    assert all(df["score"] < -10)
    # the model should not be too different from the original. Checking Fnat
    fnat = calc_fnat(
        model=Path(emscoring_module.path, "emscoring_1.pdb"),
        native=Path(GOLDEN_DATA, "protglyc_complex_1.pdb"),
    )
    assert fnat == pytest.approx(0.95, abs=0.1)

    # check the interface scoring
    expected_interface_csv = Path(emscoring_module.path, "emscoring_A_B.tsv")
    assert expected_interface_csv.exists(), f"{expected_interface_csv} does not exist"
    df_perint = pd.read_csv(expected_interface_csv, sep="\t", comment="#")
    # check that the score is equal to the global score (it's a dimer!)
    assert df_perint["score"].tolist() == df["score"].tolist()


def test_emscoring_3chains(emscoring_module):
    """Test the emscoring module with interface selection."""
    emscoring_module.previous_io = MockPreviousIO_3chains(path=emscoring_module.path)
    emscoring_module.run()

    expected_pdb1 = Path(emscoring_module.path, "emscoring_1.pdb")
    expected_csv_all = Path(emscoring_module.path, "emscoring.tsv")

    assert expected_pdb1.exists(), f"{expected_pdb1} does not exist"
    assert expected_csv_all.exists(), f"{expected_csv_all} does not exist"
    df_all = pd.read_csv(expected_csv_all, sep="\t", comment="#")

    # check the interface scoring
    expected_BH_interface_csv = Path(emscoring_module.path, "emscoring_B_H.tsv")
    assert expected_BH_interface_csv.exists(), f"{expected_BH_interface_csv} does not exist"
    expected_BL_interface_csv = Path(emscoring_module.path, "emscoring_B_L.tsv")
    assert expected_BL_interface_csv.exists(), f"{expected_BL_interface_csv} does not exist"

    # Set chain combination parameter
    emscoring_module.params["interface_combinations"] = ["B,H", "B,L"]
    emscoring_module.run()

    expected_pdb1 = Path(emscoring_module.path, "emscoring_1.pdb")
    expected_csv_combi = Path(emscoring_module.path, "emscoring.tsv")

    assert expected_pdb1.exists(), f"{expected_pdb1} does not exist"
    assert expected_csv_combi.exists(), f"{expected_csv_combi} does not exist"
    df_combi = pd.read_csv(expected_csv_combi, sep="\t", comment="#")

    # Ensure the score are different
    assert df_all["score"].tolist() != df_combi["score"].tolist()

    # Set chain combination parameter
    emscoring_module.params["interface_combinations"] = ["H,B", "L,B"]
    emscoring_module.run()

    expected_pdb2 = Path(emscoring_module.path, "emscoring_1.pdb")
    expected_csv_combi_revert = Path(emscoring_module.path, "emscoring.tsv")
    
    assert expected_pdb2.exists(), f"{expected_pdb2} does not exist"
    assert expected_csv_combi_revert.exists(), f"{expected_csv_combi_revert} does not exist"
    df_combi_revert = pd.read_csv(expected_csv_combi_revert, sep="\t", comment="#")

    # Ensure the score are the same
    combi_score = df_combi_revert["score"].tolist()[0]
    revert_score = df_combi["score"].tolist()[0]
    assert revert_score == pytest.approx(combi_score, abs=1)
