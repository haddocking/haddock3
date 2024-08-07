import tempfile
from pathlib import Path

import pytest
import shutil
import pandas as pd

from haddock.modules.scoring.emscoring import HaddockModule as EmscoringModule
from haddock.modules.scoring.emscoring import DEFAULT_CONFIG as DEFAULT_EMSCORING_CONFIG
from haddock.libs.libontology import PDBFile, TopologyFile
from integration_tests.test_alascan import MockPreviousIO

from . import CNS_EXEC, DATA_DIR, has_cns
from . import golden_data


@pytest.fixture
def emscoring_module():
    """Return a default emscoring module."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        emscoring_module = EmscoringModule(
            order=0, path=".", initial_params=DEFAULT_EMSCORING_CONFIG
        )
        # lower number of steps for faster testing
        emscoring_module.params["nemsteps"] = 5
        yield emscoring_module


class MockPreviousIO():
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(Path(golden_data, "prot.pdb"), Path(".", "prot.pdb"))
        
        # add the topology to the models
        model_list = [
            PDBFile(file_name="prot.pdb", path=".", topology=TopologyFile(Path(golden_data, "prot.psf"))),
        ]
        return model_list

    def output(self):
        return None


@has_cns
def test_emscoring_default(emscoring_module):
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
    # the model should have highly negative score
    assert all(df["score"] < -500)
