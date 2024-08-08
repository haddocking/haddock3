import tempfile
from pathlib import Path

import pytest
import shutil
import pandas as pd

from haddock.modules.scoring.mdscoring import HaddockModule as mdscoringModule
from haddock.modules.scoring.mdscoring import DEFAULT_CONFIG as DEFAULT_MDSCORING_CONFIG
from haddock.libs.libontology import PDBFile, TopologyFile

from . import CNS_EXEC, DATA_DIR, has_cns
from . import golden_data


@pytest.fixture
def mdscoring_module():
    """Return a default mdscoring module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        mdscoring_module = mdscoringModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_MDSCORING_CONFIG
        )
        # lower number of steps for faster testing
        mdscoring_module.params["nemsteps"] = 25
        mdscoring_module.params["watersteps"] = 25
        mdscoring_module.params["watercoolsteps"] = 25
        mdscoring_module.params["waterheatsteps"] = 25
        yield mdscoring_module


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(Path(golden_data, "prot.pdb"), Path(".", "prot.pdb"))
        
        # add the topology to the models
        model_list = [
            PDBFile(file_name="prot.pdb", path=".", topology=TopologyFile(Path(golden_data, "prot.psf"))),
        ]
        return model_list

    def output(self) -> None:
        return None


@has_cns
def test_mdscoring_default(mdscoring_module):
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
    assert all(df["score"] < -500)
