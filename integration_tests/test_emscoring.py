"""integration test for emscoring module"""
import tempfile
from pathlib import Path

import pytest
import shutil
import pandas as pd

from haddock.modules.scoring.emscoring import HaddockModule as EmscoringModule
from haddock.modules.scoring.emscoring import (
    DEFAULT_CONFIG as DEFAULT_EMSCORING_CONFIG)
from haddock.libs.libontology import PDBFile, TopologyFile
from haddock.modules.analysis.caprieval.capri import load_contacts

from . import has_cns
from . import golden_data


@pytest.fixture
def emscoring_module():
    """Return a default emscoring module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        emscoring_module = EmscoringModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_EMSCORING_CONFIG
        )
        # lower number of steps for faster testing
        emscoring_module.params["nemsteps"] = 5
        yield emscoring_module


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(
            Path(golden_data, "protglyc_complex_1.pdb"),
            Path(".", "protglyc_complex_1.pdb")
            )
        
        # add the topology to the models
        psf_file = Path(golden_data, "protglyc_complex_1.psf")
        model_list = [
            PDBFile(
                file_name="protglyc_complex_1.pdb",
                path=".",
                topology=TopologyFile(psf_file),
                ),
        ]
        return model_list

    def output(self) -> None:
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
    # the model should have a negative score (lower than 10)
    assert all(df["score"] < -10)
    # the model should not be too different from the original.
    # we calculated fnat with respect to the new structure, that could have new
    # weird contacts.
    start_pdb = PDBFile(Path(emscoring_module.path, "protglyc_complex_1.pdb"))
    start_contact = load_contacts(start_pdb)
    end_contact = load_contacts(expected_pdb1)
    intersection = start_contact & end_contact
    fnat = len(intersection) / float(len(end_contact))
    assert fnat > 0.95
