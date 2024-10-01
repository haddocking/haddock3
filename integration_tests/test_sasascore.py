import tempfile
from pathlib import Path

import pytest
import shutil
import pandas as pd
import numpy as np

from haddock.modules.scoring.sasascore import DEFAULT_CONFIG as DEFAULT_SASASCORE_CONFIG
from haddock.modules.scoring.sasascore import HaddockModule as SasascoreModule
from haddock.libs.libontology import PDBFile
from . import CNS_EXEC, DATA_DIR, has_cns
from tests import golden_data

@pytest.fixture
def sasascore_module():
    """Return a default sasascore module."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        sasascore = SasascoreModule(
            order=0, path=tmpdir, initial_params=DEFAULT_SASASCORE_CONFIG
        )
        # let's assume we know Val 39 should be buried upon complex formation 
        sasascore.params["resdic_buried_A"] = [39]
        # let's assume we know GLU 43 should remain accessible
        sasascore.params["resdic_accessible_A"] = [43]
        yield sasascore
    
class MockPreviousIO():
    """Mock the previous_io class."""
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(Path(golden_data, "protprot_complex_1.pdb"), Path(".", "protprot_complex_1.pdb"))
        shutil.copy(Path(golden_data, "protprot_complex_2.pdb"), Path(".", "protprot_complex_2.pdb"))
        model_list = [
            PDBFile(file_name="protprot_complex_1.pdb", path="."),
            PDBFile(file_name="protprot_complex_2.pdb", path="."),
        ]

        return model_list

    def output(self):
        return None


def test_sasascore_default(sasascore_module, mocker):
    """Test the sasascore module."""
    sasascore_module.previous_io = MockPreviousIO(path=sasascore_module.path)
    sasascore_module.run()

    expected_sasascore_csv = Path(sasascore_module.path, "sasascore.tsv")
    expected_violations_csv = Path(sasascore_module.path, "violations.tsv")
    assert expected_sasascore_csv.exists(), f"{expected_sasascore_csv} does not exist"
    assert expected_violations_csv.exists(), f"{expected_violations_csv} does not exist"

    # check sasascore.tsv
    exp_shape = (2, 4)
    df = pd.read_csv(expected_sasascore_csv, sep="\t", comment="#")
    assert df.shape == exp_shape, f"{expected_sasascore_csv} has wrong shape ({df.shape} instead of {exp_shape})"
    assert df.loc[df["structure"] == "protprot_complex_2.pdb"].iloc[0,:]["score"] == 1
    assert df.loc[df["structure"] == "protprot_complex_1.pdb"].iloc[0,:]["score"] == 0

    # check violations.tsv
    exp_shape = (2, 3)
    df = pd.read_csv(expected_violations_csv, sep="\t", comment="#")
    print(df)
    assert df.shape == exp_shape, f"{expected_violations_csv} has wrong shape ({df.shape} instead of {exp_shape})"
    assert df.loc[df["structure"] == "protprot_complex_1.pdb"].iloc[0,:]["bur_A"] == "-"
    assert df.loc[df["structure"] == "protprot_complex_2.pdb"].iloc[0,:]["bur_A"] == "39"


def test_sasascore_no_residues(sasascore_module, mocker):
    "Test the sasascore module when a non-existing chain is added."
    sasascore_module.previous_io = MockPreviousIO(path=sasascore_module.path)
    # adding a non existing chain
    sasascore_module.params["resdic_buried_C"] = [1]
    sasascore_module.run()
    expected_sasascore_csv = Path(sasascore_module.path, "sasascore.tsv")
    exp_shape = (2, 4)
    df = pd.read_csv(expected_sasascore_csv, sep="\t", comment="#")
    assert df.shape == exp_shape, f"{expected_sasascore_csv} has wrong shape ({df.shape} instead of {exp_shape})"
    assert df.loc[df["structure"] == "protprot_complex_2.pdb"].iloc[0,:]["score"] == 1
    assert df.loc[df["structure"] == "protprot_complex_1.pdb"].iloc[0,:]["score"] == 0
    # now the violations df
    expected_violations_csv = Path(sasascore_module.path, "violations.tsv")
    exp_shape = (2, 4)
    df = pd.read_csv(expected_violations_csv, sep="\t", comment="#")
    assert df.shape == exp_shape, f"{expected_violations_csv} has wrong shape ({df.shape} instead of {exp_shape})"
    assert df.columns.tolist() == ["structure", "bur_A", "bur_C", "acc_A"]
    assert df.loc[df["structure"] == "protprot_complex_1.pdb"].iloc[0,:]["bur_A"] == "-"
    assert df.loc[df["structure"] == "protprot_complex_1.pdb"].iloc[0,:]["bur_C"] == "-"
