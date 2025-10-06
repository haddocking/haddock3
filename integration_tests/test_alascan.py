import tempfile
from pathlib import Path

import pytest
import shutil
import pandas as pd

from haddock.modules.analysis.alascan import (
    DEFAULT_CONFIG as DEFAULT_ALASCAN_CONFIG,
    HaddockModule as AlascanModule,
    )
from haddock.modules.analysis.alascan.scan import RES_CODES
from haddock.libs.libio import read_from_yaml
from haddock.libs.libontology import PDBFile
from . import GOLDEN_DATA


@pytest.fixture
def alascan_module():
    """Return a default alascan module."""
    with tempfile.TemporaryDirectory() as tmpdir:
        alascan = AlascanModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_ALASCAN_CONFIG
        )
        alascan.params["int_cutoff"] = 3.5
        alascan.params["output_mutants"] = True
        yield alascan


@pytest.fixture
def alascan_module_protlig(alascan_module):
    """Return a default alascan module."""
    # Copy parameters and toplogy of the ligand
    shutil.copy(
        Path(GOLDEN_DATA, "ligand.top"),
        Path(alascan_module.path, "ligand.top"),
        )
    shutil.copy(
        Path(GOLDEN_DATA, "ligand.param"),
        Path(alascan_module.path, "ligand.param"),
        )
    # Set the parameters to point the file
    alascan_module.params["ligand_param_fname"] = Path(alascan_module.path, "ligand.param")
    alascan_module.params["ligand_top_fname"] = Path(alascan_module.path, "ligand.top")
    yield alascan_module

class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(
            Path(GOLDEN_DATA, "protprot_complex_1.pdb"),
            Path(self.path, "protprot_complex_1.pdb"),
        )
        shutil.copy(
            Path(GOLDEN_DATA, "protprot_complex_2.pdb"),
            Path(self.path, "protprot_complex_2.pdb"),
        )
        model_list = [
            PDBFile(file_name="protprot_complex_1.pdb", path=self.path),
            PDBFile(file_name="protprot_complex_2.pdb", path=self.path),
        ]

        return model_list

    def output(self):
        return None
    

class MockPreviousIO_single_model:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(
            Path(GOLDEN_DATA, "2oob.pdb"),
            Path(self.path, "2oob.pdb"),
        )
        model_list = [
            PDBFile(file_name="2oob.pdb", path=self.path),
        ]

        return model_list

    def output(self):
        return None


class MockPreviousIO_protlig:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize: bool = False):
        shutil.copy(
            Path(GOLDEN_DATA, "protlig_complex_1.pdb"),
            Path(self.path, "protlig_complex_1.pdb"),
        )
        model_list = [
            PDBFile(file_name="protlig_complex_1.pdb", path=self.path),
        ]

        return model_list

    def output(self):
        return None

def test_alascan_default(alascan_module, mocker):
    """Test the alascan module."""
    alascan_module.previous_io = MockPreviousIO(path=alascan_module.path)
    alascan_module.run()

    expected_csv1 = Path(alascan_module.path, "scan_protprot_complex_1.tsv")
    expected_csv2 = Path(alascan_module.path, "scan_protprot_complex_2.tsv")
    expected_clt_csv = Path(alascan_module.path, "scan_clt_unclustered.tsv")

    assert expected_csv1.exists(), f"{expected_csv1} does not exist"
    assert expected_csv2.exists(), f"{expected_csv2} does not exist"
    assert expected_clt_csv.exists(), f"{expected_clt_csv} does not exist"

    # check single complex csv
    df = pd.read_csv(expected_csv1, sep="\t", comment="#")
    assert df.shape == (10, 15), f"{expected_csv1} has wrong shape"
    # ARG 17 B should have a negative delta_score
    assert df.loc[df["ori_resname"] == "ARG"].iloc[0, :]["delta_score"] < 0.0

    # check cluster csv
    df_clt = pd.read_csv(expected_clt_csv, sep="\t", comment="#")
    assert df_clt.shape == (18, 16), f"{expected_clt_csv} has wrong shape"
    # average delta score of A-38-ASP should be negative
    assert (
        df_clt.loc[df_clt["full_resname"] == "A-38-ASP"].iloc[0, :]["delta_score"] < 0.0
    )


def test_alascan_single_model(alascan_module, mocker):
    """Test the alascan module with only one model (saving mutants)."""
    # use lysine as the scan residue
    alascan_module.params["scan_residue"] = "LYS"
    alascan_module.previous_io = MockPreviousIO_single_model(path=alascan_module.path)
    alascan_module.run()

    expected_csv = Path(alascan_module.path, "scan_2oob.tsv")
    expected_clt_csv = Path(alascan_module.path, "scan_clt_unclustered.tsv")

    assert expected_csv.exists(), f"{expected_csv} does not exist"
    assert expected_clt_csv.exists(), f"{expected_clt_csv} does not exist"

    # check single complex csv
    df = pd.read_csv(expected_csv, sep="\t", comment="#")
    assert df.shape == (12, 15), f"{expected_csv} has wrong shape"
    
    # there should be several mutants saved to file
    # for each mutation in df, check that the corresponding file exists
    from haddock.libs.libalign import PROT_SIDE_CHAINS_DICT

    for _, row in df.iterrows():
        ch = row["chain"]
        resi = str(row["res"])
        mut_file_identifier = f"{ch}_{RES_CODES[row['ori_resname']]}{resi}K"
        mut_file = Path(alascan_module.path, f"2oob-{mut_file_identifier}.pdb")
        assert mut_file.exists(), f"{mut_file} does not exist"
        # now let's open the file and check that the mutation is correct
        
        heavy_atoms = []
        with open(mut_file, "r") as f:
            for ln in f:
                if ln.startswith("ATOM") and ln[21] == ch and ln[22:26].strip() == resi:
                    atom_name = ln[12:16].strip()
                    if not atom_name.startswith("H"):
                        heavy_atoms.append(atom_name)
        # heavy_atoms should be = PROT_SIDE_CHAINS_DICT["LYS"] (order may vary)
        assert set(heavy_atoms) == set(
            PROT_SIDE_CHAINS_DICT["LYS"]
        ), (
            f"Heavy atoms for {mut_file_identifier} are not correct: {heavy_atoms}"
        )


def test_alascan_mutation_resiudes():
    """Test making sure list of residues in defaults and script are same."""
    default_config = read_from_yaml(DEFAULT_ALASCAN_CONFIG)
    config_allowed_resiudes = set(default_config["scan_residue"]["choices"])
    script_allowed_resiudes = set(list(RES_CODES.keys()))
    assert config_allowed_resiudes == script_allowed_resiudes


def test_alascan_with_ligand_topar(alascan_module_protlig):
    """Test the use of alascan in presence of a ligand."""
    alascan_module_protlig.previous_io = MockPreviousIO_protlig(path=alascan_module_protlig.path)
    alascan_module_protlig.run()

    expected_csv = Path(alascan_module_protlig.path, "scan_protlig_complex_1.tsv")
    expected_clt_csv = Path(alascan_module_protlig.path, "scan_clt_unclustered.tsv")

    assert expected_csv.exists(), f"{expected_csv} does not exist"
    assert expected_clt_csv.exists(), f"{expected_clt_csv} does not exist"

    # List mutated files
    mutated_filepaths = list(Path(alascan_module_protlig.path).glob("protlig_complex_1-*.pdb"))
    assert len(mutated_filepaths) > 1

    # Loop over files
    for mutated_fpath in mutated_filepaths:
        # Make sure the ligand is in it
        file_content = mutated_fpath.read_text()
        assert file_content.count("G39") > 20


def test_alascan_without_ligand_topar(alascan_module):
    """Test the use of alascan in presence of a ligand."""
    alascan_module.previous_io = MockPreviousIO_protlig(path=alascan_module.path)
    alascan_module.run()

    expected_csv = Path(alascan_module.path, "scan_protlig_complex_1.tsv")
    expected_clt_csv = Path(alascan_module.path, "scan_clt_unclustered.tsv")

    assert not expected_csv.exists(), f"{expected_csv} does exist"
    assert not expected_clt_csv.exists(), f"{expected_clt_csv} does exist"

    # List mutated files
    mutated_filepaths = list(Path(alascan_module.path).glob("protlig_complex_1-*.pdb"))
    assert len(mutated_filepaths) == 0
