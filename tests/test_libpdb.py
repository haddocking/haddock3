"""Test lib PDB."""
from pathlib import Path
import pytest


from haddock.libs import libpdb
from haddock.libs.libio import PDBFile

from . import golden_data

chainC = [
    'ATOM      3  CA  ARG C   4      37.080  43.455  -3.421  1.00  0.00      C    C  ',  # noqa: E501
    'ATOM      3  CA  GLU C   6      33.861  45.127  -2.233  1.00  0.00      C    C  ',  # noqa: E501
    'ATOM      3  CA  ALA C   7      35.081  45.036   1.305  1.00  0.00      C    C  ',  # noqa: E501
    ]


@pytest.mark.parametrize(
    "lines,expected",
    (
        (chainC, ["C", "C", "C"]),
        ),
    )
def test_read_chain_ids(lines, expected):
    result = libpdb.read_chainids(lines)
    assert result == expected


@pytest.mark.parametrize(
    "lines,expected",
    (
        (chainC, ["C", "C", "C"]),
        ),
    )
def test_read_seg_ids(lines, expected):
    result = libpdb.read_segids(lines)
    assert result == expected


@pytest.fixture(name="wrong_rigid_molecules")
def fixture_wrong_rigidbody_molecules():
    """fixture for wrong rigidbody input molecules."""
    receptor = PDBFile(Path(golden_data, "protprot_complex_1.pdb"))
    ligand = PDBFile(Path(golden_data, "protprot_complex_2.pdb"))
    return [receptor, ligand]

@pytest.fixture(name="good_rigid_molecules")
def fixture_good_rigidbody_molecules():
    """fixture for good rigidbody input molecules."""
    receptor = PDBFile(Path(golden_data, "e2aP_1F3G_haddock.pdb"))
    ligand = PDBFile(Path(golden_data, "hpr_ensemble_1_haddock.pdb"))
    return [receptor, ligand]

def test_check_combination_chains(good_rigid_molecules, wrong_rigid_molecules):
    """Test check_combination_chains."""
    exp_chains = ["A", "B"]
    obs_chains = libpdb.check_combination_chains(good_rigid_molecules)
    assert obs_chains == exp_chains
    # when input molecules share chains there should be a ValueError
    with pytest.raises(ValueError):
        libpdb.check_combination_chains(wrong_rigid_molecules)
    
