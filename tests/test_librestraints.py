"""Test the librestraints library."""
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.libs.librestraints import (
    get_unique_resids,
    parse_tbl,
    validate_ambig_fname,
    )

from . import golden_data


@pytest.fixture
def protprot_onechain_pdb():
    """Protein-Protein complex with a single chain ID."""
    return [
        PDBFile(Path(golden_data, "protprot_onechain.pdb"), path=golden_data)
        ]


@pytest.fixture
def protdna_pdb():
    """Protein-Protein dna complex."""
    return [
        PDBFile(Path(golden_data, "protdna_complex_2.pdb"), path=golden_data)
        ]


@pytest.fixture
def ligand_pdb():
    """Protein-Protein complex with a single chain ID."""
    return [PDBFile(Path(golden_data, "ligand.pdb"), path=golden_data)]


@pytest.fixture
def tbl_file():
    """Protein-Protein complex with a single chain ID."""
    return Path(golden_data, "example_ambig_1.tbl")


def test_parse_tbl(tbl_file):
    """Test the parse_tbl function."""
    obs_restraints = parse_tbl(tbl_file)
    exp_restraints = [
        (("A", 1), ("B", 1)),
        (("A", 1), ("B", 2)),
        (("A", 2), ("B", 1)),
        (("A", 2), ("B", 2)),
        (("B", 1), ("A", 1)),
        (("B", 1), ("A", 2)),
        (("B", 2), ("A", 1)),
        (("B", 2), ("A", 2)),
        ]
    assert obs_restraints == exp_restraints


def test_get_unique_resids(ligand_pdb):
    """Test the get_unique_resids function."""
    exp_resid_dict = {"B": [500]}
    obs_resid_dict = get_unique_resids(ligand_pdb)
    assert exp_resid_dict == obs_resid_dict


def test_validate_ambig_fname(protdna_pdb, tbl_file):
    """Test the validate_ambig_fname function."""
    assert validate_ambig_fname(tbl_file, protdna_pdb) is True


def test_validate_ambig_fname_error(protprot_onechain_pdb, tbl_file):
    """Test the validate_ambig_fname function in case of invalid restraints."""
    with pytest.raises(Exception):
        validate_ambig_fname(tbl_file, protprot_onechain_pdb)
