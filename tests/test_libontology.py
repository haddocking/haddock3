"""Test functions and methods in haddock.libs.libontology."""
import pytest
from pathlib import Path

from haddock.libs.libontology import (
    Molecule,
    PDBFile,
    )

from . import golden_data


@pytest.fixture
def molecule():
    return Molecule(None)


@pytest.fixture
def protein():
    return Path(golden_data, "protein.pdb")


@pytest.fixture
def ensemble_header_w_md5():
    return Path(golden_data, "ens_header.pdb")


def test_get_md5(molecule, ensemble_header_w_md5, protein):
    """Test get_md5 method."""
    observed_md5_dic = molecule.get_md5(ensemble_header_w_md5)
    expected_md5_dic = {
        1: '71098743056e0b95fbfafff690703761',
        2: 'f7ab0b7c751adf44de0f25f53cfee50b',
        3: '41e028d8d28b8d97148dc5e548672142',
        4: '761cb5da81d83971c2aae2f0b857ca1e',
        5: '6c438f941cec7c6dc092c8e48e5b1c10',
        }

    assert observed_md5_dic == expected_md5_dic
    observed_md5_dic = molecule.get_md5(protein)
    assert observed_md5_dic == {}


def test_get_ensemble_origin(molecule, ensemble_header_w_md5, protein):
    """Test get_ensemble_origin method."""
    expected_origin_dic = {
        1: 'T161-hybrid-fit-C2-NCS_complex_100w',
        2: 'T161-hybrid-fit-C2-NCS_complex_101w',
        3: 'T161-hybrid-fit-C2-NCS_complex_102w',
        4: 'T161-hybrid-fit-C2-NCS_complex_103w',
        5: 'T161-hybrid-fit-C2-NCS_complex_104w',
        }
    observed_origin = molecule.get_ensemble_origin(ensemble_header_w_md5)
    assert observed_origin == expected_origin_dic
    observed_origin = molecule.get_ensemble_origin(protein)
    assert observed_origin == {}


def test_load_single_pdb(molecule, protein):
    """Test casting into PDBFile."""
    # Re-initialize with a actual protein
    molecule.__init__(protein)
    assert isinstance(molecule.pdb_files, dict)
    for pdbfile in molecule.pdb_files.values():
        assert isinstance(pdbfile, PDBFile)
    assert len(molecule) == 1


def test_load_single_pdb(molecule, protein):
    """Test casting into PDBFile."""
    # Re-initialize with a actual protein
    molecule.__init__(protein)
    assert isinstance(molecule.pdb_files, dict)
    for pdbfile in molecule.pdb_files.values():
        assert isinstance(pdbfile, PDBFile)
    assert len(molecule) == 1
