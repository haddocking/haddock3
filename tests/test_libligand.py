"""Tests for libligand."""

import pytest

from pathlib import Path
import tempfile
import shutil
from haddock.libs.libligand import identify_unknown_hetatms, run_prodrg
from haddock.core.supported_molecules import supported_HETATM
from . import golden_data as GOLDEN_DATA


@pytest.fixture
def ligand_pdb():
    src = Path(GOLDEN_DATA, "ligand.pdb")
    with tempfile.TemporaryDirectory() as tmpdir:
        dst = Path(tmpdir, "ligand.top")
        shutil.copy(src, dst)
        yield dst


@pytest.fixture
def protlig_complex_pdb():
    src = Path(GOLDEN_DATA, "protlig_complex_1.pdb")
    with tempfile.TemporaryDirectory() as tmpdir:
        dst = Path(tmpdir, "protlig_complex_1.pdb")
        shutil.copy(src, dst)
        yield dst


@pytest.fixture
def protein_pdb():
    src = Path(GOLDEN_DATA, "protein.pdb")
    with tempfile.TemporaryDirectory() as tmpdir:
        dst = Path(tmpdir, "protein.pdb")
        shutil.copy(src, dst)
        yield dst


def test_identify_unknown_hetatms_returns_unknown(ligand_pdb):
    """Unknown HETATM residue names are returned."""
    result = identify_unknown_hetatms(ligand_pdb)
    assert "G39" in result


def test_identify_unknown_hetatms_empty_for_protein(protein_pdb):
    """Standard protein PDB yields no unknown residues."""
    result = identify_unknown_hetatms(protein_pdb)
    assert result == []


def test_identify_unknown_hetatms_in_protlig_complex(protlig_complex_pdb):
    """Find HETATMS in a prot/ligand complex"""
    result = identify_unknown_hetatms(protlig_complex_pdb)
    assert "G39" in result


def test_run_prodrg(ligand_pdb):
    """run_prodrg returns CNS topology and parameter contents from a ligand PDB."""
    top, par = run_prodrg(ligand_pdb)
    assert "MASS" in top
    assert "BOND" in par
    assert "NBONds" not in par


def test_identify_unknown_hetatms_skips_known(tmp_path):
    """Residues already in supported_HETATM are not returned."""
    # Find a known HETATM
    known = next(iter(supported_HETATM))
    # Create a dummy pdb with this known HETATM
    content = f"HETATM    1  XX  {known} A 500       0.000   0.000   0.000  1.00  0.00          ZN\nEND\n"
    p = tmp_path / "known.pdb"
    p.write_text(content)
    result = identify_unknown_hetatms(p)
    assert known not in result
