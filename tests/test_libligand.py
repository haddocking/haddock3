"""Tests for libligand."""

import pytest

from pathlib import Path
import tempfile
import shutil
from haddock.libs.libligand import (
    identify_unknown_hetatms,
    run_prodrg,
    _remove_nbonds,
    _sanitize_atom_names,
)
from haddock.core.supported_molecules import supported_HETATM
from . import golden_data as GOLDEN_DATA
from . import is_linux_x86_64


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


@is_linux_x86_64
def test_run_prodrg(ligand_pdb, tmp_path):
    """run_prodrg writes named .top/.param files and returns their paths."""
    top, par = run_prodrg(ligand_pdb, tmp_path)
    assert top.exists()
    assert par.exists()
    assert top.suffix == ".top"
    assert par.suffix == ".param"
    assert "MASS" in top.read_text()
    assert "BOND" in par.read_text()
    assert "NBONds" not in par.read_text()


@is_linux_x86_64
def test_run_prodrg_fails_on_complex(protlig_complex_pdb, tmp_path):
    """run_prodrg raises RuntimeError when given a protein-ligand complex.

    prodrg only accepts small molecule inputs; passing a full complex
    causes it to abort with 'Too many lines of text in input'.
    """
    with pytest.raises(RuntimeError):
        run_prodrg(protlig_complex_pdb, tmp_path)


def test_remove_nbonds_removes_block():
    """NBONds...END block is stripped from parameter content."""
    content = "BOND CP1 OP1 100.0\nNBONds\n  tolerance 0.5\nEND\nANGLE CP1 CP2\n"
    result = _remove_nbonds(content)
    assert "NBONds" not in result
    assert "BOND CP1 OP1 100.0" in result
    assert "ANGLE CP1 CP2" in result


def test_remove_nbonds_no_block():
    """Content without an NBONds block is returned unchanged."""
    content = "BOND CP1 OP1 100.0\nANGLE CP1 CP2\n"
    assert _remove_nbonds(content) == content


def test_sanitize_atom_names_removes_colons():
    """Colons in atom type names are removed."""
    content = "MASS HT:A   1.0080\n  ATOM H:A  TYPE=HT:A CHARGE=-0.019 END\n"
    result = _sanitize_atom_names(content)
    assert ":" not in result
    assert "MASS HTA" in result
    assert "ATOM HA  TYPE=HTA" in result


def test_sanitize_atom_names_preserves_comment_colons():
    """Colons inside comment lines are not removed."""
    content = "! cite: Author et al.\nMASS HT:A   1.0080\n"
    result = _sanitize_atom_names(content)
    assert "cite: Author et al." in result
    assert "HTA" in result


def test_sanitize_atom_names_no_colons():
    """Content without colons is returned unchanged."""
    content = "MASS HTA   1.0080\n  ATOM HA  TYPE=HTA CHARGE=-0.019 END\n"
    assert _sanitize_atom_names(content) == content


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
