import os
import pytest
import gzip
from unittest.mock import patch, mock_open, MagicMock
from haddock.libs import libnotebooks
from . import has_notebook


# Sample minimal PDB string for tests
@pytest.fixture
def sample_pdb():
    pdb_string = """ATOM      1  CA  MET A   1      11.104  13.207   2.527  1.00  0.00           CA 
ATOM      2  CA  MET A   2      16.104  13.207   2.527  1.00  0.00           CA 
ATOM      3  CA  MET A   3      21.104  13.207   2.527  1.00  0.00           CA 
END
"""
    return pdb_string


@pytest.fixture
def sample_pdb_as_file(tmp_path, sample_pdb):
    fname = tmp_path / "test.pdb"
    with open(fname, "w") as f:
        f.write(sample_pdb)

    yield fname


@pytest.fixture
def sample_pdb_as_gz(tmp_path, sample_pdb):
    fname = tmp_path / "test.pdb.gz"
    with gzip.open(fname, "wt") as f:
        f.write(sample_pdb)

    yield fname


def test_load_pdb_file_regular(sample_pdb_as_file, sample_pdb):
    result = libnotebooks.load_pdb_file(str(sample_pdb_as_file))
    assert result == sample_pdb


def test_load_pdb_file_gz(sample_pdb_as_gz, sample_pdb):
    result = libnotebooks.load_pdb_file(str(sample_pdb_as_gz))
    assert result == sample_pdb


def test_load_pdb_file_not_found():
    result = libnotebooks.load_pdb_file("non_existent_file.pdb")
    assert result is None


def test_pdb_string_to_structure_and_structure_to_pdb_string(sample_pdb):
    # Convert PDB string to structure and back
    structure = libnotebooks.pdb_string_to_structure(sample_pdb, "test")
    assert structure.id == "test"

    out_pdb = libnotebooks.structure_to_pdb_string(structure)
    assert "ATOM" in out_pdb and "END" in out_pdb


@has_notebook
def test_align_full_success(tmp_path, sample_pdb):
    # Patch py3Dmol and create two similar minimal pdb files
    pdb1 = tmp_path / "model1.pdb"
    pdb2 = tmp_path / "model2.pdb"
    pdb1.write_text(sample_pdb)
    pdb2.write_text(sample_pdb)

    view = libnotebooks.align_full(
        str(pdb1), str(pdb2), chains=["A"], atom_types=["CA"], show_labels=False
    )
    # Check this is a view
    # NOTE: Yes this import needs to be here and not on top
    import py3Dmol

    assert isinstance(view, py3Dmol.view)


@has_notebook
def test_align_full_file_not_found():

    result = libnotebooks.align_full("nofile1.pdb", "nofile2.pdb")
    assert result[1] is None

    # NOTE: Yes this import needs to be here and not on top
    import py3Dmol

    assert isinstance(result[0], py3Dmol.view)
