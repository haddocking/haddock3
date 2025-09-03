"""Test the libnotebooks library."""
from unittest.mock import patch, MagicMock
from haddock.libs import libnotebooks

# Sample minimal PDB string for tests
SAMPLE_PDB = """ATOM      1  CA  MET A   1      11.104  13.207   2.527  1.00  0.00           CA 
ATOM      2  CA  MET A   2      16.104  13.207   2.527  1.00  0.00           CA 
ATOM      3  CA  MET A   3      21.104  13.207   2.527  1.00  0.00           CA 
END
"""

def test_load_pdb_file_regular(tmp_path):
    # Create a test pdb file
    filepath = tmp_path / "test.pdb"
    filepath.write_text(SAMPLE_PDB)
    result = libnotebooks.load_pdb_file(str(filepath))
    assert result == SAMPLE_PDB

def test_load_pdb_file_gz(tmp_path):
    import gzip
    filepath = tmp_path / "test.pdb.gz"
    with gzip.open(filepath, 'wt') as f:
        f.write(SAMPLE_PDB)
    result = libnotebooks.load_pdb_file(str(filepath))
    assert result == SAMPLE_PDB

def test_load_pdb_file_not_found():
    result = libnotebooks.load_pdb_file("non_existent_file.pdb")
    assert result is None

def test_pdb_string_to_structure_and_structure_to_pdb_string():
    # Convert PDB string to structure and back
    structure = libnotebooks.pdb_string_to_structure(SAMPLE_PDB, "test")
    assert structure.id == "test"
    out_pdb = libnotebooks.structure_to_pdb_string(structure)
    assert "ATOM" in out_pdb and "END" in out_pdb

@patch("haddock.libs.libnotebooks.py3Dmol")
def test_align_full_success(mock_py3Dmol, tmp_path):
    # Patch py3Dmol and create two similar minimal pdb files
    pdb1 = tmp_path / "model1.pdb"
    pdb2 = tmp_path / "model2.pdb"
    pdb1.write_text(SAMPLE_PDB)
    pdb2.write_text(SAMPLE_PDB)
    mock_view = MagicMock()
    mock_py3Dmol.view.return_value = mock_view

    view = libnotebooks.align_full(str(pdb1), str(pdb2), chains=['A'], atom_types=['CA'], show_labels=False)
    # Should return the mocked view
    assert view == mock_view

@patch("haddock.libs.libnotebooks.py3Dmol")
def test_align_full_file_not_found(mock_py3Dmol):
    mock_view = MagicMock()
    mock_py3Dmol.view.return_value = mock_view

    result = libnotebooks.align_full("nofile1.pdb", "nofile2.pdb")
    assert result[0] == mock_view
    assert result[1] is None

