"""Test the haddock3-restraints CLI."""
from pathlib import Path
import pytest
import tempfile

from haddock.restraints.active_passive_to_ambig import actpass_to_ambig, parse_actpass_file
from haddock.restraints.validate_tbl import validate_tbl
from haddock.restraints.passive_from_active import passive_from_active
from haddock.restraints.restrain_bodies import restrain_bodies

from . import golden_data
from tests.test_module_caprieval import protdna_input_list


@pytest.fixture
def example_actpass_file():
    """Provide example actpass filename."""
    return Path(golden_data, "example.act-pass")


@pytest.fixture
def example_tbl_file():
    """Provide example tbl filename."""
    return Path(golden_data, "example_ambig_1.tbl")


@pytest.fixture
def example_pdb_file():
    """Provide example pdb filename."""
    return Path(golden_data, "protein.pdb")


@pytest.fixture
def example_protdna_file():
    """Provide example pdb filename."""
    return Path(golden_data, "protein.pdb")


def test_parse_actpass_file(example_actpass_file):
    exp_active = [97, 98, 115, 116, 117]
    exp_passive = [91, 93, 95, 118, 119, 120]
    obs_active, obs_passive = parse_actpass_file(example_actpass_file)
    assert obs_active == exp_active
    assert obs_passive == exp_passive


def test_actpass_to_ambig(capsys):
    """Test actpass_to_ambig function."""
    # create temp file
    with tempfile.NamedTemporaryFile(dir=".") as tmp:
        # write something to it
        tmp.write(b"1\n2")
        # close it
        tmp.flush()
        
        # capture stdout and stderr
        actpass_to_ambig(tmp.name, tmp.name, "A", "B")
        captured = capsys.readouterr()
        assert captured.err == ""
        out_lines = captured.out.split("\n")
        assert out_lines[0] == "assign (resi 1 and segid A)"
        assert out_lines[2] == "       (resi 1 and segid B)"
        assert out_lines[5] == ") 2.0 2.0 0.0"
        assert out_lines[7] == "assign (resi 1 and segid B)"


def test_validate_tbl(example_tbl_file, capsys):
    """Test validate_tbl function."""
    validate_tbl(example_tbl_file, pcs=False, quick=True, silent=True)
    captured = capsys.readouterr()
    assert captured.err == ""
    assert captured.out == ""
    # now let's test the silent=False option
    validate_tbl(example_tbl_file, pcs=False, quick=True, silent=False)
    captured = capsys.readouterr()
    assert captured.err == ""
    assert captured.out.startswith("!\nassign")


def test_validate_tbl_error(example_tbl_file, capsys):
    """Test validate_tbl function in case of malformed tbl."""
    lines = open(example_tbl_file, "r").readlines()
    with tempfile.NamedTemporaryFile(dir=".") as tmp:
        # let's say I forget some lines
        for ln in lines[3:]:
            tmp.write(ln.encode())
        tmp.flush()
        with pytest.raises(Exception, match=r"Invalid TBL file: Unknown statement *"):
            validate_tbl(tmp.name, pcs=False, quick=True, silent=True)


def test_passive_from_active(example_pdb_file, capsys):
    """Test passive_from_active function."""
    active_residues = "1"
    passive_from_active(example_pdb_file, active_residues)
    captured = capsys.readouterr()
    assert captured.out == "2 3\n"


def test_restrain_bodies(protdna_input_list, capsys):
    """Test restrain_bodies function."""
    restrain_bodies(protdna_input_list[0].rel_path)
    captured = capsys.readouterr()
    out_lines = captured.out.split("\n")
    assert out_lines[0] == "assign (segid A and resi 10 and name CA) (segid B and resi 7 and name P) 26.542 0.0 0.0"


def test_restrain_bodies_empty(example_pdb_file, capsys):
    """Test restrain_bodies function."""
    restrain_bodies(example_pdb_file)
    captured = capsys.readouterr()
    assert captured.out == ""


def test_restrain_bodies_exclude(protdna_input_list, capsys):
    """Test restrain_bodies function."""
    restrain_bodies(protdna_input_list[0].rel_path, exclude="A")
    captured = capsys.readouterr()
    out_lines = captured.out.split("\n")
    assert out_lines[0] == "assign (segid B and resi 3 and name P) (segid B and resi 29 and name P) 31.170 0.0 0.0"
