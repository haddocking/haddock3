"""Test the haddock3-restraints CLI."""

import logging
import tempfile
from pathlib import Path

import pytest

from haddock.clis.restraints.active_passive_to_ambig import (
    actpass_to_ambig,
    parse_actpass_file,
    )
from haddock.clis.restraints.calc_accessibility import (
    REL_ASA,
    calc_accessibility,
    )
from haddock.clis.restraints.passive_from_active import passive_from_active
from haddock.clis.restraints.restrain_bodies import restrain_bodies
from haddock.clis.restraints.validate_tbl import validate_tbl
from tests.test_module_caprieval import protdna_input_list  # noqa : F401

from . import golden_data


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


def test_parse_actpass_file(example_actpass_file):
    exp_active = [97, 98, 115, 116, 117]
    exp_passive = [91, 93, 95, 118, 119, 120]
    obs_active, obs_passive = parse_actpass_file(example_actpass_file)
    assert obs_active == exp_active
    assert obs_passive == exp_passive


def test_actpass_to_ambig(capsys):
    """Test actpass_to_ambig function."""
    # create temp file
    with tempfile.NamedTemporaryFile() as tmp:
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
    with tempfile.NamedTemporaryFile() as tmp:
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


def test_restrain_bodies(protdna_input_list, capsys):  # noqa : F811
    """Test restrain_bodies function."""
    restrain_bodies(protdna_input_list[0].rel_path)
    captured = capsys.readouterr()
    out_lines = captured.out.split("\n")
    assert (
        out_lines[0]
        == "assign (segid A and resi 10 and name CA) (segid B and resi 7 and name P) 26.542 0.0 0.0"
    )  # noqa : E501


def test_restrain_bodies_empty(example_pdb_file, capsys):
    """Test restrain_bodies function."""
    restrain_bodies(example_pdb_file)
    captured = capsys.readouterr()
    assert captured.out == ""


def test_restrain_bodies_exclude(protdna_input_list, capsys):  # noqa : F811
    """Test restrain_bodies function."""
    restrain_bodies(protdna_input_list[0].rel_path, exclude="A")
    captured = capsys.readouterr()
    out_lines = captured.out.split("\n")
    assert (
        out_lines[0]
        == "assign (segid B and resi 6 and name P) (segid B and resi 35 and name P) 15.187 0.0 0.0"
    )  # noqa : E501


def test_calc_accessibility_rel_asa_data():
    """Test content matching in REL_ASA."""
    all_entries = set([k for d in REL_ASA.values() for k in d.keys()])
    for dic in REL_ASA.values():
        assert all_entries == set(dic.keys())


def test_calc_accessibility(protdna_input_list, caplog):  # noqa : F811
    """Test calc_accessibility function."""
    caplog.set_level(logging.INFO)
    calc_accessibility(str(protdna_input_list[0].rel_path))
    assert caplog.messages[-2].endswith(
        "Chain A - -1,0,1,7,8,11,12,13,14,16,18,19,22,23,25,27,36,37,38,40,41,43,46,47,50,53,55,57,61"
    )  # noqa : E501
    assert caplog.messages[-1].endswith(
        "Chain B - 2,3,4,5,6,7,8,9,10,11,28,29,30,31,32,33,34,35,36,37,38"
    )  # noqa : E501
    # now let's change the cutoff
    caplog.clear()
    calc_accessibility(str(protdna_input_list[0].rel_path), cutoff=0.9)
    assert caplog.messages[-2].endswith("Chain A - 14,25,53")
    # DNA accessibility does not change with a higher cutoff!
    assert caplog.messages[-1].endswith(
        "Chain B - 2,3,4,5,6,7,8,9,10,11,28,29,30,31,32,33,34,35,36,37,38"
    )  # noqa : E501

