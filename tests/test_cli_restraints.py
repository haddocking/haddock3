"""Test the haddock3-restraints CLI."""

import logging
import tempfile
from pathlib import Path
import tarfile
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
from haddock.clis.restraints.random_removal import random_removal
from haddock.clis.restraints.restrain_ligand import restrain_ligand
from haddock.clis.restraints.z_surface_restraints import (
    compute_barycenter,
    get_z_coords,
    load_selected_resiudes_coords,
    load_selections,
    shape_bead,
    step_coords,
    )
from haddock.libs.libpdb import (
    slc_chainid,
    slc_name,
    slc_resname,
    slc_resseq,
    slc_serial,
    slc_temp,
    slc_x,
    slc_y,
    slc_z,
    )
from haddock.libs.librestraints import extract_restraint_entries

from . import golden_data


@pytest.fixture
def example_actpass_file():
    """Provide example actpass filename."""
    return Path(golden_data, "example.act-pass")


@pytest.fixture(name="example_tbl_file")
def fixture_example_tbl_file():
    """Provide example tbl filename."""
    return Path(golden_data, "example_ambig_1.tbl")


@pytest.fixture
def example_pdb_file():
    """Provide example pdb filename."""
    return Path(golden_data, "protein.pdb")


@pytest.fixture
def example_liganded_pdbfile():
    """Provide example pdb file containing a ligand."""
    return Path(golden_data, "protlig_complex_1.pdb")


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


def test_actpass_to_ambig_no_active(capsys):
    """Test actpass_to_ambig function with no active residues."""
    with tempfile.NamedTemporaryFile() as tmp:
        tmp.write(b"\n2")
        tmp.flush()
        # running actpass to ambig should raise a ValueError
        with pytest.raises(ValueError, match=r"No active residues found"):
            actpass_to_ambig(tmp.name, tmp.name, "A", "B")


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

@pytest.mark.skip
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

@pytest.mark.skip
def test_restrain_bodies_exclude(protdna_input_list, capsys):  # noqa : F811
    """Test restrain_bodies function."""
    restrain_bodies(protdna_input_list[0].rel_path, exclude="A")
    captured = capsys.readouterr()
    out_lines = captured.out.split("\n")
    assert (
        out_lines[0]
        == "assign (segid B and resi 3 and name P) (segid B and resi 29 and name P) 31.170 0.0 0.0"
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


def test_step_coords():
    """Test small size versus spacing."""
    coords = [coord for coord in step_coords(20, 2)]
    ref_coords = [coord for coord in range(-10, 12, 2)]
    assert coords == ref_coords


def test_step_coords_wrong_spacing():
    """Test small size versus spacing."""
    with pytest.raises(ValueError):
        step_coords(2, 20)


def test_shape_bead_valid_pdb_format():
    """Test that the generated shape beads are well formated."""
    bead = shape_bead(1, 1, 1, 1, chain="Z", atindex=1234, bfactor=100)
    assert bead[slc_serial].strip() == "1234"
    assert bead[slc_name].strip() == "SHA"
    assert bead[slc_chainid].strip() == "Z"
    assert bead[slc_resseq].strip() == "1"
    assert bead[slc_resname].strip() == "SHA"
    assert bead[slc_x].strip() == "1.000"
    assert bead[slc_y].strip() == "1.000"
    assert bead[slc_z].strip() == "1.000"
    assert bead[slc_temp].strip() == "100.00"


def test_load_selections():
    """Test residue selection loading."""
    loaded = load_selections(["1,2,3", "4,5,6"])
    assert list(loaded.values()) == [[1, 2, 3], [4, 5, 6]]
    loaded = load_selections(["1,2,3", "4,5,a"])
    assert list(loaded.values()) == [[1, 2, 3], [4, 5]]


def test_compute_barycenter():
    """Test center of mass computation."""
    barycenter = compute_barycenter([[0, 0, 0], [1, 1, 1]])
    assert barycenter == (0.5, 0.5, 0.5)
    barycenter = compute_barycenter([[0, 0, 0], [1, 1, 1], [2, 2, 2]])
    assert barycenter == (1, 1, 1)


def test_load_selected_resiudes_coords(example_pdb_file):
    """Test loading of selected resiudes from pdb file."""
    select_coords, _chainids, _atnames = load_selected_resiudes_coords(
        example_pdb_file,
        {
            "select_1": [1],
            "select_2": [2, 3],
        },
    )
    assert select_coords["select_1"][0] == (3.439, 7.910, -11.913)
    assert select_coords["select_2"][0] == (1.091, 9.158, -9.167)
    assert select_coords["select_2"][1] == (2.045, 9.187, -5.471)


def test_get_z_coords():
    """Test return z positions for two selections."""
    z_coords = get_z_coords(
        {
            "select_1": [(0, 0, 0)],
            "select_2": [(0, 0, 0)],
        },
        padding=0,
    )
    assert any(
        [
            z_coords == {"select_1": 0.0, "select_2": -0.0},
            z_coords == {"select_1": -0.0, "select_2": 0.0},
        ]
    )

    z_coords = get_z_coords(
        {
            "select_1": [(0, 0, 0)],
            "select_2": [(0, 0, -3)],
        },
        padding=1,
    )
    assert any(
        [
            z_coords == {"select_1": 2, "select_2": -2},
            z_coords == {"select_1": -2, "select_2": 2},
        ]
    )


def test_get_z_coords_empty_select():
    """Test return of position 0 if no selection."""
    z_coords = get_z_coords({})
    assert list(z_coords.values()) == [0]


def test_get_z_coords_1_select():
    """Test return of position 0 if only one selection."""
    z_coords = get_z_coords({"select_1": [(1, 2, 3)]})
    assert list(z_coords.values()) == [0]


def test_extract_restraint_entries(example_tbl_file):
    """Test read and extract restraints from tbl file."""
    list_of_restraints = extract_restraint_entries(example_tbl_file)
    assert len(list_of_restraints) == 4


def test_random_removal_nbfiles(example_tbl_file):
    """Test random removal number of files."""
    output_filepath = random_removal(example_tbl_file, 0.5, nb_tbl=5)
    assert output_filepath.exists()
    with tarfile.open(output_filepath, "r:gz") as tar:
        members = tar.getmembers()
        assert len(members) == 5
    output_filepath.unlink()


def test_random_removal_ratio(example_tbl_file):
    """Test random removal ratio."""
    output_filepath = random_removal(example_tbl_file, 0.5, nb_tbl=5)
    assert output_filepath.exists()
    with tarfile.open(output_filepath, "r:gz") as tar:
        members = tar.getmembers()
        assert len(members) == 5
        # get an example
        ex_0_5 = tar.extractfile(members[0]).read().decode("utf-8")
        assert ex_0_5.lower().count("assi") == 2
    output_filepath.unlink()
    output_filepath2 = random_removal(example_tbl_file, 0.2, nb_tbl=3)
    assert output_filepath2.exists()
    with tarfile.open(output_filepath2, "r:gz") as tar:
        members2 = tar.getmembers()
        assert len(members2) == 3
        # get an example
        ex_0_2 = tar.extractfile(members2[0]).read().decode("utf-8")
    # Checking that removing 50 % leads to lower number of assign statements than 20%
    assert ex_0_5.lower().count("assi") < ex_0_2.lower().count("assi")
    output_filepath2.unlink()

    # Check for errors in ratio selection
    # Ratio removed == 0
    with pytest.raises(ValueError):
        output_error = random_removal(example_tbl_file, 0, nb_tbl=1)
        assert output_error is None
    # Ratio removed == 1
    with pytest.raises(ValueError):
        output_error2 = random_removal(example_tbl_file, 1, nb_tbl=1)
        assert output_error2 is None


def test_random_removal_seed(example_tbl_file):
    """Test random removal seed."""
    output_filepath = random_removal(example_tbl_file, 0.2, nb_tbl=1, seed=1)
    assert output_filepath.exists()
    with tarfile.open(output_filepath, "r:gz") as tar:
        members = tar.getmembers()
        # get an example
        rd_rm_tbl1 = tar.extractfile(members[0]).read().decode("utf-8")
    output_filepath.unlink()
    output_filepath2 = random_removal(example_tbl_file, 0.2, nb_tbl=1, seed=1)
    with tarfile.open(output_filepath2, "r:gz") as tar2:
        members2 = tar2.getmembers()
        # get an example
        rd_rm_tbl2 = tar2.extractfile(members2[0]).read().decode("utf-8")
    output_filepath2.unlink()
    output_filepath3 = random_removal(example_tbl_file, 0.2, nb_tbl=1, seed=5)
    with tarfile.open(output_filepath2, "r:gz") as tar3:
        members3 = tar3.getmembers()
        # get an example
        rd_rm_tbl3 = tar3.extractfile(members3[0]).read().decode("utf-8")
    output_filepath3.unlink()

    # Check that 1 and 2 are equal !
    assert rd_rm_tbl1 == rd_rm_tbl2
    # Check that 2 and 3 are different !
    assert rd_rm_tbl2 != rd_rm_tbl3
    assert rd_rm_tbl1 != rd_rm_tbl3


def test_restrain_ligand(example_liganded_pdbfile):
    """Tests related to the restrain_ligand subcommand."""
    ligand_restraints_default = restrain_ligand(example_liganded_pdbfile, "G39")
    assert ligand_restraints_default is not None
    assert ligand_restraints_default != ""

    # Reducing the distance deviation
    ligand_restraints_small_devi = restrain_ligand(
        example_liganded_pdbfile, "G39",
        deviation=0.5,
        )
    assert ligand_restraints_small_devi is not None
    assert ligand_restraints_small_devi != ""
    assert ligand_restraints_default != ligand_restraints_small_devi

    # Reducing the radius
    ligand_restraints_small_radius = restrain_ligand(
        example_liganded_pdbfile, "G39",
        radius=2.0,
        )
    assert ligand_restraints_small_radius is not None
    assert ligand_restraints_small_radius != ""
    assert ligand_restraints_default != ligand_restraints_small_radius
    assert ligand_restraints_small_radius.count("assi") <= ligand_restraints_default.count("assi")

    # Reducing the number of restaints
    ligand_restraints_low_nb_restraints = restrain_ligand(
        example_liganded_pdbfile, "G39",
        max_restraints=2,
        )
    assert ligand_restraints_low_nb_restraints is not None
    assert ligand_restraints_low_nb_restraints != ""
    assert ligand_restraints_default != ligand_restraints_low_nb_restraints
    assert ligand_restraints_low_nb_restraints.count("assi") <= ligand_restraints_default.count("assi")
    assert ligand_restraints_low_nb_restraints.count("assi") == 2
