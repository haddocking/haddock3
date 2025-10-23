"""Test lib PDB."""

import shutil
import pytest
import tempfile

from pathlib import Path

from . import data_folder, golden_data
from haddock.libs.libio import PDBFile
from haddock.libs.libpdb import (
    add_TER_on_chain_breaks,
    check_combination_chains,
    read_chainids,
    read_segids,
    sanitize,
)


chainC = [
    "ATOM      3  CA  ARG C   4      37.080  43.455  -3.421  1.00  0.00      C    C  ",  # noqa: E501
    "ATOM      3  CA  GLU C   6      33.861  45.127  -2.233  1.00  0.00      C    C  ",  # noqa: E501
    "ATOM      3  CA  ALA C   7      35.081  45.036   1.305  1.00  0.00      C    C  ",  # noqa: E501
]


@pytest.mark.parametrize(
    "lines,expected",
    ((chainC, ["C", "C", "C"]),),
)
def test_read_chain_ids(lines, expected):
    result = read_chainids(lines)
    assert result == expected


@pytest.mark.parametrize(
    "lines,expected",
    ((chainC, ["C", "C", "C"]),),
)
def test_read_seg_ids(lines, expected):
    result = read_segids(lines)
    assert result == expected


def test_add_TER_on_chain_breaks(monkeypatch):
    """Test proper detection of chain breaks."""
    with tempfile.TemporaryDirectory(".") as tdir:
        output_fpath = Path(tdir, "test_withTER.pdb")
        add_TER_on_chain_breaks(
            Path(data_folder, "noTER.pdb"),
            output_fpath,
        )
        # Make sure file was created
        assert output_fpath.exists()
        assert output_fpath.stat().st_size != 0
        # Make sure the file resemble the reference one
        with (
            open(output_fpath, "r") as testfile,
            open(Path(data_folder, "withTER.pdb")) as reffile,
        ):
            # Line by line comparisons
            for test_, ref_ in zip(testfile, reffile):
                assert test_.rstrip() == ref_.rstrip()


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


@pytest.fixture
def protlig_complex_pdb():
    src = Path(golden_data, "protlig_complex_1.pdb")
    with tempfile.TemporaryDirectory() as tmpdir:
        dst = Path(tmpdir, "protlig_complex_1.pdb")
        shutil.copy(src, dst)
        yield dst


@pytest.fixture
def ligand_topology_file():
    src = Path(golden_data, "ligand.top")
    with tempfile.TemporaryDirectory() as tmpdir:
        dst = Path(tmpdir, "ligand.top")
        shutil.copy(src, dst)
        yield dst


def test_check_combination_chains(good_rigid_molecules, wrong_rigid_molecules):
    """Test check_combination_chains."""
    exp_chains = ["A", "B"]
    obs_chains = check_combination_chains(good_rigid_molecules)
    assert obs_chains == exp_chains
    # when input molecules share chains there should be a ValueError
    with pytest.raises(ValueError):
        check_combination_chains(wrong_rigid_molecules)


def test_sanitize(protlig_complex_pdb, monkeypatch):

    monkeypatch.chdir(protlig_complex_pdb.parent)

    output = sanitize(
        pdb_file_path=protlig_complex_pdb, overwrite=False, custom_topology=None
    )

    assert len(output.read_text()) == 476608


def test_sanitize_with_custom_topology(
    protlig_complex_pdb, ligand_topology_file, monkeypatch
):
    monkeypatch.chdir(protlig_complex_pdb.parent)

    output = sanitize(
        pdb_file_path=protlig_complex_pdb,
        overwrite=False,
        custom_topology=ligand_topology_file,
    )

    assert len(output.read_text()) == 478552


def test_sanitize_not_changing_global(
    protlig_complex_pdb, ligand_topology_file, monkeypatch
):
    monkeypatch.chdir(protlig_complex_pdb.parent)

    from haddock.libs.libpdb import _to_keep

    initial_global_to_keep = len(_to_keep)

    _ = sanitize(
        pdb_file_path=protlig_complex_pdb,
        overwrite=False,
        custom_topology=ligand_topology_file,
    )

    current_global_to_keep = len(_to_keep)

    assert current_global_to_keep == initial_global_to_keep
