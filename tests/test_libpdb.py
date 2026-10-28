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
    handle_input_reference,
    read_chainids,
    read_segids,
    rename_residues,
    sanitize,
    split_ensemble,
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


@pytest.fixture
def atom_lines():
    return (
        "ATOM      1  N   MET A   1      -0.497   0.778  -7.571  1.00  0.00      A    N  ",
        "ATOM      2  CA  MET A   1       0.814   0.227  -7.319  1.00  0.00      A    C  ",
    )


def test_split_ensemble_single_model(tmp_path, monkeypatch, atom_lines):
    monkeypatch.chdir(tmp_path)
    content = "\n".join([*atom_lines, "END", ""])
    dst = tmp_path / "single.pdb"
    dst.write_text(content)

    result = split_ensemble(dst, dest=tmp_path)

    assert len(result) == 1
    assert result[0].name == "single.pdb"
    assert result[0] == dst


def test_split_ensemble_proper_ensemble(tmp_path, monkeypatch, atom_lines):
    monkeypatch.chdir(tmp_path)
    content = "\n".join(
        [
            "MODEL        1",
            *atom_lines,
            "ENDMDL",
            "MODEL        2",
            *atom_lines,
            "ENDMDL",
            "END",
            "",
        ]
    )
    dst = tmp_path / "ensemble.pdb"
    dst.write_text(content)

    result = split_ensemble(dst, dest=tmp_path)

    assert len(result) == 2
    assert [r.name for r in result] == ["ensemble_1.pdb", "ensemble_2.pdb"]
    assert all(r.stat().st_size > 0 for r in result)


def test_split_ensemble_model_without_endmdl(tmp_path, monkeypatch, atom_lines):
    monkeypatch.chdir(tmp_path)
    content = "\n".join(["MODEL        1", *atom_lines, "END", ""])
    dst = tmp_path / "model_without_endmdl.pdb"
    dst.write_text(content)

    result = split_ensemble(dst, dest=tmp_path)

    assert len(result) == 1
    assert result[0].name == "model_without_endmdl_1.pdb"
    assert result[0].stat().st_size > 0


def test_split_ensemble_multiple_models_without_endmdl(
    tmp_path, monkeypatch, atom_lines
):
    monkeypatch.chdir(tmp_path)
    content = "\n".join(
        [
            "MODEL        1",
            *atom_lines,
            "MODEL        2",
            *atom_lines,
            "END",
            "",
        ]
    )
    dst = tmp_path / "two_models_no_endmdl.pdb"
    dst.write_text(content)

    result = split_ensemble(dst, dest=tmp_path)

    assert len(result) == 2
    assert [r.name for r in result] == [
        "two_models_no_endmdl_1.pdb",
        "two_models_no_endmdl_2.pdb",
    ]
    assert all(r.stat().st_size > 0 for r in result)


def test_split_ensemble_non_sequential_model_numbers(tmp_path, monkeypatch, atom_lines):
    """pdb_tidy is used in `split_ensemble` and it renumbers MODEL records sequentially,
    so non-sequential IDs (ex: 1, 5) produce sequential output files (_1.pdb, _2.pdb)
    rather than the original numbers (_1.pdb, _5.pdb)."""
    monkeypatch.chdir(tmp_path)
    content = "\n".join(
        [
            "MODEL        1",
            *atom_lines,
            "ENDMDL",
            "MODEL        5",
            *atom_lines,
            "ENDMDL",
            "END",
            "",
        ]
    )
    dst = tmp_path / "non_sequential.pdb"
    dst.write_text(content)

    result = split_ensemble(dst, dest=tmp_path)

    assert len(result) == 2
    assert [r.name for r in result] == ["non_sequential_1.pdb", "non_sequential_2.pdb"]
    assert all(r.stat().st_size > 0 for r in result)


# test handle_input_reference for capri

# Minimal ATOM line shared across fixtures
_ATOM_LINE = (
    "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00  0.00           C"
)


def test_handle_single_model_reference(tmp_path):
    ref = tmp_path / "single.pdb"
    ref.write_text(f"{_ATOM_LINE}\nEND\n")
    result = handle_input_reference(ref)
    assert len(result) == 1
    assert result[0] == ref


def test_handle_ensemble_reference(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    ref = tmp_path / "ensemble.pdb"
    ref.write_text(
        f"MODEL        1\n{_ATOM_LINE}\nENDMDL\n"
        f"MODEL        2\n{_ATOM_LINE}\nENDMDL\n"
        "END\n"
    )
    result = handle_input_reference(ref)
    assert len(result) == 2
    assert [p.name for p in result] == [
        "reference_model_1.pdb",
        "reference_model_2.pdb",
    ]
    assert all(p.exists() for p in result)


def test_handle_empty_file_raises(tmp_path):
    ref = tmp_path / "empty.pdb"
    ref.write_text("")
    with pytest.raises(ValueError, match="No atoms found"):
        handle_input_reference(ref)


def test_handle_remark_only_raises(tmp_path):
    ref = tmp_path / "remark_only.pdb"
    ref.write_text("REMARK this file has no coordinates\nEND\n")
    with pytest.raises(ValueError, match="No atoms found"):
        handle_input_reference(ref)


def test_rename_residues_aliases(tmp_path):
    """PDB residue names are renamed to the HADDOCK3 ones (e.g. TPO to TOP)."""
    pdb = Path(tmp_path, "aliases.pdb")
    pdb.write_text(
        "ATOM      1  P   TPO A   1       0.000   0.000   0.000  1.00  0.00           P\n"
        "HETATM    2 FE   HEM A   2       0.000   0.000   0.000  1.00  0.00          FE\n"
        "HETATM    3 AS   ARS A   3       0.000   0.000   0.000  1.00  0.00          AS\n"
        "ATOM      4  CA  HSD A   4       0.000   0.000   0.000  1.00  0.00           C\n"
        "END\n"
    )
    rename_residues(pdb)
    resnames = [line[17:20] for line in pdb.read_text().splitlines()[:-1]]
    assert resnames == ["TOP", "HEB", " AS", "HIS"]
    # column layout is unchanged
    assert all(len(line) == 78 for line in pdb.read_text().splitlines()[:-1])
    # sanitize keeps the renamed residues
    sanitize(pdb, overwrite=True)
    assert "TPO" not in pdb.read_text()
    assert " TOP A   1" in pdb.read_text()


def test_rename_residues_selenocysteine(tmp_path):
    """PDB selenocysteines become cysteines; HADDOCK3 SEC (serine) is kept."""
    pdb = Path(tmp_path, "sec.pdb")
    pdb.write_text(
        "ATOM      1  CB  SEC A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2 SE   SEC A   1       1.900   0.000   0.000  1.00  0.00          SE\n"
        "ATOM      3  CB  SEC A   2       5.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      4  OG  SEC A   2       6.400   0.000   0.000  1.00  0.00           O\n"
        "END\n"
    )
    converted = rename_residues(pdb)
    assert converted == ["chain A, residue 1"]
    lines = pdb.read_text().splitlines()
    assert lines[0][17:20] == "CYS"
    assert lines[1] == (
        "ATOM      2  SG  CYS A   1       1.900   0.000   0.000  1.00  0.00           S"
    )
    # the HADDOCK3 SEC residue (serine for covalent docking) is not changed
    assert [line[17:20] for line in lines[2:4]] == ["SEC", "SEC"]
