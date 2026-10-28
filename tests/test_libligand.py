"""Tests for libligand."""

import pytest

from pathlib import Path
import tempfile
import shutil
from haddock.libs.libligand import (
    identify_unknown_hetatms,
    extract_ligand,
    run_prodrg,
    _remove_nbonds,
    _sanitize_atom_names,
    _used_cofactor_separators,
    _unsupported_metal_symbols,
    _demetalise_atom,
    _prodrg_messages,
    find_chain_residues,
    find_linked_sugars,
    ion_naming_hint,
)
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
def azs_complex_pdb():
    src = Path(GOLDEN_DATA, "1AZS_l_u.pdb")
    with tempfile.TemporaryDirectory() as tmpdir:
        dst = Path(tmpdir, "1AZS_l_u.pdb")
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


def test_run_prodrg_with_relative_output_dir(ligand_pdb, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    # Simulate topoaa passing Path(".") while CWD == the module output dir
    top, par = run_prodrg(ligand_pdb, Path("."))
    assert top.parent == tmp_path
    assert par.parent == tmp_path
    assert top.exists()
    assert par.exists()


def test_run_prodrg_fails_on_complex(azs_complex_pdb, tmp_path):
    """run_prodrg raises RuntimeError when given a protein-ligand complex.

    prodrg only accepts small molecule inputs; passing a full complex
    causes it to abort with 'Too many lines of text in input'.
    """
    with pytest.raises(RuntimeError):
        run_prodrg(azs_complex_pdb, tmp_path)


def test_extract_ligand_keeps_only_selected_residues(azs_complex_pdb, tmp_path):
    """extract_ligand strips everything except the requested residues."""
    dest = tmp_path / "ligand_only.pdb"
    out = extract_ligand(azs_complex_pdb, ["GSP"], dest)
    assert out == dest
    lines = dest.read_text().splitlines()
    atom_lines = [ln for ln in lines if ln.startswith(("ATOM  ", "HETATM"))]
    assert atom_lines  # at least one ligand atom was kept
    assert all(ln[17:20].strip() == "GSP" for ln in atom_lines)
    # no leftover protein residues
    assert identify_unknown_hetatms(dest) == ["GSP"]


def test_unsupported_metal_symbols():
    """Metal symbols are flagged; prodrg-supported halides are not."""
    metals = _unsupported_metal_symbols()
    assert "PB" in metals  # lead -> mislabelling to fix
    assert "CL" not in metals  # chlorine -> prodrg supports it
    assert "BR" not in metals  # bromine -> prodrg supports it


def test_demetalise_atom_fixes_mislabeled_metal():
    """A phosphate PB mislabelled as lead is collapsed back to phosphorus."""
    metals = _unsupported_metal_symbols()
    # CNS writes the beta-phosphorus with element "PB" and the name shifted to
    # the two-letter-element column.
    line = (
        "ATOM   3437 PB   GSP B 395      "
        "-0.953   3.951   3.655  1.00 15.00      B   PB  \n"
    )
    assert line[76:78] == "PB"  # precondition: mislabelled
    fixed = _demetalise_atom(line, metals)
    assert fixed[76:78] == " P"  # element collapsed to first character
    assert fixed[12:16] == " PB "  # name re-justified to 1-letter-element column
    assert fixed[17:20] == "GSP"  # residue name preserved


def test_demetalise_atom_keeps_halogen():
    """A genuine chlorine atom is left untouched (prodrg supports Cl)."""
    metals = _unsupported_metal_symbols()
    line = (
        "ATOM   3437 CL   LIG B 395      "
        "-0.953   3.951   3.655  1.00 15.00      B   CL  \n"
    )
    assert _demetalise_atom(line, metals) == line


def test_extract_ligand_demetalises_mislabeled_atom(tmp_path):
    """extract_ligand corrects a metal-mislabelled ligand atom."""
    pdb = tmp_path / "mislabeled.pdb"
    pdb.write_text(
        "ATOM   3437 PB   GSP B 395      "
        "-0.953   3.951   3.655  1.00 15.00      B   PB  \n"
        "END\n"
    )
    dest = tmp_path / "ligand_only.pdb"
    extract_ligand(pdb, ["GSP"], dest)
    atom_lines = [
        ln
        for ln in dest.read_text().splitlines()
        if ln.startswith(("ATOM  ", "HETATM"))
    ]
    assert len(atom_lines) == 1
    assert atom_lines[0][76:78] == " P"
    assert atom_lines[0][12:16] == " PB "


def test_extract_ligand_keeps_single_copy(tmp_path):
    """Only the first copy of a repeated ligand residue is kept."""
    pdb = tmp_path / "multi_copy.pdb"
    pdb.write_text(
        "ATOM      1  C1  LIG A 284       0.000   0.000   0.000  1.00  0.00      A    C\n"
        "ATOM      2  C2  LIG A 284       1.000   0.000   0.000  1.00  0.00      A    C\n"
        "ATOM      3  C1  LIG A 285       0.000   0.000   0.000  1.00  0.00      A    C\n"
        "ATOM      4  C2  LIG A 285       1.000   0.000   0.000  1.00  0.00      A    C\n"
        "END\n"
    )
    dest = tmp_path / "ligand_only.pdb"
    extract_ligand(pdb, ["LIG"], dest)
    atom_lines = [
        ln
        for ln in dest.read_text().splitlines()
        if ln.startswith(("ATOM  ", "HETATM"))
    ]
    assert len(atom_lines) == 2
    assert all(ln[21:27] == "A 284 " for ln in atom_lines)


def test_run_prodrg_on_complex_with_ligand_resnames(azs_complex_pdb, tmp_path):
    """run_prodrg succeeds on a complex when the ligand residues are given.

    The ligand is extracted from the surrounding protein before prodrg runs,
    which is the fix for feeding a full system to prodrg.
    """
    top, par = run_prodrg(azs_complex_pdb, tmp_path, ligand_resnames=["GSP"])
    assert top.exists()
    assert par.exists()
    assert "MASS" in top.read_text()
    assert "BOND" in par.read_text()
    assert "NBONds" not in par.read_text()


def test_run_prodrg_concatenates_multiple_ligands(tmp_path, monkeypatch):
    """Each distinct ligand runs prodrg once; the outputs are concatenated."""
    pdb = tmp_path / "two_ligands.pdb"
    pdb.write_text(
        "ATOM      1  C1  LGA A 284       0.000   0.000   0.000  1.00  0.00      A    C\n"
        "ATOM      2  C1  LGB A 285       0.000   0.000   0.000  1.00  0.00      A    C\n"
        "END\n"
    )

    calls: list[tuple[str, object]] = []

    def fake_single(ligand_pdb, cnssep=None):
        resname = ligand_pdb.stem
        calls.append((resname, cnssep))
        return (f"TOP {resname}", f"PAR {resname}")

    monkeypatch.setattr("haddock.libs.libligand._run_prodrg_single", fake_single)

    top, par = run_prodrg(pdb, tmp_path, ligand_resnames=["LGA", "LGB"])

    # prodrg is invoked once per distinct ligand
    assert [c[0] for c in calls] == ["LGA", "LGB"]
    # each ligand gets a distinct CNSSEP character that avoids the cofactor ones
    seps = [c[1] for c in calls]
    assert len(set(seps)) == len(seps)
    assert not (set(seps) & _used_cofactor_separators())
    # both ligand topologies/parameters are concatenated
    assert "TOP LGA" in top.read_text() and "TOP LGB" in top.read_text()
    assert "PAR LGA" in par.read_text() and "PAR LGB" in par.read_text()


def test_run_prodrg_single_ligand_uses_cofactor_safe_cnssep(tmp_path, monkeypatch):
    """A single ligand also gets a cofactor-safe CNSSEP separator."""
    pdb = tmp_path / "one.pdb"
    pdb.write_text(
        "ATOM      1  C1  LGA A 284       0.000   0.000   0.000  1.00  0.00      A    C\n"
        "END\n"
    )

    seps: list[object] = []

    def fake_single(ligand_pdb, cnssep=None):
        seps.append(cnssep)
        return ("TOP", "PAR")

    monkeypatch.setattr("haddock.libs.libligand._run_prodrg_single", fake_single)

    run_prodrg(pdb, tmp_path, ligand_resnames=["LGA"])
    assert len(seps) == 1
    assert seps[0] is not None
    assert seps[0] not in _used_cofactor_separators()


def test_run_prodrg_deduplicates_ligand_resnames(tmp_path, monkeypatch):
    """A repeated residue name only triggers a single prodrg run."""
    pdb = tmp_path / "dup.pdb"
    pdb.write_text(
        "ATOM      1  C1  LGA A 284       0.000   0.000   0.000  1.00  0.00      A    C\n"
        "END\n"
    )

    calls: list[str] = []

    def fake_single(ligand_pdb, cnssep=None):
        calls.append(ligand_pdb.stem)
        return ("TOP", "PAR")

    monkeypatch.setattr("haddock.libs.libligand._run_prodrg_single", fake_single)

    run_prodrg(pdb, tmp_path, ligand_resnames=["LGA", "LGA"])
    assert calls == ["LGA"]


def test_run_prodrg_failure_names_residue(tmp_path):
    """A prodrg failure names the residue and gives prodrg's own reason."""
    pdb = tmp_path / "bad.pdb"
    pdb.write_text(
        "HETATM    1  C1  XYZ A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "HETATM    2  C2  XYZ A   1       0.000   1.500   0.000  1.00  0.00           C\n"
        "HETATM    3  U1  XYZ A   1       1.500   0.000   0.000  1.00  0.00           U\n"
        "END\n"
    )
    with pytest.raises(RuntimeError) as err:
        run_prodrg(pdb, tmp_path, ligand_resnames=["XYZ"])
    message = str(err.value)
    assert message.startswith("residue XYZ:")
    assert "unsupported atom type  U" in message
    assert "ERRDRG>" in message


def test_run_prodrg_rejects_single_atoms(tmp_path):
    """Ions and other one- or two-atom residues are not passed to prodrg."""
    pdb = tmp_path / "zn.pdb"
    pdb.write_text(
        "HETATM    1 ZN    ZN A   1       0.000   0.000   0.000  1.00  0.00          ZN\n"
        "END\n"
    )
    with pytest.raises(RuntimeError, match=r"residue ZN: .*\(found 1\)"):
        run_prodrg(pdb, tmp_path, ligand_resnames=["ZN"])


def test_prodrg_messages():
    """Only the error lines of the prodrg output are kept."""
    output = (
        "PRODRG> banner\nPRODRG>\nPRODRG> Found unsupported atom type  U\n"
        "ERRDRG> Currently only N, C, O, S, P, Cl, Br, F, I are supported.\n"
        "PRODRG> Program terminated unsuccessfully, sorry!\n"
    )
    assert _prodrg_messages(output).splitlines() == [
        "PRODRG> Found unsupported atom type  U",
        "ERRDRG> Currently only N, C, O, S, P, Cl, Br, F, I are supported.",
        "PRODRG> Program terminated unsuccessfully, sorry!",
    ]
    # without an error line, the last lines are returned
    assert _prodrg_messages("a\nb\nc", nlines=2) == "b\nc"


def test_ion_naming_hint(tmp_path):
    """Ion-like residues get the supported residue and signed atom names."""
    pdb = tmp_path / "ions.pdb"
    pdb.write_text(
        "HETATM    1 CU   CU3 A 101       0.000   0.000   0.000  1.00  0.00          CU\n"
        "HETATM    2 CL   CL- A 102       5.000   0.000   0.000  1.00  0.00          CL\n"
        "HETATM    3 GD   GD  A 103      10.000   0.000   0.000  1.00  0.00          GD\n"
        "HETATM    4  C1  XYZ A   1       0.000   5.000   0.000  1.00  0.00           C\n"
        "HETATM    5  O1  XYZ A   1       1.500   5.000   0.000  1.00  0.00           O\n"
        "END\n"
    )
    hint = ion_naming_hint(pdb, ["CU3", "CL-", "GD", "XYZ"])
    assert "CU3: if this is an ion, use residue CU1 with atom CU+1" in hint
    assert "residue CU2 with atom CU+2" in hint
    assert "CL-: if this is an ion, use residue CL1 with atom CL-1" in hint
    assert "GD: single atom, not a supported ion" in hint
    assert "XYZ" not in hint
    assert ion_naming_hint(pdb, ["XYZ"]) == ""


def test_used_cofactor_separators_excludes_known_chars():
    """The cofactor separator set reflects the atom types in cofactors.top."""
    seps = _used_cofactor_separators()
    # cofactors.top atom types use these separators (2nd character)
    assert {"A", "C", "F", "G", "N"} <= seps


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


def _atom(serial, name, resname, chain, resid, x, element):
    """Format a PDB ATOM line."""
    return (
        f"ATOM  {serial:5d} {name:<4s} {resname:>3s} {chain}{resid:4d}    "
        f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00          {element:>2s}\n"
    )


def test_find_chain_residues(tmp_path):
    """Unknown residues bonded to their chain neighbours are found."""
    pdb = tmp_path / "chain.pdb"
    pdb.write_text(
        # protein chain A: ALA 1 - XXX 2 - ALA 3, peptide bonds of 1.33 A
        _atom(1, "N", "ALA", "A", 1, 0.0, "N")
        + _atom(2, "C", "ALA", "A", 1, 1.0, "C")
        + _atom(3, "N", "XXX", "A", 2, 2.33, "N")
        + _atom(4, "C", "XXX", "A", 2, 3.33, "C")
        + _atom(5, "N", "ALA", "A", 3, 4.66, "N")
        # C-terminal unknown residue of chain B, bonded to the previous one only
        + _atom(6, "C", "ALA", "B", 1, 10.0, "C")
        + _atom(7, "N", "YYY", "B", 2, 11.33, "N")
        + _atom(8, "C", "YYY", "B", 2, 12.33, "C")
        # free amino-acid-like ligand: same atoms, far from any chain
        + _atom(9, "N", "LIG", "B", 3, 30.0, "N")
        + _atom(10, "C", "LIG", "B", 3, 31.0, "C")
        # unknown residue whose neighbour is in another chain
        + _atom(11, "N", "ZZZ", "C", 1, 13.66, "N")
        # nucleic acid: O3' of DA 1 bonded to P of NUC 2
        + _atom(12, "O3'", "DA", "D", 1, 50.0, "O")
        + _atom(13, "P", "NUC", "D", 2, 51.6, "P")
        + "END\n"
    )
    linked = find_chain_residues(pdb, ["XXX", "YYY", "LIG", "ZZZ", "NUC"])
    assert linked == [
        "XXX (chain A, residue 2)",
        "YYY (chain B, residue 2)",
        "NUC (chain D, residue 2)",
    ]


def _sugar(serial, resname, chain, resid, x, ring, anomeric=None, anomeric_x=0.0):
    """Format the ring atoms of a sugar, the ``anomeric`` carbon at ``anomeric_x``."""
    lines = ""
    for i, name in enumerate(sorted(ring)):
        element = name[0]
        pos = anomeric_x if name == anomeric else x + 3 * i
        lines += _atom(serial + i, name, resname, chain, resid, pos, element)
    return lines


def test_find_linked_sugars(tmp_path):
    """Unknown sugars bonded to a protein or another sugar are found."""
    pyranose = {"C1", "C2", "C3", "C4", "C5", "O5"}
    ketose = {"C2", "C3", "C4", "C5", "C6", "O6"}
    pdb = tmp_path / "glycan.pdb"
    pdb.write_text(
        # O-linked sugar: C1 of XGA bonded to OG of SER (1.43 A)
        _atom(1, "OG", "SER", "A", 1, 0.0, "O")
        + _sugar(2, "XGA", "B", 1, 100.0, pyranose, "C1", 1.43)
        # sialic-acid-like ketose: C2 bonded to O3 of a known sugar
        + _atom(20, "O3", "GAL", "C", 1, 200.0, "O")
        + _sugar(21, "XIA", "C", 2, 300.0, ketose, "C2", 201.43)
        # free sugar, not bonded to anything
        + _sugar(40, "FRE", "D", 1, 500.0, pyranose)
        + "END\n"
    )
    linked = find_linked_sugars(pdb, ["XGA", "XIA", "FRE"])
    assert linked == [
        "XGA (chain B, residue 1) linked to SER 1 OG",
        "XIA (chain C, residue 2) linked to GAL 1 O3",
    ]
