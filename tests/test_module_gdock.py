"""Test the gdock module."""

import tempfile
from pathlib import Path

from haddock.libs import libpdb
from haddock.libs.libio import working_directory
from haddock.modules.sampling.gdock.gdock import (
    GdockWrapper,
    extract_pairs_from_tbl,
    parse_restraints,
)

from . import golden_data, has_gdock


def test_parse_restraints_receptor_to_ligand():
    """Anchor residue in the receptor chain, partners in the ligand chain."""
    tbl_content = (
        "assign (resi 933 and segid A)\n"
        "(\n"
        "       (resi 6 and segid B)\n"
        "        or\n"
        "       (resi 8 and segid B)\n"
        ") 2.0 2.0 0.0\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False) as tbl_file:
        tbl_file.write(tbl_content)
        tbl_path = tbl_file.name

    pairs = parse_restraints(tbl_path, receptor_chain="A", ligand_chain="B")

    assert pairs == [(933, 6), (933, 8)]


def test_parse_restraints_resid_keyword():
    """The `resid` keyword (used by HADDOCK AIR files) is also supported."""
    tbl_content = (
        "assign (resid 933 and segid A)\n"
        "(\n"
        "       (resid 6 and segid B)\n"
        ") 2.0 2.0 0.0\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False) as tbl_file:
        tbl_file.write(tbl_content)
        tbl_path = tbl_file.name

    pairs = parse_restraints(tbl_path, receptor_chain="A", ligand_chain="B")

    assert pairs == [(933, 6)]


def test_parse_restraints_ligand_to_receptor():
    """Anchor residue in the ligand chain, partner in the receptor chain."""
    tbl_content = (
        "assign (resi 6 and segid B)\n(\n       (resi 933 and segid A)\n) 2.0 2.0 0.0\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False) as tbl_file:
        tbl_file.write(tbl_content)
        tbl_path = tbl_file.name

    pairs = parse_restraints(tbl_path, receptor_chain="A", ligand_chain="B")

    assert pairs == [(933, 6)]


def test_parse_restraints_unrelated_chains_skipped():
    """Restraints that do not span receptor/ligand chains are ignored."""
    tbl_content = (
        "assign (resi 1 and segid A)\n(\n       (resi 2 and segid A)\n) 2.0 2.0 0.0\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False) as tbl_file:
        tbl_file.write(tbl_content)
        tbl_path = tbl_file.name

    pairs = parse_restraints(tbl_path, receptor_chain="A", ligand_chain="B")

    assert pairs == []


def test_extract_pairs_active_active_only():
    """Only the residues that appear as anchors on both sides form pairs.

    993-A and 6-B are both anchors (active). 8-B and 42-B appear only in OR
    groups (passive). 936-A and 55-A appear only in OR groups (passive on A
    side). The only active-active pair is (933, 6).
    """
    tbl_content = (
        "assign (resi 933 and segid A)\n"
        "(\n"
        "       (resi 6 and segid B)\n"
        "        or\n"
        "       (resi 8 and segid B)\n"
        "        or\n"
        "       (resi 42 and segid B)\n"
        ") 2.0 2.0 0.0\n"
        "\n"
        "assign (resi 6 and segid B)\n"
        "(\n"
        "       (resi 933 and segid A)\n"
        "        or\n"
        "       (resi 936 and segid A)\n"
        "        or\n"
        "       (resi 55 and segid A)\n"
        ") 2.0 2.0 0.0\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False) as f:
        f.write(tbl_content)
        tbl_path = f.name

    pairs = extract_pairs_from_tbl(tbl_path, receptor_chain="A", ligand_chain="B")

    assert pairs == [(933, 6)]


def test_extract_pairs_multiple_actives_both_sides():
    """All active-active combinations are returned when both sides have multiple actives.

    active A: {933, 936}  active B: {6, 8}  passive B: {42}  passive A: {55}
    Expected cross-product: (933,6), (933,8), (936,6), (936,8)
    """
    tbl_content = (
        "assign (resi 933 and segid A)\n"
        "(\n"
        "       (resi 6 and segid B)\n"
        "        or\n"
        "       (resi 8 and segid B)\n"
        "        or\n"
        "       (resi 42 and segid B)\n"
        ") 2.0 2.0 0.0\n"
        "\n"
        "assign (resi 936 and segid A)\n"
        "(\n"
        "       (resi 6 and segid B)\n"
        "        or\n"
        "       (resi 8 and segid B)\n"
        "        or\n"
        "       (resi 42 and segid B)\n"
        ") 2.0 2.0 0.0\n"
        "\n"
        "assign (resi 6 and segid B)\n"
        "(\n"
        "       (resi 933 and segid A)\n"
        "        or\n"
        "       (resi 936 and segid A)\n"
        "        or\n"
        "       (resi 55 and segid A)\n"
        ") 2.0 2.0 0.0\n"
        "\n"
        "assign (resi 8 and segid B)\n"
        "(\n"
        "       (resi 933 and segid A)\n"
        "        or\n"
        "       (resi 936 and segid A)\n"
        "        or\n"
        "       (resi 55 and segid A)\n"
        ") 2.0 2.0 0.0\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False) as f:
        f.write(tbl_content)
        tbl_path = f.name

    pairs = extract_pairs_from_tbl(tbl_path, receptor_chain="A", ligand_chain="B")

    assert pairs == [(933, 6), (933, 8), (936, 6), (936, 8)]


def test_extract_pairs_no_active_partners_returns_empty():
    """If OR-group contains only passive residues, no pair is emitted.

    Only 933-A is an anchor (active). The OR group on B side has 42-B and
    44-B, neither of which appears as an anchor — both are passive. Result: [].
    """
    tbl_content = (
        "assign (resi 933 and segid A)\n"
        "(\n"
        "       (resi 42 and segid B)\n"
        "        or\n"
        "       (resi 44 and segid B)\n"
        ") 2.0 2.0 0.0\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False) as f:
        f.write(tbl_content)
        tbl_path = f.name

    pairs = extract_pairs_from_tbl(tbl_path, receptor_chain="A", ligand_chain="B")

    assert pairs == []


def test_extract_pairs_same_resi_different_segid_not_confused():
    """Residue number 6 on segid A and segid B are treated as distinct.

    6-A is an active anchor; 6-B appears only in an OR group (passive on B
    side). They must not be conflated just because they share residue number 6.
    """
    tbl_content = (
        "assign (resi 6 and segid A)\n"
        "(\n"
        "       (resi 6 and segid B)\n"
        "        or\n"
        "       (resi 99 and segid B)\n"
        ") 2.0 2.0 0.0\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False) as f:
        f.write(tbl_content)
        tbl_path = f.name

    pairs = extract_pairs_from_tbl(tbl_path, receptor_chain="A", ligand_chain="B")

    # 6-B and 99-B never appear as anchors, so no active partner exists → []
    assert pairs == []


@has_gdock
def test_gdock_wrapper_run_and_save_models():
    """Test running gdock and saving the resulting models to disk."""
    with tempfile.TemporaryDirectory() as tempdir:
        tmp_complex = Path(tempdir, "complex.pdb")
        tmp_complex.write_bytes(
            Path(golden_data, "protprot_complex_1.pdb").read_bytes()
        )

        with working_directory(tempdir):
            new_models = [Path(m).resolve() for m in libpdb.split_by_chain(tmp_complex)]
        receptor_pdb_file = next(m for m in new_models if m.stem.endswith("_A"))
        ligand_pdb_file = next(m for m in new_models if m.stem.endswith("_B"))

        wrapper = GdockWrapper(
            receptor_pdb_file=receptor_pdb_file,
            ligand_pdb_file=ligand_pdb_file,
            max_generations=2,
            seed=1,
            sampling=2,
        )
        wrapper.run()
        models = wrapper.save_models(tempdir)

        assert len(models) <= 2
        for model in models:
            model_path = Path(tempdir, model["file_name"])
            assert model_path.exists()
            content = model_path.read_text()
            assert "ATOM" in content
            for key in ("fitness", "vdw", "elec", "desolv", "air"):
                assert isinstance(model[key], float)
