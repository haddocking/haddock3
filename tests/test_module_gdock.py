"""Test the gdock module."""

import tempfile
from pathlib import Path

from haddock.libs import libpdb
from haddock.libs.libio import working_directory
from haddock.modules.sampling.gdock.gdock import GdockWrapper, parse_restraints

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


def test_parse_restraints_ligand_to_receptor():
    """Anchor residue in the ligand chain, partner in the receptor chain."""
    tbl_content = (
        "assign (resi 6 and segid B)\n"
        "(\n"
        "       (resi 933 and segid A)\n"
        ") 2.0 2.0 0.0\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False) as tbl_file:
        tbl_file.write(tbl_content)
        tbl_path = tbl_file.name

    pairs = parse_restraints(tbl_path, receptor_chain="A", ligand_chain="B")

    assert pairs == [(933, 6)]


def test_parse_restraints_unrelated_chains_skipped():
    """Restraints that do not span receptor/ligand chains are ignored."""
    tbl_content = (
        "assign (resi 1 and segid A)\n"
        "(\n"
        "       (resi 2 and segid A)\n"
        ") 2.0 2.0 0.0\n"
    )
    with tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False) as tbl_file:
        tbl_file.write(tbl_content)
        tbl_path = tbl_file.name

    pairs = parse_restraints(tbl_path, receptor_chain="A", ligand_chain="B")

    assert pairs == []


@has_gdock
def test_gdock_wrapper_run_and_save_poses():
    """Test running gdock and saving the resulting poses to disk."""
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
        )
        wrapper.run()
        poses = wrapper.save_poses(tempdir, top=2)

        assert len(poses) <= 2
        for pose in poses:
            pose_path = Path(tempdir, pose["file_name"])
            assert pose_path.exists()
            content = pose_path.read_text()
            assert "ATOM" in content
            for key in ("fitness", "vdw", "elec", "desolv", "air"):
                assert isinstance(pose[key], float)
