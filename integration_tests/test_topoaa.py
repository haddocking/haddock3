import tempfile
from pathlib import Path

import pytest

from haddock.modules.topology.topoaa import \
    DEFAULT_CONFIG as DEFAULT_TOPOAA_CONFIG
from haddock.modules.topology.topoaa import HaddockModule as TopoaaModule

from . import CNS_EXEC, DATA_DIR, GOLDEN_DATA


@pytest.fixture(name="topoaa_module")
def fixture_topoaa_module():
    """Topoaa module with protein-protein input"""
    with tempfile.TemporaryDirectory() as tmpdir:
        topoaa = TopoaaModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_TOPOAA_CONFIG
        )
        topoaa.params["molecules"] = [
            Path(DATA_DIR, "docking-protein-protein/data/e2aP_1F3G.pdb"),
            Path(DATA_DIR, "docking-protein-protein/data/hpr_ensemble.pdb"),
        ]
        topoaa.params["mol1"] = {"prot_segid": "A"}
        topoaa.params["mol2"] = {"prot_segid": "B"}

        topoaa.params["cns_exec"] = CNS_EXEC

        yield topoaa


@pytest.fixture(name="topoaa_module_with_ligand")
def fixture_topoaa_module_with_peptide():
    """Topoaa fixture with peptide as input"""
    with tempfile.TemporaryDirectory() as tmpdir:
        topoaa = TopoaaModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_TOPOAA_CONFIG
        )
        topoaa.params["molecules"] = [
            Path(GOLDEN_DATA, "cyclic-peptide.pdb"),
        ]

        topoaa.params["cns_exec"] = CNS_EXEC

        yield topoaa


@pytest.fixture(name="topoaa_module_with_peptide")
def fixture_topoaa_module_with_ligand():
    """Topoaa fixture with ligand as input"""
    with tempfile.TemporaryDirectory() as tmpdir:
        topoaa = TopoaaModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_TOPOAA_CONFIG
        )
        topoaa.params["molecules"] = [
            Path(GOLDEN_DATA, "oseltamivir.pdb"),
        ]

        topoaa.params["cns_exec"] = CNS_EXEC

        yield topoaa


@pytest.fixture(name="ligand_toppar")
def fixture_ligand_toppar():
    """Fixture of a small molecule param/top"""
    return (Path(GOLDEN_DATA, "ligand.param"), Path(GOLDEN_DATA, "ligand.top"))


def test_topoaa_default(topoaa_module):
    """Test the topoaa module with default values"""

    topoaa_module.run()

    expected_inp = Path(topoaa_module.path, "e2aP_1F3G.inp")
    expected_psf = Path(topoaa_module.path, "e2aP_1F3G_haddock.psf")
    expected_pdb = Path(topoaa_module.path, "e2aP_1F3G_haddock.pdb")
    expected_gz = Path(topoaa_module.path, "e2aP_1F3G.out.gz")

    assert expected_inp.exists(), f"{expected_inp} does not exist"
    assert expected_psf.exists(), f"{expected_psf} does not exist"
    assert expected_gz.exists(), f"{expected_gz} does not exist"
    assert expected_pdb.exists(), f"{expected_pdb} does not exist"

    assert expected_inp.stat().st_size > 0, f"{expected_inp} is empty"

    for i in range(1, 10 + 1):
        expected_inp = Path(topoaa_module.path, f"hpr_ensemble_{i}.inp")
        expected_psf = Path(topoaa_module.path, f"hpr_ensemble_{i}_haddock.psf")
        expected_pdb = Path(topoaa_module.path, f"hpr_ensemble_{i}_haddock.pdb")
        expected_gz = Path(topoaa_module.path, f"hpr_ensemble_{i}.out.gz")

        assert expected_inp.exists(), f"{expected_inp} does not exist"
        assert expected_psf.exists(), f"{expected_psf} does not exist"
        assert expected_pdb.exists(), f"{expected_pdb} does not exist"
        assert expected_gz.exists(), f"{expected_gz} does not exist"

        assert expected_inp.stat().st_size > 0, f"{expected_inp} is empty"
        assert expected_psf.stat().st_size > 0, f"{expected_psf} is empty"
        assert expected_pdb.stat().st_size > 0, f"{expected_pdb} is empty"
        assert expected_gz.stat().st_size > 0, f"{expected_gz} is empty"


def test_topoaa_cyclic(topoaa_module_with_peptide):
    """Test the topoaa module to generate a cyclic peptide."""

    topoaa_module_with_peptide.params["cyclicpept_dist"] = 3.5
    topoaa_module_with_peptide.params["disulphide_dist"] = 4.0
    topoaa_module_with_peptide.params["mol1"] = {"cyclicpept": True}

    topoaa_module_with_peptide.run()

    expected_inp = Path(topoaa_module_with_peptide.path, "cyclic-peptide.inp")
    expected_psf = Path(topoaa_module_with_peptide.path, "cyclic-peptide_haddock.psf")
    expected_pdb = Path(topoaa_module_with_peptide.path, "cyclic-peptide_haddock.pdb")
    expected_gz = Path(topoaa_module_with_peptide.path, "cyclic-peptide.out.gz")

    assert expected_inp.exists(), f"{expected_inp} does not exist"
    assert expected_psf.exists(), f"{expected_psf} does not exist"
    assert expected_gz.exists(), f"{expected_gz} does not exist"
    assert expected_pdb.exists(), f"{expected_pdb} does not exist"

    with open(expected_psf, encoding="utf-8", mode="r") as f:
        file_content = f.read()

    assert "detected" in file_content
    assert "disulphide" in file_content


def test_topoaa_ligand(topoaa_module_with_ligand, ligand_toppar):
    """Test the topoaa module to generate topologies of a ligand"""

    ligand_param_f, ligand_top_f = ligand_toppar

    topoaa_module_with_ligand.params["ligand_param_fname"] = ligand_param_f
    topoaa_module_with_ligand.params["ligand_top_fname"] = ligand_top_f
    topoaa_module_with_ligand.params["delenph"] = False

    topoaa_module_with_ligand.run()

    expected_inp = Path(topoaa_module_with_ligand.path, "oseltamivir.inp")
    expected_psf = Path(topoaa_module_with_ligand.path, "oseltamivir_haddock.psf")
    expected_pdb = Path(topoaa_module_with_ligand.path, "oseltamivir_haddock.pdb")
    expected_gz = Path(topoaa_module_with_ligand.path, "oseltamivir.out.gz")

    assert expected_inp.exists()
    assert expected_psf.exists()
    assert expected_pdb.exists()
    assert expected_gz.exists()
