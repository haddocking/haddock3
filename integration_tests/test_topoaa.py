import tempfile
from pathlib import Path

import pytest

from haddock.modules.topology.topoaa import \
    DEFAULT_CONFIG as DEFAULT_TOPOAA_CONFIG
from haddock.modules.topology.topoaa import HaddockModule as TopoaaModule

from . import CNS_EXEC
from integration_tests import GOLDEN_DATA


@pytest.fixture
def topoaa_module():
    with tempfile.TemporaryDirectory() as tmpdir:
        topoaa = TopoaaModule(
            order=0, path=Path(tmpdir), initial_params=DEFAULT_TOPOAA_CONFIG
        )
        yield topoaa


def test_topoaa_module_protein(topoaa_module):
    """Topoaa module with protein-protein input"""
    topoaa_module.params["molecules"] = [
        Path(GOLDEN_DATA, "e2aP_1F3G.pdb"),
        Path(GOLDEN_DATA, "hpr_ensemble.pdb"),
    ]
    topoaa_module.params["mol1"] = {"prot_segid": "A"}
    topoaa_module.params["mol2"] = {"prot_segid": "B"}
    topoaa_module.params["cns_exec"] = CNS_EXEC
    topoaa_module.params["debug"] = True

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


def test_topoaa_module_ligand(topoaa_module):
    """Topoaa with ligand as input"""
    topoaa_module.params["molecules"] = [
        Path(GOLDEN_DATA, "oseltamivir.pdb"),
    ]
    topoaa_module.params["ligand_top_fname"] = Path(GOLDEN_DATA, "ligand.top")
    topoaa_module.params["ligand_param_fname"] = Path(GOLDEN_DATA, "ligand.param")
    topoaa_module.params["delenph"] = False
    topoaa_module.params["preprocess"] = False
    topoaa_module.params["cns_exec"] = CNS_EXEC
    topoaa_module.params["debug"] = True

    topoaa_module.run()

    expected_inp = Path(topoaa_module.path, "oseltamivir.inp")
    expected_psf = Path(topoaa_module.path, "oseltamivir_haddock.psf")
    expected_pdb = Path(topoaa_module.path, "oseltamivir_haddock.pdb")
    expected_gz = Path(topoaa_module.path, "oseltamivir.out.gz")

    assert expected_inp.exists()
    assert expected_psf.exists()
    assert expected_pdb.exists()
    assert expected_gz.exists()


def test_topoaa_cyclic(topoaa_module):
    """Test the topoaa module to generate a cyclic peptide."""
    topoaa_module.params["molecules"] = [
        Path(GOLDEN_DATA, "cyclic-peptide.pdb"),
    ]
    topoaa_module.params["cyclicpept_dist"] = 3.5
    topoaa_module.params["disulphide_dist"] = 4.0
    topoaa_module.params["mol1"] = {"cyclicpept": True}
    topoaa_module.params["cns_exec"] = CNS_EXEC
    topoaa_module.params["debug"] = True

    topoaa_module.run()

    expected_inp = Path(topoaa_module.path, "cyclic-peptide.inp")
    expected_psf = Path(topoaa_module.path, "cyclic-peptide_haddock.psf")
    expected_pdb = Path(topoaa_module.path, "cyclic-peptide_haddock.pdb")
    expected_gz = Path(topoaa_module.path, "cyclic-peptide.out.gz")

    assert expected_inp.exists(), f"{expected_inp} does not exist"
    assert expected_psf.exists(), f"{expected_psf} does not exist"
    assert expected_gz.exists(), f"{expected_gz} does not exist"
    assert expected_pdb.exists(), f"{expected_pdb} does not exist"

    with open(expected_psf, encoding="utf-8", mode="r") as f:
        file_content = f.read()

    assert "detected" in file_content
    assert "disulphide" in file_content
