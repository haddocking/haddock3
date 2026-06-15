"""Integration test for the gdock sampling module."""

import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.sampling.gdock import DEFAULT_CONFIG as GDOCK_CONF
from haddock.modules.sampling.gdock import HaddockModule as GdockModule

from integration_tests import GOLDEN_DATA, has_gdock


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize=False):
        src = GOLDEN_DATA / "protprot_complex_1.pdb"
        dst = Path(self.path, src.name)
        shutil.copy(src, dst)
        return [PDBFile(file_name=dst.name, path=self.path)]


@pytest.fixture
def gdock_module():
    with tempfile.TemporaryDirectory() as tmpdir:
        module = GdockModule(
            order=0,
            path=Path(tmpdir),
            initial_params=GDOCK_CONF,
        )
        module.params["receptor_chains"] = "A"
        module.params["ligand_chains"] = "B"
        module.params["max_generations"] = 2
        module.params["seed"] = 1
        module.params["sampling"] = 2
        yield module


@has_gdock
def test_gdock_run(gdock_module, mocker):
    gdock_module.previous_io = MockPreviousIO(path=gdock_module.path)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )

    gdock_module.run()

    assert len(gdock_module.output_models) == gdock_module.params["sampling"]

    for model in gdock_module.output_models:
        assert model.score is not None
        assert isinstance(model.score, float)
        assert model.unw_energies is not None
        for key in ("vdw", "elec", "desolv", "air"):
            assert key in model.unw_energies

        model_path = Path(model.path, model.file_name)
        assert model_path.exists()
        content = model_path.read_text()
        assert "ATOM" in content
