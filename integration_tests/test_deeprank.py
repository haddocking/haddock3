"""Integration test for the deeprank scoring module."""

import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.scoring.deeprank import DEFAULT_CONFIG as DEEPRANK_CONF
from haddock.modules.scoring.deeprank import HaddockModule as DeeprankModule

from integration_tests import GOLDEN_DATA, has_deeprank


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize=False):
        target_models = ["protprot_complex_1.pdb", "protprot_complex_2.pdb"]
        model_list = []
        for pdb in target_models:
            src = GOLDEN_DATA / pdb
            dst = Path(self.path, src.name)
            shutil.copy(src, dst)
            model_list.append(PDBFile(file_name=str(dst), path=self.path))

        return model_list


@pytest.fixture
def deeprank_module():
    with tempfile.TemporaryDirectory() as tmpdir:
        module = DeeprankModule(
            order=0,
            path=Path(tmpdir),
            init_params=DEEPRANK_CONF,
        )
        module.params["ncores"] = 1
        module.params["chain_i"] = "A"
        module.params["chain_j"] = "B"
        yield module


@has_deeprank
def test_deeprank_run(deeprank_module, mocker):
    deeprank_module.previous_io = MockPreviousIO(path=deeprank_module.path)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )

    deeprank_module.run()

    assert len(deeprank_module.output_models) == 2
    model1 = deeprank_module.output_models[0]
    assert model1.score is not None
    assert isinstance(model1.score, float)

    model2 = deeprank_module.output_models[1]
    assert model2.score is not None
    assert isinstance(model2.score, float)

    # FIXME: Make sure the results are consistent
    # assert model1.score == pytest.approx(0.119)
    # assert model2.score == pytest.approx(0.081)
