"""Test related to the [seletop] module."""

import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.seletop import DEFAULT_CONFIG
from haddock.modules.analysis.seletop import HaddockModule as SeletopModule

from . import golden_data


class MockPreviousIO:
    """Mocking class to retrieve models."""

    def __init__(self, models):
        self.models = models
        self.output = models

    def retrieve_models(self, individualize: bool = False):
        """Provide a set of models."""
        return self.models


@pytest.fixture(name="scored_models")
def fixture_scored_models() -> list[PDBFile]:
    """Set of PDBfiles to be selected."""
    models: list[PDBFile] = []
    # Create 5 models with 5 different scores
    for i in range(5, 0, -1):  # 5, 4, 3, 2, 1
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_f:
            dst = temp_f.name
            src = Path(golden_data, "protprot_complex_1.pdb")
            shutil.copy(src, dst)
            pdbfile = PDBFile(file_name=dst, path=Path(dst).parent)
            pdbfile.score = i
            models.append(pdbfile)
    return models


@pytest.fixture(name="seletop")
def fixture_seletop(monkeypatch):
    """Test module __init__()."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        yield SeletopModule(
            order=1,
            path=Path("."),
            initial_params=DEFAULT_CONFIG,
        )


def test_confirm_installation(seletop):
    """Test confirm install."""
    assert seletop.confirm_installation() is None


def test_seletop_default_sort(seletop, mocker, scored_models):
    """By default the models with the lowest score are selected first."""
    seletop.params["select"] = 3
    seletop.previous_io = MockPreviousIO(scored_models)
    mocker.patch(
        "haddock.modules.analysis.seletop.HaddockModule.previous_path",
        return_value="0_rigidbody",
    )
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )
    seletop.run()
    # Lowest scores selected: 1, 2, 3
    assert [m.score for m in seletop.output_models] == [1, 2, 3]


def test_seletop_reversed_sort(seletop, mocker, scored_models):
    """With sort_ascending disabled the highest scores are selected first."""
    seletop.params["select"] = 3
    seletop.params["sort_ascending"] = False
    seletop.previous_io = MockPreviousIO(scored_models)
    mocker.patch(
        "haddock.modules.analysis.seletop.HaddockModule.previous_path",
        return_value="0_rigidbody",
    )
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )
    seletop.run()
    # Highest scores selected: 5, 4, 3
    assert [m.score for m in seletop.output_models] == [5, 4, 3]
