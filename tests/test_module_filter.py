"""Test related to the [filter] module."""

import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.analysis.filter import DEFAULT_CONFIG
from haddock.modules.analysis.filter import HaddockModule as FilterModule

from . import golden_data


class MockPreviousIO:
    """A mocking class holding specific methods."""

    def __init__(self, models):
        self.models = models
        self.output = models

    # In the mocked method, add the arguments that are called by the original
    #  method that is being tested
    def retrieve_models(self, individualize: bool = False):
        """Provide a set of models."""
        return self.models


@pytest.fixture(name="scored_models")
def fixture_scored_models() -> list[PDBFile]:
    """Set of clustered PDBfiles."""
    models: list[PDBFile] = []
    # Create 5 models with 5 different scores
    for i in range(-10, -15, -1):  # -10, -11, -12, -13, -14
        with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_f:
            dst = temp_f.name
            src = Path(golden_data, "protprot_complex_1.pdb")
            shutil.copy(src, dst)
            pdbfile = PDBFile(file_name=dst, path=Path(dst).parent)
            # Add attributs
            pdbfile.score = i
            models.append(pdbfile)
    return models


@pytest.fixture(name="filter")
def fixture_filter(monkeypatch):
    """Test module __init__()."""
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        yield FilterModule(
            order=1,
            path=Path("."),
            initial_params=DEFAULT_CONFIG,
        )


def test_confirm_installation(filter):
    """Test confirm install."""
    assert filter.confirm_installation() is None


def test_init(filter):
    """Test __init__ function."""
    filter.__init__(
        order=42,
        path=Path("0_anything"),
        initial_params=DEFAULT_CONFIG,
    )

    # Once a module is initialized, it should have the following attributes
    assert filter.path == Path("0_anything")
    assert filter._origignal_config_file == DEFAULT_CONFIG
    assert type(filter.params) == dict
    assert len(filter.params.keys()) != 0


def test_filter_expected_behavior(filter, mocker, scored_models):
    """Test finish_with_error due to wrong cluster parameter type."""
    # Change parameter
    cutoff = -12
    filter.params["by"] = "score"
    filter.params["cutoff"] = cutoff
    # Mock some functions
    filter.previous_io = MockPreviousIO(scored_models)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )
    filter.run()
    assert len(filter.output_models) == 3
    for m in filter.output_models:
        assert m.score <= cutoff


def test_filter_wrong_by(filter, mocker, scored_models):
    """Test finish_with_error due to wrong nb models parameter value."""
    # Change parameter
    filter.params["by"] = "HopeItsNotThereAndWillNeverBe"
    # Mock some functions
    filter.previous_io = MockPreviousIO(scored_models)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.finish_with_error",
        side_effect=Exception("mocked error"),
    )
    with pytest.raises(Exception) as moked_finish_with_error:
        # run main module _run() function
        filter.run()
    assert moked_finish_with_error.value.__str__() == "mocked error"


def test_filter_cutoff_too_strigent(filter, mocker, scored_models):
    """Test finish_with_error due to wrong cluster parameter type."""
    # Change parameter
    filter.params["by"] = "score"
    filter.params["cutoff"] = -1000
    # Mock some functions
    filter.previous_io = MockPreviousIO(scored_models)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.finish_with_error",
        side_effect=Exception("mocked error"),
    )
    with pytest.raises(Exception) as moked_finish_with_error:
        # run main module _run() function
        filter.run()
    assert moked_finish_with_error.value.__str__() == "mocked error"
