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
    """Mocking class to retrieve models."""

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
    """Set of PDBfiles to be filtered."""
    models: list[PDBFile] = []
    # Create 5 models with 5 different scores
    for i in range(5, -5, -1):  # 5, 4, 3, 2, 1, 0, -1, -2, -3, -4
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


def test_filter_default_parameters(filter, mocker, scored_models):
    """Test normal behavior of the filtering module."""
    # Change parameter
    threshold = 0.0
    filter.params["threshold"] = threshold
    # Mock some functions
    filter.previous_io = MockPreviousIO(scored_models)
    mocker.patch(
        "haddock.modules.BaseHaddockModule.export_io_models",
        return_value=None,
    )
    filter.run()
    # Check that models with positive values were filtered out
    assert len(filter.output_models) == 5
    # Check that their values is under the threshold
    for m in filter.output_models:
        assert m.score <= threshold


def test_filter_cutoff_too_strigent(filter, mocker, scored_models):
    """Test expected error when filtering cutoff is too strigent.
    
    Here we expect that if no models are passed to the next module,
    we must terminate the run as nothing can be performed afterwards.
    """
    # Change parameter
    filter.params["threshold"] = -1000
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
