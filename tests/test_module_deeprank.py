"""Tests for the deeprank scoring module wrapper."""

import tempfile
import shutil
from pathlib import Path

import pytest

from haddock.modules.scoring.deeprank.deeprank import DeeprankWrapper

from . import golden_data as GOLDEN_DATA
from .conftest import has_deeprank


@pytest.fixture
def deeprank_wrapper():
    with tempfile.TemporaryDirectory() as temp_dir:
        src = GOLDEN_DATA / "protprot_complex_1.pdb"
        dst = Path(temp_dir, src.name)
        shutil.copy(src, dst)
        yield DeeprankWrapper(
            models=[dst],
            ncores=1,
            chain_i="A",
            chain_j="B",
        )


@pytest.fixture
def deeprank_wrapper_ensemble():
    with tempfile.TemporaryDirectory() as temp_dir:
        models = []
        for pdb in ("protprot_complex_1.pdb", "protprot_complex_2.pdb"):
            src = GOLDEN_DATA / pdb
            dst = Path(temp_dir, src.name)
            shutil.copy(src, dst)
            models.append(dst)
        yield DeeprankWrapper(
            models=models,
            ncores=2,
            chain_i="A",
            chain_j="B",
        )


@has_deeprank
def test_run(deeprank_wrapper):
    """`run()` must return the score for the single input model."""
    scores = deeprank_wrapper.run()

    model = deeprank_wrapper.models[0]
    assert str(model) in scores
    assert isinstance(scores[str(model)], float)


@has_deeprank
def test_run_leaves_no_files_behind(deeprank_wrapper):
    """Deeprank's intermediate files must not leak outside the temp workspace."""
    model = deeprank_wrapper.models[0]
    model_dir = Path(model).parent

    deeprank_wrapper.run()

    assert list(model_dir.iterdir()) == [model]


@has_deeprank
def test_run_ensemble(deeprank_wrapper_ensemble):
    """Scores from a multi-model run must map back to the correct model."""
    scores = deeprank_wrapper_ensemble.run()

    assert len(scores) == 2
    model1, model2 = deeprank_wrapper_ensemble.models
    assert scores[str(model1)] == pytest.approx(0.147, abs=1e-3)
    assert scores[str(model2)] == pytest.approx(0.102, abs=1e-3)
