"""Tests for the deeprank scoring module wrapper."""

import builtins
import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.modules.scoring.deeprank.deeprank import (
    DeeprankWrapper,
    deeprank_is_available,
)

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
            score_type="mean",
            chain_i="A",
            chain_j="B",
        )


@pytest.fixture
def deeprank_wrapper_two_models():
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
            score_type="mean",
            chain_i="A",
            chain_j="B",
        )


def test_deeprank_is_available_requires_sqlite3(monkeypatch):
    """Must raise a clear error when the interpreter lacks sqlite3 support."""
    real_import = builtins.__import__

    # mock not being able to find `sqlite3`
    def fake_import(name, *args, **kwargs):
        if name == "sqlite3":
            raise ImportError("No module named '_sqlite3'")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)

    with pytest.raises(ImportError, match="sqlite3"):
        deeprank_is_available()


@has_deeprank
def test_run(deeprank_wrapper, tmp_path, monkeypatch):
    """`run()` must score the model and write only deeprank's csv files."""
    model = deeprank_wrapper.models[0]
    model_dir = Path(model).parent
    monkeypatch.chdir(tmp_path)

    scores = deeprank_wrapper.run()

    assert isinstance(scores[str(model)], float)
    assert list(model_dir.iterdir()) == [model]
    assert sorted(p.name for p in tmp_path.iterdir()) == [
        DeeprankWrapper.INTERFACES_CSV,
        DeeprankWrapper.SUMMARY_CSV,
    ]


@has_deeprank
def test_run_two_models(deeprank_wrapper_two_models, tmp_path, monkeypatch):
    """Scores from a multi-model run must map back to the correct model."""
    monkeypatch.chdir(tmp_path)
    scores = deeprank_wrapper_two_models.run()

    assert len(scores) == 2
    model1, model2 = deeprank_wrapper_two_models.models
    assert scores[str(model1)] == pytest.approx(0.147, abs=1e-3)
    assert scores[str(model2)] == pytest.approx(0.102, abs=1e-3)


@pytest.mark.parametrize(
    "score_type,chain_i,chain_j,expected",
    [
        ("mean", "A", "B", {"model_1.pdb": 0.25, "model_2.pdb": 0.0}),
        ("min", "A", "B", {"model_1.pdb": 0.2, "model_2.pdb": 0.0}),
        ("pair", "C", "A", {"model_1.pdb": 0.3, "model_2.pdb": 0.0}),
    ],
)
def test_retrieve_scores(tmp_path, score_type, chain_i, chain_j, expected):
    """Each `score_type` must read its fnat; no interface (`NA`) gives 0."""
    Path(tmp_path, DeeprankWrapper.INTERFACES_CSV).write_text(
        "pdb_id,chain_i,chain_j,predicted_fnat\n"
        "model_1,A,B,0.200\n"
        "model_1,A,C,0.300\n"
        "model_2,A,B,NA\n"
        "model_2,A,C,NA\n"
    )
    Path(tmp_path, DeeprankWrapper.SUMMARY_CSV).write_text(
        "pdb_id,scored_pairs,mean_fnat,std_fnat,min_fnat\n"
        "model_1,2,0.250,0.050,0.200\n"
        "model_2,0,NA,NA,NA\n"
    )
    wrapper = DeeprankWrapper(
        models=["model_1.pdb", "model_2.pdb"],
        ncores=1,
        score_type=score_type,
        chain_i=chain_i,
        chain_j=chain_j,
    )

    assert wrapper._retrieve_scores(tmp_path) == expected
