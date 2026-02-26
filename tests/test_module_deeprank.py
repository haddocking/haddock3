"""Tests for the deeprank scoring module wrapper."""

import sys
import tempfile
import shutil
from unittest.mock import MagicMock, patch
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.scoring.deeprank.deeprank import DeeprankWraper

from . import golden_data as GOLDEN_DATA, has_deeprank


@pytest.fixture
def deeprank_wrapper():
    with tempfile.TemporaryDirectory() as temp_dir:
        src = GOLDEN_DATA / "protprot_complex_1.pdb"
        dst = Path(temp_dir, src.name)
        shutil.copy(src, dst)
        yield DeeprankWraper(
            models=[dst],
            ncores=1,
            chain_i="A",
            chain_j="B",
            path=temp_dir,
        )


@has_deeprank
def test_run(deeprank_wrapper):
    """Test the execution method of the wrapper."""
    deeprank_wrapper.run()

    model = deeprank_wrapper.models[0]
    expected_csv = (
        Path(deeprank_wrapper.path)
        / f"{Path(model).stem}-gnn_esm_pred_{deeprank_wrapper.chain_i}_{deeprank_wrapper.chain_j}"
        / "GNN_esm_prediction.csv"
    )

    # Check if the results folders were created
    assert expected_csv.exists()


@has_deeprank
def test_retrieve_scores(deeprank_wrapper):
    """Check the method that retrieves the scores."""
    deeprank_wrapper.run()
    scores = deeprank_wrapper.retrieve_scores()

    assert len(scores) == len(deeprank_wrapper.models)
    for model in deeprank_wrapper.models:
        assert str(model) in scores
        assert isinstance(scores[str(model)], float)
