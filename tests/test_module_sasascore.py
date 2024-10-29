"""Test the sasascore module."""
import tempfile
import os
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.scoring.sasascore.sasascore import (
    AccScore,
    calc_acc_score,
    )

from . import golden_data


@pytest.fixture
def scoring_models():
    """Prot-DNA models using for emscoring output."""
    return [
        PDBFile(
            Path(golden_data, "protdna_complex_1.pdb"),
            path=golden_data,
            score=42.0
            ),
        PDBFile(
            Path(golden_data, "protdna_complex_2.pdb"),
            path=golden_data,
            score=-28.0
            )]


@pytest.fixture
def buried_resdic():
    """Provide example buried_resdic."""
    return {
        "A": [17, 29],
        }


@pytest.fixture
def acc_resdic():
    """Provide example acc_resdic."""
    return {
        "A": [43],
        "B": [35, 36]  # dna will always be accessible
        }

@pytest.fixture
def acc_resdic_more_chains():
    """Provide example acc_resdic with chains that are not in the model."""
    return {
        "A": [43],
        "B": [35, 36],
        "C": [1, 2, 3],
        "D": [4, 5, 6],
    }


def test_accscore(scoring_models, buried_resdic, acc_resdic):
    """Test accscore expected output."""
    with tempfile.TemporaryDirectory() as tmpdir:
        os.chdir(tmpdir)
        accscore = AccScore(
            model=scoring_models[0],
            path=Path("."),
            buried_resdic=buried_resdic,
            acc_resdic=acc_resdic,
            cutoff=0.1,
            probe_radius=1.4,
            )
        accscore.run()
        exp_score = 1
        obs_score = accscore.data[-1]
        # now the violations
        exp_viol_data = [scoring_models[0].file_name, '29', None, None]
        obs_viol_data = accscore.violations_data
        assert exp_score == obs_score
        assert exp_viol_data == obs_viol_data


def test_get_calc_acc_score(buried_resdic, acc_resdic):
    """Test calc_acc_score."""
    result_dict = {
        'A':
            [17, 18, 19],
        'B':
            [35, 36, 37],
        }
    acc_score, b_viols, a_viols = calc_acc_score(result_dict,
                                                 buried_resdic,
                                                 acc_resdic)
    assert acc_score == 2
    assert b_viols == {"A": set([17])}
    assert a_viols == {"A": set([43]), "B": set()}
