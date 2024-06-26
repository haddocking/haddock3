"""Test the sasascore module."""
from haddock.libs.libontology import PDBFile
from haddock.modules.scoring.sasascore.sasascore import (
    AccScore,
    calc_acc_score
    )
import pytest

import tempfile
import pandas as pd
from pathlib import Path

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
    """Provide example buried_resdic."""
    return {
        "A": [43],
        "B": [35, 36]  # dna will always be accessible
        }


def test_accscore(scoring_models, buried_resdic, acc_resdic):
    """Test emscoring expected output."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        output_name = Path(tmpdir, "sasascore")
        viol_output_name = Path(tmpdir, "violations")
        accscore = AccScore(
            model_list=scoring_models,
            output_name=output_name,
            core=0,
            path=Path("."),
            buried_resdic=buried_resdic,
            acc_resdic=acc_resdic,
            cutoff=0.15,
            viol_output_name=viol_output_name,
            probe_radius=1.4,
            )
        accscore.run()
        accscore.output()
        assert output_name.exists()
        exp_scores = pd.Series([1, 2], name="score")
        obs_out_df = pd.read_csv(output_name, sep="\t", header=None)
        assert obs_out_df.iloc[:, -1].equals(exp_scores)
        # now the violations file
        assert viol_output_name.exists()
        exp_viol = pd.Series(["17", "17,29"], name="viols")
        obs_viol_df = pd.read_csv(viol_output_name, sep="\t", header=None)
        assert obs_viol_df.iloc[:, 1].equals(exp_viol)


def test_get_calc_acc_score(buried_resdic, acc_resdic):
    """Test calc_acc_score."""
    result_dict = {
        'A':
        {
            17: {'side_chain_rel': 10.98, 'main_chain_rel': 0.61},
            29: {'side_chain_rel': 6.08, 'main_chain_rel': 3.5},
            43: {'side_chain_rel': 69.79, 'main_chain_rel': 73.01}
            },
        'B':
            {
                35: {'side_chain_rel': 18095.95, 'main_chain_rel': 663},
                36: {'side_chain_rel': 16510.00, 'main_chain_rel': 679}
            }
        }
    acc_score, b_viols, a_viols = calc_acc_score(result_dict,
                                                 buried_resdic,
                                                 acc_resdic)
    assert acc_score == 0
    assert b_viols == {"A": set()}
    assert a_viols == {"A": set(), "B": set()}
