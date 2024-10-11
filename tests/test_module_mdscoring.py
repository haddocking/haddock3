"""Test the mdscoring module."""

import os
import tempfile
from pathlib import Path

import pytest

from haddock.modules.scoring.mdscoring import \
    DEFAULT_CONFIG as MDSCORING_DEFAULT_PARAMS
from haddock.modules.scoring.mdscoring import HaddockModule as MDScoring


@pytest.fixture(name="mdscoring")
def fixture_mdscoring():
    """mdscoring module fixture"""

    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        yield MDScoring(
            order=1, path=Path("."), initial_params=MDSCORING_DEFAULT_PARAMS
        )


@pytest.fixture(name="mdscoring_dna")
def fixture_mdscoring_dna(mdscoring, protdna_input_list):
    """mdscoring module fixture protdna"""

    protdna_input_list[0].score = -28.0
    protdna_input_list[0].ori_name = "original_name_1.pdb"

    protdna_input_list[1].score = 42.0
    protdna_input_list[1].ori_name = "original_name_2.pdb"

    mdscoring.output_models = protdna_input_list

    yield mdscoring


def test_mdscoring_output(mdscoring_dna, protdna_input_list):
    """Test mdscoring expected output."""

    output_fname = Path("mdscoring.tsv")
    mdscoring_dna.output(output_fname)

    with open(output_fname, "r", encoding="utf-8") as fh:
        observed_outf_l = [l.split() for l in fh.readlines() if not l.startswith("#")]

    # expected output
    expected_outf_l = [
        ["structure", "original_name", "md5", "score"],
        [
            Path(protdna_input_list[0].rel_path).name,
            Path(protdna_input_list[0].ori_name).name,
            "None",
            str(protdna_input_list[0].score),
        ],
        [
            Path(protdna_input_list[1].rel_path).name,
            Path(protdna_input_list[1].ori_name).name,
            "None",
            str(protdna_input_list[1].score),
        ],
    ]

    assert len(observed_outf_l) == len(expected_outf_l)

    for i, j in zip(observed_outf_l, expected_outf_l):
        assert i == j

    assert observed_outf_l == expected_outf_l

    os.unlink(output_fname)
