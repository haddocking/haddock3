"""Test the emscoring module."""

import os
import tempfile
from pathlib import Path

import pytest

from haddock.modules.scoring.emscoring import \
    DEFAULT_CONFIG as DEFAULT_EMSCORING_PARAMS
from haddock.modules.scoring.emscoring import HaddockModule as EMScoring


@pytest.fixture(name="emscoring")
def fixture_emscoring():
    """???"""

    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        yield EMScoring(
            order=1, path=Path("."), initial_params=DEFAULT_EMSCORING_PARAMS
        )


@pytest.fixture(name="emscoring_dna")
def fixture_emscoring_dna(emscoring, protdna_input_list):
    """???"""

    protdna_input_list[0].score = -28.0
    protdna_input_list[0].ori_name = "original_name_1.pdb"

    protdna_input_list[1].score = 42.0
    protdna_input_list[1].ori_name = "original_name_2.pdb"

    emscoring.output_models = protdna_input_list

    yield emscoring


def test_emscoring_output(emscoring_dna, protdna_input_list):
    """Test emscoring expected output."""

    output_fname = Path("emscoring.tsv")
    emscoring_dna.output(output_fname)

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
