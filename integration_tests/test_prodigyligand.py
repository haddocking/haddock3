"""Integration tests for module prodigyligand."""

import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.scoring.prodigyligand import (
    DEFAULT_CONFIG as DEFAULT_prodigylig_CONFIG,
    HaddockModule as ProdigyModule,
    )

from integration_tests import GOLDEN_DATA


@pytest.fixture(name="prodigylig")
def fixture_prodigylig():
    with tempfile.TemporaryDirectory() as tmpdir:
        prodigylig = ProdigyModule(
            order=0,
            path=Path(tmpdir),
            initial_params=DEFAULT_prodigylig_CONFIG,
        )
        # Fix parameters for chains and ligand resname
        prodigylig.params["receptor_chain"] = "A"
        prodigylig.params["ligand_chain"] = "B"
        prodigylig.params["ligand_resname"] = "NDG"
        yield prodigylig


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize = True):
        shutil.copy(
            Path(GOLDEN_DATA, "protglyc_complex_1.pdb"),
            Path(self.path, "protglyc_complex_1.pdb"),
        )
        shutil.copy(
            Path(GOLDEN_DATA, "protglyc_complex_2.pdb"),
            Path(self.path, "protglyc_complex_2.pdb"),
        )

        models = [
            PDBFile(
                file_name="protglyc_complex_1.pdb",
                path=self.path,
            ),
            PDBFile(
                file_name="protglyc_complex_2.pdb",
                path=self.path,
            ),
        ]

        return models

    def output(self):
        return None


def read_scoring_tsv(fpath) -> list[list[str]]:
    """Load tsv as string skipping header."""
    with open(fpath, "r") as fin:
        return [
            [e.strip() for e in _.split("\t")]
            for i, _ in enumerate(fin)
            if i != 0
            ]


def test_prodigyligand_defaults(prodigylig):
    """Integration test with default parameters."""
    prodigylig.previous_io = MockPreviousIO(path=prodigylig.path)
    prodigylig.run()
    # Build expected output file path
    expected_output_filepath = Path(prodigylig.path, "prodigyligand.tsv")
    # Check that the file was created
    assert expected_output_filepath.exists()
    assert expected_output_filepath.stat().st_size != 0
    # Check that the scores are correct
    scores = read_scoring_tsv(expected_output_filepath)
    assert len(scores) == 2
    assert float(scores[0][-1]) == pytest.approx(-4.62, abs=0.1)
    assert float(scores[1][-1]) == pytest.approx(-4.42, abs=0.1)


def test_prodigyligand_modified_temperature(prodigylig):
    """Integration test with modified temperature parameter."""
    prodigylig.previous_io = MockPreviousIO(path=prodigylig.path)
    prodigylig.params["temperature"] = 37.0
    prodigylig.run()
    # Build expected output file path
    expected_output_filepath = Path(prodigylig.path, "prodigyligand.tsv")
    # Check that the file was created
    assert expected_output_filepath.exists()
    assert expected_output_filepath.stat().st_size != 0
    # Check that the scores are correct
    scores = read_scoring_tsv(expected_output_filepath)
    assert len(scores) == 2
    assert float(scores[0][-1]) == pytest.approx(-4.44, abs=0.1)
    assert float(scores[1][-1]) == pytest.approx(-4.25, abs=0.1)


def test_prodigyligand_deltaG(prodigylig):
    """Integration test with deltaG output."""
    prodigylig.previous_io = MockPreviousIO(path=prodigylig.path)
    prodigylig.params["to_pkd"] = False
    prodigylig.run()
    # Build expected output file path
    expected_output_filepath = Path(prodigylig.path, "prodigyligand.tsv")
    # Check that the file was created
    assert expected_output_filepath.exists()
    assert expected_output_filepath.stat().st_size != 0
    # Check that the scores are correct
    scores = read_scoring_tsv(expected_output_filepath)
    assert len(scores) == 2
    assert float(scores[0][-1]) == pytest.approx(-6.296, abs=0.1)
    assert float(scores[1][-1]) == pytest.approx(-6.027, abs=0.1)


def test_prodigyligand_modified_distance_cutoff(prodigylig):
    """Integration test with modified distance_cutoff parameter."""
    prodigylig.previous_io = MockPreviousIO(path=prodigylig.path)
    prodigylig.params["distance_cutoff"] = 7.0
    prodigylig.run()
    # Build expected output file path
    expected_output_filepath = Path(prodigylig.path, "prodigyligand.tsv")
    # Check that the file was created
    assert expected_output_filepath.exists()
    assert expected_output_filepath.stat().st_size != 0
    # Check that the scores are correct
    scores = read_scoring_tsv(expected_output_filepath)
    assert len(scores) == 2
    assert float(scores[0][-1]) == pytest.approx(-3.96, abs=0.1)
    assert float(scores[1][-1]) == pytest.approx(-3.9, abs=0.1)


def test_prodigyligand_no_electrostatics(prodigylig):
    """Integration test without the use of electrostatics information."""
    prodigylig.previous_io = MockPreviousIO(path=prodigylig.path)
    prodigylig.params["electrostatics"] = False
    prodigylig.run()
    # Build expected output file path
    expected_output_filepath = Path(prodigylig.path, "prodigyligand.tsv")
    # Check that the file was created
    assert expected_output_filepath.exists()
    assert expected_output_filepath.stat().st_size != 0
    # Check that the scores are correct
    scores = read_scoring_tsv(expected_output_filepath)
    assert len(scores) == 2
    assert float(scores[0][-1]) == pytest.approx(-5.0, abs=0.1)
    assert float(scores[1][-1]) == pytest.approx(-4.63, abs=0.1)
