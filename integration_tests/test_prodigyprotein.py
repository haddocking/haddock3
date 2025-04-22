"""Integration tests for module prodigyprotein."""

import shutil
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import PDBFile
from haddock.modules.scoring.prodigyprotein import (
    DEFAULT_CONFIG as DEFAULT_PRODIGYPROT_CONFIG,
    HaddockModule as ProdigyModule,
    )

from integration_tests import GOLDEN_DATA


@pytest.fixture(name="prodigyprot")
def fixture_prodigyprot():
    with tempfile.TemporaryDirectory() as tmpdir:
        prodigyprot = ProdigyModule(
            order=0,
            path=Path(tmpdir),
            initial_params=DEFAULT_PRODIGYPROT_CONFIG,
        )
        yield prodigyprot


class MockPreviousIO:
    def __init__(self, path):
        self.path = path

    def retrieve_models(self, individualize = True):
        shutil.copy(
            Path(GOLDEN_DATA, "2oob.pdb"),
            Path(self.path, "2oob.pdb"),
        )

        model = PDBFile(
            file_name="2oob.pdb",
            path=self.path,
        )

        return [model]

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


def test_prodigyprotein_defaults(prodigyprot):
    """Integration test with default parameters."""
    prodigyprot.previous_io = MockPreviousIO(path=prodigyprot.path)
    prodigyprot.run()
    # Build expected output file path
    expected_output_filepath = Path(prodigyprot.path, "prodigyprotein.tsv")
    # Check that the file was created
    assert expected_output_filepath.exists()
    assert expected_output_filepath.stat().st_size != 0
    # Check that the scores are correct
    scores = read_scoring_tsv(expected_output_filepath)
    assert float(scores[0][-1]) == pytest.approx(-4.75, abs=0.1)


def test_prodigyprotein_modified_temperature(prodigyprot):
    """Integration test with modified temperature parameter."""
    prodigyprot.previous_io = MockPreviousIO(path=prodigyprot.path)
    prodigyprot.params["temperature"] = 37.0
    prodigyprot.run()
    # Build expected output file path
    expected_output_filepath = Path(prodigyprot.path, "prodigyprotein.tsv")
    # Check that the file was created
    assert expected_output_filepath.exists()
    assert expected_output_filepath.stat().st_size != 0
    # Check that the scores are correct
    scores = read_scoring_tsv(expected_output_filepath)
    assert float(scores[0][-1]) == pytest.approx(-4.57, abs=0.1)


def test_prodigyprotein_deltaG(prodigyprot):
    """Integration test with deltaG output."""
    prodigyprot.previous_io = MockPreviousIO(path=prodigyprot.path)
    prodigyprot.params["to_pkd"] = False
    prodigyprot.run()
    # Build expected output file path
    expected_output_filepath = Path(prodigyprot.path, "prodigyprotein.tsv")
    # Check that the file was created
    assert expected_output_filepath.exists()
    assert expected_output_filepath.stat().st_size != 0
    # Check that the scores are correct
    scores = read_scoring_tsv(expected_output_filepath)
    assert float(scores[0][-1]) == pytest.approx(-6.482, abs=0.1)


def test_prodigyprotein_modified_distance_cutoff(prodigyprot):
    """Integration test with modified distance_cutoff parameter."""
    prodigyprot.previous_io = MockPreviousIO(path=prodigyprot.path)
    prodigyprot.params["distance_cutoff"] = 4.0
    prodigyprot.run()
    # Build expected output file path
    expected_output_filepath = Path(prodigyprot.path, "prodigyprotein.tsv")
    # Check that the file was created
    assert expected_output_filepath.exists()
    assert expected_output_filepath.stat().st_size != 0
    # Check that the scores are correct
    scores = read_scoring_tsv(expected_output_filepath)
    assert float(scores[0][-1]) == pytest.approx(-4.08, abs=0.1)


def test_prodigyprotein_modified_accessibility(prodigyprot):
    """Integration test with modified accessibility parameter."""
    prodigyprot.previous_io = MockPreviousIO(path=prodigyprot.path)
    prodigyprot.params["accessibility_cutoff"] = 0.1
    prodigyprot.run()
    # Build expected output file path
    expected_output_filepath = Path(prodigyprot.path, "prodigyprotein.tsv")
    # Check that the file was created
    assert expected_output_filepath.exists()
    assert expected_output_filepath.stat().st_size != 0
    # Check that the scores are correct
    scores = read_scoring_tsv(expected_output_filepath)
    assert float(scores[0][-1]) == pytest.approx(-5.08, abs=0.1)
