"""Test module prodigy-protein."""

import pytest

from pathlib import Path

from haddock.core.typing import Generator
from haddock.libs.libontology import PDBFile
from haddock.modules.scoring.prodigyprotein import ProdigyJob

from . import golden_data


@pytest.fixture(name="prodigyprot")
def fixture_prodigyprot() -> Generator[ProdigyJob, None, None]:
    yield ProdigyJob(
        model=PDBFile(
            Path(golden_data, "protprot_1bkd_1.pdb"),
            path=golden_data,
            ),
        params={
            "chains": ["R", "S"],
            "distance_cutoff": 5.5,
            "accessibility_cutoff": 0.05,
            "temperature": 25.0,
            "to_pkd": True,
            },
        index=1,
        )


def test_prodigyprot_output_pKd(prodigyprot):
    """Test prodigy pKd values."""
    output = prodigyprot.run()
    assert not prodigyprot.score.error
    assert prodigyprot.score.index == 1
    assert prodigyprot.score.score == pytest.approx(-10.39, 0.01)
    assert not output.error
    assert output.index == 1
    assert output.score == pytest.approx(-10.39, 0.01)


def test_prodigyprot_output_temperature(prodigyprot):
    """Test prodigy pKd values when modifying the temperature."""
    prodigyprot.worker.temperature = 36
    output = prodigyprot.run()
    assert not prodigyprot.score.error
    assert prodigyprot.score.index == 1
    assert prodigyprot.score.score == pytest.approx(-10.02, 0.01)
    assert not output.error
    assert output.index == 1
    assert output.score == pytest.approx(-10.02, 0.01)


def test_prodigyprot_output_DeltaG(prodigyprot):
    """Test prodigy output as DeltaG."""
    prodigyprot.worker.topKd = False
    output = prodigyprot.run()
    assert not prodigyprot.score.error
    assert prodigyprot.score.index == 1
    assert prodigyprot.score.score == pytest.approx(-14.169, 0.01)
    assert not output.error
    assert output.index == 1
    assert output.score == pytest.approx(-14.169, 0.01)


def test_prodigyprot_output_wrong_chains(prodigyprot):
    """Test prodigy output when wrong chains are selected."""
    prodigyprot.worker.chains = ["A", "B"]
    output = prodigyprot.run()
    assert prodigyprot.score.index == 1
    assert not prodigyprot.score.score
    assert prodigyprot.score.error
    assert output.index == 1
    assert not output.score
    assert output.error