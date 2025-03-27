"""Test module prodigy-protein."""

import pytest
import tempfile

from pathlib import Path

from haddock.core.typing import Generator
from haddock.libs.libontology import PDBFile
from haddock.modules.scoring.prodigyligand import ProdigyJob

from . import golden_data


@pytest.fixture(name="prodigyligand")
def fixture_prodigyligand(monkeypatch) -> Generator[ProdigyJob, None, None]:
    with tempfile.TemporaryDirectory() as tmpdir:
        monkeypatch.chdir(tmpdir)
        yield ProdigyJob(
            model=PDBFile(
                Path(golden_data, "protlig_complex_1.pdb"),
                path=golden_data,
                ),
            params={
                "receptor_chain": "A",
                "ligand_chain": "B",
                "ligand_resname": "G39",
                "distance_cutoff": 10.5,
                "temperature": 25.0,
                "to_pkd": True,
                },
            index=1,
            )


def test_prodigyligand_output_pKd(prodigyligand):
    """Test prodigy-lig pKd values."""
    output = prodigyligand.run()
    assert not prodigyligand.score.error
    assert prodigyligand.score.index == 1
    assert prodigyligand.score.score == pytest.approx(-5.12, 0.01)
    assert not output.error
    assert output.index == 1
    assert output.score == pytest.approx(-5.12, 0.01)


def test_prodigyligand_output_temperature(prodigyligand):
    """Test prodigy-lig pKd values when modifying the temperature."""
    prodigyligand.worker.temperature = 36
    output = prodigyligand.run()
    assert not prodigyligand.score.error
    assert prodigyligand.score.index == 1
    assert prodigyligand.score.score == pytest.approx(-4.94, 0.01)
    assert not output.error
    assert output.index == 1
    assert output.score == pytest.approx(-4.94, 0.01)


def test_prodigyligand_output_DeltaG(prodigyligand):
    """Test prodigy-lig output as DeltaG."""
    prodigyligand.worker.topKd = False
    output = prodigyligand.run()
    assert not prodigyligand.score.error
    assert prodigyligand.score.index == 1
    assert prodigyligand.score.score == pytest.approx(-6.986, 0.01)
    assert not output.error
    assert output.index == 1
    assert output.score == pytest.approx(-6.986, 0.01)


def test_prodigyligand_output_wrong_receptor_chain(prodigyligand):
    """Test prodigy-lig output when wrong receptor chain."""
    # Set wrong receptor chain
    prodigyligand.worker.receptor_chain = "F"
    output = prodigyligand.run()
    assert prodigyligand.score.index == 1
    assert not prodigyligand.score.score
    assert prodigyligand.score.error
    assert output.index == 1
    assert not output.score
    assert output.error


def test_prodigyligand_output_wrong_ligand_chain(prodigyligand):
    """Test prodigy-lig output when wrong ligand chain."""
    # Set wrong ligand chain
    prodigyligand.worker.lig_chain = "F"
    output = prodigyligand.run()
    assert prodigyligand.score.index == 1
    assert not prodigyligand.score.score
    assert prodigyligand.score.error
    assert output.index == 1
    assert not output.score
    assert output.error


def test_prodigyligand_output_wrong_ligand_resname(prodigyligand):
    """Test prodigy-lig output when wrong ligand chain."""
    # Set wrong ligand chain
    prodigyligand.worker.lig_resname = "WRG"
    output = prodigyligand.run()
    assert prodigyligand.score.index == 1
    assert not prodigyligand.score.score
    assert prodigyligand.score.error
    assert output.index == 1
    assert not output.score
    assert output.error
