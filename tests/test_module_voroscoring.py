"""Test the voroscoring module."""
import os
import pytest
import pytest_mock  # noqa : F401
import tempfile
import subprocess
import shutil

from numpy import isnan
from pathlib import Path

from haddock.libs.libontology import NaN, PDBFile
from haddock.modules.scoring.voroscoring import (
    DEFAULT_CONFIG as params,
    HaddockModule as VoroScoringModule,
    )
from haddock.modules.scoring.voroscoring.voroscoring import (
    VoroMQA,
    update_models_with_scores,
    )

from . import golden_data


@pytest.fixture
def output_models():
    """Prot-DNA models using for emscoring output."""
    return [
        PDBFile(
            Path(golden_data, "protdna_complex_1.pdb"),
            path=golden_data,
            score=-0.28,
            ),
        PDBFile(
            Path(golden_data, "protdna_complex_2.pdb"),
            path=golden_data,
            score=-0.42,
            ),
        PDBFile(
            Path(golden_data, "protdna_complex_3.pdb"),
            path=golden_data,
            score=NaN,
            ),
        ]


@pytest.fixture
def voromqa(output_models):
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        voromqa_object = VoroMQA(
            output_models,
            tmpdir,
            params,
            Path("raw_voromqa_scores.tsv"),
            )
        yield voromqa_object


def test_voroscoring_output(output_models):
    """Test voroscoring expected output."""
    voro_module = VoroScoringModule(
        order=1,
        path=Path("1_voroscoring"),
        initial_params=params
        )
    # original names
    voro_module.output_models = output_models
    for mod in range(len(output_models)):
        ori_name = "original_name_" + str(mod) + ".pdb"
        voro_module.output_models[mod].ori_name = ori_name
    # creating output
    output_fname = Path("voroscoring.tsv")
    voro_module.output(output_fname)
    observed_outf_l = [
        e.split()
        for e in open(output_fname).readlines()
        if not e.startswith('#')
        ]
    # expected output
    expected_outf_l = [
        ["structure", "original_name", "md5", "score"],
        ["protdna_complex_2.pdb", "original_name_1.pdb", "None", "-0.420"],
        ["protdna_complex_1.pdb", "original_name_0.pdb", "None", "-0.280"],
        ["protdna_complex_3.pdb", "original_name_2.pdb", "None", "None"],
        ]
        
    assert observed_outf_l == expected_outf_l
    output_fname.unlink()


def test_wait_for_termination(voromqa):
    """Test waiting for results function behavior in voromqa."""
    nested_batch_dir = Path(voromqa.workdir, "batch_1")
    os.mkdir(nested_batch_dir)
    expected_ssv = Path(nested_batch_dir, "voro_scores.ssv")
    # Trick to fake the generation of a file
    delay_scriptpath = Path(nested_batch_dir, "delay.sh")
    delay_scriptpath.write_text(
        "\n".join(["sleep 0.1", f'echo "haddock3" > {expected_ssv}'])
        )
    assert delay_scriptpath.exists()
    os.system(f"chmod u+x {delay_scriptpath}")
    os.system(f"./{delay_scriptpath} &")
    assert not expected_ssv.exists()
    # The actual test of the function
    batches_ssv = voromqa.wait_for_termination(wait_time=0.1)
    assert expected_ssv.exists()
    assert batches_ssv[0] == expected_ssv
    shutil.rmtree(nested_batch_dir)


def test_batched(voromqa):
    """Test batched function behavior in voromqa."""
    for batch in voromqa.batched(list(range(10)), size=2):
        assert len(batch) == 2
    batches = list(voromqa.batched(list(range(100)), size=99))
    assert len(batches[0]) == 99
    assert len(batches[1]) == 1


def test_update_models_with_scores(output_models):
    """Test to update PDBFiles with scores from voromqa tsv."""
    # Generate fake voro output file
    output_fname = Path("fake_voro.tsv")
    output_fname.write_text(
        """ID\tjury_score\tfake_energy
protdna_complex_2.pdb\t0.5256\t-2
protdna_complex_1.pdb\t0.1234\t-1
"""
        )
    updated_models = update_models_with_scores(
        output_fname,
        output_models,
        metric="jury_score",
        )
    assert updated_models[0].score == -0.1234
    assert updated_models[0].rank == 2
    assert updated_models[1].score == -0.5256
    assert updated_models[1].rank == 1
    assert isnan(updated_models[2].score)
    assert isnan(updated_models[2].rank)

    updated_models = update_models_with_scores(
        output_fname,
        output_models,
        metric="fake_energy",
        )
    assert updated_models[0].score == -1
    assert updated_models[0].rank == 2
    assert updated_models[1].score == -2
    assert updated_models[1].rank == 1
    assert isnan(updated_models[2].score)
    assert isnan(updated_models[2].rank)

    # Test error raising
    with pytest.raises(ValueError):
        updated_models2 = update_models_with_scores(
            output_fname,
            output_models,
            metric="wrong",
            )
        assert updated_models2 is None
    output_fname.unlink()
