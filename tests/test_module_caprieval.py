"""Test the CAPRI module."""

import random
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libcapri import (
    CAPRI,
    capri_cluster_analysis,
    rank_according_to_score,
)
from haddock.modules.analysis.caprieval.capri import get_previous_cns_step


def read_capri_file(fname):
    """Helper function that reads a capri tsv file."""
    file_content = [
        e.split()[1:] for e in open(fname).readlines() if not e.startswith("#")
    ]
    return file_content


def test_capri_run(mocker, monkeypatch):
    """Test that CAPRI.run() calls each metric method exactly once."""

    mock_get_align_func = mocker.Mock(
        return_value={"numbering": {1: 1, 2: 2}, "chain_dict": {1: 2}}
    )
    mocker.patch(
        "haddock.libs.libcapri.get_align",
        return_value=mock_get_align_func,
    )

    mocker.patch.object(CAPRI, "_load_atoms", return_value=None)
    with tempfile.TemporaryDirectory() as tempdir:
        monkeypatch.chdir(tempdir)
        capri = CAPRI(
            identificator="test",
            model=Path("some-file"),
            path=Path("."),
            reference=Path("some-file"),
            params={
                "fnat": True,
                "irmsd": True,
                "lrmsd": True,
                "ilrmsd": True,
                "dockq": True,
                "global_rmsd": True,
                "allatoms": True,
                "receptor_chain": "A",
                "ligand_chains": ["B"],
                "alignment_method": None,
                "lovoalign_exec": None,
                "fnat_cutoff": 10,
                "irmsd_cutoff": 10,
                "keep_hetatm": False,
            },
        )
        rand_fnat = random.random()
        mocker.patch.object(
            capri,
            "calc_fnat",
            side_effect=lambda cutoff: setattr(capri, "fnat", rand_fnat),
        )
        rand_irmsd = random.random()
        mocker.patch.object(
            capri,
            "calc_irmsd",
            side_effect=lambda cutoff: setattr(capri, "irmsd", rand_irmsd),
        )
        rand_lrmsd = random.random()
        mocker.patch.object(
            capri, "calc_lrmsd", side_effect=lambda: setattr(capri, "lrmsd", rand_lrmsd)
        )
        rand_ilrmsd = random.random()
        mocker.patch.object(
            capri,
            "calc_ilrmsd",
            side_effect=lambda cutoff: setattr(capri, "ilrmsd", rand_ilrmsd),
        )
        rand_dockq = random.random()
        mocker.patch.object(
            capri, "calc_dockq", side_effect=lambda: setattr(capri, "dockq", rand_dockq)
        )

        rand_global_rmsd = random.random()
        mocker.patch.object(
            capri,
            "calc_global_rmsd",
            side_effect=lambda: setattr(capri, "rmsd", rand_global_rmsd),
        )

        capri.run()

        # The only logic to be tested is if the methods are called
        assert capri.fnat == pytest.approx(rand_fnat)
        assert capri.irmsd == pytest.approx(rand_irmsd)
        assert capri.lrmsd == pytest.approx(rand_lrmsd)
        assert capri.ilrmsd == pytest.approx(rand_ilrmsd)
        assert capri.dockq == pytest.approx(rand_dockq)
        assert capri.rmsd == pytest.approx(rand_global_rmsd)


def test_capri_cluster_analysis(protprot_caprimodule, protprot_input_list, monkeypatch):
    """Test the cluster analysis."""
    model1, model2 = protprot_input_list[0], protprot_input_list[1]
    model1.clt_rank, model2.clt_rank = 1, 2
    model1.clt_id, model2.clt_id = 1, 2
    model1.score, model2.score = 42.0, 50.0
    protprot_caprimodule.irmsd = 0.1
    protprot_caprimodule.fnat = 1.0
    protprot_caprimodule.lrmsd = 1.2
    protprot_caprimodule.ilrmsd = 4.3
    protprot_caprimodule.rmsd = (0.01,)
    with tempfile.TemporaryDirectory() as tempdir:
        monkeypatch.chdir(tempdir)
        capri_cluster_analysis(
            capri_list=[protprot_caprimodule, protprot_caprimodule],
            model_list=[model1, model2],
            output_fname="capri_clt.txt",
            clt_threshold=5,
            sort_key="score",
            sort_ascending=True,
            path=Path("."),
        )

        assert Path("capri_clt.txt").stat().st_size != 0

        observed_outf_l = read_capri_file("capri_clt.txt")
        expected_outf_l = [
            [
                "cluster_id",
                "n",
                "under_eval",
                "score",
                "score_std",
                "irmsd",
                "irmsd_std",
                "fnat",
                "fnat_std",
                "lrmsd",
                "lrmsd_std",
                "dockq",
                "dockq_std",
                "ilrmsd",
                "ilrmsd_std",
                "rmsd",
                "rmsd_std",
                "caprieval_rank",
            ],
            [
                "1",
                "1",
                "yes",
                "42.000",
                "0.000",
                "0.100",
                "0.000",
                "1.000",
                "0.000",
                "1.200",
                "0.000",
                "nan",
                "nan",
                "4.300",
                "0.000",
                "0.010",
                "0.000",
                "1",
            ],
            [
                "2",
                "1",
                "yes",
                "50.000",
                "0.000",
                "0.100",
                "0.000",
                "1.000",
                "0.000",
                "1.200",
                "0.000",
                "nan",
                "nan",
                "4.300",
                "0.000",
                "0.010",
                "0.000",
                "2",
            ],
        ]
        assert observed_outf_l == expected_outf_l

        # test sorting
        capri_cluster_analysis(
            capri_list=[protprot_caprimodule, protprot_caprimodule],
            model_list=[model1, model2],
            output_fname="capri_clt.txt",
            clt_threshold=5,
            sort_key="cluster_rank",
            sort_ascending=False,
            path=Path("."),
        )

        # With sort_ascending=False the caprieval_rank is assigned by score in
        # descending order, so the highest-scoring cluster (id 2, score 50) is
        # ranked first.
        observed_outf_l = read_capri_file("capri_clt.txt")
        expected_outf_l = [
            [
                "cluster_id",
                "n",
                "under_eval",
                "score",
                "score_std",
                "irmsd",
                "irmsd_std",
                "fnat",
                "fnat_std",
                "lrmsd",
                "lrmsd_std",
                "dockq",
                "dockq_std",
                "ilrmsd",
                "ilrmsd_std",
                "rmsd",
                "rmsd_std",
                "caprieval_rank",
            ],
            [
                "2",
                "1",
                "yes",
                "50.000",
                "0.000",
                "0.100",
                "0.000",
                "1.000",
                "0.000",
                "1.200",
                "0.000",
                "nan",
                "nan",
                "4.300",
                "0.000",
                "0.010",
                "0.000",
                "1",
            ],
            [
                "1",
                "1",
                "yes",
                "42.000",
                "0.000",
                "0.100",
                "0.000",
                "1.000",
                "0.000",
                "1.200",
                "0.000",
                "nan",
                "nan",
                "4.300",
                "0.000",
                "0.010",
                "0.000",
                "2",
            ],
        ]
        assert observed_outf_l == expected_outf_l


def test_rank_according_to_score_ascending():
    """Test that the lowest score is ranked first when sort_ascending is True."""
    data = {
        0: {"score": 42.0, "irmsd": 1.0},
        1: {"score": 50.0, "irmsd": 2.0},
        2: {"score": 10.0, "irmsd": 3.0},
    }
    ranked = rank_according_to_score(data, sort_key="score", sort_ascending=True)

    # Returned dict is keyed by rank (1..n) in ascending score order
    assert [v["score"] for v in ranked.values()] == [10.0, 42.0, 50.0]
    # caprieval_rank follows the score: lowest score gets rank 1
    assert ranked[1]["score"] == 10.0
    assert ranked[1]["caprieval_rank"] == 1
    assert ranked[2]["caprieval_rank"] == 2
    assert ranked[3]["score"] == 50.0
    assert ranked[3]["caprieval_rank"] == 3


def test_rank_according_to_score_descending():
    """Test that the highest score is ranked first when sort_ascending is False."""
    data = {
        0: {"score": 42.0, "irmsd": 1.0},
        1: {"score": 50.0, "irmsd": 2.0},
        2: {"score": 10.0, "irmsd": 3.0},
    }
    ranked = rank_according_to_score(data, sort_key="score", sort_ascending=False)

    # Returned dict is keyed by rank (1..n) in descending score order
    assert [v["score"] for v in ranked.values()] == [50.0, 42.0, 10.0]
    # caprieval_rank follows the score: highest score gets rank 1
    assert ranked[1]["score"] == 50.0
    assert ranked[1]["caprieval_rank"] == 1
    assert ranked[2]["caprieval_rank"] == 2
    assert ranked[3]["score"] == 10.0
    assert ranked[3]["caprieval_rank"] == 3


def test_get_previous_cns_step():
    "Test getting the previous CNS step."
    mock_steps = ["0_topoaa", "1_emscoring", "1_emscoring_fake", "2_clustfcc"]
    assert get_previous_cns_step(mock_steps, 2) == "1_emscoring"
    mock_steps_2 = ["bla", "test"]
    assert get_previous_cns_step(mock_steps_2, 1) is None
