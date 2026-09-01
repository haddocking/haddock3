"""Test the rigidbody module."""

import random
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.libs.libsubprocess import CNSJob
from haddock.modules.sampling.rigidbody import (
    DEFAULT_CONFIG as DEFAULT_RIGIDBODY_PARAMS,
)
from haddock.modules.sampling.rigidbody import HaddockModule as RigidbodyModule
from haddock.modules.sampling.rigidbody import (
    _chainids_for_sampled_combinations,
    _repeats_of_sampled_combinations,
)


@pytest.fixture(name="rigidbody_module")
def fixture_rigidbody_module(monkeypatch):
    """???"""
    with (
        tempfile.TemporaryDirectory() as tempdir,
    ):
        monkeypatch.chdir(tempdir)

        yield RigidbodyModule(
            order=1,
            path=Path("."),
            initial_params=DEFAULT_RIGIDBODY_PARAMS,
        )


def test_prev_fnames(monkeypatch):
    """Tests the correct retrieval of ambiguous restraints information."""
    with tempfile.TemporaryDirectory() as tempdir:
        monkeypatch.chdir(tempdir)
        rigidbody = RigidbodyModule(
            order=1,
            path=Path("1_rigidbody"),
            initial_params=DEFAULT_RIGIDBODY_PARAMS,
        )
        prev_ambig_fnames = [None for md in range(rigidbody.params["sampling"])]
        diff_ambig_fnames = rigidbody.get_ambig_fnames(prev_ambig_fnames)  # type: ignore
        assert diff_ambig_fnames is None


def test_sample_models_to_dock_is_prefix_stable(rigidbody_module):
    """Increasing sampling must only append scheduled combinations."""
    combinations = [[object()], [object()], [object()]]

    for sampling in range(1, 9):
        scheduled = rigidbody_module._sample_models_to_dock(combinations, sampling)
        assert len(scheduled) == sampling
        for prefix_length in range(1, sampling + 1):
            assert scheduled[:prefix_length] == rigidbody_module._sample_models_to_dock(
                combinations, prefix_length
            )


def test_repeats_count_per_combination_not_per_job():
    """A job's repeat index counts its own combination, not the schedule.

    This is what keeps a seed out of the schedule's numbering: combination 2
    is on its first repeat at job 2 of a three-combination run and at job 2
    of a four-combination one alike, so both jobs are the same computation
    and receive the same seed.
    """
    first, second, third = [object()], [object()], [object()]
    schedule = [first, second, third, first, second, third, first]

    assert _repeats_of_sampled_combinations(schedule) == [0, 0, 0, 1, 1, 1, 2]
    assert _repeats_of_sampled_combinations(schedule[:4]) == [0, 0, 0, 1]


def test_chainids_are_checked_once_per_distinct_combination(mocker):
    """Repeated samples reuse the source combination's checked chain IDs."""
    combinations = [[object()], [object()]]
    sampled = [combinations[0], combinations[1], combinations[0]]
    check = mocker.patch(
        "haddock.modules.sampling.rigidbody.check_combination_chains",
        side_effect=[["A"], ["B"]],
    )

    observed = _chainids_for_sampled_combinations(combinations, sampled)

    assert observed == [["A"], ["B"], ["A"]]
    assert check.call_count == len(combinations)


def test_rigidbody_make_cns_jobs(rigidbody_module):
    "???"

    rigidbody_module.output_models = []
    rigidbody_module.envvars = {}

    cns_input_file = Path("rigidbody.inp")
    cns_output_file = "rigidbody_1.out"
    output_pdb_name = "rigidbody_1.pdb"
    ambig_fname = "ambig.tbl"
    seed = random.randint(0, 1000)
    topology = Persistent(file_name="topology.psf", path=".", file_type=Format.TOPOLOGY)
    input_pdb_1 = PDBFile(
        Path("model1.pdb"), path=".", restr_fname="ambig1.tbl", topology=topology
    )
    input_pdb_2 = PDBFile(
        Path("model2.pdb"), path=".", restr_fname="ambig2.tbl", topology=topology
    )
    inp_list = [
        (
            (
                input_pdb_1,
                input_pdb_2,
            ),
            cns_input_file,
            ambig_fname,
            seed,
        )
    ]

    observed_jobs = rigidbody_module.make_cns_jobs(inp_list)

    assert isinstance(observed_jobs[0], CNSJob)
    assert observed_jobs[0].input_file == cns_input_file
    assert observed_jobs[0].output_file == str(cns_output_file)
    assert rigidbody_module.output_models[0].restr_fname == ambig_fname
    assert rigidbody_module.output_models[0].file_name == output_pdb_name
    assert rigidbody_module.output_models[0].topology == [topology, topology]
    assert rigidbody_module.output_models[0].seed == seed
    assert rigidbody_module.output_models[0].ligand_top_fname is None
    assert rigidbody_module.output_models[0].ligand_param_fname is None


def test_rigidbody_make_cns_jobs_with_toppar(rigidbody_module):
    """Test rigidbody module's make_cns_jobs with topology and parameter files.

    Verifies that the rigidbody module correctly creates CNS jobs when provided
    with input PDB files that have associated topology files and ligand parameters.
    Checks that the output contains proper CNSJob objects and that the output models
    have the correct topology, restraint files, and ligand parameter references.
    """

    rigidbody_module.output_models = []
    rigidbody_module.envvars = {}

    cns_input_file = Path("rigidbody.inp")
    cns_output_file = "rigidbody_1.out"
    output_pdb_name = "rigidbody_1.pdb"
    ambig_fname = "ambig.tbl"
    seed = random.randint(0, 1000)
    topology = Persistent(file_name="topology.psf", path=".", file_type=Format.TOPOLOGY)
    input_pdb_1 = PDBFile(
        Path("model1.pdb"), path=".", restr_fname="ambig1.tbl", topology=topology
    )
    input_pdb_2 = PDBFile(
        Path("model2.pdb"),
        path=".",
        restr_fname="ambig2.tbl",
        topology=topology,
        ligand_top_fname="top",
        ligand_param_fname="param",
    )
    inp_list = [
        (
            (
                input_pdb_1,
                input_pdb_2,
            ),
            cns_input_file,
            ambig_fname,
            seed,
        )
    ]

    observed_jobs = rigidbody_module.make_cns_jobs(inp_list)

    assert isinstance(observed_jobs[0], CNSJob)
    assert observed_jobs[0].input_file == cns_input_file
    assert observed_jobs[0].output_file == str(cns_output_file)
    assert rigidbody_module.output_models[0].restr_fname == ambig_fname
    assert rigidbody_module.output_models[0].file_name == output_pdb_name
    assert rigidbody_module.output_models[0].topology == [topology, topology]
    assert rigidbody_module.output_models[0].seed == seed
    assert rigidbody_module.output_models[0].ligand_top_fname == "top"
    assert rigidbody_module.output_models[0].ligand_param_fname == "param"


def _models_on_disk() -> tuple[PDBFile, PDBFile]:
    """Two models whose bytes exist, because a job's seed is read from them."""
    topology = Persistent(file_name="topology.psf", path=".", file_type=Format.TOPOLOGY)
    models = []
    for index, restraints in enumerate(("ambig1.tbl", "ambig2.tbl"), start=1):
        name = f"model{index}.pdb"
        Path(name).write_text(f"ATOM  model {index}\n", encoding="utf-8")
        models.append(
            PDBFile(Path(name), path=".", restr_fname=restraints, topology=topology)
        )
    return models[0], models[1]


def test_prepare_cns_input_sequential(mocker, rigidbody_module):
    """???"""

    mocker.patch(
        "haddock.modules.sampling.rigidbody.prepare_cns_input",
        return_value="cns_input",
    )

    input_pdb_1, input_pdb_2 = _models_on_disk()
    observed_cns_input_list = rigidbody_module.prepare_cns_input_sequential(
        models_to_dock=[
            [
                input_pdb_1,
                input_pdb_2,
            ]
        ],
        ambig_fnames=["ambig1.tbl"],
    )

    assert observed_cns_input_list[0][0] == [input_pdb_1, input_pdb_2]
    assert observed_cns_input_list[0][1] == "cns_input"
    assert observed_cns_input_list[0][2] == "ambig1.tbl"


def test_prepare_cns_input_parallel(mocker, rigidbody_module):
    """???"""
    mock_engine_cls = mocker.Mock()
    mocker.patch(
        "haddock.modules.sampling.rigidbody.get_engine", return_value=mock_engine_cls
    )
    mock_prepare_engine = mock_engine_cls.return_value
    mock_prepare_engine.run.return_value = None
    mock_prepare_engine.results = ["cns_input"]
    mocker.patch(
        "haddock.modules.sampling.rigidbody.prepare_cns_input",
        return_value="cns_input",
    )

    input_pdb_1, input_pdb_2 = _models_on_disk()
    observed_cns_input_list = rigidbody_module.prepare_cns_input_parallel(
        models_to_dock=[
            [
                input_pdb_1,
                input_pdb_2,
            ]
        ],
        ambig_fnames=["ambig1.tbl"],
    )

    assert observed_cns_input_list[0][0] == [input_pdb_1, input_pdb_2]
    assert observed_cns_input_list[0][1] == "cns_input"
    assert observed_cns_input_list[0][2] == "ambig1.tbl"
