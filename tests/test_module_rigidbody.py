"""Test the rigidbody module."""

import os
import random
import tempfile
from pathlib import Path

import pytest

from haddock.libs.libontology import Format, PDBFile, Persistent
from haddock.libs.libsubprocess import CNSJob
from haddock.modules.sampling.rigidbody import \
    DEFAULT_CONFIG as DEFAULT_RIGIDBODY_PARAMS
from haddock.modules.sampling.rigidbody import HaddockModule as RigidbodyModule


@pytest.fixture(name="rigidbody_module")
def fixture_rigidbody_module():
    """???"""
    with (tempfile.TemporaryDirectory() as tempdir,):
        os.chdir(tempdir)

        yield RigidbodyModule(
            order=1,
            path=Path("."),
            initial_params=DEFAULT_RIGIDBODY_PARAMS,
        )


def test_prev_fnames():
    """Tests the correct retrieval of ambiguous restraints information."""
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        rigidbody = RigidbodyModule(
            order=1,
            path=Path("1_rigidbody"),
            initial_params=DEFAULT_RIGIDBODY_PARAMS,
        )
        prev_ambig_fnames = [None for md in range(rigidbody.params["sampling"])]
        diff_ambig_fnames = rigidbody.get_ambig_fnames(prev_ambig_fnames)  # type: ignore
        assert diff_ambig_fnames is None


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
            [
                input_pdb_1,
                input_pdb_2,
            ],
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


def test_prepare_cns_input_sequential(mocker, rigidbody_module):
    """???"""

    mocker.patch(
        "haddock.modules.sampling.rigidbody.prepare_cns_input",
        return_value="cns_input",
    )

    topology = Persistent(file_name="topology.psf", path=".", file_type=Format.TOPOLOGY)
    input_pdb_1 = PDBFile(
        Path("model1.pdb"), path=".", restr_fname="ambig1.tbl", topology=topology
    )
    input_pdb_2 = PDBFile(
        Path("model2.pdb"), path=".", restr_fname="ambig2.tbl", topology=topology
    )
    observed_cns_input_list = rigidbody_module.prepare_cns_input_sequential(
        models_to_dock=[
            [
                input_pdb_1,
                input_pdb_2,
            ]
        ],
        sampling_factor=1,
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
    mocker.patch(
        "haddock.modules.sampling.rigidbody.prepare_cns_input",
        return_value="cns_input",
    )

    topology = Persistent(file_name="topology.psf", path=".", file_type=Format.TOPOLOGY)
    input_pdb_1 = PDBFile(
        Path("model1.pdb"), path=".", restr_fname="ambig1.tbl", topology=topology
    )
    input_pdb_2 = PDBFile(
        Path("model2.pdb"), path=".", restr_fname="ambig2.tbl", topology=topology
    )
    observed_cns_input_list = rigidbody_module.prepare_cns_input_sequential(
        models_to_dock=[
            [
                input_pdb_1,
                input_pdb_2,
            ]
        ],
        sampling_factor=1,
        ambig_fnames=["ambig1.tbl"],
    )

    assert observed_cns_input_list[0][0] == [input_pdb_1, input_pdb_2]
    assert observed_cns_input_list[0][1] == "cns_input"
    assert observed_cns_input_list[0][2] == "ambig1.tbl"
