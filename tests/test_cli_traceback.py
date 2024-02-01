"""Test haddock3-traceback client."""

import os
import shutil
from pathlib import Path
import tempfile

import pandas as pd
import pytest

from haddock.clis.cli_traceback import main, get_steps_without_pdbs

from . import golden_data


@pytest.fixture
def rigid_json():
    """Provide example rigidbody io.json file."""
    return Path(golden_data, "io_rigid.json")


@pytest.fixture
def flexref_json():
    """Provide example flexref io.json file."""
    return Path(golden_data, "io_flexref.json")


def test_main(rigid_json, flexref_json):
    """Test haddock3-traceback client."""
    # build fake run_dir
    run_dir = "example_dir"
    step_dirs = [os.path.join(run_dir, "1_rigidbody"),
                 os.path.join(run_dir, "4_flexref")]
    
    if os.path.isdir(run_dir):
        shutil.rmtree(run_dir)
    # Loop over directories to be created
    for d in [run_dir, *step_dirs]:
        os.mkdir(d)
    shutil.copy(rigid_json, os.path.join(step_dirs[0], "io.json"))
    shutil.copy(flexref_json, os.path.join(step_dirs[1], "io.json"))

    # create fake, empty pdb files
    for i in range(1, 5):
        open(os.path.join(step_dirs[0], f"rigidbody_{i}.pdb"), "w").close()
    for i in range(1, 3):
        open(os.path.join(step_dirs[1], f"flexref_{i}.pdb"), "w").close()
    
    # run haddock3-traceback
    main(run_dir)

    # check traceback folder exists
    assert os.path.isdir(os.path.join(run_dir, "traceback"))

    # check traceback files exist
    tr_file = os.path.join(run_dir, "traceback", "traceback.tsv")
    assert os.path.isfile(tr_file)

    obs_tr = pd.read_csv(tr_file, sep="\t", dtype=str)
    exp_tr = [["00_topo1", "00_topo2", "1_rigidbody", "1_rigidbody_rank", "4_flexref", "4_flexref_rank"],  # noqa: E501
              ["4G6K_fv_haddock.psf", "4I1B-matched_haddock.psf", "rigidbody_3.pdb", "1", "flexref_1.pdb", "1"],  # noqa: E501
              ["4G6K_fv_haddock.psf", "4I1B-matched_haddock.psf", "rigidbody_1.pdb", "2", "flexref_2.pdb", "2"],  # noqa: E501
              ["4G6K_fv_haddock.psf", "4I1B-matched_haddock.psf", "rigidbody_2.pdb", "4", "-", "-"],  # noqa: E501
              ["4G6K_fv_haddock.psf", "4I1B-matched_haddock.psf", "rigidbody_4.pdb", "3", "-", "-"]]  # noqa: E501
    exp_tr_df = pd.DataFrame(exp_tr[1:], columns=exp_tr[0])

    assert obs_tr.columns.tolist() == exp_tr_df.columns.tolist()
    assert obs_tr.equals(exp_tr_df)

    # clean up
    shutil.rmtree(run_dir)


def test_analysis():
    """Test traceback on a pure analysis run."""
    # build fake run_dir
    run_dir = "example_dir"
    step_dirs = [os.path.join(run_dir, "0_topoaa"),
                 os.path.join(run_dir, "1_caprieval")]
    
    if os.path.isdir(run_dir):
        shutil.rmtree(run_dir)
    # Loop over directories to be created
    for d in [run_dir, *step_dirs]:
        os.mkdir(d)

    # run haddock3-traceback
    main(run_dir)

    # check traceback folder does not exist
    assert not os.path.isdir(os.path.join(run_dir, "traceback"))


def test_get_steps_without_pdbs():
    """Test get_steps_without_pdbs."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        # build fake run_dir
        run_dir = Path(tmpdir)
        steps = ["0_topoaa", "1_rigidbody", "2_caprieval"]
        
        for st in steps:
            os.mkdir(Path(run_dir, st))

        # create fake, empty pdb files
        fake_pdb_path = Path(run_dir, steps[1] , f"rigidbody_1.pdb")
        open(fake_pdb_path, "w").close()

        # get steps without pdbs
        obs_steps = get_steps_without_pdbs(run_dir, steps)

        # check steps without pdbs
        exp_steps = ["0_topoaa", "2_caprieval"]
        assert obs_steps == exp_steps

        # now we remove the fake pdb file
        os.remove(fake_pdb_path)

        # get steps without pdbs
        obs_steps = get_steps_without_pdbs(run_dir, steps)
        exp_steps = ["0_topoaa", "1_rigidbody", "2_caprieval"]
        assert obs_steps == exp_steps
