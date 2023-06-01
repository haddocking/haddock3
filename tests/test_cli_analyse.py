"""Test haddock3-analyse client."""
import os
import shutil
from pathlib import Path

import pytest

from haddock.clis.cli_analyse import (
    get_cluster_ranking,
    main,
    update_capri_dict,
    )
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.modules.analysis.caprieval import \
    DEFAULT_CONFIG as caprieval_params

from . import golden_data


@pytest.fixture
def default_capri():
    """Prot-prot input."""
    return read_from_yaml_config(caprieval_params)


@pytest.fixture
def example_capri_clt():
    """Provide example capri_clt.tsv filename."""
    return Path(golden_data, "capri_clt_example.tsv")


@pytest.fixture
def example_capri_ss():
    """Provide example capri_ss.tsv filename."""
    return Path(golden_data, "capri_ss_example.tsv")


def test_update_capri_dict(default_capri):
    """Test update_capri_dict."""
    ext_params = {"reference_fname": "example.pdb"}
    wrong_ext_params = {"wrong_parameter_name": 32}
    obs_capri_dict = update_capri_dict(default_capri, ext_params)
    assert obs_capri_dict["reference_fname"].name == "example.pdb"
    # when an invalid parameter is provided, it should exit
    with pytest.raises(SystemExit):
        update_capri_dict(default_capri, wrong_ext_params)


def test_get_cluster_ranking(example_capri_clt):
    """Test get_cluster_ranking."""
    obs_cl_ranking = get_cluster_ranking(example_capri_clt, 5)
    exp_cl_ranking = {16: 1,
                      1: 2,
                      13: 3,
                      4: 4,
                      5: 5}
    assert exp_cl_ranking == obs_cl_ranking


def test_main(example_capri_ss, example_capri_clt):
    """Test cli_analyse main."""
    # build fake run_dir
    run_dir = "example_dir"
    if os.path.isdir(run_dir):
        shutil.rmtree(run_dir)
    step_name = "2_caprieval"
    step_dir = Path(run_dir, step_name)
    os.mkdir(run_dir)
    os.mkdir(step_dir)
    shutil.copy(example_capri_ss, Path(step_dir, "capri_ss.tsv"))
    shutil.copy(example_capri_clt, Path(step_dir, "capri_clt.tsv"))

    # run haddock3-analyse
    main(run_dir, [2], 5, format=None, scale=None)

    # check analysis directory exists
    ana_dir = Path(run_dir, "analysis/")
    assert os.path.isdir(ana_dir) is True

    # check whether there are some html files
    ana_subdir = Path(ana_dir, f"{step_name}_analysis")
    html_files = [el for el in os.listdir(ana_subdir) if el.endswith(".html")]
    assert len(html_files) > 0

    shutil.rmtree(run_dir)
