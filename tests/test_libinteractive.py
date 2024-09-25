import tempfile
from haddock.libs.libplots import read_capri_table 
from haddock.libs.libinteractive import handle_ss_file, look_for_capri, handle_clt_file
import pytest
from pathlib import Path
import numpy as np

from . import golden_data


@pytest.fixture
def example_capri_df_ss():
    """Provide example capri_ss.tsv filename."""
    return read_capri_table(Path(golden_data, "capri_ss_example.tsv"))


def test_handle_ss_file(example_capri_df_ss):
    """Test handle_ss_file function."""
    prev_cluster_ranking = example_capri_df_ss['cluster_ranking']
    df_ss = example_capri_df_ss
    df_ss['cluster_ranking'] = df_ss['cluster_id']
    df_ss, clt_ranks_dict = handle_ss_file(df_ss)
    assert (df_ss['cluster_ranking'].values == prev_cluster_ranking.values).all()

    assert len(clt_ranks_dict) == 39
    assert clt_ranks_dict[16] == 1
    assert clt_ranks_dict[1] == 2
    assert clt_ranks_dict[34] == 39


def test_handle_ss_file_unclustered(example_capri_df_ss):
    """Test handle_ss_file function with unclustered data."""
    df_ss = example_capri_df_ss.copy()
    df_ss['cluster_id'] = '-'
    df_ss['cluster_ranking'] = "-"
    df_ss['model-cluster_ranking'] = "-"
    df_ss, clt_ranks_dict = handle_ss_file(df_ss)
    assert clt_ranks_dict == {"-": 1}
    assert df_ss['cluster_ranking'].values[0] == "-"


def test_handle_clt_file(example_capri_df_ss):
    """Test handle_clt_file function."""
    clt_ranks_dict = {
        16: 1,
        1: 2,
        13: 3,
        4: 4,
        5: 5,
        23: 6,
        18: 7,
        17: 8,
        12: 9,
        7: 10,
        8: 11,
        28: 12,
        31: 13,
        11: 14,
        20: 15,
        2: 16,
        6: 17,
        3: 18,
        24: 19,
        15: 20,
        32: 21,
        26: 22,
        25: 23,
        9: 24,
        30: 25,
        19: 26,
        38: 27,
        22: 28,
        10: 29,
        33: 30,
        21: 31,
        36: 32,
        29: 33,
        14: 34,
        27: 35,
        35: 36,
        37: 37,
        39: 38,
        34: 39}
    df_clt = handle_clt_file(example_capri_df_ss, clt_ranks_dict)
    assert df_clt.shape == (39, 27)
    # first row
    first_row = df_clt.iloc[0]
    assert first_row['cluster_rank'] == 1
    assert first_row['cluster_id'] == 16
    assert first_row['n'] == 14
    assert np.isclose(first_row['score'], -36.37725)
    # last row
    last_row = df_clt.iloc[-1]
    assert last_row['cluster_rank'] == 39
    assert last_row['cluster_id'] == 34
    assert last_row['n'] == 5
    assert np.isclose(last_row['elec'], -2.79250)


def test_look_for_capri():
    """Test look_for_capri function."""
    with tempfile.TemporaryDirectory() as tmpdirname:
        run_dir = Path(tmpdirname)
        # touch some files
        Path(run_dir, "0_topoaa").mkdir()
        Path(run_dir, "1_caprieval").mkdir()
        Path(run_dir, "2_clustrsmd").mkdir()
        capri_folder = look_for_capri(run_dir, 2)
        exp_capri_folder = Path(run_dir, "1_caprieval")
        assert capri_folder == exp_capri_folder
        