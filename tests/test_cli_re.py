"""Test the haddock3-re CLI."""
from haddock.clis.cli_re import maincli as cli_re
import pytest
import tempfile
import subprocess
from pathlib import Path
from . import golden_data
import json
import shutil
import os
from numpy import isclose
from haddock.libs.libplots import read_capri_table


@pytest.fixture
def weights_dict():
    """Provide example rigidbody parameters."""
    return {
        "w_elec": 1.0,
        "w_vdw": 0.01,
        "w_desolv": 1.0,
        "w_bsa": -0.01,
        "w_air": 0.01,
        }


def test_cli_re_empty():
    """Test haddock3-re with no arguments."""
    with pytest.raises(SystemExit):
        cli_re()


def test_cli_rescore(weights_dict):
    """Test haddock3-re rescore subcommand."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        with tempfile.TemporaryDirectory(dir=tmpdir) as nested_tmpdir:
            # weights json file
            weights_json = Path(nested_tmpdir, "weights_params.json")
            weights_json.write_text(json.dumps(weights_dict))

            # capri ss and clt file
            capri_ss = Path(golden_data, "capri_ss_example.tsv")
            capri_clt = Path(golden_data, "capri_clt_example.tsv")
            shutil.copy(capri_ss, Path(nested_tmpdir, "capri_ss.tsv"))
            shutil.copy(capri_clt, Path(nested_tmpdir, "capri_clt.tsv"))
            subprocess.run(["haddock3-re", "score", nested_tmpdir, "-e", "0.1"])

            # check if the files are created
            interactive_folder = [
                el for el in os.listdir(tmpdir)
                if el.endswith("interactive")
                ]
            assert len(interactive_folder) == 1
            interactive_folder = Path(tmpdir, interactive_folder[0])

            # check ss file
            capritable_path = Path(interactive_folder, "capri_ss.tsv")
            rescored_ss = read_capri_table(capritable_path)
            first_model = rescored_ss.iloc[0]
            assert isclose(first_model["dockq"], 0.455, atol=0.001)
            assert first_model["model"] == "../1_rigidbody/rigidbody_645.pdb"
            assert isclose(first_model["irmsd"], 2.621, atol=0.001)
            assert isclose(first_model["score"], -28.743, atol=0.001)

            # check clt file
            capriclt_path = Path(interactive_folder, "capri_clt.tsv")
            rescored_clt = read_capri_table(capriclt_path)
            first_cluster = rescored_clt.iloc[0]
            assert first_cluster["cluster_id"] == 16
            assert first_cluster["cluster_rank"] == 1
            assert isclose(first_cluster["score"], -27.537, atol=0.001)


def test_cli_reclustfcc():
    """Test haddock3-re clustfcc subcommand."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        nested_tmpdir = Path(tmpdir, "03_clustfcc")
        os.mkdir(nested_tmpdir)
        # json file
        flexref_json = Path(golden_data, "io_flexref.json")
        shutil.copy(flexref_json, Path(nested_tmpdir, "io.json"))
        # params.cfg
        clustfcc_params_cfg = Path(golden_data, "params_clustfcc.cfg")
        shutil.copy(clustfcc_params_cfg, Path(nested_tmpdir, "params.cfg"))
        # fcc matrix
        fcc_matrix = Path(golden_data, "example_fcc.matrix")
        shutil.copy(fcc_matrix, Path(nested_tmpdir, "fcc.matrix"))
        subprocess.run([
            "haddock3-re", "clustfcc", nested_tmpdir,
            "-f", "0.65",
            "-p"  # shortcut to --plot_matrix
            ])

        # check if the interactive folders is created
        interactive_folder = [
            el for el in os.listdir(tmpdir)
            if el.endswith("interactive")
            ]
        assert len(interactive_folder) == 1
        # check that the clustfcc.tsv file is correctly created
        interactive_folder = Path(tmpdir, interactive_folder[0])
        clustfcc_tsv = Path(interactive_folder, "clustfcc.tsv")
        assert clustfcc_tsv.exists()
        lines = clustfcc_tsv.read_text().splitlines()
        assert len(lines) == 4
        assert lines[1] == "1\trigidbody_3.pdb\t4.53\t2"

        # clustfcc.txt
        clustfcc_txt = Path(interactive_folder, "clustfcc.txt")
        assert clustfcc_txt.exists()
        lines = clustfcc_txt.read_text().splitlines()
        assert lines[4] == "> clust_cutoff=0.65"

        # Test generation of plot
        clustfcc_html_matrix = Path(interactive_folder, "fcc_matrix.html")
        assert clustfcc_html_matrix.exists()
        assert clustfcc_html_matrix.stat().st_size != 0
    

def test_cli_reclustrmsd():
    """Test haddock3-re clustrmsd subcommand."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        # fake ilrmsdmatrix module files
        nested_tmpdir_previousstep = Path(tmpdir, "1_ilrmsdmatrix")
        os.mkdir(nested_tmpdir_previousstep)
        # rmsdmatrix.json
        rmsdmatrix_json = Path(golden_data, "example_rmsd_matrix.json")
        shutil.copy(
            rmsdmatrix_json,
            Path(nested_tmpdir_previousstep, "rmsd_matrix.json"),
            )

        # Fake clustrmsd module files
        nested_tmpdir = Path(tmpdir, "2_clustrmsd")
        os.mkdir(nested_tmpdir)
        # json file
        flexref_json = Path(golden_data, "io_clustrmsd.json")
        shutil.copy(flexref_json, Path(nested_tmpdir, "io.json"))
        # params.cfg
        clustrmsd_params_cfg = Path(golden_data, "params_clustrmsd.cfg")
        shutil.copy(clustrmsd_params_cfg, Path(nested_tmpdir, "params.cfg"))
        # dendrogram
        dendrogram = Path(golden_data, "example_dendrogram.txt")
        shutil.copy(dendrogram, Path(nested_tmpdir, "dendrogram.txt"))
        subprocess.run([
            "haddock3-re", "clustrmsd", nested_tmpdir,
            "-n", "2",
            '-p'  # shortcut to --plot_matrix
            ])
        # check if the interactive folders is created
        interactive_folder = [
            el for el in os.listdir(tmpdir)
            if el.endswith("interactive")
            ]
        assert len(interactive_folder) == 1

        # check that the clustrmsd.tsv file is correctly created
        interactive_folder = Path(tmpdir, interactive_folder[0])
        clustrmsd_tsv = Path(interactive_folder, "clustrmsd.tsv")
        assert clustrmsd_tsv.exists()
        lines = clustrmsd_tsv.read_text().splitlines()
        assert lines[1] == "1\tensemble_4G6M_1_haddock.pdb\tnan\t1"

        # clustrmsd.txt
        clustrmsd_txt = Path(interactive_folder, "clustrmsd.txt")
        assert clustrmsd_txt.exists()
        lines = clustrmsd_txt.read_text().splitlines()
        assert lines[4] == "> criterion=maxclust"
        assert lines[5] == "> n_clusters=2"

        # cluster.out
        cluster_out = Path(interactive_folder, "cluster.out")
        assert cluster_out.exists()
        lines = cluster_out.read_text().splitlines()
        assert lines[0] == "Cluster 1 -> 1 2 3 4 5 7 8 9"
        assert lines[1] == "Cluster 2 -> 6 10"

        # Test generation of plot
        clustrmsd_html_matrix = Path(interactive_folder, "rmsd_matrix.html")
        assert clustrmsd_html_matrix.exists()
        assert clustrmsd_html_matrix.stat().st_size != 0
