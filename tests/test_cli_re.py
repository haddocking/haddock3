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
    print(f"weights_dict: {weights_dict}")
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
            interactive_folder = [el for el in os.listdir(tmpdir) if el.endswith("interactive")]
            assert len(interactive_folder) == 1
            interactive_folder = Path(tmpdir, interactive_folder[0])
            rescored_ss = read_capri_table(Path(interactive_folder, "capri_ss.tsv"))
            first_model = rescored_ss.iloc[0]
            assert first_model["dockq"] == 0.455
            assert first_model["model"] == "../1_rigidbody/rigidbody_645.pdb"
            assert first_model["irmsd"] == 2.621
            assert first_model["score"] == -28.743

def test_cli_reclustfcc():
    """Test haddock3-re clustfcc subcommand."""
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        os.mkdir(Path(tmpdir, "03_clustfcc"))
        nested_tmpdir = Path(tmpdir, "03_clustfcc")
        # json file
        flexref_json = Path(golden_data, "io_flexref.json")
        shutil.copy(flexref_json, Path(nested_tmpdir, "io.json"))
        # params.cfg
        clustfcc_params_cfg = Path(golden_data, "params_clustfcc.cfg")
        shutil.copy(clustfcc_params_cfg, Path(nested_tmpdir, "params.cfg"))
        # fcc matrix
        fcc_matrix = Path(golden_data, "example_fcc.matrix")
        shutil.copy(fcc_matrix, Path(nested_tmpdir, "fcc.matrix"))
        subprocess.run(["haddock3-re", "clustfcc", nested_tmpdir])

        # check if the interactive folders is created
        interactive_folder = [el for el in os.listdir(tmpdir) if el.endswith("interactive")]
        assert len(interactive_folder) == 1
        # check that the clustfcc.tsv file is correctly created
        interactive_folder = Path(tmpdir, interactive_folder[0])
        clustfcc_tsv = Path(interactive_folder, "clustfcc.tsv")
        assert clustfcc_tsv.exists()
        lines = clustfcc_tsv.read_text().splitlines()
        assert len(lines) == 4
        assert lines[1] == "1\tcluster_1_model_1.pdb\t6.20\t2"

        


            

    
