"""Fitness calculation module."""
import subprocess
import uuid
from pathlib import Path

from haddock import toppar_path as TOPPAR_DIR
from haddock.core.defaults import cns_exec as CNS_EXEC
from haddock.gear.haddockmodel import HaddockModel
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.libs.libcns import prepare_cns_input
from haddock.libs.libontology import PDBFile, TopologyFile
from haddock.modules.scoring.emscoring import \
    DEFAULT_CONFIG as SCORING_DEFAULT_CONFIG
from haddock.modules.scoring.emscoring import \
    RECIPE_PATH as emscoring_module_folder
from haddock.modules.topology.topoaa import \
    DEFAULT_CONFIG as TOPO_DEFAULT_CONFIG
from haddock.modules.topology.topoaa import RECIPE_PATH as topoaa_module_folder
from haddock.modules.topology.topoaa import generate_topology


def run_topoaa(pdb):
    """Generate the topology."""
    main_topoaa_cns_script_as_string = Path(
        topoaa_module_folder, "cns/generate-topology.cns"
        ).read_text()
    param_dict = read_from_yaml_config(TOPO_DEFAULT_CONFIG)
    mol_dict = param_dict.pop("mol1")

    topo_inp_file_name = generate_topology(
        Path(pdb.name),
        main_topoaa_cns_script_as_string,
        param_dict,
        mol_dict,
        TOPPAR_DIR,
        )
    topo_out_file_name = topo_inp_file_name.with_suffix(".out")
    topo_err_file_name = topo_inp_file_name.with_suffix(".err")

    env = {
        "MODULE": Path(topoaa_module_folder, "cns").resolve(),
        "MODDIR": ".",
        "TOPPAR": TOPPAR_DIR,
        }

    inp = open(topo_inp_file_name, "r")
    out = open(topo_out_file_name, "w")
    err = open(topo_err_file_name, "w")
    _ = subprocess.run(CNS_EXEC, env=env, stdin=inp, stdout=out, stderr=err)
    topo_inp_file_name.unlink()
    topo_out_file_name.unlink()
    topo_err_file_name.unlink()

    topoaa_output_pdb = Path(str(pdb).replace(".pdb", "_haddock.pdb"))
    topoaa_output_psf = Path(str(pdb).replace(".pdb", "_haddock.psf"))

    if not topoaa_output_pdb.exists():
        return

    # IMPORTANT: pass `file_name=Path.name` or CNS will fail
    topology = TopologyFile(file_name=topoaa_output_psf.name, path=".")
    pdb_object = PDBFile(file_name=topoaa_output_pdb.name,
                         topology=topology, path=".")

    return pdb_object


def run_emscoring(pdb_object):
    """Run emscoring."""
    main_emscoring_cns_script_as_string = Path(
        emscoring_module_folder, "cns/emscoring.cns"
        ).read_text()

    model_num = 1
    root = str(uuid.uuid4())[:8]
    scoring_inp_file_name = prepare_cns_input(
        model_num,
        pdb_object,
        ".",
        main_emscoring_cns_script_as_string,
        read_from_yaml_config(SCORING_DEFAULT_CONFIG),
        root,
        native_segid=True,
        )
    env = {
        "MODULE": Path(emscoring_module_folder, "cns").resolve(),
        "MODDIR": ".",
        "TOPPAR": TOPPAR_DIR,
        }

    scoring_out_file_name = scoring_inp_file_name.with_suffix(".out")
    scoring_err_file_name = scoring_inp_file_name.with_suffix(".err")

    inp = open(scoring_inp_file_name, "r")
    out = open(scoring_out_file_name, "w")
    err = open(scoring_err_file_name, "w")
    _ = subprocess.run(CNS_EXEC, env=env, stdin=inp, stdout=out, stderr=err)
    scoring_out_file_name.unlink()
    scoring_err_file_name.unlink()
    scoring_inp_file_name.unlink()

    pdb_object.rel_path.unlink()

    emscoring_pdb = Path(f"{root}_{model_num}.pdb")
    pdb_object.topology.rel_path.rename(f"{root}_{model_num}.psf")

    return emscoring_pdb


def calc_haddockscore(pdb):
    """Calculate the haddock-score."""
    topoaa_pdb_obj = run_topoaa(pdb)
    emscoring_pdb_f = run_emscoring(topoaa_pdb_obj)
    haddock_obj = HaddockModel(emscoring_pdb_f)
    weights = {"w_vdw": 1.0, "w_elec": 0.2, "w_desolv": 1.0, "w_air": 0.1}
    haddock_score = haddock_obj.calc_haddock_score(**weights)
    return haddock_score, emscoring_pdb_f
