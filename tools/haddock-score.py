"""
A simple tool to calculate the HADDOCK-score of a complex.

Usage:
    python haddock-score complex.pdb
"""

import argparse
import subprocess
from pathlib import Path

from pdbtools.pdb_tidy import tidy_pdbfile

from haddock import toppar_path as TOPPAR_DIR
from haddock.core.defaults import cns_exec as CNS_EXEC
from haddock.gear.haddockmodel import HaddockModel
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.libs.libcns import prepare_cns_input
from haddock.libs.libontology import PDBFile, TopologyFile
from haddock.modules.scoring.emscoring import DEFAULT_CONFIG as SCORING_DEFAULT_CONFIG
from haddock.modules.scoring.emscoring import RECIPE_PATH as emscoring_module_folder
from haddock.modules.topology.topoaa import DEFAULT_CONFIG as TOPO_DEFAULT_CONFIG
from haddock.modules.topology.topoaa import RECIPE_PATH as topoaa_module_folder
from haddock.modules.topology.topoaa import generate_topology


def topo_wrapper(pdb):
    """
    Wrapper for the topology generation.

    Parameters
    ----------
    pdb : str
        Path of the PDB file.

    Returns
    -------
    pdb_object : PDBFile
        The PDBFile object.
    """
    model_path = Path(pdb).resolve()

    tidy_model_path = Path(f"{model_path.stem}_tidy.pdb").resolve()
    with open(model_path, "r") as inp_fh:
        with open(tidy_model_path, "w") as out_fh:
            for line in tidy_pdbfile(inp_fh):
                out_fh.write(line)

    main_topoaa_cns_script_as_string = Path(
        topoaa_module_folder, "cns/generate-topology.cns"
        ).read_text()
    param_dict = read_from_yaml_config(TOPO_DEFAULT_CONFIG)
    mol_dict = param_dict.pop("mol1")

    topo_inp_file_name = generate_topology(
        tidy_model_path,
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

    output_pdb = Path(
        Path.cwd(),
        f"{tidy_model_path.stem}_haddock.pdb").resolve()
    output_psf = Path(
        Path.cwd(),
        f"{tidy_model_path.stem}_haddock.psf").resolve()

    pdb_obj = PDBFile(
        file_name=output_pdb,
        topology=TopologyFile(
            output_psf,
            path="."),
        path=".")

    topo_inp_file_name.unlink()
    topo_out_file_name.unlink()
    topo_err_file_name.unlink()
    tidy_model_path.unlink()

    return pdb_obj


def emscoring_wrapper(pdb_object, root="calc-hs", model_num=1):
    """
    Wrapper for the EMScoring.

    Parameters
    ----------
    pdb_object : PDBFile
        The PDBFile object.

    root : str
        The root of the output PDB name.

    model_num : int
        The model number.

    Returns
    -------
    component_dic : dict
        The dictionary with the energetic components.
    """

    main_emscoring_cns_script_as_string = Path(
        emscoring_module_folder, "cns/emscoring.cns"
        ).read_text()

    topology = pdb_object.topology

    scoring_inp_file_name = prepare_cns_input(
        1,
        pdb_object,
        Path.cwd(),
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

    output_pdb = Path(Path.cwd(), f"{root}_{model_num}.pdb").resolve()

    pdb_obj = PDBFile(file_name=output_pdb, topology=topology, path=".")

    # Clean everything except the minimized pdb
    scoring_out_file_name.unlink()
    scoring_err_file_name.unlink()
    scoring_inp_file_name.unlink()

    return pdb_obj


def main():

    parser = argparse.ArgumentParser(
        prog="haddock3-score",
        description=__doc__,
        )

    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument(
        "--full", action="store_true", help="Print all energy components"
        )
    parser.add_argument(
        "--outputpdb", action="store_true",
        help="Save the output PDB file (minimized structure)")
    parser.add_argument(
        "--outputpsf", action="store_true",
        help="Save the output PSF file (topology)")

    args = parser.parse_args()

    # Create the topology and return a PDBFile object
    pdb_object = topo_wrapper(args.pdb_file)

    # Run emscoring and return the energy components
    emscoring_pdb = emscoring_wrapper(pdb_object)

    pdb_object.rel_path.unlink()

    haddock_score_component_dic = HaddockModel(emscoring_pdb.rel_path).energies

    vdw = haddock_score_component_dic["vdw"]
    elec = haddock_score_component_dic["elec"]
    desolv = haddock_score_component_dic["desolv"]
    air = haddock_score_component_dic["air"]
    bsa = haddock_score_component_dic["bsa"]

    # emscoring is equivalent to itw
    haddock_score_itw = (
        1.0 * vdw) + (0.2 * elec) + (1.0 * desolv) + (0.1 * air)

    print("> HADDOCK-score = (1.0 * vdw) + (0.2 * elec) + (1.0 * desolv) + (0.1 * air)")
    if args.full:
        print(f"> vdw={vdw},elec={elec},desolv={desolv},air={air},bsa={bsa}")

    original_name = Path(args.pdb_file).stem
    if args.outputpdb:
        # haddock_model
        output_pdb = Path(f"{original_name}_hs.pdb")
        emscoring_pdb.rel_path.rename(output_pdb)
    else:
        emscoring_pdb.rel_path.unlink()

    if args.outputpsf:
        output_psf = Path(f"{original_name}_hs.psf")
        emscoring_pdb.topology.rel_path.rename(output_psf)
    else:
        emscoring_pdb.topology.rel_path.unlink()

    print(f"# HADDOCK-score (emscoring) = {haddock_score_itw:.4f}")


if __name__ == "__main__":
    main()
