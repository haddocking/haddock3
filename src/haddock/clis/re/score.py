"""haddock3-re score subcommand."""
import json
from pathlib import Path
import sys

from haddock import log
from haddock.core.defaults import INTERACTIVE_RE_SUFFIX
from haddock.core.typing import Union
from haddock.libs.libinteractive import handle_ss_file, handle_clt_file
from haddock.libs.libplots import read_capri_table


def add_rescore_arguments(rescore_subcommand):
    """Add arguments to the score subcommand."""
    rescore_subcommand.add_argument(
        "capri_dir",
        help="The caprieval directory to rescore.",
        )

    rescore_subcommand.add_argument(
        "-e",
        "--w_elec",
        help="weight of the electrostatic component.",
        required=False,
        type=float,
        )

    rescore_subcommand.add_argument(
        "-w",
        "--w_vdw",
        help="weight of the van-der-Waals component.",
        required=False,
        type=float,
        )

    rescore_subcommand.add_argument(
        "-d",
        "--w_desolv",
        help="weight of the desolvation component.",
        required=False,
        type=float,
        )

    rescore_subcommand.add_argument(
        "-b",
        "--w_bsa",
        help="weight of the BSA component.",
        required=False,
        type=float,
        )

    rescore_subcommand.add_argument(
        "-a",
        "--w_air",
        help="weight of the AIR component.",
        required=False,
        type=float,
        )

    return rescore_subcommand


ANA_FOLDER = "interactive"  # name of the analysis folder


def rescore(
        capri_dir: Union[str, Path],
        w_elec: Union[bool, float] = None,
        w_vdw: Union[bool, float] = None,
        w_desolv: Union[bool, float] = None,
        w_bsa: Union[bool, float] = None,
        w_air: Union[bool, float] = None,
        ) -> Path:
    """Rescore the CAPRI models."""
    log.info(f"Rescoring {capri_dir}")
    # load the scoring pars via json
    weights_filename = Path(capri_dir, "weights_params.json")
    if not weights_filename.exists():
        log.error(f"weights file {weights_filename} not found. Exiting.")
        sys.exit(1)
    scoring_pars = json.load(open(weights_filename, "r"))
    log.info(f"Previous scoring parameters: {scoring_pars}")
    if w_elec is not None:
        scoring_pars.update({"w_elec": w_elec})
    if w_vdw is not None:
        scoring_pars.update({"w_vdw": w_vdw})
    if w_desolv is not None:
        scoring_pars.update({"w_desolv": w_desolv})
    if w_bsa is not None:
        scoring_pars.update({"w_bsa": w_bsa})
    if w_air is not None:
        scoring_pars.update({"w_air": w_air})

    log.info(f"Rescoring parameters: {scoring_pars}")

    capri_ss = Path(capri_dir, "capri_ss.tsv")
    capri_clt = Path(capri_dir, "capri_clt.tsv")
    if not capri_ss.exists() or not capri_clt.exists():
        log.error("capri_ss.tsv or capri_clt.tsv not found. Exiting.")
        sys.exit(1)
    # ss file
    df_ss = read_capri_table(capri_ss)
    # now we want to rewrite the score parameter
    new_scores = scoring_pars["w_vdw"] * df_ss["vdw"] + \
        scoring_pars["w_elec"] * df_ss["elec"] + \
        scoring_pars["w_bsa"] * df_ss["bsa"] + \
        scoring_pars["w_desolv"] * df_ss["desolv"] + \
        scoring_pars["w_air"] * df_ss["air"]

    df_ss["score"] = new_scores
    
    # now we want to write the new capri_ss.tsv file
    run_dir = Path(capri_dir).parent
    capri_name = Path(capri_dir).name
    
    # create the interactive folder
    outdir = Path(run_dir, f"{capri_name}_{INTERACTIVE_RE_SUFFIX}")
    outdir.mkdir(exist_ok=True)

    # handle ss file first
    df_ss, clt_ranks_dict = handle_ss_file(df_ss)
    capri_ss_file = Path(outdir, "capri_ss.tsv")
    log.info(f"Saving capri_ss file to {capri_ss_file}")
    df_ss.to_csv(capri_ss_file, sep="\t", index=False, float_format="%.3f")

    df_clt = handle_clt_file(df_ss, clt_ranks_dict)
    capri_clt_file = Path(outdir, "capri_clt.tsv")
    df_clt.to_csv(
            capri_clt_file,
            sep="\t",
            index=False,
            float_format="%.3f",
            )

    # Write the latest parameters file
    # define output fname
    rescoring_params_fname = Path(outdir, "weights_params.json")
    # write json file
    with open(rescoring_params_fname, 'w', encoding='utf-8') as jsonf:
        json.dump(
            scoring_pars,
            jsonf,
            indent=4,
            )
    log.info(f"new rescoring parameters written in: {rescoring_params_fname}")

    # return path to the new interactive folder
    return outdir
