"""
Rescore a step of a run.

Considering the example run::

    run1/
        0_topoaa/
        1_rigidbody/
        2_seletop/
        3_flexref/
        (etc...)


USAGE::

    haddock3-int_rescore -r <run_dir> -m <num_modules>
    haddock3-analyse -r run1 -m 1


Where, ``-m 1 3`` means that the analysis will be performed on ``1_rigidbody``.

IMPORTANT: only one module can be selected (as the interactive analysis will
be performed only on one module at a time).
"""
import argparse
import sys
from pathlib import Path

import numpy as np

from haddock import log
from haddock.gear.config import load as read_config
from haddock.libs.libplots import read_capri_table
from haddock.modules import get_module_steps_folders


# TODO: see if it's possible to avoid hard-coding
CNS_MODULES = ["rigidbody",
               "flexref",
               "emscoring",
               "mdscoring",
               "mdref",
               "emref"]
WEIGHTS = ["w_elec", "w_vdw", "w_desolv", "w_bsa", "w_air"]
ANA_FOLDER = "interactive"  # name of the analysis folder

# Command line interface parser
ap = argparse.ArgumentParser(
    prog="haddock3-int_rescore",
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

ap.add_argument(
    "-r",
    "--run-dir",
    help="The input run directory.",
    required=True,
    )

ap.add_argument(
    "-m",
    "--module",
    help="The number of the step to rescore.",
    required=True,
    type=int,
    )

ap.add_argument(
    "-e",
    "--w_elec",
    help="weight of the electrostatic component.",
    required=False,
    type=float,
    )

ap.add_argument(
    "-w",
    "--w_vdw",
    help="weight of the van-der-Waals component.",
    required=False,
    type=float,
    )

ap.add_argument(
    "-d",
    "--w_desolv",
    help="weight of the desolvation component.",
    required=False,
    type=float,
    )

ap.add_argument(
    "-b",
    "--w_bsa",
    help="weight of the BSA component.",
    required=False,
    type=float,
    )

ap.add_argument(
    "-a",
    "--w_air",
    help="weight of the AIR component.",
    required=False,
    type=float,
    )


def _ap():
    return ap


def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = vars(load_args(ap))
    main(**cmd)


def maincli():
    """Execute main client."""
    cli(ap, main)


def main(run_dir, module, w_elec, w_vdw, w_desolv, w_bsa, w_air, **kwargs):
    """
    Analyse CLI.

    Parameters
    ----------
    run_dir : str or Path
        Path to the original run directory.

    module : int
        List of the integer prefix of the modules to copy.
    """
    log.level = 20
    log.info(f"Running haddock3-int_rescore on {run_dir}, module {module}, "
             f"with w_elec = {w_elec}")

    # Reading steps
    log.info("Reading input run directory")
    # get the module folders from the run_dir input
    module_list = list(range(1, module + 1))
    sel_steps = get_module_steps_folders(Path(run_dir), module_list)

    mod = len(sel_steps) - 2
    while mod > -1:
        st_name = sel_steps[mod].split("_")[1]
        if st_name in CNS_MODULES:
            cns_step = sel_steps[mod]
            break
        mod -= 1
    
    log.info(f"selected step is {sel_steps[-1]}, last CNS step is {cns_step}")

    # now we read the params.cfg file from the CNS module
    cns_params = read_config(Path(run_dir, cns_step, "params.cfg"))
    key = list(cns_params['final_cfg'].keys())[0]
    scoring_pars = {kv: cns_params['final_cfg'][key][kv] for kv in WEIGHTS}

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

    log.info(f"rescoring_pars: {scoring_pars}")
    capri_ss = Path(run_dir, sel_steps[-1], "capri_ss.tsv")
    capri_clt = Path(run_dir, sel_steps[-1], "capri_clt.tsv")
    # ss file
    df_ss = read_capri_table(capri_ss)
    # now we want to rewrite the score parameter
    new_scores = scoring_pars["w_vdw"] * df_ss["vdw"] + \
        scoring_pars["w_elec"] * df_ss["elec"] + \
        scoring_pars["w_bsa"] * df_ss["bsa"] + \
        scoring_pars["w_desolv"] * df_ss["desolv"] + \
        scoring_pars["w_air"] * df_ss["air"]

    df_ss["score"] = new_scores
    # sort the dataframe by score
    df_ss.sort_values(by=["score"], inplace=True)
    # assign to the column caprieval_rank the index of the dataframe
    df_ss.index = range(1, len(df_ss) + 1)
    df_ss["caprieval_rank"] = df_ss.index
    
    # now we want to write the new capri_ss.tsv file
    outdir = Path(run_dir, f"{sel_steps[-1]}_interactive")
    outdir.mkdir(exist_ok=True)
    df_ss.to_csv(Path(outdir, "capri_ss.tsv"), sep="\t", index=False)

    # now we want to calculate mean and std dev of the scores on df_ss
    # first groupby score
    df_ss_grouped = df_ss.groupby("cluster-ranking")
    # calculate the mean and standard deviation of the first 4 elements
    # of each group
    new_values = np.zeros((len(df_ss_grouped), 2))
    # loop over df_ss_grouped with enumerate
    for i, clt_id in enumerate(df_ss_grouped):
        ave_score = np.mean(clt_id[1]["score"].iloc[:4])
        std_score = np.std(clt_id[1]["score"].iloc[:4])
        new_values[i] = [ave_score, std_score]
    
    # get the index that sorts the array by the first column
    clt_ranks = np.argsort(new_values[:, 0])

    # CLT file
    df_clt = read_capri_table(capri_clt)
    # it may not be ordered by cluster_rank
    df_clt.sort_values(by=["cluster_rank"], inplace=True)
    df_clt["score"] = new_values[:, 0]
    df_clt["score_std"] = new_values[:, 1]
    if df_clt["cluster_rank"].iloc[0] != "-":
        df_clt["cluster_rank"] = clt_ranks + 1
        df_clt["caprieval_rank"] = clt_ranks + 1
    df_clt.to_csv(Path(outdir, "capri_clt.tsv"),
                  sep="\t",
                  index=False,
                  float_format="%.3f")
    return


if __name__ == "__main__":
    sys.exit(maincli())
