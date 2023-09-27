"""haddock3-re score subcommand."""
import json
from pathlib import Path

import numpy as np

from haddock import log
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


def rescore(capri_dir, w_elec=None, w_vdw=None, w_desolv=None, w_bsa=None, w_air=None):  # noqa: E501
    """Rescore the CAPRI models."""
    log.info(f"Rescoring {capri_dir}")
    # load the scoring pars via json
    scoring_pars = json.load(open(Path(capri_dir, "weights_params.json"), "r"))
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
    run_dir = Path(capri_dir).parent
    capri_name = Path(capri_dir).name
    
    # create the interactive folder
    outdir = Path(run_dir, f"{capri_name}_interactive")
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
