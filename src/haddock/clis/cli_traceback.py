"""
Traces back PDB files from a HADDOCK run directory.

Given an input run directory, haddock3-traceback traces back each model to the
initial input molecules used, providing the rank of each intermediate model.

USAGE::

    haddock3-traceback -r <run_dir>

"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from haddock import log
from haddock.core.typing import Any, FilePath
from haddock.libs import libcli
from haddock.libs.libontology import ModuleIO, PDBFile
from haddock.libs.libplots import make_traceback_plot
from haddock.modules import get_module_steps_folders


TRACK_FOLDER = "traceback"  # name of the traceback folder


def get_steps_without_pdbs(run_dir, all_steps):
    """
    Get the modules that do not produce PDB files.

    Parameters
    ----------
    run_dir : str or pathlib.Path
        Path to the run directory.
    
    all_steps : list
        List of all the steps in the run directory.
    
    Returns
    -------
    steps_without_pdbs : list
        List of steps that did not produce PDB files.
    """
    steps_without_pdbs = []
    for step in all_steps:
        if step.endswith("topoaa"):
            steps_without_pdbs.append(step)
        else:
            step_dir = Path(run_dir, step)
            if step_dir.is_dir():
                pdbs = list(step_dir.glob("*.pdb*"))
                if len(pdbs) == 0:
                    steps_without_pdbs.append(step)
    return steps_without_pdbs


def get_ori_names(n: int, pdbfile: PDBFile, max_topo_len: int) -> tuple[list, int]:
    """
    Get the original name(s) of the PDB file.

    Parameters
    ----------
    n : int
        Step number.
    pdbfile : PDBFile
        PDBFile object.
    max_topo_len : int
        Maximum length of the topologies found so far.

    Returns
    -------
    ori_names : list
        List of original names.
    max_topo_len : int
        Maximum length of the topologies found so far.
    """
    if n != 0:  # not the first step, ori_name should be defined
        ori_names = [pdbfile.ori_name]
    else:  # first step, we get topology files instead of ori_name
        # topology can either be a list of topologies or a single
        # topology
        if isinstance(pdbfile.topology, list):
            ori_names = [el.file_name for el in pdbfile.topology]
            if len(pdbfile.topology) > max_topo_len:
                max_topo_len = len(pdbfile.topology)
        else:
            ori_names = [pdbfile.topology.file_name]  # type: ignore
            max_topo_len = 1
    return ori_names, max_topo_len


def traceback_dataframe(
    data_dict: dict, rank_dict: dict, sel_step: list, max_topo_len: int
) -> pd.DataFrame:
    """
    Create traceback dataframe by combining together ranks and data.

    Parameters
    ----------
    data_dict : dict
        Dictionary containing the data to be traced back.
    rank_dict : dict
        Dictionary containing the ranks of the data to be traced back.
    sel_step : list
        List of selected steps.
    max_topo_len : int
        Maximum length of the topologies.

    Returns
    -------
    df_ord : pandas.DataFrame
        Dataframe containing the traceback data.
    """
    # get last step of the workflow
    last_step = sel_step[-1]
    # data dict to dataframe
    df_data = pd.DataFrame.from_dict(data_dict, orient="index")
    df_data.reset_index(inplace=True)
    # assign columns
    data_cols = [el for el in reversed(sel_step)]
    data_cols.extend([f"00_topo{i+1}" for i in range(max_topo_len)])
    df_data.columns = data_cols

    # same for the rank_dict
    df_ranks = pd.DataFrame.from_dict(rank_dict, orient="index")
    df_ranks.reset_index(inplace=True)
    ranks_col = [last_step]  # the key to merge the dataframes
    ranks_col.extend([f"{el}_rank" for el in reversed(sel_step)])
    df_ranks.columns = ranks_col

    # merging the data and ranks dataframes
    df_merged = pd.merge(df_data, df_ranks, on=last_step)
    ordered_cols = sorted(df_merged.columns)
    df_ord = df_merged[ordered_cols]
    # last thing: substituting unk records with - in the last step
    unk_records = df_ord[f"{last_step}"].str.startswith("unk")
    df_ord.loc[unk_records, last_step] = "-"
    return df_ord


def order_traceback_df(df_output, sel_step):
    """
    Order the traceback dataframe. Each step is ordered by rank.

    Parameters
    ----------
    df_output : pandas.DataFrame
        Dataframe containing the traceback data.
    
    sel_step : list
        List of selected steps.
    
    Returns
    -------
    df_output : pandas.DataFrame
        Dataframe containing the ordered traceback data.
    """
    # loop over sel_step in reverse order
    sorted_list = []
    indexes = []
    for n in range(len(sel_step) - 1, -1, -1):
        rank_col = sel_step[n] + "_rank"
        # take only models with a rank
        df_last = df_output[df_output[rank_col] != "-"]
        # remove from df_last the indexes that are already in the dataframe
        df_last = df_last[~df_last.index.isin(indexes)]
        # sorting the dataframe by rank
        sorted_df_last = df_last.sort_values(by=rank_col)
        sorted_list.append(sorted_df_last)
        # concat the current indexes with the previous ones
        indexes = sorted_df_last.index.tolist() + indexes
    df_output = pd.concat(sorted_list)
    return df_output


def subset_traceback(traceback_df: pd.DataFrame, cons_filename: Path) -> pd.DataFrame:
    """
    Generate a subset the traceback dataframe with the top 40 models.

    Parameters
    ----------
    traceback_df : pandas.DataFrame
        Dataframe containing the traceback data.
    
    cons_filename : pathlib.Path
        name of the consensus file.
    
    Returns
    -------
    rank_data_subset : pandas.DataFrame
        Dataframe containing the subset of the traceback data.
    """
    red_traceback_df = traceback_df.head(40)
    rank_data = red_traceback_df.filter(regex='rank$')
    # get the last column: this will define the name of the models
    last_column = red_traceback_df.columns[-2]
    last_column_data = red_traceback_df[last_column]
    # concat the ranks with the model name
    rank_data = pd.concat([last_column_data, rank_data], axis=1)
    # the rank of the last column must be defined
    last_column_rank = last_column + "_rank"
    # copy avoids SettingWithCopyWarning
    rank_data_subset = rank_data[rank_data[last_column_rank] != '-'].copy()
    rank_columns = rank_data_subset.columns[1:].tolist()
    rank_data_subset[rank_columns] = rank_data_subset[rank_columns].apply(pd.to_numeric, errors='coerce')
    rank_data_subset = rank_data_subset.sort_values(by=last_column_rank, ascending=True)
    
    # rename the columns
    rank_data_subset.columns = ["Model"] + rank_columns
    # sum the ranks
    rank_data_subset["Sum-of-Ranks"] = rank_data_subset[rank_columns].sum(axis=1)
    rank_data_subset.to_csv(cons_filename, index=False, sep="\t")
    return rank_data_subset


# Command line interface parser
ap = argparse.ArgumentParser(
    prog="haddock3-traceback",
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

libcli.add_rundir_arg(ap)


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


def main(run_dir: FilePath, offline: bool = False) -> None:
    """
    Traceback CLI.

    Parameters
    ----------
    run_dir : str or Path
        Path to the original run directory.
    """
    log.level = 20
    log.info(f"Running haddock3-traceback on {run_dir}")

    # Reading steps
    log.info("Reading input run directory")
    # get the module folders from the run_dir input
    all_steps = get_module_steps_folders(Path(run_dir))
    log.info(f"All_steps: {', '.join(all_steps)}")
    ana_modules = get_steps_without_pdbs(run_dir, all_steps)
    log.info(f"Modules not to be analysed: {', '.join(ana_modules)}")
    sel_step = [st for st in all_steps if st not in ana_modules]
    # check if there are steps to traceback
    if len(sel_step) == 0:
        log.info("No steps to trace back. Exiting.")
        return
    else:
        log.info(f"Steps to trace back: {', '.join(sel_step)}")

    # creating traceback folder
    outdir = Path(run_dir, TRACK_FOLDER)
    try:
        outdir.mkdir(exist_ok=False)
        log.info(f"Created directory: {str(outdir.resolve())}")
    except FileExistsError:
        log.warning(f"Directory {str(outdir.resolve())} already exists.")

    data_dict: dict[Any, Any] = {}
    rank_dict: dict[Any, Any] = {}
    unk_idx, max_topo_len = 0, 0
    # this cycle goes through the steps in reverse order
    for n in range(len(sel_step) - 1, -1, -1):
        log.info(f"Tracing back step {sel_step[n]}")
        # correcting names in the dictionary. The ori_name must be complemented
        # with the step folder name
        for key in data_dict.keys():
            if data_dict[key][-1] != "-":
                data_dict[key][-1] = f"../{sel_step[n]}/{data_dict[key][-1]}"

        delta = len(sel_step) - n - 1  # how many steps have we gone back?
        # loading the .json file
        json_path = Path(run_dir, sel_step[n], "io.json")
        io = ModuleIO()
        io.load(json_path)
        # list all the values in the data_dict
        ls_values = [x for val in data_dict.values() for x in val]
        # getting and sorting the ranks for the current step folder
        ranks = [pdbfile.score for pdbfile in io.output]
        ranks_argsort = np.argsort(ranks)

        # iterating through the pdbfiles to fill data_dict and rank_dict
        for i, pdbfile in enumerate(io.output):
            rank = np.where(ranks_argsort == i)[0][0] + 1
            # getting the original names
            ori_names, max_topo_len = get_ori_names(n, pdbfile, max_topo_len)
            if n != len(sel_step) - 1:
                if str(pdbfile.rel_path) not in ls_values:
                    # this is the first step in which the pdbfile appears.
                    # This means that it was discarded for the subsequent steps
                    # We need to add the pdbfile to the data_dict
                    keys = [f"unk{unk_idx}"]
                    data_dict[keys[0]] = ["-" for el in range(delta - 1)]
                    data_dict[keys[0]].append(str(pdbfile.rel_path))
                    rank_dict[keys[0]] = ["-" for el in range(delta)]
                    unk_idx += 1
                else:
                    # we've already seen this pdb before.
                    idxs = [i for i, el in enumerate(ls_values) if el==str(pdbfile.rel_path)]
                    keys = [list(data_dict.keys())[idx // delta] for idx in idxs]

                # assignment
                for el in ori_names:
                    for key in keys:
                        data_dict[key].append(el)
                for key in keys:
                    rank_dict[key].append(rank)
            else:  # last step of the workflow
                data_dict[str(pdbfile.rel_path)] = [on for on in ori_names]
                rank_dict[str(pdbfile.rel_path)] = [rank]

        # print(f"rank_dict {rank_dict}")
        # print(f"data_dict {data_dict}, maxtopo {max_topo_len}")

    # stripping away relative paths
    final_data_dict = {}
    for key in data_dict.keys():
        new_key = key.split("/")[-1]
        final_data_dict[new_key] = [el.split("/")[-1] for el in data_dict[key]]
    final_rank_dict = {}
    for key in rank_dict.keys():
        new_key = key.split("/")[-1]
        final_rank_dict[new_key] = rank_dict[key]
    # dumping the data into a dataframe
    df_output = traceback_dataframe(
        final_data_dict, final_rank_dict, sel_step, max_topo_len
    )

    # ordering the dataframe
    df_output = order_traceback_df(df_output, sel_step)
    # dumping the dataframe
    track_filename = Path(run_dir, TRACK_FOLDER, "traceback.tsv")
    log.info(
        f"Output dataframe {track_filename} " f"created with shape {df_output.shape}"
    )
    df_output.to_csv(track_filename, sep="\t", index=False)

    # taking (and writing) a subset of the dataframe
    consensus_filename = Path(run_dir, TRACK_FOLDER, "consensus.tsv")
    rank_data_subset = subset_traceback(df_output, consensus_filename)

    # plotting the traceback dataframe
    plot_filename = Path(run_dir, TRACK_FOLDER, "traceback.html")
    make_traceback_plot(rank_data_subset, plot_filename, offline=offline)
    return


if __name__ == "__main__":
    sys.exit(maincli())
