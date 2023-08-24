"""
Traces back PDB files from a HADDOCK run directory.

Given an input run directory, haddock3-traceback traces back each model to the
initial input molecules used, providing the rank of each intermediate model.

USAGE::

    haddock3-traceback -r <run_dir>

For more information please check the tutorial at
docs/tutorials/analysing_runs.md
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from haddock import log
from haddock.libs import libcli
from haddock.libs.libontology import ModuleIO, PDBFile
from haddock.modules import get_module_steps_folders


TRACK_FOLDER = "traceback"  # name of the traceback folder

ANA_MODULES = ["caprieval",
               "seletop",
               "topoaa",
               "rmsdmatrix",
               "clustrmsd",
               "clustfcc"]


def get_ori_names(n: int, pdbfile: PDBFile, max_topo_len: int):
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
            ori_names = [pdbfile.topology.file_name]
            max_topo_len = 1
    return ori_names, max_topo_len


def traceback_dataframe(data_dict: dict,
                        rank_dict: dict,
                        sel_step: list,
                        max_topo_len: int):
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
    unk_records = df_ord[f'{last_step}'].str.startswith('unk')
    df_ord.loc[unk_records, last_step] = "-"
    return df_ord


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


def main(run_dir):
    """
    Analyse CLI.

    Parameters
    ----------
    run_dir : str or Path
        Path to the original run directory.
    """
    log.level = 20
    log.info(f"Running haddock3-traceback on {run_dir}")

    outdir = Path(run_dir, TRACK_FOLDER)
    try:
        outdir.mkdir(exist_ok=False)
        log.info(f"Created directory: {str(outdir.resolve())}")
    except FileExistsError:
        log.warning(f"Directory {str(outdir.resolve())} already exists.")

    # Reading steps
    log.info("Reading input run directory")
    # get the module folders from the run_dir input
    all_steps = get_module_steps_folders(Path(run_dir))
    log.info(f"All_steps: {', '.join(all_steps)}")
    sel_step = [st for st in all_steps if st.split("_")[1] not in ANA_MODULES]
    log.info(f"Steps to trace back: {', '.join(sel_step)}")
    
    data_dict, rank_dict = {}, {}
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
                    key = f"unk{unk_idx}"
                    data_dict[key] = ["-" for el in range(delta - 1)]
                    data_dict[key].append(str(pdbfile.rel_path))
                    rank_dict[key] = ["-" for el in range(delta)]
                    unk_idx += 1
                else:
                    # we've already seen this pdb before.
                    idx = ls_values.index(str(pdbfile.rel_path))
                    key = list(data_dict.keys())[idx // delta]
                 
                # assignment
                for el in ori_names:
                    data_dict[key].append(el)
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
    df_output = traceback_dataframe(final_data_dict,
                                    final_rank_dict,
                                    sel_step,
                                    max_topo_len)
    # dumping the dataframe
    track_filename = Path(run_dir, TRACK_FOLDER, "traceback.tsv")
    log.info(f"Output dataframe {track_filename} "
             f"created with shape {df_output.shape}")
    df_output.to_csv(track_filename, sep="\t", index=False)
    return


if __name__ == "__main__":
    sys.exit(maincli())
