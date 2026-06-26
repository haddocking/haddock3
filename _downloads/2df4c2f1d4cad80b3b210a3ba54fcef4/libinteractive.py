"""Set of functions related to haddock3 interactive rescoring `haddock3-re`."""

from haddock import log
from haddock.clis.cli_traceback import get_steps_without_pdbs
from haddock.core.typing import Union
import pandas as pd
import numpy as np
from pathlib import Path
from haddock.libs.libplots import read_capri_table
from haddock.modules import get_module_steps_folders


def handle_ss_file(
        df_ss: pd.DataFrame,
        ) -> tuple[pd.DataFrame, dict]:
    """
    Manage a caprieval capri_ss file focusing on 4 first elements
    per cluster.

    Parameters
    ----------
    df_ss : pd.DataFrame
        The caprieval ss data.

    Returns
    -------
    df_ss : pd.DataFrame
        The input dataframe
    clt_ranks_dict : dict
        Dictionary with the cluster ranks
    """
    # now we want to calculate mean and std dev of the scores on df_ss
    # first sort the dataframe by score
    df_ss.sort_values(by=["score", "caprieval_rank"], inplace=True)
    # groupby cluster_id
    df_ss_grouped = df_ss.groupby("cluster_id")
    # calculate the mean and standard deviation of the first 4 elements
    # of each group
    new_values = []
    # loop over df_ss_grouped with enumerate
    for clt_id in df_ss_grouped:
        ave_score = np.mean(clt_id[1]["score"].iloc[:4])
        std_score = np.std(clt_id[1]["score"].iloc[:4])
        new_values.append([ave_score, std_score, clt_id[0]])
    # get the index that sorts the array by the first column
    new_values_arr = np.array(new_values)
    clt_ranks = np.argsort(new_values_arr[:, 0])
    # the ranked clusters are the third column of the new_values array
    clt_sorted = new_values_arr[clt_ranks, 2]
    clt_ranks_dict = {clt_sorted[i]: i + 1 for i in range(len(clt_sorted))}
    # adjust clustering values if there are clusters
    if list(np.unique(df_ss["cluster_id"])) != ["-"]:
        df_ss['model-cluster_ranking'] = df_ss.groupby('cluster_id')['score'].rank(ascending=True).astype(int)  # noqa : E501
        # assign to the values of cluster_ranking the corresponding clt_ranks
        df_ss["cluster_ranking"] = df_ss["cluster_id"].apply(lambda x: clt_ranks_dict[x])  # noqa : E501
    # assign to the column caprieval_rank the index of the dataframe
    df_ss.index = range(1, len(df_ss) + 1)
    df_ss["caprieval_rank"] = df_ss.index
    return df_ss, clt_ranks_dict


def rewrite_capri_tables(
        caprieval_folder: str,
        clt_dic: dict,
        outdir: str,
        ) -> None:
    """Rewrite the capri tables with new values.

    Parameters
    ----------
    caprieval_folder : str
        Path to the capriveal folder to be changed
    clt_dic : dict
        Data for each cluster
    outdir : str
        Output directory
    """
    capri_ss = Path(caprieval_folder, "capri_ss.tsv")
    capri_clt = Path(caprieval_folder, "capri_clt.tsv")
    if not capri_ss.exists() or not capri_clt.exists():
        # raise warning and exit
        log.warning("Capri evaluation files not found. Skipping...")
        return
    # ss file
    df_ss = read_capri_table(capri_ss)
    for cl in clt_dic:
        models = [
            f"../{model.path.split('/')[-1]}/{model.file_name}"
            for model in clt_dic[cl]
            ]
        # all the models should now have the cluster_id field
        df_ss.loc[df_ss['model'].isin(models), 'cluster_id'] = cl
    
    # delete all the models that are not in the clusters
    df_ss = df_ss[df_ss['cluster_id'] != "-"]
    # assign cluster_ranking to cluster_id (aka random assignment)
    df_ss['cluster_ranking'] = df_ss['cluster_id']
    # handle ss file
    df_ss, clt_ranks_dict = handle_ss_file(df_ss)
    
    # save capri_ss file
    capri_ss_file = Path(outdir, "capri_ss.tsv")
    log.info(f"Saving capri_ss file to {capri_ss_file}")
    df_ss.to_csv(capri_ss_file, sep="\t", index=False)

    # retrieve df_clt object
    df_clt = handle_clt_file(df_ss, clt_ranks_dict)
    # save capri_clt file
    capri_clt_file = Path(outdir, "capri_clt.tsv")
    log.info(f"Saving capri_clt file to {capri_clt_file}")
    df_clt.to_csv(capri_clt_file, sep="\t", index=False, float_format='%.3f')
    return


def look_for_capri(run_dir: str, module_id: int) -> Union[Path, None]:
    """Look for capri evaluation files previous to clustfcc_dir.

    Parameters
    ----------
    run_dir : str
        Path to the haddock3 run directory
    module_id : int
        Id of the module.
    
    Returns
    -------
    capri_eval : Path
        Path to the capri evaluation file
    """
    prev_modules_id = range(1, module_id)
    prev_modules = get_module_steps_folders(run_dir, prev_modules_id)
    # remove possible interactive modules
    prev_modules = [
        mod for mod in prev_modules
        if not mod.endswith("interactive")
        ]
    # analysis modules
    ana_modules = get_steps_without_pdbs(run_dir, prev_modules)
    # loop over the reversed list of previous modules
    capri_folder = None
    for prev_module in reversed(prev_modules):
        log.info(f"prev_module {prev_module}")
        if prev_module.endswith("caprieval"):
            # caprieval module found before any module that generates models
            capri_folder = Path(run_dir, prev_module)
            break
        elif prev_module not in ana_modules:
            break
        else:
            continue
    log.info(f"capri_folder {capri_folder}")
    return capri_folder


def handle_clt_file(df_ss, clt_ranks_dict):
    """Handle capri.clt file.

    Parameters
    ----------
    df_ss : pd.DataFrame
        The caprieval ss data.
    clt_ranks_dict : dict
        New cluster ranking dictionary.
    
    Reclustering modifies the cluster data, so capri_clt.tsv must be updated.
    """
    capri_keys = ["score", "irmsd", "fnat", "lrmsd", "dockq"]
    model_keys = ["air", "bsa", "desolv", "elec", "total", "vdw"]
    df_ss_grouped = df_ss.groupby("cluster_id")
    # loop over df_ss_grouped
    cl_data = []
    for clt_id in df_ss_grouped:
        cl_rank = clt_ranks_dict[clt_id[0]] # the rank of the cluster is the key
        data = [cl_rank, clt_id[0], clt_id[1].shape[0], "-", ]
        # updating capri quantities
        for column in capri_keys:
            ave_score = np.mean(clt_id[1][column].iloc[:4])
            std_score = np.std(clt_id[1][column].iloc[:4])
            data.extend([ave_score, std_score])
        # updating model quantities
        for column in model_keys:
            ave_score = np.mean(clt_id[1][column].iloc[:4])
            std_score = np.std(clt_id[1][column].iloc[:4])
            data.extend([ave_score, std_score])
        cl_data.append(data)
    
    # create the dataframe
    capri_clt_columns = ["cluster_rank", "cluster_id", "n", "under_eval"]
    for column in capri_keys:
        capri_clt_columns.append(f"{column}")
        capri_clt_columns.append(f"{column}_std")
    for column in model_keys:
        capri_clt_columns.append(f"{column}")
        capri_clt_columns.append(f"{column}_std")
    # capri_clt_columns.append("caprieval_rank")
    df_clt = pd.DataFrame(cl_data, columns=capri_clt_columns)
    df_clt.sort_values(by="score", inplace=True)
    df_clt.index = range(1, len(df_clt) + 1)
    df_clt["caprieval_rank"] = df_clt.index
    #
    df_clt.sort_values(by="caprieval_rank", inplace=True)
    return df_clt
