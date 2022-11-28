"""
Analyse a set of steps of a run.

Considering the example run::

    run1/
        0_topoaa/
        1_rigidbody/
        2_seletop/
        3_flexref/
        (etc...)


USAGE::

    haddock3-analyse -r <run_dir> -m <num_modules>
    haddock3-analyse -r run1 -m 1 3


Where, ``-m 1 3`` means that the analysis will be performed on ``1_rigidbody``
 and ``3_flexref``.
"""
import argparse
import os
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.colors as px_colors
import plotly.express as px
import plotly.graph_objects as go

from haddock import log
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.libs.libcli import _ParamsToDict
from haddock.libs.libontology import ModuleIO
from haddock.modules import get_module_steps_folders, modules_names
from haddock.modules.analysis.caprieval import \
    DEFAULT_CONFIG as caprieval_params
from haddock.modules.analysis.caprieval import HaddockModule


ANA_FOLDER = "analysis"  # name of the analysis folder

SCATTER_PAIRS = [
    ("irmsd", "score"),
    ("irmsd", "desolv"),
    ("irmsd", "vdw"),
    ("irmsd", "elec"),
    ("irmsd", "air"),
    ("dockq", "score"),
    ("dockq", "desolv"),
    ("dockq", "vdw"),
    ("dockq", "elec"),
    ("dockq", "air"),
    ("lrmsd", "score"),
    ("lrmsd", "desolv"),
    ("lrmsd", "vdw"),
    ("lrmsd", "elec"),
    ("lrmsd", "air"),
    ("ilrmsd", "score"),
    ("ilrmsd", "desolv"),
    ("ilrmsd", "vdw"),
    ("ilrmsd", "elec"),
    ("ilrmsd", "air"),
    ]

TITLE_NAMES = {
    "score": "HADDOCK score",
    "irmsd": "i-RMSD",
    "lrmsd": "l-RMSD",
    "ilrmsd": "il-RMSD",
    "dockq": "DOCKQ",
    "desolv": "Edesolv",
    "vdw": "Evdw",
    "elec": "Eelec",
    "air": "Eair",
    "fnat": "FCC"
    }

AXIS_NAMES = {
    "score": "HADDOCK score [a.u.]",
    "irmsd": "interface RMSD [A]",
    "lrmsd": "ligand RMSD [A]",
    "ilrmsd": "interface-ligand RMSD [A]",
    "desolv": "Desolvation Energy",
    "vdw": "Van der Waals Energy",
    "elec": "Electrostatic Energy",
    "air": "Restraints Energy",
    "fnat": "Fraction of Common Contacts",
    "dockq": "DOCKQ",
    }


def read_capri_table(capri_filename, comment="#"):
    """
    Read capri table with pandas.
    
    Parameters
    ----------
    capri_filename : str or Path
        capri single structure filename
    comment : str
        the string used to denote a commented line in capri tables

    Returns
    -------
    capri_df : pandas DataFrame
        dataframe of capri values
    """
    capri_df = pd.read_csv(
        capri_filename,
        sep="\t",
        comment=comment)
    return capri_df


def get_cluster_ranking(capri_clt_filename, top_cluster):
    """
    Get capri cluster ranking.

    Parameters
    ----------
    capri_clt_filename : str or Path
        capri cluster filename
    top_cluster : int
        number of clusters to consider

    Returns
    -------
    cl_ranking : dict
        {cluster_id : cluster_rank} dictionary
    """
    cl_ranking = {}
    dfcl = read_capri_table(capri_clt_filename)
    for n in range(min(top_cluster, dfcl.shape[0])):
        cl_ranking[dfcl["cluster_id"].iloc[n]] = dfcl["caprieval_rank"].iloc[n]
    return cl_ranking


def box_plots(capri_filename, cl_ranking):
    """
    Create box plots.

    The idea is that for each of the top X-ranked clusters we create a box plot
    showing how the basic statistics are distributed within each model.
    
    Parameters
    ----------
    capri_filename : str or Path
        capri single structure filename
    cl_ranking : dict
        {cluster_id : cluster_rank} dictionary
    """
    # generating the correct dataframe
    capri_df = read_capri_table(capri_filename, comment="#")
    gb_cluster = capri_df.groupby("cluster-id")
    gb_other = pd.DataFrame([])
    gb_good = pd.DataFrame([])
    for cl_id, cl_df in gb_cluster:
        if cl_id not in cl_ranking.keys():
            gb_other = pd.concat([gb_other, cl_df])
        else:
            cl_df["capri_rank"] = cl_ranking[cl_id]
            gb_good = pd.concat([gb_good, cl_df])
    gb_other["cluster-id"] = "Other"
    gb_other["capri_rank"] = len(cl_ranking.keys()) + 1
    gb_cluster = pd.concat([gb_good, gb_other])
    
    for x_ax in AXIS_NAMES.keys():
        if x_ax not in capri_df.columns:
            log.info(f"x axis quantity {x_ax} not present in capri table")
            continue
        fig = px.box(gb_cluster,
                     x="capri_rank",
                     y=f"{x_ax}",
                     color="cluster-id",
                     boxmode="overlay",
                     points="outliers"
                     )
        # layout
        px_fname = f"{x_ax}_clt.html"
        fig.update_layout(
            xaxis=dict(
                title="Cluster rank",
                tickfont_size=14,
                titlefont_size=16,
                ),
            yaxis=dict(
                title=AXIS_NAMES[x_ax],
                titlefont_size=16,
                tickfont_size=14,
                ),
            legend=dict(x=1.01, y=1.0, font_family="Helvetica", font_size=16),
            hoverlabel=dict(font_size=16, font_family="Helvetica"),
            )

        fig.write_html(px_fname, full_html=False, include_plotlyjs='cdn')
        

def scatter_plots(capri_filename, cl_ranking):
    """
    Create scatter plots.

    The idea is that for each pair of variables of interest (SCATTER_PAIRS,
     declared as global) we create a scatter plot.
    If available, each scatter plot containts cluster information.

    Parameters
    ----------
    capri_filename : str or Path
        capri single structure filename
    cl_ranking : dict
        {cluster_id : cluster_rank} dictionary
    """
    capri_df = read_capri_table(capri_filename, comment="#")
    for x_ax, y_ax in SCATTER_PAIRS:
        log.debug(f"x_ax, y_ax {x_ax}, {y_ax}")
        if x_ax not in capri_df.columns:
            log.info(f"x axis quantity {x_ax} not present in capri table")
            continue
        if y_ax not in capri_df.columns:
            log.info(f"y axis quantity {y_ax} not present in capri table")
            continue
        
        gb_cluster = capri_df.groupby("cluster-id")
        gb_other = pd.DataFrame([])
        fig = go.Figure(layout={"width": 1000, "height": 800})
        traces = []
        for cl_id, cl_df in gb_cluster:
            if cl_id not in cl_ranking.keys():
                gb_other = pd.concat([gb_other, cl_df])
            else:
                if cl_id == "-":
                    cl_name = "Unclustered"
                else:
                    cl_name = f"Cluster {cl_id}"
                x_mean = np.mean(cl_df[x_ax])
                y_mean = np.mean(cl_df[y_ax])
                text_list = [f"Model: {cl_df['model'].iloc[n].split('/')[-1]}<br>Score: {cl_df['score'].iloc[n]}" for n in range(cl_df.shape[0])]  # noqa:E501
                colors = px_colors.qualitative.Safe
                traces.append(
                    go.Scatter(
                        x=cl_df[x_ax],
                        y=cl_df[y_ax],
                        name=cl_name,
                        mode="markers",
                        text=text_list,
                        legendgroup=cl_name,
                        marker_color=colors[cl_ranking[cl_id] - 1],
                        hoverlabel=dict(
                            bgcolor=colors[cl_ranking[cl_id] - 1],
                            font_size=16,
                            font_family="Helvetica"
                            )
                        )
                    )
                clt_text = f"{cl_name}<br>"
                if 'score' not in [x_ax, y_ax]:
                    clt_text += f"Score: {np.mean(cl_df['score']):.3f}<br>"
                clt_text += f"{x_ax}: {x_mean:.3f}<br>{y_ax}: {y_mean:.3f}"
                clt_text_list = [clt_text]
                traces.append(
                    go.Scatter(
                        x=[x_mean],
                        y=[y_mean],
                        # error bars
                        error_x=dict(
                            type='data',
                            array=[np.std(cl_df[x_ax])],
                            visible=True),
                        error_y=dict(
                            type='data',
                            array=[np.std(cl_df[y_ax])],
                            visible=True),
                        # color and text
                        marker_color=colors[cl_ranking[cl_id] - 1],
                        text=clt_text_list,
                        legendgroup=cl_name,
                        showlegend=False,
                        mode="markers",
                        marker=dict(
                            size=10,
                            symbol="square-dot"),
                        hovertemplate=f'<b>{clt_text}</b><extra></extra>',
                        hoverlabel=dict(
                            bgcolor=colors[cl_ranking[cl_id] - 1],
                            font_size=16,
                            font_family="Helvetica"
                            )
                        )
                    )
        # append trace other
        text_list_other = [f"Model: {gb_other['model'].iloc[n].split('/')[-1]}<br>Score: {gb_other['score'].iloc[n]}" for n in range(gb_other.shape[0])]  # noqa:E501
        traces.append(
            go.Scatter(
                x=gb_other[x_ax],
                y=gb_other[y_ax],
                name="Other",
                mode="markers",
                text=text_list_other,
                legendgroup="Other",
                marker=dict(
                    color="white",
                    line=dict(
                        width=2,
                        color='DarkSlateGrey')
                    ),
                hoverlabel=dict(
                    bgcolor="white",
                    font_size=16,
                    font_family="Helvetica"
                    )
                )
            )
        for trace in traces:
            fig.add_trace(trace)
        px_fname = f"{x_ax}_{y_ax}.html"

        fig.update_layout(
            title=f"{TITLE_NAMES[x_ax]} vs {TITLE_NAMES[y_ax]}",
            xaxis=dict(
                title=AXIS_NAMES[x_ax],
                tickfont_size=14,
                titlefont_size=16,
                ),
            yaxis=dict(
                title=AXIS_NAMES[y_ax],
                titlefont_size=16,
                tickfont_size=14,
                ),
            legend=dict(x=1.01, y=1.0, font_family="Helvetica", font_size=16),
            hoverlabel=dict(font_size=16, font_family="Helvetica"),
            )

        fig.write_html(px_fname, full_html=False, include_plotlyjs='cdn')


# Command line interface parser
ap = argparse.ArgumentParser(
    prog="haddock3-analyse",
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
    "--modules",
    nargs="+",
    help="The number of the steps to copy.",
    required=True,
    type=int,
    )

ap.add_argument(
    "-t",
    "--top_cluster",
    help="The number of clusters to show.",
    required=False,
    type=int,
    default=10
    )

ap.add_argument(
    "-p",
    "--other-params",
    dest="other_params",
    help=(
        "Any other parameter of the `caprieval` module."
        "For example: -p reference_fname target.pdb."
        "You can give any number of parameters."
        ),
    action=_ParamsToDict,
    default={},
    nargs="*",
    )


def _ap():
    return ap


def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = vars(load_args(ap))
    kwargs = cmd.pop("other_params")
    main(**cmd, **kwargs)


def maincli():
    """Execute main client."""
    cli(ap, main)


def run_capri_analysis(step, run_dir, capri_dict, target_path):
    """
    Run the caprieval analysis.

    Parameters
    ----------
    step : str
        step name
    run_dir : str or Path
        path to run directory
    capri_dict : dict
        capri dictionary of parameters
    target_path : Path
        path to the output folder
    """
    new_capri_dict = capri_dict.copy()
    for key in new_capri_dict:
        if key.endswith("fname") and new_capri_dict[key] not in ['', None]:
            try:
                ref_path = Path(target_path, "reference.pdb")
                shutil.copy(new_capri_dict[key], ref_path)
                new_capri_dict[key] = Path("reference.pdb")
            except FileNotFoundError:
                sys.exit(f'file not found {new_capri_dict[key]}')
    os.chdir(target_path)
    # retrieve json file with all information
    io = ModuleIO()
    filename = Path("..", f"{step}/io.json")
    io.load(filename)
    # create capri
    caprieval_module = HaddockModule(
        order=1,
        path=Path(run_dir),
        initial_params=caprieval_params,
        )
    caprieval_module.update_params(**new_capri_dict)
    # update model info
    caprieval_module.previous_io = io
    # run capri module
    caprieval_module._run()


def get_steps(run_dir, modules):
    """
    Get steps to be analysed.

    Parameters
    ----------
    run_dir : str or Path
        path to run directory
    modules : list
        list of integers (i.e. the step IDs)

    Returns
    -------
    selected_steps : list
        list of selected step names
    """
    steps = get_module_steps_folders(run_dir)
    log.info(f"Available steps {steps}")
    selected_steps = []
    for st in steps:
        splt_stepname = st.split("_")
        try:
            st_num = int(splt_stepname[0])
        except Exception as e:  # noqa:F841
            log.info(f"found badly formatted step {st}, skipping..")
            continue
        st_name = "".join(splt_stepname[1:])
        if st_num in modules and st_name in modules_names:
            selected_steps.append(st)
    return selected_steps


def analyse_step(step, run_dir, capri_dict, target_path, top_cluster):
    """
    Analyse a step.

    If the step is a caprieval step, use the available capri files.
    Otherwise, launch a capri analysis.

    Parameters
    ----------
    step : str
        step name
    run_dir : str or Path
        path to run directory
    capri_dict : dict
        capri dictionary of parameters
    target_path : Path
        path to the output folder
    """
    log.info(f"Analysing step {step}")
    
    target_path.mkdir(parents=True, exist_ok=False)
    if step.split("_")[1] != "caprieval":
        run_capri_analysis(step, run_dir, capri_dict, target_path)
    else:
        log.info(f"step {step} is caprieval, files should be already available")
        ss_fname = Path(run_dir, f"{step}/capri_ss.tsv")
        clt_fname = Path(run_dir, f"{step}/capri_clt.tsv")
        shutil.copy(ss_fname, target_path)
        shutil.copy(clt_fname, target_path)
        os.chdir(target_path)

    # plotting
    ss_file = Path("capri_ss.tsv")
    clt_file = Path("capri_clt.tsv")
    if clt_file.exists():
        cluster_ranking = get_cluster_ranking(clt_file, top_cluster)
    else:
        raise Exception(f"clustering file {clt_file} does not exist")
    if ss_file.exists():
        log.info("Plotting results..")
        scatter_plots(ss_file, cluster_ranking)
        box_plots(ss_file, cluster_ranking)


def main(run_dir, modules, top_cluster, **kwargs):
    """
    Analyse CLI.

    Parameters
    ----------
    run_dir : str or Path
        Path to the original run directory.

    modules : list of ints
        List of the integer prefix of the modules to copy.
    """
    log.info(f"Running haddock3-analyse on {run_dir}, modules {modules}, "
             f"with top_cluster = {top_cluster}")

    # Create analysis folder
    ori_cwd = os.getcwd()
    outdir = Path(run_dir, ANA_FOLDER)
    try:
        outdir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        log.error(f"Directory {str(outdir.resolve())} already exists.")
        sys.exit(1)
    log.info(f"Created directory: {str(outdir.resolve())}")

    # Reading steps
    log.info("Reading input run directory")
    # get the module folders from the run_dir input
    selected_steps = get_steps(run_dir, modules)
    log.info(f"selected steps: {', '.join(selected_steps)}")

    # modifying the parameters
    default_capri = read_from_yaml_config(caprieval_params)
    capri_dict = default_capri.copy()
    for param in kwargs:
        if param not in default_capri:
            sys.exit(f'* ERROR * Parameter {param!r} is not a valid `caprieval` parameter')  # noqa:E501
        else:
            capri_dict[param] = kwargs[param]
            log.info(f"setting {param} to {kwargs[param]}")
    
    # analysis
    good_folder_paths, bad_folder_paths = [], []
    for step in selected_steps:
        target_path = Path(run_dir, f"{step}_analysis")
        error = False
        try:
            analyse_step(step, run_dir, capri_dict, target_path, top_cluster)
        except Exception as e:
            error = True
            log.warning(
                f"""Could not execute the analysis for step {step}.
                The following error occurred {e}"""
                )
        if error:
            bad_folder_paths.append(target_path)
        else:
            good_folder_paths.append(target_path)
        
        # going back
        os.chdir(ori_cwd)

    # moving files into analysis folder
    log.info("moving files to analysis folder")
    for directory in good_folder_paths:
        shutil.move(directory, outdir)

    if bad_folder_paths != []:
        log.info("cancelling unsuccesful analysis folders")
        for directory in bad_folder_paths:
            if directory.exists():
                shutil.rmtree(directory)

    return


if __name__ == "__main__":
    sys.exit(maincli())
