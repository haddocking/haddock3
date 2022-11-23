#!/usr/bin/env python3
"""
Analyse a set of steps of a run.

In HADDOCK3 you can analyze the successful steps of a run.

Considering the example::

    run1/
        0_topoaa/
        1_rigidbody/
        2_seletop/
        3_flexref/
        (etc...)

You can use `4_flexref` step folder as a starting point for a new run.

USAGE::

    haddock3-analyse -r <run_dir> -m <num_modules>
    haddock3-analyse -r run1 -m 1 3

Where, ``-m 1 3`` means that the analysis will be performed on ``1_rigidbody`` and ``3_flexref``.
"""
import argparse
import os
import sys
import shutil
import ast
from haddock import log
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

SCATTER_PAIRS = [("score", "irmsd")]

NAMES_TO_PLOT ={
    "score" : "HADDOCK score",
    "irmsd" : "i-RMSD",
}

class _ParamsToDict(argparse.Action):
    """
    Convert command-line parameters in an argument to a dictionary.

    Example
    -------

    Where ``-x`` is an optional argument of the command-line client
    interface.

        >>> par1 1 par2 'my name' par3 [1,2,3] par4 True
        >>> {'par1': 1, 'par2': 'my name', 'par3': [1, 2, 3]}

    """

    def __call__(self, parser, namespace, ivalues, option_string=None):
        """Execute."""
        params = ivalues[::2]
        values = ivalues[1::2]

        if len(params) != len(values):
            raise parser.error(
                "The parameters and value pairs "
                "do not match for argument `-p`"
                )

        param_dict = {}
        for k, v in zip(params, values):
            print(f"k {k} v {v}")
            try:
                param_dict[k] = v
            except (ValueError, TypeError, SyntaxError):
                raise parser.error(f"Parameter {k} with invalid value {v}")

        setattr(namespace, self.dest, param_dict)

def read_capri_table(capri_filename, rows_to_skip=1):
    capri_df = pd.read_csv(capri_filename, sep="\t", header=rows_to_skip)
    return capri_df

                

def ss_plots(ss_filename):
    capri_df = read_capri_table(ss_filename, rows_to_skip=0)
    log.info(f"capri_df {capri_df}")
    for x_ax, y_ax in SCATTER_PAIRS:
        log.info(f"x_ax, y_ax {x_ax}, {y_ax}")
        if x_ax not in capri_df.columns:
            log.info(f"x axis quantity {x_ax} not present in capri table")
            continue
        if y_ax not in capri_df.columns:
            log.info(f"y axis quantity {y_ax} not present in capri table")
            continue
        plt_fname = f"{x_ax}_{y_ax}.png"
        plt.scatter(capri_df[x_ax], capri_df[y_ax])
        plt.savefig(plt_fname, dpi=200)
        plt.close()


# Command line interface parser
ap = argparse.ArgumentParser(
    prog="haddock3-analyse",
    description=__doc__,
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
    "-p",
    "--other-params",
    dest="other_params",
    help=(
        "Any other parameter of the `emscoring` module."
        "For example: -p nemsteps 1000. You can give any number of "
        "parameters."
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


def main(run_dir, modules, **kwargs):
    """
    

    Parameters
    ----------
    run_dir : str or Path
        Path to the original run directory.

    modules : list of ints
        List of the integer prefix of the modules to copy.
    """
    from pathlib import Path

    from haddock.gear.zerofill import zero_fill
    from haddock.modules import get_module_steps_folders
    from haddock.modules.analysis.caprieval import HaddockModule
    from haddock.modules.analysis.caprieval import DEFAULT_CONFIG as caprieval_params
    from haddock.libs.libontology import ModuleIO

    from haddock.gear.yaml2cfg import read_from_yaml_config
    from haddock.libs.libio import working_directory
    from haddock.libs.libworkflow import WorkflowManager

    # create analysis folder
    ana_folder = "analysis" # name of the analysis folder
    ori_cwd = os.getcwd()
    outdir = Path(run_dir, ana_folder)
    try:
        outdir.mkdir(parents=True, exist_ok=False)
    except FileExistsError:
        log.error(f"Directory {str(outdir.resolve())} already exists.")
        sys.exit(1)
    log.info(f"Created directory: {str(outdir.resolve())}")

    # reading steps
    log.info("Reading input run directory")
    # get the module folders from the run_dir input
    steps = get_module_steps_folders(run_dir)
    selected_steps = [steps[i] for i in range(len(steps)) if i in modules]
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
    log.info(f"capri_dict is {capri_dict}")

    # calling the analysis
    folder_paths = []
    for step in selected_steps:
        print(f"step {step}")
        # posticipate this
        target_path = Path(run_dir, f"{step}_analysis")
        target_path.mkdir(parents=True, exist_ok=False)
        folder_paths.append(target_path)

        new_capri_dict = capri_dict.copy()
        for key in new_capri_dict:
            if key.endswith("fname") and key is not None:
                try:
                    shutil.copy(new_capri_dict[key], Path(target_path, "reference.pdb"))
                    new_capri_dict[key] = Path("reference.pdb")
                except FileNotFoundError:
                    sys.exit(f'file not found {new_capri_dict[key]}')
        os.chdir(target_path)
        io = ModuleIO()
        filename = Path(f"/trinity/login/mgiulini/haddock3_dev/postprocessing/run1-test/{step}/io.json")
        io.load(filename)
        caprieval_module = HaddockModule(
            order=1,
            #path=target_path,
            path=Path(run_dir),
            initial_params=caprieval_params,
            )
        caprieval_module.update_params(**new_capri_dict)
        # update model info
        caprieval_module.previous_io = io
        # run capri module
        caprieval_module._run()
        # plotting
        print(f"cwd pre plot {os.getcwd()}")
        ss_file = Path("capri_ss.tsv")
        if ss_file.exists():
            log.info("ss plotting")
            ss_plots(ss_file)

        # going back
        os.chdir(ori_cwd)
        print(f"cwd {os.getcwd()}")


    
    # moving files into analysis folder
    #os.chdir(Path(run_dir))
    log.info("moving files to analysis folder")
    for directory in folder_paths:
        print(f"directory {directory}")
        shutil.move(directory, outdir)

    
    return


if __name__ == "__main__":
    sys.exit(maincli())
