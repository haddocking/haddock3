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

from haddock import log
from haddock.clis.cli_unpack import main as haddock3_unpack
from haddock.clis.cli_clean import main as haddock3_clean
from haddock.core.defaults import INTERACTIVE_RE_SUFFIX
from haddock.core.typing import (
    Any,
    ArgumentParser,
    Callable,
    FilePath,
    ImgFormat,
    Namespace,
    Optional,
    ParamDict,
    ParamMap,
)
from haddock.gear.yaml2cfg import read_from_yaml_config
from haddock.gear.clean_steps import _unpack_gz
from haddock.libs.libcli import _ParamsToDict
from haddock.libs.libio import archive_files_ext
from haddock.libs.libontology import ModuleIO
from haddock.libs.libplots import (
    ClRank,
    box_plot_handler,
    clt_table_handler,
    read_capri_table,
    report_generator,
    scatter_plot_handler,
)
from haddock.modules import get_module_steps_folders
from haddock.modules.analysis.caprieval import DEFAULT_CONFIG as caprieval_params
from haddock import modules_defaults_path
from haddock.modules.analysis.caprieval import HaddockModule


ANA_FOLDER = "analysis"  # name of the analysis folder
INTER_STR = INTERACTIVE_RE_SUFFIX  # suffix of interactive analysis folders


def get_cluster_ranking(capri_clt_filename: FilePath, top_cluster: int) -> ClRank:
    """
    Get capri cluster ranking.

    Parameters
    ----------
    capri_clt_filename : str or Path
        capri cluster filename
    top_cluster : int
        Number of clusters to be considered

    Returns
    -------
    cl_ranking : dict
        {cluster_id : cluster_rank} dictionary
    """
    cl_ranking: ClRank = {}
    dfcl = read_capri_table(capri_clt_filename)
    for n in range(min(top_cluster, dfcl.shape[0])):
        cl_ranking[dfcl["cluster_id"].iloc[n]] = dfcl["caprieval_rank"].iloc[n]
    return cl_ranking


def update_paths(
    capri_ss_filename: FilePath, toch: str = "../", toadd: str = "../../"
) -> None:
    """
    Update paths in capri_ss_filename.

    Parameters
    ----------
    capri_ss_filename : str or Path
        capri ss filename
    toch : str
        string to be replaced
    toadd : str
        string to be added
    """
    new_lines: list[str] = []
    with open(capri_ss_filename, "r") as rfile:
        for ln in rfile:
            new_ln = ln.replace(toch, toadd)
            new_lines.append(new_ln)

    with open(capri_ss_filename, "w") as wfile:
        for ln in new_lines:
            wfile.write(ln)
    return


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
    default=10,
)

ap.add_argument(
    "--format",
    help="produce images in the desired format",
    required=False,
    type=str,
    default=None,
    choices=["png", "pdf", "svg", "jpeg", "webp"],
)

ap.add_argument(
    "--scale", help="scale for images", required=False, type=float, default=1.0
)

ap.add_argument(
    "--inter",
    help="interactive analysis",
    required=False,
    type=bool,
    default=False
    )

ap.add_argument(
    "--is_cleaned",
    help="is the directory going to be cleaned?",
    required=False,
    type=bool,
    default=False
)

ap.add_argument(
    "--offline",
    help="Should plots js functions be self-contained?",
    required=False,
    type=bool,
    default=False
)

ap.add_argument(
    "--mode",
    help="mode of execution",
    required=False,
    type=str,
    default="local",
    choices=["local", "batch", "mpi"],
)

ap.add_argument(
    "--ncores",
    help="number of cores to use",
    required=False,
    type=int,
    default=1,
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


def _ap() -> ArgumentParser:
    return ap


def load_args(ap: ArgumentParser) -> Namespace:
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap: ArgumentParser, main: Callable[..., None]) -> None:
    """Command-line interface entry point."""
    cmd = vars(load_args(ap))
    kwargs = cmd.pop("other_params")
    main(**cmd, **kwargs)


def maincli() -> None:
    """Execute main client."""
    cli(ap, main)


def run_capri_analysis(step: str, run_dir: FilePath, capri_dict: ParamMap, is_cleaned: bool, mode: str, ncores: int) -> None:
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
    """
    # retrieve json file with all information
    io = ModuleIO()
    filename = Path("..", f"{step}/io.json")
    io.load(filename)
    # unpack the files if they are compressed
    if is_cleaned:
        default_general_params = read_from_yaml_config(modules_defaults_path)
        path_to_unpack = io.output[0].path
        haddock3_unpack(path_to_unpack, ncores=ncores)
    # define step_order. We add one to it, as the caprieval module will
    # interpret itself as being after the selected step
    step_order = int(step.split("_")[0]) + 1
    # create capri
    caprieval_module = HaddockModule(
        order=step_order,
        path=Path(run_dir),
        initial_params=caprieval_params,
    )
    caprieval_module.update_params(**capri_dict)
    # updating mode and ncores
    caprieval_module.params["mode"] = mode
    caprieval_module.params["ncores"] = ncores
    # update model info
    caprieval_module.previous_io = io
    # run capri module
    caprieval_module._run()
    # compress files if they should be compressed
    if is_cleaned:
        haddock3_clean(path_to_unpack, ncores=ncores)


def update_capri_dict(default_capri: ParamDict, kwargs: ParamMap) -> ParamDict:
    """
    Update capri dictionary.

    Parameters
    ----------
    default_capri : dict
        default capri dictionary of parameters
    kwargs : dict
        dictionary of input elements

    Returns
    -------
    capri_dict : dict
        updated capri dictionary of parameters
    """
    capri_dict = default_capri.copy()
    for param in kwargs:
        if param not in default_capri:
            sys.exit(
                f"* ERROR * Parameter {param!r} is not a valid `caprieval` parameter"
            )  # noqa:E501
        else:
            if param.endswith("fname"):  # using full path for files
                rel_path = Path(kwargs[param])
                _param = rel_path.resolve()
                kwargs[param] = _param
            capri_dict[param] = kwargs[param]
            log.info(f"setting {param} to {kwargs[param]}")

    return capri_dict


def update_paths_in_capri_dict(
    capri_dict: ParamDict, target_path: FilePath
) -> ParamDict:
    """
    Make capri_dict specific to target_path.

    Parameters
    ----------
    capri_dict : dict
        capri dictionary of parameters
    target_path : Path
        path to the output folder

    Returns
    -------
    new_capri_dict : dict
        target_path-specific capri dictionary of parameters
    """
    new_capri_dict = capri_dict.copy()
    for key in new_capri_dict:
        if key.endswith("fname") and new_capri_dict[key] not in ["", None]:
            try:
                ref_path = Path(target_path, "reference.pdb")
                shutil.copy(new_capri_dict[key], ref_path)
                new_capri_dict[key] = Path("reference.pdb")
            except FileNotFoundError:
                sys.exit(f"file not found {new_capri_dict[key]}")
    return new_capri_dict


def zip_top_ranked(capri_filename: FilePath, cluster_ranking: ClRank, summary_name: FilePath) -> None:
    """
    Zip the top ranked structures.

    Parameters
    ----------
    cluster_ranking : dict
        {cluster_id : cluster_rank} dictionary
    ss_file : str or Path
        capri ss filename

    Returns
    -------
    output_zipfile : str or Path
        path to the zipped file
    """
    capri_df = read_capri_table(capri_filename, comment="#")
    gb_cluster = capri_df.groupby("cluster_id")
    for cl_id, cl_df in gb_cluster:
        if cl_id in cluster_ranking.keys():
            if cl_id != "-":
                structs = cl_df.loc[cl_df["model-cluster_ranking"] <= 4][["model", "model-cluster_ranking"]]
            else:
                structs = cl_df.loc[cl_df["caprieval_rank"] <= 10][["model", "caprieval_rank"]]
            structs.columns = ["model", "rank"]
            # iterate over the structures
            for _, row in structs.iterrows():
                struct = Path(row["model"])
                struct_gz = Path(f"{struct}.gz")
                rank = row["rank"]
                # set target name
                if cl_id != "-":
                    target_name = f"cluster_{cluster_ranking[cl_id]}_model_{rank}.pdb"
                else:
                    target_name = f"model_{rank}.pdb"
                # copy the structure
                if Path(struct).exists():
                    shutil.copy(struct, Path(target_name))
                elif struct_gz.exists():
                    shutil.copy(struct_gz, ".")
                    # unpack the file
                    _unpack_gz(Path(".", struct_gz.name))
                    shutil.move(struct.name, Path(target_name))
                else:
                    log.warning(f"structure {struct} not found")

    # now make the archive and delete the pdb files
    archive_files_ext(".", "pdb")
    for file in Path(".").glob("*.pdb"):
        file.unlink()
    # move archive to summary
    expected_archive = Path(".", "pdb.tgz")
    if expected_archive.exists():
        shutil.move("pdb.tgz", summary_name)
        log.info(f"Summary archive {summary_name} created!")
    else:
        log.warning(f"Summary archive {summary_name} not created!")


def analyse_step(
    step: str,
    run_dir: FilePath,
    capri_dict: ParamDict,
    target_path: Path,
    top_cluster: int,
    format: Optional[ImgFormat],
    scale: Optional[float],
    is_cleaned: Optional[bool],
    offline: bool = False,
    mode: str = "local",
    ncores: int = 4,
) -> None:
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
    top_cluster : int
        Number of clusters to be considered
    format : str
        Produce images in the selected format.
    scale : int
        scale for images.
    """
    log.info(f"Analysing step {step}")

    target_path.mkdir(parents=True, exist_ok=False)
    step_name = step.split("_")[1]
    ss_fname = Path(run_dir, f"{step}/capri_ss.tsv")
    clt_fname = Path(run_dir, f"{step}/capri_clt.tsv")
    if step_name != "caprieval":
        if ss_fname.exists() and clt_fname.exists():
            log.info(f"step {step} has caprieval data, files are available")
            run_capri = False
        else:
            capri_dict = update_paths_in_capri_dict(capri_dict, target_path)
            run_capri = True
    else:
        log.info(f"step {step} is caprieval, files should be already available")
        run_capri = False

    if run_capri == False:
        shutil.copy(ss_fname, target_path)
        shutil.copy(clt_fname, target_path)

    # Go to directory where to write all the analysis figures / report
    os.chdir(target_path)
    # if the step is not caprieval, caprieval must be run
    if run_capri == True:
        run_capri_analysis(step, run_dir, capri_dict, is_cleaned, mode, ncores)

    log.info("CAPRI files identified")
    # plotting
    ss_file = Path("capri_ss.tsv")
    clt_file = Path("capri_clt.tsv")
    if clt_file.exists():
        cluster_ranking = get_cluster_ranking(clt_file, top_cluster)
    else:
        raise Exception(f"clustering file {clt_file} does not exist")
    if ss_file.exists():
        log.info("Plotting results..")
        scatters = scatter_plot_handler(
            ss_file,
            cluster_ranking,
            format,
            scale,
            offline=offline,
            )
        boxes = box_plot_handler(
            ss_file,
            cluster_ranking,
            format,
            scale,
            offline=offline,
            )
        tables = clt_table_handler(clt_file, ss_file, is_cleaned)
        report_generator(boxes, scatters, tables, step, '.', offline)
        # provide a zipped archive of the top ranked structures
        zip_top_ranked(ss_file, cluster_ranking, Path("summary.tgz"))


def main(
    run_dir: FilePath,
    modules: list[int],
    top_cluster: int,
    format: Optional[ImgFormat],
    scale: Optional[float],
    inter: Optional[bool],
    is_cleaned: Optional[bool],
    offline: bool = False,
    mode: Optional[str] = None,
    ncores: Optional[int] = None,
    **kwargs: Any,
) -> None:
    """
    Analyse CLI.

    Parameters
    ----------
    run_dir : str or Path
        Path to the original run directory.

    modules : list of ints
        List of the integer prefix of the modules to copy.

    top_cluster : int
        Number of clusters to be considered.

    format : str
        Produce images in the selected format.

    scale : int
        scale for images.

    inter: bool
        analyse only steps labelled as 'interactive'

    is_cleaned: bool
        is the directory going to be cleaned?
    
    offline: bool
        Should plots js functions be self-contained?
    
    mode: str
        mode of execution
    
    ncores: int
        number of cores to use
    """
    log.level = 20
    log.info(
        f"Running haddock3-analyse on {run_dir}, modules {modules}, "
        f"with top_cluster = {top_cluster}"
    )
    ori_cwd = os.getcwd()
    # modifying the parameters
    default_capri = read_from_yaml_config(caprieval_params)
    capri_dict = update_capri_dict(default_capri, kwargs)

    os.chdir(run_dir)
    # Create analysis folder
    rundir_cwd = os.getcwd()
    outdir = Path(ANA_FOLDER)
    try:
        outdir.mkdir(exist_ok=False)
        log.info(f"Created directory: {str(outdir.resolve())}")
    except FileExistsError:
        log.warning(f"Directory {str(outdir.resolve())} already exists.")

    # Reading steps
    log.info("Reading input run directory")
    # get the module folders from the run_dir input
    sel_steps = get_module_steps_folders(Path("./"), modules)
    if inter:
        sel_steps = [st for st in sel_steps if st.endswith(INTER_STR)]
    else:
        sel_steps = [st for st in sel_steps if not st.endswith(INTER_STR)]
    log.info(f"selected steps: {', '.join(sel_steps)}")

    # analysis
    good_folder_paths: list[Path] = []
    bad_folder_paths: list[Path] = []
    for step in sel_steps:
        subfolder_name = f"{step}_analysis"
        target_path = Path(Path("./"), subfolder_name)

        # check if subfolder is already present
        dest_path = Path(ANA_FOLDER, subfolder_name)
        if dest_path.exists():
            if len(os.listdir(dest_path)) != 0 and not inter:
                log.warning(
                    f"{dest_path} exists and is not empty. "
                    "Skipping analysis..."
                    )
                continue
            else:  # subfolder is empty or is interactive, remove it.
                log.info(f"Removing folder {dest_path}.")
                shutil.rmtree(dest_path)

        # run the analysis
        error = False
        try:
            analyse_step(
                step,
                Path("./"),
                capri_dict,
                target_path,
                top_cluster,
                format,
                scale,
                is_cleaned,
                offline=offline,
                mode=mode,
                ncores=ncores,
            )
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
        os.chdir(rundir_cwd)

    # moving files into analysis folder
    if good_folder_paths != []:
        log.info("moving files to analysis folder")
        for directory in good_folder_paths:
            shutil.move(directory, outdir)

    if bad_folder_paths != []:
        log.info("cancelling unsuccesful analysis folders")
        for directory in bad_folder_paths:
            if directory.exists():
                shutil.rmtree(directory)

    # substituting the correct paths in the capri_ss files
    # after moving files into analysis folder
    for directory in good_folder_paths:
        ss_file = Path(outdir, directory, "capri_ss.tsv")
        if ss_file.exists():
            log.info(f"updating paths in {ss_file}")
            update_paths(ss_file, "../", "../../")
        report_file = Path(outdir, directory, "report.html")
        log.info(f"View the results in {report_file}")
        info_msg = (
            "To view structures or download the structure files, "
            f"in a terminal run the command "
            f"`python -m http.server --directory {rundir_cwd}`. "
            "By default, http server runs on `http://0.0.0.0:8000/`. "
            f"Open the link http://0.0.0.0:8000/{report_file} "
            "in a web browser."
        )
        log.info(info_msg)
    os.chdir(ori_cwd)
    return


if __name__ == "__main__":
    sys.exit(maincli())  # type: ignore
