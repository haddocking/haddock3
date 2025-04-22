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
    SUPPORTED_OUTPUT_FORMATS,
    )
from haddock.modules import get_module_steps_folders
from haddock.modules.analysis.caprieval import (
    DEFAULT_CONFIG as caprieval_params,
    )
from haddock.modules.analysis.caprieval import HaddockModule


ANA_FOLDER = "analysis"  # name of the analysis folder
INTER_STR = INTERACTIVE_RE_SUFFIX  # suffix of interactive analysis folders


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
    "--top_clusters",
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
    choices=list(SUPPORTED_OUTPUT_FORMATS),
)

ap.add_argument(
    "--scale", help="scale for images", required=False, type=float, default=1.0
)

ap.add_argument(
    "--inter",
    help="interactive analysis",
    required=False,
    type=bool,
    default=False,
    )

ap.add_argument(
    "--is_cleaned",
    help="is the directory going to be cleaned?",
    required=False,
    type=bool,
    default=False,
)

ap.add_argument(
    "--offline",
    help="Should plots js functions be self-contained?",
    required=False,
    type=bool,
    default=False,
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
    "--self-contained",
    help="If self-contained is set, models will be copied locally in the analysis directory, allowing to visualize structures outside of the haddock3 run.",
    required=False,
    default=False,
    type=bool,
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
    cli(_ap(), main)


def get_cluster_ranking(
        capri_clt_filename: FilePath,
        top_clusters: int,
        ) -> ClRank:
    """
    Get capri cluster ranking.

    Parameters
    ----------
    capri_clt_filename : str or Path
        capri cluster filename
    top_clusters : int
        Number of clusters to be considered

    Returns
    -------
    cl_ranking : dict
        {cluster_id : cluster_rank} dictionary
    """
    cl_ranking: ClRank = {}
    dfcl = read_capri_table(capri_clt_filename)
    for n in range(min(top_clusters, dfcl.shape[0])):
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


def run_capri_analysis(
        step: str,
        run_dir: FilePath,
        capri_dict: ParamMap,
        is_cleaned: bool,
        mode: str,
        ncores: int,
        ) -> None:
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
        path_to_unpack = io.output[0].rel_path.parent
        haddock3_unpack(path_to_unpack, ncores=ncores)
    # define step_order. We add one to it, as the caprieval module will
    # interpret itself as being after the selected step
    step_order = int(step.split("_")[0]) + 1
    # create capri
    caprieval_module = HaddockModule(
        order=step_order,
        path=Path(run_dir),
        init_params=caprieval_params,
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
                f"* ERROR * Parameter {param!r} is not "
                "a valid `caprieval` parameter"
                )
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


def get_top_ranked_mapping(
        capri_filename: FilePath,
        cluster_ranking: ClRank,
        clustered_topX: int = 4,
        unclustered_topX: int = 10,
        ) -> dict[Path, Path]:
    """Obtain mapping of top ranked files to their future paths.

    Parameters
    ----------
    capri_filename : FilePath
        capri_ss.tsv filepath
    cluster_ranking : ClRank
        Cluster ranking, if the form of {cluster_id: cluster_rank} dictionary
    clustered_topX : int, optional
        Number of models to access per cluster. Default is 4.
    unclustered_topX : int, optional
        Number of unclustered models to access. Default is 10.

    Returns
    -------
    top_ranked_mapping : dict[Path, Path]
        Mapping between original filepath and new (future) one
    """
    # Set mapping of generated files
    top_ranked_mapping: dict[Path, Path] = {}

    # Read table
    capri_df = read_capri_table(capri_filename, comment="#")
    # Group by clusters
    gb_cluster = capri_df.groupby("cluster_id")

    # Loop over clusters
    for cl_id, cl_df in gb_cluster:
        # Filter only top clusters
        if cl_id in cluster_ranking.keys():
            # If clustered structure
            if cl_id != "-":
                # Retrieve only top 4 models per cluster
                structs = cl_df.loc[cl_df["model-cluster_ranking"] <= clustered_topX][["model", "model-cluster_ranking"]]  # noqa : E501
            # If un-clustered structures
            else:
                # Retrieve top 10
                structs = cl_df.loc[cl_df["caprieval_rank"] <= unclustered_topX][["model", "caprieval_rank"]]  # noqa : E501
            # Rename columns to access them using same keywords
            structs.columns = ["model", "rank"]
            # iterate over the structures
            for _, row in structs.iterrows():
                # Point rank
                rank = row["rank"]
                # set target name
                if cl_id != "-":
                    # Give it its cluster name
                    target_name = (
                        f"cluster_{cluster_ranking[cl_id]}"
                        f"_model_{rank}.pdb"
                        )
                else:
                    # Give it its rank name
                    target_name = f"model_{rank}.pdb"

                # Generate structure path
                struct = Path(row["model"])
                struct_gz = Path(f"{struct}.gz")
                # copy the structure
                if Path(struct).exists():
                    top_ranked_mapping[struct] = Path(target_name)
                elif struct_gz.exists():
                    top_ranked_mapping[struct_gz] = Path(target_name)
                else:
                    log.warning(f"structure {struct} not found")
    return top_ranked_mapping

def zip_top_ranked(
        top_ranked_mapping: dict[Path, Path],
        summary_name: str,
        gen_archive: bool,
        ) -> Optional[Path]:
    """
    Zip the top ranked structures.

    Parameters
    ----------
    capri_filename : str or Path
        capri ss filename
    cluster_ranking : dict
        {cluster_id : cluster_rank} dictionary
    summary_name: str
        Base name of the archive to be generated
    gen_archive: bool
        Should the archive be generated?
    clustered_topX: int
        Number of models to access per cluster. Default is 4.
    unclustered_topX: int
        Number of models to access when no clusters. Default is 10.

    Return
    ------
    output_fname : Optional[Path]
        Path to the generated output. Can be a .tgz archive or a directory.
    """
    for ori_fpath, new_name in top_ranked_mapping.items():
        # If already compressed
        if ori_fpath.suffix == ".gz":
            copied_fpath = Path(shutil.copy(ori_fpath, "."))
            # unpack the file
            _unpack_gz(copied_fpath)
            # Rename it
            shutil.move(copied_fpath.name.replace(".gz", ""), new_name)
        else:
            shutil.copy(ori_fpath, new_name)

    # Compress pdb files
    if gen_archive:
        archive_was_created = archive_files_ext(".", "pdb")
        # Delete the pdb files
        for file_ in top_ranked_mapping.values():
            Path(file_).unlink()
        output_fname = Path(f"{summary_name}.tgz")
        if archive_was_created:
            # move archive to summary
            shutil.move("pdb.tgz", output_fname)
            log.info(f"Top structures summary archive {output_fname} created!")
            return 
        else:
            log.warning(f"Summary archive {output_fname} not created!")
            return None
    # Generate a directory holding all the structures
    else:
        output_fname = Path(summary_name)
        output_fname.mkdir(parents=True, exist_ok=True)
        for ori_fpath, new_name in top_ranked_mapping.items():
            # Create new path
            next_filepath = Path(output_fname, str(new_name))
            # Hold it in mapping dict
            top_ranked_mapping[ori_fpath] = Path(next_filepath)
            # Displace file
            shutil.move(new_name, top_ranked_mapping[ori_fpath])
        log.info(f"Top structures copied into {output_fname}!")
        return output_fname


def analyse_step(
    step: str,
    run_dir: FilePath,
    capri_dict: ParamDict,
    target_path: Path,
    top_clusters: int,
    format: Optional[ImgFormat],
    scale: Optional[float],
    is_cleaned: Optional[bool],
    offline: bool = False,
    mode: str = "local",
    ncores: int = 4,
    self_contained: bool = False,
    clustered_topX: int = 4,
    unclustered_topX: int = 10,
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
    top_clusters : int
        Number of clusters to be considered
    format : str
        Produce images in the selected format.
    scale : int
        scale for images.
    is_cleaned: bool
        is the directory going to be cleaned?
    offline: bool
        Should plots js functions be self-contained?
    mode: str
        mode of execution
    ncores: int
        number of cores to use
    self_contained : bool
        Should the analysis directory contain the models?
    clustered_topX: int
        Number of models to access per cluster. Default is 4.
    unclustered_topX: int
        Number of models to access when no clusters. Default is 10.
    """
    log.info(f"Analysing step {step}")
    # Create directory
    target_path.mkdir(parents=True, exist_ok=False)
    # Build caprieval output file names/paths
    ss_filename = Path("capri_ss.tsv")
    clt_filename = Path("capri_clt.tsv")
    step_name = step.split("_")[1]
    ss_fname = Path(run_dir, f"{step}/{ss_filename}")
    clt_fname = Path(run_dir, f"{step}/{clt_filename}")
    # Search for caprieval output files
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

    # If caprieval data available, just copy them
    if not run_capri:
        shutil.copy(ss_fname, target_path)
        shutil.copy(clt_fname, target_path)

    # Go to directory where to write all the analysis figures / report
    os.chdir(target_path)
    # if the step is not caprieval, caprieval must be run
    if run_capri:
        run_capri_analysis(step, run_dir, capri_dict, is_cleaned, mode, ncores)

    log.info("CAPRI files identified")
    # plotting

    if clt_filename.exists():
        cluster_ranking = get_cluster_ranking(clt_filename, top_clusters)
    else:
        raise Exception(f"clustering file {clt_filename} does not exist")
    if ss_filename.exists():
        # Generate file mapping for top ranked structures
        top_ranked_mapping = get_top_ranked_mapping(
            ss_filename,
            cluster_ranking,
            clustered_topX=clustered_topX,
            unclustered_topX=unclustered_topX,
            )
        # provide a zipped archive of the top ranked structures
        zip_top_ranked(
            top_ranked_mapping,
            "summary",
            not self_contained,
            )
        log.info("Plotting results..")
        scatters = scatter_plot_handler(
            ss_filename,
            cluster_ranking,
            format,
            scale,
            offline=offline,
            )
        boxes = box_plot_handler(
            ss_filename,
            cluster_ranking,
            format,
            scale,
            offline=offline,
            )
        tables = clt_table_handler(
            clt_filename,
            ss_filename,
            is_cleaned,
            topX_clusters=top_clusters,
            clustered_topX=clustered_topX,
            unclustered_topX=unclustered_topX,
            top_ranked_mapping=top_ranked_mapping if self_contained else None,
            )
        report_generator(boxes, scatters, tables, step, ".", offline)


def validate_format(_format: Optional[ImgFormat]) -> Optional[ImgFormat]:
    """Validate the optional argument `format`.

    Parameters
    ----------
    _format : Optional[ImgFormat]
        A optential output format set by the user.

    Returns
    -------
    Optional[ImgFormat]
        A valid output format for the figures to be generated.

    Raises
    ------
    ImportError
        When the kaleido package is not installed.
    ValueError
        When the export format is not supported by plotly / kaleido.
    """
    # Optional, so if not defined, return _format (None)
    if not _format:
        return _format
    # Make sure it is lower case
    format = _format.lower()
    # Check if part of supported output formats
    if format in SUPPORTED_OUTPUT_FORMATS:
        # Make sure the `kaleido` package is installed
        try:
            import kaleido  # noqa : F401
        except ImportError:
            raise ImportError(
                f"Exporting with format {format} requires the use of `kaleido`"
                f" package, that seems not to be installed. {os.linesep}"
                "Please install it with: `pip install kaleido==0.2.*`"
                )
        # At this stage, everything should go smooth with the export format
        else:
            # Return supported format
            return format
    # Raise error if unsupported format
    raise ValueError(
        f"Format `{format}` is not supported.{os.linesep}"
        f"Supported formats: {SUPPORTED_OUTPUT_FORMATS}."
        )


def main(
    run_dir: FilePath,
    modules: list[int],
    top_clusters: int,
    format: Optional[ImgFormat],
    scale: Optional[float],
    inter: Optional[bool],
    is_cleaned: Optional[bool],
    offline: bool = False,
    mode: Optional[str] = None,
    ncores: Optional[int] = None,
    self_contained: bool = False,
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
    top_clusters : int
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
    self_contained : bool
        Should the analysis directory contain the models?
    """
    log.level = 20
    log.info(
        f"Running haddock3-analyse on {run_dir}, modules {modules}, "
        f"with top_clusters = {top_clusters}"
    )
    # Validate output format
    format = validate_format(format)

    # Obtain starting working directory
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
        target_path = Path(subfolder_name)

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
        try:
            analyse_step(
                step,
                Path("./"),
                capri_dict,
                target_path,
                top_clusters,
                format,
                scale,
                is_cleaned,
                offline=offline,
                mode=mode,
                ncores=ncores,
                self_contained=self_contained,
                )
        except Exception as e:
            log.warning(
                f"Could not execute the analysis for step {step}. "
                f"The following error occurred: {e}"
                )
            bad_folder_paths.append(target_path)
        else:
            good_folder_paths.append(target_path)

        # going back
        os.chdir(rundir_cwd)

    # moving files into analysis folder
    if good_folder_paths != []:
        log.info("moving files to analysis folder")
        urls: list[str] = []
        for directory in good_folder_paths:
            shutil.move(directory, outdir)
            url = f"- [{Path(directory, 'report.html')}](http://0.0.0.0:8000/{Path(directory, 'report.html')}) "  # noqa : E501
            urls.append(url)
        
        # Adding instructions on how to setup the server
        readme_fpath = Path(outdir, "README.md")
        readme_fpath.write_text(
            f"# Usage{os.linesep}{os.linesep}"
            "To view structures or download the structure files, "
            f"in a terminal run the command:{os.linesep}```bash{os.linesep}"
            f"python -m http.server --directory .{os.linesep}```{os.linesep}"
            f"And open the link following links in a web browser:{os.linesep}"
            f"{os.linesep.join(urls)}{os.linesep}"
            )
        assert readme_fpath.exists()

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
            f"in a terminal run the command: {os.linesep}"
            f">python -m http.server --directory {rundir_cwd}{os.linesep}"
            # "By default, http server runs on `http://0.0.0.0:8000/`. "
            f"And open the link http://0.0.0.0:8000/{report_file} "
            "in a web browser."
            )
        log.info(info_msg)
    os.chdir(ori_cwd)
    return


if __name__ == "__main__":
    sys.exit(maincli())  # type: ignore
