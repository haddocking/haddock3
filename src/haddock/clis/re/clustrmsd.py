"""haddock3-re clustrmsd subcommand."""
from pathlib import Path
import numpy as np

from haddock import log
from haddock.core.defaults import INTERACTIVE_RE_SUFFIX
from haddock.core.typing import Union, Optional
from haddock.gear.config import load as read_config
from haddock.gear.config import save as save_config
from haddock.modules import get_module_steps_folders
from haddock.libs.libclust import (
    add_cluster_info,
    clustrmsd_tolerance_params,
    get_cluster_matrix_plot_clt_dt,
    plot_cluster_matrix,
    rank_clusters,
    write_structure_list,
    )
from haddock.libs.libinteractive import look_for_capri, rewrite_capri_tables
from haddock.libs.libontology import ModuleIO
from haddock.modules.analysis.clustrmsd.clustrmsd import (
    get_clusters,
    get_matrix_path,
    iterate_min_population,
    order_clusters,
    write_clusters,
    write_clustrmsd_file,
    )


def add_clustrmsd_arguments(clustrmsd_subcommand):
    """Add arguments to the clustrmsd subcommand."""
    clustrmsd_subcommand.add_argument(
        "clustrmsd_dir",
        help="The clustrmsd directory to recluster.",
        )
    
    clustrmsd_subcommand.add_argument(
        "-n",
        "--n_clusters",
        help="number of clusters to generate.",
        required=False,
        type=int,
        )
    
    clustrmsd_subcommand.add_argument(
        "-d",
        "--clust_cutoff",
        help="clustering cutoff distance.",
        required=False,
        type=float,
        )
    
    clustrmsd_subcommand.add_argument(
        "-t",
        "--min_population",
        help="minimum cluster population.",
        required=False,
        type=int,
        )
    
    clustrmsd_subcommand.add_argument(
        "-p",
        "--plot_matrix",
        help="Generate the matrix plot with the clusters.",
        required=False,
        default=False,
        action='store_true',
        )

    return clustrmsd_subcommand


def reclustrmsd(
        clustrmsd_dir: str,
        n_clusters: Union[bool, int] = None,
        clust_cutoff: Union[bool, float] = None,
        min_population: Union[bool, int] = None,
        plot_matrix: bool = True,
        ) -> Path:
    """
    Recluster the models in the clustrmsd directory.
    
    Parameters
    ----------
    clustrmsd_dir : str
        Path to the clustrmsd directory.
    
    n_clusters : Union[bool, int]
        Number of clusters to generate.
    
    clust_cutoff : Union[bool, float]
        Clustering cutoff distance.
    
    min_population : Union[bool, int]
        Cluster population min_population.
    
    plot_matrix : bool
        Should the corresponding matrix plot be generated.
    
    Returns
    -------
    outdir : Path
        Path to the interactive directory.
    """
    log.info(f"Reclustering {clustrmsd_dir}")

    run_dir = Path(clustrmsd_dir).parent
    clustrmsd_name = Path(clustrmsd_dir).name
    # create the interactive folder
    outdir = Path(run_dir, f"{clustrmsd_name}_{INTERACTIVE_RE_SUFFIX}")
    outdir.mkdir(exist_ok=True)

    # create an io object
    io = ModuleIO()
    filename = Path(clustrmsd_dir, "io.json")
    io.load(filename)
    models = io.input

    # load the original clustering parameters via json
    clustrmsd_params = read_config(Path(clustrmsd_dir, "params.cfg"))
    key = list(clustrmsd_params['final_cfg'].keys())[0]
    clustrmsd_params = clustrmsd_params['final_cfg'][key]
    log.info(f"Previous clustering parameters: {clustrmsd_params}")

    # setting previous tolerance, just in case no new parameters are given
    tolerance_param_name, tolerance = clustrmsd_tolerance_params(
        clustrmsd_params,
        )

    # adjust the parameters
    if n_clusters is not None:
        clustrmsd_params["n_clusters"] = n_clusters
        clustrmsd_params["criterion"] = "maxclust"
        tolerance = n_clusters
    else:
        if clust_cutoff is not None:
            clustrmsd_params["clust_cutoff"] = clust_cutoff
            clustrmsd_params["criterion"] = "distance"
            tolerance = clust_cutoff
    
    if min_population is not None:
        clustrmsd_params["min_population"] = min_population

    clustrmsd_params["plot_matrix"] = plot_matrix

    log.info(
        f"Clustering with {tolerance_param_name} = {tolerance}, "
        f"and criterion {clustrmsd_params['criterion']}"
        )
    
    # load the clustering dendrogram
    dendrogram = np.loadtxt(Path(clustrmsd_dir, "dendrogram.txt"))

    # get the clusters
    cluster_arr = get_clusters(
        dendrogram,
        tolerance,
        clustrmsd_params["criterion"],
        )
    log.info(f"clusters {cluster_arr}")

    if clustrmsd_params['criterion'] == "distance":
        cluster_arr, min_population = iterate_min_population(
            cluster_arr,
            clustrmsd_params['min_population']
            )
        clustrmsd_params['min_population'] = min_population
    log.info(f"Updated clustering parameters = {clustrmsd_params}")
    
    # processing the clusters
    clusters, cluster_arr = order_clusters(cluster_arr)
    log.info(f"clusters = {clusters}")
    log.info(f"cluster_arr = {cluster_arr}")

    clt_dic, cluster_centers = write_clusters(
        clusters,
        cluster_arr,
        models,
        out_filename=Path(outdir, "cluster.out"),
        rmsd_matrix=None,
        centers=False
        )

    score_dic, sorted_score_dic = rank_clusters(
        clt_dic,
        clustrmsd_params["min_population"]
        )
    
    output_models = add_cluster_info(sorted_score_dic, clt_dic)
    
    write_structure_list(
        models,
        output_models,
        out_fname=Path(outdir, "clustrmsd.tsv"),
        )
    
    write_clustrmsd_file(
        clusters,
        clt_dic,
        cluster_centers,
        score_dic,
        sorted_score_dic,
        clustrmsd_params,
        output_fname=Path(outdir, "clustrmsd.txt"),
        )
    
    # Draw the matrix
    if clustrmsd_params["plot_matrix"]:
        if not (matrix_json_path := search_previousstep_matrix(clustrmsd_dir)):
            log.warn(
                "Could not find the rmsd matrix in previous step."
                " Unable to produce a graph out of it!"
                )
        else:
            log.info("Generating graphical representation of the clusters.")
            matrix_io = ModuleIO()
            matrix_io.load(matrix_json_path)
            # Obtain final models indices
            final_order_idx, labels, cluster_ids = [], [], []
            for pdb in output_models:
                final_order_idx.append(models.index(pdb))
                labels.append(pdb.file_name.replace('.pdb', ''))
                cluster_ids.append(pdb.clt_id)
            # Get custom cluster data
            matrix_cluster_dt, cluster_limits = get_cluster_matrix_plot_clt_dt(
                cluster_ids
                )
            # Define output filename
            html_matrix_basepath = Path(outdir, 'rmsd_matrix')
            # Plot matrix
            html_matrixpath = plot_cluster_matrix(
                get_matrix_path(matrix_io.input[0]),
                final_order_idx,
                labels,
                dttype='RMSD(Å)',
                reverse=True,
                diag_fill=0,
                output_fname=html_matrix_basepath,
                matrix_cluster_dt=matrix_cluster_dt,
                cluster_limits=cluster_limits,
                )
            log.info(f"Plotting matrix in {html_matrixpath}")

    # save the io.json file
    io.save(outdir)

    # save the updated parameters in a json file
    save_config(clustrmsd_params, Path(outdir, "params.cfg"))

    # analysis
    clustrmsd_id = int(clustrmsd_name.split("_")[0])
    caprieval_folder = look_for_capri(run_dir, clustrmsd_id)
    if caprieval_folder:
        log.info("Rewriting capri tables")
        rewrite_capri_tables(caprieval_folder, clt_dic, outdir)

    return outdir


def search_previousstep_matrix(clustrmsd_dir: str) -> Optional[Path]:
    """Retrieve the path of the previous step matrix_json file.

    Parameters
    ----------
    clustrmsd_dir : str
        Path to the clustrmsd directory.

    Returns
    -------
    matrix_json : Optional[Path]
        Path to the matrix_json file.
    """
    # Compute previous step index
    previous_step_ind = int(str(Path(clustrmsd_dir).name).split('_')[0]) - 1
    workflow_dir = Path(clustrmsd_dir).parent
    # Try to get previous step directory name
    try:
        previous_steps = get_module_steps_folders(
            workflow_dir,
            [previous_step_ind],
            )
        previous_step = previous_steps[0]
    except IndexError:
        return None
    else:
        matrix_json = Path(workflow_dir, previous_step, "rmsd_matrix.json")
        if matrix_json.exists():
            return matrix_json
