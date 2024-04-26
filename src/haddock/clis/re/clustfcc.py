"""haddock3-re clustfcc subcommand."""

from pathlib import Path
import shutil
from fcc.scripts import cluster_fcc


from haddock import log
from haddock.core.defaults import INTERACTIVE_RE_SUFFIX
from haddock.core.typing import Union
from haddock.gear.config import load as read_config
from haddock.gear.config import save as save_config
from haddock.libs.libclust import (
    add_cluster_info,
    rank_clusters,
    write_structure_list,
    )
from haddock.libs.libontology import ModuleIO
from haddock.libs.libinteractive import look_for_capri, rewrite_capri_tables
from haddock.modules.analysis.clustfcc.clustfcc import (
    get_cluster_centers,
    iterate_clustering,
    write_clusters,
    write_clustfcc_file,
    )


def add_clustfcc_arguments(clustfcc_subcommand):
    """Add arguments to the clustfcc subcommand."""
    clustfcc_subcommand.add_argument(
        "clustfcc_dir",
        help="The clustfcc directory to recluster.",
        )

    clustfcc_subcommand.add_argument(
        "-f",
        "--clust_cutoff",
        help="Minimum fraction of common contacts to be considered in a cluster.",  # noqa: E501
        required=False,
        type=float,
        )

    clustfcc_subcommand.add_argument(
        "-s",
        "--strictness",
        help="Strictness factor.",
        required=False,
        type=float,
        )
    
    clustfcc_subcommand.add_argument(
        "-t",
        "--min_population",
        help="Clustering population threshold.",
        required=False,
        type=int,
        )
    
    return clustfcc_subcommand


def reclustfcc(
        clustfcc_dir: str,
        clust_cutoff: Union[bool, float] = None,
        strictness: Union[bool, float] = None,
        min_population: Union[bool, int] = None,
        ) -> Path:
    """
    Recluster the models in the clustfcc directory.
    
    Parameters
    ----------
    clustfcc_dir : str
        Path to the clustfcc directory.
    
    clust_cutoff : Union[bool, float]
        Fraction of common contacts to not be considered a singleton model.
    
    strictness : Union[bool, float]
        Fraction of common contacts to be considered to be part of the same
         cluster.
    
    min_population : Union[bool, int]
        Minimum cluster population.
    
    Returns
    -------
    outdir : Path
        Path to the interactive directory.
    """
    log.info(f"Reclustering {clustfcc_dir}")

    # create the interactive folder
    run_dir = Path(clustfcc_dir).parent
    clustfcc_name = Path(clustfcc_dir).name
    outdir = Path(run_dir, f"{clustfcc_name}_{INTERACTIVE_RE_SUFFIX}")
    outdir.mkdir(exist_ok=True)

    # create an io object
    io = ModuleIO()
    filename = Path(clustfcc_dir, "io.json")
    io.load(filename)
    models = io.input
    # copying io.json to the new directory
    shutil.copy(filename, Path(outdir, "io.json"))

    # load the original clustering parameters via json
    clustfcc_params = read_config(Path(clustfcc_dir, "params.cfg"))
    key = list(clustfcc_params['final_cfg'].keys())[0]
    clustfcc_params = clustfcc_params['final_cfg'][key]
    log.info(f"Previous clustering parameters: {clustfcc_params}")

    # adjust the parameters
    if clust_cutoff is not None:
        clustfcc_params["clust_cutoff"] = clust_cutoff
    if strictness is not None:
        clustfcc_params["strictness"] = strictness
    if min_population is not None:
        clustfcc_params["min_population"] = min_population
    
    # load the fcc matrix
    pool = cluster_fcc.read_matrix(
        Path(clustfcc_dir, "fcc.matrix"),
        clustfcc_params['clust_cutoff'],
        clustfcc_params['strictness'],
        )
    
    # iterate clustering until at least one cluster is found
    clusters, min_population = iterate_clustering(
        pool,
        clustfcc_params['min_population']
        )
    clustfcc_params['min_population'] = min_population
    log.info(f"Updated clustering parameters: {clustfcc_params}")

    # Prepare output and read the elements
    clt_dic = {}
    if clusters:
        write_clusters(clusters, out_filename=Path(outdir, "cluster.out"))

        # Get the cluster centers
        clt_dic, clt_centers = get_cluster_centers(clusters, models)

        _score_dic, sorted_score_dic = rank_clusters(clt_dic, min_population)

        output_models = add_cluster_info(sorted_score_dic, clt_dic)
        
        # Write unclustered structures
        write_structure_list(models,
                             output_models,
                             out_fname=Path(outdir, "clustfcc.tsv"))
        
        write_clustfcc_file(
            clusters,
            clt_centers,
            clt_dic,
            clustfcc_params,
            sorted_score_dic,
            output_fname=Path(outdir, 'clustfcc.txt')
            )

        save_config(clustfcc_params, Path(outdir, "params.cfg"))

        # analysis
        clustfcc_id = int(clustfcc_name.split("_")[0])
        caprieval_folder = look_for_capri(run_dir, clustfcc_id)
        if caprieval_folder:
            log.info("Rewriting capri tables")
            rewrite_capri_tables(caprieval_folder, clt_dic, outdir)
            
    return outdir
