"""optint wrapper."""
import itertools

from haddock import log
from haddock.modules.design.optint.modules.ga import GA


def run(pdb_list, params):
    """
    `optint` wrapper.

    Parameters
    ----------
    pdb_list : list
        List of PDB files to be used as input.
    params : dict
        Parameters for the module.

    Returns
    -------
    list
        Nested list containing the output of the module
        [(pdb, psf, score), ... ]
    """
    log.info("Executing Interface Optimization Module")
    output_list = []

    for i, pdb in enumerate(pdb_list, start=1):
        if isinstance(pdb, tuple):
            pdb = pdb[0]

        log.info(f"Optimizing {pdb.rel_path.name} ID: {i}")
        ga_worker = GA(i, pdb, params)
        ga_worker.run()
        output_models = ga_worker.retrieve_models()
        output_list.append(output_models)

    return list(itertools.chain(*output_list))
