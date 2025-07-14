"""
HADDOCK3 module for alanine scan.

This module is responsible for the alanine (or any other residue) scan analysis
of the model(s) generated in the previous step of the workflow.
For each model, this module will mutate the interface residues and calculate
the haddock score differences between the wild type and the mutant, thus providing
a measure of the impact of such mutation. Such difference (delta_score) is always
calculated as:

    delta_score = score_wildtype - score_mutant

Therefore, a _positive_ delta_score indicates that the mutation is destabilizing
while a _negative_ delta_score indicates that the mutation is stabilizing.

If cluster information is available, the module will also calculate the
average haddock score difference for each cluster of models. For each amino acid,
a Z score is calculated as:

    Z = (delta_score - mean) / std

where mean and std are the mean and standard deviation of the delta_score over 
all the amino acids.

The module will also generate plots of the alanine scan data, showing the
distribution of the delta_score (and every component) for each amino acid at the
interface.

You can use the parameters below to customize the behavior of the module:

    * `chains`: list of chains to be considered for the alanine scan. In some
      cases you may want to limit the analysis to a single chain.
    * `output_mutants`: if True, the module will output the models with the
      mutations applied (only possible if there is only one model)
    * `output_bfactor`: if True, the module will output the non-mutated models
      with the rescaled delta_score in the B-factor column
    * `plot`: if True, the module will generate plots of the alanine scan data
    * `scan_residue`: the residue to scan (default is 'ALA')
    * `resdic`: list of residues to be used for the scanning. An example is:

    >>> resdic_A = [1,2,3,4]
    >>> resdic_B = [2,3,4]
"""
import os
from pathlib import Path

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.libs.libparallel import get_index_list
from haddock.libs.libparallel import AlascanScheduler
from haddock.libs.libutil import parse_ncores
from haddock.modules import BaseHaddockModule
#from haddock.modules import get_engine
#from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.alascan.scan import (
    Scan,
    ScanJob,
    alascan_cluster_analysis,
    create_alascan_plots,
    generate_alascan_output,
    calculate_core_allocation,
    )


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for alanine scan."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, *ignore, init_params=DEFAULT_CONFIG,
                 **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm if module is installed."""
        return

    def _run(self):
        """Execute module."""
        # Get the models generated in previous step
        try:
            models = self.previous_io.retrieve_models(individualize=True)
        except Exception as e:
            self.finish_with_error(e)
        # Parallelization : optimal dispatching of models
        nmodels = len(models)
        # Outputting mutants is possible if there is just one model
        if self.params["output_mutants"]:
            if nmodels != 1:
                log.warning(
                    "output_mutants is set to True, but more than one model "
                    "was found. Setting output_mutants to False."
                    )
                self.params["output_mutants"] = False

        # Core maths: Is there enough cores to parallelize per-residue after 
        # giving enough cores per model ?
        ncores = parse_ncores(n=self.params['ncores'])
        model_cores, residue_cores_per_model = calculate_core_allocation(nmodels, ncores)
        index_list = get_index_list(nmodels, model_cores)
        
        if residue_cores_per_model == 1:
            log.info(f"Using {model_cores} core(s) to process {nmodels} model(s)")
        else:
            log.info(f"Processing {model_cores} model(s), using {residue_cores_per_model} cores per model")
        
        alascan_jobs = []
        for core in range(model_cores):
            models_for_core = models[index_list[core]:index_list[core + 1]]
            scan_obj = Scan(
                model_list=models_for_core,
                core=core,
                residue_ncores=residue_cores_per_model,
                path=Path("."),
                params=self.params,
            )
            # running now the ScanJob
            # init ScanJob
            job = ScanJob(
                self.params,
                scan_obj,
            )
            alascan_jobs.append(job)

        log.info(f"Created {len(alascan_jobs)} scan jobs")
        # next libutil and libparallel will log info about per-model cores/tasks.
        # This is misleading, if per-residue parallelization is present.
        # This log makes log look more coherent, in a way.
       
        ## Actually. This can cause semaphore leaks
        # log.info(f"Model-level parallelization:")

        # exec_mode = get_analysis_exec_mode(self.params["mode"])

        # Engine = get_engine(exec_mode, self.params)
        # engine = Engine(alascan_jobs)
        # engine.run()
        ## So replacing it with AlascanScheduler

        log.info("Using AlascanScheduler for model-level parallelization")
        #try:
        with AlascanScheduler(alascan_jobs, ncores=model_cores) as scheduler:
            scheduler.run()

        #except Exception as e:
        #    self.finish_with_error(f"Alascan model-level parallelization failed: {e}")

        # cluster-based analysis
        clt_alascan = alascan_cluster_analysis(models)
        # now plot the data
        if self.params["plot"] is True:
            create_alascan_plots(
                clt_alascan,
                self.params["scan_residue"],
                offline=self.params["offline"],
            )
        # if output_bfactor is true, write the models and export them
        if self.params["output_bfactor"] is True:
            models_to_export = generate_alascan_output(models, self.path)
            self.output_models = models_to_export
        else:
            # Send models to the next step, no operation is done on them
            self.output_models = models
        self.export_io_models()