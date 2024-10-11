"""HADDOCK3 module for alanine scan.

This module is responsible for the alanine scan analysis of the models
generated in the previous step of the workflow. For each model, the module
will mutate the interface residues and calculate the energy differences
between the wild type and the mutant, thus providing a measure of the impact
of such mutation.

If cluster information is available, the module will also calculate the
average energy difference for each cluster of models.
"""
from pathlib import Path

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.libs.libparallel import get_index_list
from haddock.libs.libutil import parse_ncores
from haddock.modules import BaseHaddockModule
from haddock.modules import get_engine
from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.alascan.scan import (
    Scan,
    ScanJob,
    alascan_cluster_analysis,
    create_alascan_plots,
    generate_alascan_output,
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
        # Parallelisation : optimal dispatching of models
        nmodels = len(models)
        ncores = parse_ncores(n=self.params['ncores'], njobs=nmodels)

        log.info(f"Running on {ncores} cores")

        index_list = get_index_list(nmodels, ncores)

        alascan_jobs = []
        for core in range(ncores):
            output_name = "alascan_" + str(core) + ".scan"
            scan_obj = Scan(
                model_list=models[index_list[core]:index_list[core + 1]],
                output_name=output_name,
                core=core,
                path=Path("."),
                params=self.params,
                )
            # running now the ScanJob
            job_f = Path(output_name)
            # init ScanJob
            job = ScanJob(
                job_f,
                self.params,
                scan_obj,
                )
            alascan_jobs.append(job)

        exec_mode = get_analysis_exec_mode(self.params["mode"])

        Engine = get_engine(exec_mode, self.params)
        engine = Engine(alascan_jobs)
        engine.run()
        
        # cluster-based analysis
        clt_alascan = alascan_cluster_analysis(models)
        # now plot the data
        if self.params["plot"] is True:
            create_alascan_plots(
                clt_alascan,
                self.params["scan_residue"],
                offline=self.params["offline"],
                )
        # if output is true, write the models and export them
        if self.params["output"] is True:
            models_to_export = generate_alascan_output(models, self.path)
            self.output_models = models_to_export
        else:
            # Send models to the next step,
            #  no operation is done on them
            self.output_models = models
        self.export_io_models()
