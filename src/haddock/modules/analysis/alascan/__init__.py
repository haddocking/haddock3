"""
HADDOCK3 module for alanine scan.

This module is responsible for the alanine (or any other residue) scan analysis
of the model(s) generated in the previous step of the workflow.
For each model, this module will mutate the interface residues and calculate 
the differences in the haddock score and its individual components between the
wild type and the mutant, thus providing a measure of the impact of such mutation.
Such difference (delta_score) is always calculated as:

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

from pathlib import Path

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.libs.libparallel import GenericTask
from haddock.modules import BaseHaddockModule, get_engine
from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.alascan.scan import (
    AddDeltaBFactor,
    ClusterOutputer,
    group_scan_by_cluster,
    InterfaceScanner,
    write_scan_out,
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
    
    def validate_ouput_mutant_parameter(self, nmodels: int) -> None:
        """Validate the output mutant parameter.
        
        This parameter can be set to True if only one input model is provided,
        otherewise we risk to generate too many PDB files.

        Parameters
        ----------
        nmodels: int
            Number of input models.
        """
        if self.params["output_mutants"]:
            # output mutants is only possible if there is only one model
            if nmodels > 1:
                log.warning(
                    "'output_mutants' parameter is set to True, "
                    "but more than one model was found. "
                    "Setting 'output_mutant' parameter to False."
                    )
                self.params["output_mutants"] = False

    def _run(self):
        """Execute module."""
        # Get the models generated in previous step
        try:
            models = self.previous_io.retrieve_models(individualize=True)
        except Exception as e:
            self.finish_with_error(e)
        
        # Compute number of input model
        nmodels = len(models)
        
        # Validate `output_mutant` parameter
        self.validate_ouput_mutant_parameter(nmodels)
        
        # Step1: "get mutations" i.e. get target interface residues per input model
        # 1 scan_obj per input model, merged into scan_objects to give to Engine  
        scan_objects = [
            InterfaceScanner(
                mutation_res=self.params["scan_residue"],
                model=model,
                params=self.params,
                library_mode = False
                )
            for model in models
            ]

        log.info(f"Scanning {nmodels} models for possible mutations")
        exec_mode = get_analysis_exec_mode(self.params["mode"])
        Engine = get_engine(exec_mode, self.params)
        engine = Engine(scan_objects)
        engine.run()

        # Step2: perform mutations
        # Collect all point mutations to be performed
        mutation_objects = []
        for mutations_to_perform in engine.results:
            if mutations_to_perform:
                mutation_objects.extend(mutations_to_perform) 

        total_mutations = len(mutation_objects)
        log.info(f"Found {total_mutations} mutations")

        if mutation_objects:
            # let engine take care of parallelization  
            engine = Engine(mutation_objects)
            engine.run()

            # Organize engine output by model
            results_by_model = {}
            for result in engine.results:
                if result and result.success:
                    if result.model_id not in results_by_model:
                        results_by_model[result.model_id] = []
                    results_by_model[result.model_id].append(result)

            # Save to .tsv
            scan_writter_jobs = [
                GenericTask(write_scan_out, results, model_id)
                for model_id, results in results_by_model.items()
            ]
            engine = Engine(scan_writter_jobs)
            engine.run()

            # Generate output models with bfactors if requested  
            if self.params["output_bfactor"]:
                update_with_bfactor_jobs = []
                for model in models:
                    model_id = model.file_name.removesuffix(".pdb")
                    try:
                        model_results = results_by_model[model_id]
                    except KeyError:
                        # Case when no data computed for this model
                        model_results = []
                    update_with_bfactor_jobs.append(
                        AddDeltaBFactor(model, self.path, model_results)
                        )
                engine = Engine(update_with_bfactor_jobs)
                engine.run()
                models_to_export = engine.results
                self.output_models = models_to_export
            else:
                # Send models to the next step, no operation is done on them 
                self.output_models = models

            # Cluster-based analysis
            clt_scan, clt_pops = group_scan_by_cluster(models, results_by_model)
            alascan_cluster_jobs = [
                ClusterOutputer(
                    clt_data,
                    clt_id,
                    clt_pops[clt_id],
                    scan_residue=self.params["scan_residue"],
                    generate_plot=self.params["plot"],
                    offline=self.params["offline"],
                )
                for clt_id, clt_data in clt_scan.items()
            ]
            engine = Engine(alascan_cluster_jobs)
            engine.run()

        else:
            log.info("No interface residues found - skipping mutation analysis")
            # Send models to the next step, no operation is done on them
            self.output_models = models

        self.export_io_models()