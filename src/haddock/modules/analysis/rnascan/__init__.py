"""
HADDOCK3 module for RNA base scan.

This module performs a mutagenesis scan of RNA bases for the model(s)
generated in the previous step of the workflow. For each model, every selected
interface nucleic acid residue is mutated into each of the four RNA bases
(A, C, G and U) and the differences in the haddock score and its individual
components between the wild type and each mutant are calculated, thus providing
a measure of the impact of such mutations. By default, all four possible bases
are tested for every selected interface nucleotide (the wild-type base is
skipped). Such difference (delta_score) is always calculated as:

    delta_score = score_wildtype - score_mutant

As the score are typically negative with lower values being better,
a _positive_ delta_score indicates that the mutation is stabilizing
while a _negative_ delta_score indicates that the mutation is destabilizing.


If cluster information is available, the module will also calculate the
average haddock score difference for each cluster of models. For each mutation,
a Z score is calculated as:

    Z = (delta_score - mean) / std

where mean and std are the mean and standard deviation of the delta_score over
all the mutations.

The module will also generate plots of the RNA scan data, showing the
distribution of the delta_score (and every component) for each mutation at the
interface.

You can use the parameters below to customize the behavior of the module:

    * `scan_bases`: list of target bases to be tested for each nucleotide.
      By default all four RNA bases (A, C, G, U) are tested.
    * `chains`: list of chains to be considered for the RNA scan. In some
      cases you may want to limit the analysis to a single chain.
    * `output_mutants`: if True, the module will output the models with the
      mutations applied (only possible if there is only one model)
    * `output_bfactor`: if True, the module will output the non-mutated models
      with the rescaled delta_score in the B-factor column
    * `plot`: if True, the module will generate plots of the RNA scan data
    * `splitplot`: if True, the scan plot shows one panel per energy component;
      if False (default) all components are overlaid in a single panel
    * `resdic`: list of residues to be used for the scanning. An example is:

    >>> resdic_A = [1,2,3,4]
    >>> resdic_B = [2,3,4]

Only nucleic acid (RNA) residues at the interface are mutated; protein and
other residues are ignored.

When applying a mutation, the ribose-phosphate backbone is always kept. For
mutations between bases of the same ring type (purine <-> purine, i.e.
A <-> G, or pyrimidine <-> pyrimidine, i.e. C <-> U), the shared base ring
atoms are also kept so that the base orientation is preserved and CNS only
rebuilds the differing substituents. For cross-type mutations (purine <->
pyrimidine), the three glycosidic-region anchor atoms are kept and renamed to
their counterpart in the target ring system (pyrimidine N1/C2/C6 <-> purine
N9/C4/C8), so that the base orientation is preserved and CNS rebuilds the rest
of the new base.
"""

from pathlib import Path

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.libs.libparallel import GenericTask
from haddock.modules import BaseHaddockModule, get_engine
from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.rnascan.scan import (
    AddDeltaBFactor,
    ClusterOutputer,
    group_scan_by_cluster,
    InterfaceScanner,
    validate_scan_bases,
    write_scan_out,
)


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for RNA base scan."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, *ignore, init_params=DEFAULT_CONFIG, **everything):
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
        # Validate and normalise the requested target RNA bases
        try:
            self.params["scan_bases"] = validate_scan_bases(self.params["scan_bases"])
        except Exception as e:
            self.finish_with_error(e)

        # Get the models generated in previous step
        try:
            models = self.previous_io.retrieve_models(individualize=True)
        except Exception as e:
            self.finish_with_error(e)

        # Compute number of input model
        nmodels = len(models)

        # Validate `output_mutant` parameter
        self.validate_ouput_mutant_parameter(nmodels)

        # Step1: "get mutations" i.e. get target interface nucleotides per input model
        # 1 scan_obj per input model, merged into scan_objects to give to Engine
        scan_objects = [
            InterfaceScanner(
                scan_bases=self.params["scan_bases"],
                model=model,
                params=self.params,
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

        if not mutation_objects:
            log.info("No interface nucleotides found - skipping mutation analysis")
            # Send models to the next step, no operation is done on them
            self.output_models = models
            self.export_io_models()
            return

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
        rnascan_cluster_jobs = [
            ClusterOutputer(
                clt_data,
                clt_id,
                clt_pops[clt_id],
                scan_residue="RNA base",
                generate_plot=self.params["plot"],
                offline=self.params["offline"],
                splitplot=self.params["splitplot"],
            )
            for clt_id, clt_data in clt_scan.items()
        ]
        engine = Engine(rnascan_cluster_jobs)
        engine.run()

        self.export_io_models()
