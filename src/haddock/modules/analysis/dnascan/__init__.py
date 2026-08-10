"""
HADDOCK3 module for DNA base-pair scan.

This module performs a mutagenesis scan of DNA base pairs for the model(s)
generated in the previous step of the workflow. DNA is double stranded and its
bases are engaged in Watson-Crick base pairs (A:T and G:C). Mutating a single
base in isolation would break the base pair and is not physically meaningful,
so every mutation performed by this module is a *double* mutation: for each
selected interface nucleotide, the nucleotide is mutated into a target base and
its base-pairing partner on the complementary strand is simultaneously mutated
into the complementary base, so that a valid Watson-Crick pair is preserved
(e.g. A:T -> G:C). The differences in the haddock score and its individual
components between the wild type and each mutant base pair are calculated, thus
providing a measure of the impact of such mutations. By default, all four
possible bases are tested for each selected interface nucleotide (the wild-type
base is skipped). Such difference (delta_score) is always calculated as:

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

The module will also generate plots of the DNA scan data, showing the
distribution of the delta_score (and every component) for each base-pair
mutation at the interface.

You can use the parameters below to customize the behavior of the module:

    * `scan_bases`: list of target bases to be tested for each nucleotide.
      By default all four DNA bases (DA, DC, DG, DT) are tested.
    * `bp_cutoff`: distance cutoff (Å) used to detect Watson-Crick base pairs
      from the distance between the hydrogen-bonding ring nitrogens.
    * `chains`: list of chains to be considered for the DNA scan. In some
      cases you may want to limit the analysis to a single chain.
    * `output_mutants`: if True, the module will output the models with the
      mutations applied (only possible if there is only one model)
    * `output_bfactor`: if True, the module will output the non-mutated models
      with the rescaled delta_score in the B-factor column
    * `plot`: if True, the module will generate plots of the DNA scan data
    * `splitplot`: if True, the scan plot shows one panel per energy component;
      if False (default) all components are overlaid in a single panel
    * `resdic`: list of residues to be used for the scanning. An example is:

    >>> resdic_A = [1,2,3,4]
    >>> resdic_B = [2,3,4]

Only nucleic acid (DNA) residues at the interface are mutated; protein and
other residues are ignored. Interface nucleotides for which no base-pairing
partner can be found are skipped, since a double (base-pair) mutation cannot
be constructed for them.

When applying a mutation, the deoxyribose-phosphate backbone is always kept for
both nucleotides of the pair. For mutations between bases of the same ring type
(purine <-> purine, i.e. DA <-> DG, or pyrimidine <-> pyrimidine, i.e.
DC <-> DT), the shared base ring atoms are also kept so that the base
orientation is preserved and CNS only rebuilds the differing substituents; the
whole base pair is then scored in a single CNS call. For cross-type mutations,
the two nucleotides swap ring type (one purine -> pyrimidine and its partner
pyrimidine -> purine). Rebuilding a whole purine ring from its three anchor
atoms in the same pass as the partner mutation is unreliable, so such a base pair
is
mutated in two sequential, CNS-regularised steps: first the purine ->
pyrimidine mutation is applied and energy-minimised by CNS, then the
pyrimidine -> purine mutation is applied to that minimised intermediate and
scored by CNS. The score and energies of the second (final) CNS call are the
ones reported and plotted. In both cases the glycosidic-region anchor atoms are
kept and, for cross-type mutations, renamed to their counterpart in the target
ring system (pyrimidine N1/C2/C6 <-> purine N9/C4/C8).

To keep the delta_score comparison consistent, each mutant is compared against a
wild-type baseline that went through the same number of CNS minimisation passes:
a one-pass baseline for same-ring-type mutants and a two-pass baseline (the wild
type minimised and then re-scored) for cross-ring-type mutants.
"""

from collections import defaultdict
from pathlib import Path

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.libs.libparallel import GenericTask
from haddock.modules import BaseHaddockModule, get_engine
from haddock.modules.analysis import get_analysis_exec_mode
from haddock.modules.analysis.dnascan.dnascan import (
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
    """HADDOCK3 module for DNA base-pair scan."""

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
        # output mutants is only possible if there is a single input model
        if self.params["output_mutants"] and nmodels > 1:
            log.warning(
                "'output_mutants' parameter is set to True, "
                "but more than one model was found. "
                "Setting 'output_mutant' parameter to False."
            )
            self.params["output_mutants"] = False

    def _run(self):
        """Execute module."""
        # Validate and normalise the requested target DNA bases
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

        # Step1: "get mutations" i.e. get target interface base pairs per input model
        # 1 scan_obj per input model, merged into scan_objects to give to Engine
        scan_objects = [
            InterfaceScanner(
                scan_bases=self.params["scan_bases"],
                model=model,
                params=self.params,
            )
            for model in models
        ]

        log.info(f"Scanning {nmodels} models for possible base-pair mutations")
        exec_mode = get_analysis_exec_mode(self.params["mode"])
        Engine = get_engine(exec_mode, self.params)
        engine = Engine(scan_objects)
        engine.run()

        # Step2: perform mutations
        # Collect all base-pair mutations to be performed
        mutation_objects = []
        for mutations_to_perform in engine.results:
            if mutations_to_perform:
                mutation_objects.extend(mutations_to_perform)

        total_mutations = len(mutation_objects)
        log.info(f"Found {total_mutations} base-pair mutations")

        if not mutation_objects:
            log.info("No interface base pairs found - skipping mutation analysis")
            # Send models to the next step, no operation is done on them
            self.output_models = models
            self.export_io_models()
            return

        # let engine take care of parallelization
        engine = Engine(mutation_objects)
        engine.run()

        # Organize engine output by model
        results_by_model = defaultdict(list)
        for result in engine.results:
            if result and result.success:
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
                # empty list when no data was computed for this model
                model_results = results_by_model.get(model_id, [])
                update_with_bfactor_jobs.append(
                    AddDeltaBFactor(model, self.path, model_results)
                )
            engine = Engine(update_with_bfactor_jobs)
            engine.run()
            self.output_models = engine.results
        else:
            # Send models to the next step, no operation is done on them
            self.output_models = models

        # Cluster-based analysis
        clt_scan, clt_pops = group_scan_by_cluster(models, results_by_model)
        dnascan_cluster_jobs = [
            ClusterOutputer(
                clt_data,
                clt_id,
                clt_pops[clt_id],
                scan_residue="DNA base pair",
                generate_plot=self.params["plot"],
                offline=self.params["offline"],
                splitplot=self.params["splitplot"],
            )
            for clt_id, clt_data in clt_scan.items()
        ]
        engine = Engine(dnascan_cluster_jobs)
        engine.run()

        self.export_io_models()
