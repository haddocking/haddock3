"""Backmapping of Coarse-Grained structures into All-atom structures with CNS.

The ``[cgtoaa]`` module translate the CG conformations into AA representations, 
implemented in CNS.

For this module to be functional, it needs to be run in a workflow where ``[topocg]``
is present upstream.

For more details about this module, please `refer to the haddock3 user manual
<https://www.bonvinlab.org/haddock3-user-manual/modules/refinement.html#cgtoaa-module>`_
"""

from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath
from haddock.gear.haddockmodel import HaddockModel
from haddock.libs.libcns import prepare_cns_input, prepare_expected_pdb
from haddock.libs.libontology import PDBFile
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import get_engine
from haddock.modules.base_cns_module import BaseCNSModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseCNSModule):
    """HADDOCK3 module energy minimization refinement."""

    name = RECIPE_PATH.name

    def __init__(
        self, order: int, path: Path, initial_params: FilePath = DEFAULT_CONFIG
    ) -> None:
        """."""
        cns_script = Path(RECIPE_PATH, "cns", "cgtoaa.cns")
        super().__init__(order, path, initial_params, cns_script=cns_script)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm module is installed."""
        return

    def _run(self) -> None:
        """Execute module."""
        # Get the models generated in previous step
        try:
            models_to_refine = self.previous_io.retrieve_models(individualize=True)
        except Exception as e:
            self.finish_with_error(e)
        self.output_models = []
        sampling_factor = self.params["sampling_factor"]
        if sampling_factor == 0:
            self.log("[Warning] sampling_factor cannot be 0, setting it to 1")
            sampling_factor = 1
        elif sampling_factor > 100:
            self.log("[Warning] sampling_factor is larger than 100")

        max_nmodels = self.params["max_nmodels"]
        nmodels = len(models_to_refine) * sampling_factor
        if nmodels > max_nmodels:
            self.finish_with_error(
                f"Too many models ({nmodels}) to refine, max_nmodels ="
                f" {max_nmodels}. Please reduce the number of models or"
                " decrease the sampling_factor."
            )

        # Pool of jobs to be executed by the CNS engine
        jobs: list[CNSJob] = []
        idx = 1
        for model in models_to_refine:
            if not model.seed: model.seed = self.params["iniseed"]
            for s_ind in range(sampling_factor):
                # Prepare CNS input
                cgtoaa_input = prepare_cns_input(
                    idx,
                    model,
                    self.path,
                    self.recipe_str,
                    self.params,
                    self.name,
                    native_segid=True,
                    debug=self.params["debug"],
                    seed=(model.seed + s_ind) if isinstance(model, PDBFile) else None,
                    cgtoaa=True,
                )
                # Build CNS Job
                out_file = f"{self.name}_{idx}.out"
                err_fname = f"{self.name}_{idx}.cnserr"
                job = CNSJob(cgtoaa_input, out_file, err_fname, envvars=self.envvars)
                jobs.append(job)

                # create the expected PDBobject
                expected_pdb = prepare_expected_pdb(model, idx, ".", self.name)
                expected_pdb.topology = expected_pdb.aa_topology
                try:
                    expected_pdb.ori_name = model.file_name
                except AttributeError:
                    expected_pdb.ori_name = None
                self.output_models.append(expected_pdb)

                idx += 1

        # Run CNS Jobs
        self.log(f"Running CNS Jobs n={len(jobs)}")
        Engine = get_engine(self.params["mode"], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("CNS jobs have finished")

        # Get the weights needed for the CNS module
        _weight_keys = ("w_vdw", "w_elec", "w_desolv", "w_air", "w_bsa")
        weights = {e: self.params[e] for e in _weight_keys}

        # Loop over expected PDB outputs
        for pdb in self.output_models:
            if pdb.is_present():
                # Compute score
                haddock_model = HaddockModel(pdb.file_name)
                haddock_score = haddock_model.calc_haddock_score(**weights)
                # Hold score as attribute
                pdb.unw_energies = haddock_model.energies
                pdb.score = haddock_score

        self.export_io_models(faulty_tolerance=self.params["tolerance"])
