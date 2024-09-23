"""Energy minimization refinement with CNS.

The ``[emref]`` module refine the input complexes by energy minimization using
conjugate gradient method implemented in CNS.

Coordinates of the energy minimized structures are saved, and each
complex is then evaluated using the HADDOCK scoring function.
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
        cns_script = Path(RECIPE_PATH, "cns", "emref.cns")
        super().__init__(order, path, initial_params, cns_script=cns_script)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm module is installed."""
        return

    def _run(self) -> None:
        """Execute module."""
        # Pool of jobs to be executed by the CNS engine
        jobs: list[CNSJob] = []

        # Get the models generated in previous step
        try:
            models_to_refine = self.previous_io.retrieve_models()
        except Exception as e:
            self.finish_with_error(e)

        self.output_models = []
        sampling_factor = self.params["sampling_factor"]
        if sampling_factor > 1:
            self.log(f"sampling_factor={sampling_factor}")
        if sampling_factor == 0:
            self.log("[Warning] sampling_factor cannot be 0, setting it to 1")
            sampling_factor = 1
        if sampling_factor > 100:
            self.log("[Warning] sampling_factor is larger than 100")

        max_nmodels = self.params["max_nmodels"]
        nmodels = len(models_to_refine) * sampling_factor
        if nmodels > max_nmodels:
            self.finish_with_error(
                f"Too many models ({nmodels}) to refine, max_nmodels ="
                f" {max_nmodels}. Please reduce the number of models or"
                " decrease the sampling_factor."
            )

        # checking the ambig_fname:
        try:
            prev_ambig_fnames = [mod.restr_fname for mod in models_to_refine]
        except Exception as e:  # noqa:F841
            # cannot extract restr_fname info from tuples
            prev_ambig_fnames = [None for model in models_to_refine]

        ambig_fnames = self.get_ambig_fnames(prev_ambig_fnames)

        model_idx = 0
        idx = 1
        for model in models_to_refine:
            # assign ambig_fname
            if ambig_fnames:
                ambig_fname = ambig_fnames[model_idx]
            else:
                ambig_fname = self.params["ambig_fname"]
            model_idx += 1

            for _ in range(self.params["sampling_factor"]):
                emref_input = prepare_cns_input(
                    idx,
                    model,
                    self.path,
                    self.recipe_str,
                    self.params,
                    "emref",
                    ambig_fname=ambig_fname,
                    native_segid=True,
                    debug=self.params["debug"],
                    seed=model.seed if isinstance(model, PDBFile) else None,
                )
                out_file = f"emref_{idx}.out"
                err_fname = f"emref_{idx}.cnserr"

                # create the expected PDBobject
                expected_pdb = prepare_expected_pdb(model, idx, ".", "emref")
                expected_pdb.restr_fname = ambig_fname
                try:
                    expected_pdb.ori_name = model.file_name
                except AttributeError:
                    expected_pdb.ori_name = None
                self.output_models.append(expected_pdb)

                job = CNSJob(emref_input, out_file, err_fname, envvars=self.envvars)

                jobs.append(job)

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

        for pdb in self.output_models:
            if pdb.is_present():
                haddock_model = HaddockModel(pdb.file_name)
                pdb.unw_energies = haddock_model.energies

                haddock_score = haddock_model.calc_haddock_score(**weights)
                pdb.score = haddock_score

        self.export_io_models(faulty_tolerance=self.params["tolerance"])
