"""Explicit solvent MD refinement with CNS.

The ``[mdref]`` module (also known as ``itw`` in HADDOCK2.X series), is a small
molecular dynamics simulation, in cartesian space, using explicit solvent.

A layer of solvent (8A for water, 12.5A for DMSO) is generated around
surface residues.

The `mdref` protocol is composed of 4 sequential steps:
- Short energy minimization
- 3 stages of molecular dynamics to reach 300K (at 100, 200 and 300K)
- Molecular dynamics at 300K.
- 3 stages of molecular dynamics, to reach 100K (at 300, 200 and 100K)

Using this protocol, with default parameters, no spectacular changes are
expected, however, the scoring of the various structures might be improved.


Number of MD steps can be modified using the ``waterheatsteps``, ``watersteps``
and ``watercoolsteps``.
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
    """HADDOCK3 module for water refinement."""

    name = RECIPE_PATH.name

    def __init__(
        self, order: int, path: Path, initial_params: FilePath = DEFAULT_CONFIG
    ) -> None:
        cns_script = Path(RECIPE_PATH, "cns", "mdref.cns")
        super().__init__(order, path, initial_params, cns_script=cns_script)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if module is installed."""
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

        self.output_models: list[PDBFile] = []

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
            prev_ambig_fnames = [None for mod in models_to_refine]

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
                mdref_input = prepare_cns_input(
                    idx,
                    model,
                    self.path,
                    self.recipe_str,
                    self.params,
                    "mdref",
                    ambig_fname=ambig_fname,
                    native_segid=True,
                    debug=self.params["debug"],
                    seed=model.seed if isinstance(model, PDBFile) else None,
                )
                out_file = f"mdref_{idx}.out"
                err_fname = f"mdref_{idx}.cnserr"

                # create the expected PDBobject
                expected_pdb = prepare_expected_pdb(model, idx, ".", "mdref")
                expected_pdb.restr_fname = ambig_fname
                try:
                    expected_pdb.ori_name = model.file_name
                except AttributeError:
                    expected_pdb.ori_name = None
                self.output_models.append(expected_pdb)

                job = CNSJob(mdref_input, out_file, err_fname, envvars=self.envvars)

                jobs.append(job)

                idx += 1

        # Run CNS Jobs
        self.log(f"Running CNS Jobs n={len(jobs)}")
        Engine = get_engine(self.params["mode"], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("CNS jobs have finished")

        # Get the weights from the defaults
        _weight_keys = ("w_vdw", "w_elec", "w_desolv", "w_air", "w_bsa")
        weights = {e: self.params[e] for e in _weight_keys}

        for pdb in self.output_models:
            if pdb.is_present():
                haddock_model = HaddockModel(pdb.file_name)
                pdb.unw_energies = haddock_model.energies

                haddock_score = haddock_model.calc_haddock_score(**weights)
                pdb.score = haddock_score

        # Save module information
        self.export_io_models(faulty_tolerance=self.params["tolerance"])
