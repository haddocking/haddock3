"""
Rigid-body docking sampling module
==================================

This module does a **randomization of orientations and rigid-body
minimization.** It corresponds to the classical ``it0`` step in the HADDOCK2.x
series.

In this module, the interacting partners are treated as rigid bodies, meaning
that all geometrical parameters such as bond lengths, bond angles, and dihedral
angles are frozen. The partners are first separated in space and randomly
rotated around their respective centres of mass. Afterwards, the molecules are
brought together by rigid-body energy minimisation with rotations and
translation as the only degrees of freedom.

The driving force for this energy minimisation is the energy function, which
consists of the intermolecular van der Waals and electrostatic energy terms and
the restraints defined to guide the docking. The restraints are distance-based
and can consist of unambiguous or ambiguous interactions restraints (AIRS). In
*ab-initio* docking mode those restraints can be automatically defined in
various ways; e.g. between center of masses (CM restraints) or between randomly
selected patches on the surface (random AIRs).

The definition of those restraints is particularly important as they effectively
guide the minimisation process. For example, with a stringent set of AIRs or
unambiguous distance restraints, the solutions of the minimisation will converge
much better and the sampling can be limited. In *ab-initio* mode, however, very
diverse solutions will be obtained and the sampling should be increased to make
sure to sample enough the possible interaction space.
"""

from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath
from haddock.gear.haddockmodel import HaddockModel
from haddock.libs.libcns import prepare_cns_input
from haddock.libs.libontology import PDBFile
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import get_engine
from haddock.modules.base_cns_module import BaseCNSModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseCNSModule):
    """HADDOCK3 module for rigid body sampling."""

    name = RECIPE_PATH.name

    def __init__(
        self, order: int, path: Path, initial_params: FilePath = DEFAULT_CONFIG
    ) -> None:
        cns_script = Path(RECIPE_PATH, "cns", "rigidbody.cns")
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
            if self.params["crossdock"]:
                self.log("crossdock=true")
            models_to_dock = self.previous_io.retrieve_models(
                crossdock=self.params["crossdock"]
            )
        except Exception as e:
            self.finish_with_error(e)

        # How many times each combination should be sampled,
        #  cannot be smaller than 1
        sampling_factor = int(self.params["sampling"] / len(models_to_dock))
        if sampling_factor < 1:
            self.finish_with_error(
                "Sampling is smaller than the number"
                " of model combinations "
                f"#model_combinations={len(models_to_dock)},"
                f' sampling={self.params["sampling"]}.'
            )

        # get all the different ambig files
        prev_ambig_fnames = [None for model in range(self.params["sampling"])]
        diff_ambig_fnames = self.get_ambig_fnames(prev_ambig_fnames)  # type: ignore
        # if no files are found, we will stick to self.params["ambig_fname"]
        if diff_ambig_fnames:
            n_diffs = len(diff_ambig_fnames)
            ambig_fnames = [
                diff_ambig_fnames[n % n_diffs] for n in range(self.params["sampling"])
            ]  # noqa: E501
        else:
            ambig_fnames = None

        # Prepare the jobs
        idx = 1
        self.output_models: list[PDBFile] = []
        self.log("Preparing jobs...")
        for combination in models_to_dock:
            for _i in range(sampling_factor):
                # assign ambig_fname
                if ambig_fnames:
                    ambig_fname = ambig_fnames[idx - 1]
                else:
                    ambig_fname = self.params["ambig_fname"]

                seed = self.params["iniseed"] * idx
                # prepare cns input
                rigidbody_input = prepare_cns_input(
                    idx,
                    combination,
                    self.path,
                    self.recipe_str,
                    self.params,
                    "rigidbody",
                    ambig_fname=ambig_fname,
                    default_params_path=self.toppar_path,
                    native_segid=True,
                    less_io=self.params["less_io"],
                    seed=seed,
                )

                log_fname = f"rigidbody_{idx}.out"
                output_pdb_fname = f"rigidbody_{idx}.pdb"

                # Create a model for the expected output
                model = PDBFile(output_pdb_fname, path=".", restr_fname=ambig_fname)
                model.topology = [e.topology for e in combination]  # type: ignore
                model.seed = seed

                self.output_models.append(model)

                job = CNSJob(rigidbody_input, log_fname, envvars=self.envvars)
                jobs.append(job)

                idx += 1

        # Run CNS Jobs
        self.log(f"Running CNS Jobs n={len(jobs)}")
        Engine = get_engine(self.params["mode"], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("CNS jobs have finished")

        # Get the weights according to CNS parameters
        _weight_keys = ("w_vdw", "w_elec", "w_desolv", "w_air", "w_bsa")
        weights = {e: self.params[e] for e in _weight_keys}

        for model in self.output_models:
            if model.is_present():
                # Score the model
                haddock_model = HaddockModel(model.file_name)
                model.unw_energies = haddock_model.energies

                haddock_score = haddock_model.calc_haddock_score(**weights)
                model.score = haddock_score

        self.export_io_models(faulty_tolerance=self.params["tolerance"])
