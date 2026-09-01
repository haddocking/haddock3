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

For more details about this module, please `refer to the haddock3 user manual
<https://www.bonvinlab.org/haddock3-user-manual/modules/sampling.html#rigidbody-module>`_
"""

from pathlib import Path

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath, Optional, Sequence, Union
from haddock.gear.haddockmodel import HaddockModel
from haddock.libs.libcns import derive_seed, prepare_cns_input, prepare_expected_pdb
from haddock.libs.libontology import PDBFile
from haddock.libs.libparallel import GenericTask
from haddock.libs.libpdb import check_combination_chains
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import get_engine
from haddock.modules.base_cns_module import BaseCNSModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


def _chainids_for_sampled_combinations(
    combinations: Sequence[list[PDBFile]],
    sampled_combinations: Sequence[list[PDBFile]],
) -> list[list[str]]:
    """Resolve chain IDs once for each distinct source combination."""
    chainids_by_combination = {
        id(combination): check_combination_chains(combination)
        for combination in combinations
    }
    return [
        chainids_by_combination[id(combination)] for combination in sampled_combinations
    ]


def _repeats_of_sampled_combinations(
    sampled_combinations: Sequence[list[PDBFile]],
) -> list[int]:
    """How many times each sampled combination has already been scheduled.

    A job's seed is derived from the combination it docks and from which
    repeat of that combination it is, never from where it sits in the
    schedule, so this is counted per combination rather than taken from the
    job index.  Counting it from the schedule itself rather than computing
    ``index // n_combinations`` keeps the two in step by construction if the
    round-robin is ever replaced by another prefix-stable order.
    """
    seen: dict[int, int] = {}
    repeats = []
    for combination in sampled_combinations:
        repeat = seen.get(id(combination), 0)
        repeats.append(repeat)
        seen[id(combination)] = repeat + 1
    return repeats


class HaddockModule(BaseCNSModule):
    """HADDOCK3 module for rigid body sampling."""

    name = RECIPE_PATH.name

    def __init__(
        self, order: int, path: Path, initial_params: FilePath = DEFAULT_CONFIG
    ) -> None:
        cns_script = Path(RECIPE_PATH, "cns", f"{self.name}.cns")
        super().__init__(order, path, initial_params, cns_script=cns_script)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm module is installed."""
        return

    def make_cns_jobs(
        self,
        inp_list: Sequence[
            tuple[list[PDBFile], Union[Path, str], Union[str, None], int]
        ],
    ) -> list[CNSJob]:
        jobs: list[CNSJob] = []
        for idx, job_entries in enumerate(inp_list, start=1):
            combination, cns_input, ambig_fname, seed = job_entries

            log_fname = f"{self.name}_{idx}.out"
            err_fname = f"{self.name}_{idx}.cnserr"

            # Create a model for the expected output
            model = prepare_expected_pdb(combination, idx, ".", self.name)
            # Set additional attributes
            model.restr_fname=str(ambig_fname)
            model.seed = seed

            self.output_models.append(model)

            job = CNSJob(cns_input, log_fname, err_fname, envvars=self.envvars)
            jobs.append(job)
        return jobs

    @staticmethod
    def _sample_models_to_dock(
        models_to_dock: list[list[PDBFile]],
        sampling: int,
    ) -> list[list[PDBFile]]:
        """Return a prefix-stable model-combination schedule."""
        return [
            models_to_dock[model_idx % len(models_to_dock)]
            for model_idx in range(sampling)
        ]

    def prepare_cns_input_sequential(
        self,
        models_to_dock: list[list[PDBFile]],
        ambig_fnames: Union[list, None],
        chainid_lists: Optional[list[list[str]]] = None,
    ) -> list[tuple[list[PDBFile], Union[Path, str], Union[str, None], int]]:
        _l = []
        repeats = _repeats_of_sampled_combinations(models_to_dock)
        idx = 1
        for combination_index, combination in enumerate(models_to_dock):
            # assign ambig_fname
            if ambig_fnames:
                ambig_fname = ambig_fnames[idx - 1]
            else:
                ambig_fname = self.params["ambig_fname"]
            # prepare cns input
            seed = derive_seed(
                self.params["iniseed"],
                combination,
                repeats[combination_index],
            )
            rigidbody_input = prepare_cns_input(
                idx,
                combination,
                self.path,
                self.recipe_str,
                self.params,
                self.name,
                ambig_fname=ambig_fname,
                default_params_path=self.toppar_path,
                native_segid=True,
                debug=self.params["debug"],
                seed=seed,
                chainid_list=(
                    chainid_lists[combination_index] if chainid_lists else None
                ),
            )
            _l.append((combination, rigidbody_input, ambig_fname, seed))
            idx += 1
        return _l

    def prepare_cns_input_parallel(
        self,
        models_to_dock: list[list[PDBFile]],
        ambig_fnames: Union[list, None],
        chainid_lists: Optional[list[list[str]]] = None,
    ) -> list[tuple[list[PDBFile], Union[Path, str], Union[str, None], int]]:
        prepare_tasks = []
        _l = []
        repeats = _repeats_of_sampled_combinations(models_to_dock)
        idx: int = 1
        for combination_index, combination in enumerate(models_to_dock):
            ambig_fname = (
                ambig_fnames[idx - 1] if ambig_fnames else self.params["ambig_fname"]
            )
            seed = derive_seed(
                self.params["iniseed"],
                combination,
                repeats[combination_index],
            )
            task = GenericTask(
                function=prepare_cns_input,
                model_number=idx,
                input_element=combination,
                step_path=self.path,
                recipe_str=self.recipe_str,
                defaults=self.params,
                identifier=self.name,
                ambig_fname=ambig_fname,
                native_segid=True,
                default_params_path=self.toppar_path,
                debug=self.params["debug"],
                seed=seed,
                chainid_list=(
                    chainid_lists[combination_index] if chainid_lists else None
                ),
            )

            prepare_tasks.append(task)
            _l.append((combination, task, ambig_fname, seed))
            idx += 1
        Engine = get_engine(self.params["mode"], self.params)
        prepare_engine = Engine(prepare_tasks)
        prepare_engine.run()

        # Replace the task with the result of the task
        prepared_inputs = []
        for element, task_result in zip(_l, prepare_engine.results):
            prepared_inputs.append((element[0], task_result, element[2], element[3]))

        return prepared_inputs

    def restraints_guardrail(self, ambig_fnames: Optional[list[str]]) -> None:
        """Makes sure any restraints are available for the docking."""
        # List all types of restraints
        all_restraints = (
            ambig_fnames,
            self.params["ambig_fname"],
            self.params["unambig_fname"],
            self.params["hbond_fname"],
            self.params["cmrest"],
            self.params["ranair"],
            self.params["surfrest"],
        )
        # If not any restraints provided
        if not any(all_restraints):
            # Terminate docking
            self.finish_with_error(
                f"No restraints found for [{self.name}] module. "
                "For targeted docking, supply CNS-valid restraints file(s) "
                "using 'ambig_fname' and/or 'unambig_fname' and/or "
                "'hbond_fname' parameter(s). "
                "For ab-initio docking, set 'cmrest' or 'ranair' "
                "parameters to true."
            )

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

        # Each model combination should be sampled at least once. The global
        # sampling sequence itself is prefix-stable: lowering sampling by one
        # keeps all preceding CNS job inputs and seeds unchanged.
        if self.params["sampling"] < len(models_to_dock):
            self.finish_with_error(
                "Sampling is smaller than the number"
                " of model combinations "
                f"#model_combinations={len(models_to_dock)},"
                f" sampling={self.params['sampling']}."
            )
        sampled_models_to_dock = self._sample_models_to_dock(
            models_to_dock,
            self.params["sampling"],
        )
        sampled_chainid_lists = _chainids_for_sampled_combinations(
            models_to_dock,
            sampled_models_to_dock,
        )

        # get all the different ambig files
        prev_ambig_fnames = [None for _model in range(self.params["sampling"])]
        diff_ambig_fnames = self.get_ambig_fnames(prev_ambig_fnames)
        # if no files are found, we will stick to self.params["ambig_fname"]
        if diff_ambig_fnames:
            n_diffs = len(diff_ambig_fnames)
            ambig_fnames = [
                diff_ambig_fnames[n % n_diffs] for n in range(self.params["sampling"])
            ]
        else:
            ambig_fnames = None

        self.restraints_guardrail(ambig_fnames)

        self.log("Preparing CNS jobs...")
        if self.params["mode"] != "local":
            # Note: `batch` and (pseudo)-`mpi` mode uses files to communicate and cannot extract the information from the task object.
            cns_input = self.prepare_cns_input_sequential(
                sampled_models_to_dock,
                ambig_fnames,  # type: ignore
                sampled_chainid_lists,
            )
        else:
            cns_input = self.prepare_cns_input_parallel(
                sampled_models_to_dock,
                ambig_fnames,  # type: ignore
                sampled_chainid_lists,
            )

        self.output_models: list[PDBFile] = []
        jobs = self.make_cns_jobs(cns_input)

        # Run CNS Jobs
        self.log(f"Running CNS Jobs n={len(jobs)}")
        Engine = get_engine(self.params["mode"], self.params)
        engine = Engine(jobs)
        engine.run()
        self.log("CNS jobs have finished")

        # Get the weights according to CNS parameters
        _weight_keys = ("w_vdw", "w_elec", "w_desolv", "w_air", "w_bsa")
        weights = {e: self.params[e] for e in _weight_keys}
        # Loop over models
        for model in self.output_models:
            if model.is_present():
                # Obtain the model scores
                haddock_model = HaddockModel(model.file_name)
                haddock_score = haddock_model.calc_haddock_score(**weights)
                # Set the attributes
                model.unw_energies = haddock_model.energies
                model.score = haddock_score

        self.export_io_models(faulty_tolerance=self.params["tolerance"])
