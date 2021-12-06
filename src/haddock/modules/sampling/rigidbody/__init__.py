"""HADDOCK3 rigid-body docking module."""
from itertools import chain, product
from os import linesep
from pathlib import Path

from haddock import log
from haddock.gear.haddockmodel import HaddockModel
from haddock.libs.libcns import (
    generate_default_header,
    load_ambig,
    load_workflow_params,
    prepare_multiple_input,
    )
from haddock.libs.libontology import ModuleIO, PDBFile
from haddock.libs.libparallel import Scheduler
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


def generate_docking(
        identifier,
        input_list,
        step_path,
        recipe_str,
        defaults,
        ambig_fname=None,
        ):
    """Generate the .inp file that will run the docking."""
    # prepare the CNS header that will read the input

    # read the default parameters
    default_params = load_workflow_params(defaults)
    param, top, link, topology_protonation, \
        trans_vec, tensor, scatter, \
        axis, water_box = generate_default_header()

    pdb_list = []
    psf_list = []
    for element in input_list:
        pdb_fname = element.full_name
        pdb_list.append(pdb_fname)
        for psf in element.topology:
            psf_list.append(psf.full_name)

    input_str = prepare_multiple_input(pdb_list, psf_list)

    if ambig_fname:
        ambig_str = load_ambig(ambig_fname)
    else:
        ambig_str = ""

    output_pdb_filename = step_path / f'rigidbody_{identifier}.pdb'
    output = f"{linesep}! Output structure{linesep}"
    output += (f"eval ($output_pdb_filename="
               f" \"{output_pdb_filename}\"){linesep}")
    output += (f"eval ($count="
               f"{identifier}){linesep}")
    inp = default_params + param + top + input_str + output \
        + topology_protonation + ambig_str + recipe_str

    inp_file = step_path / f'rigidbody_{identifier}.inp'
    with open(inp_file, 'w') as fh:
        fh.write(inp)

    return inp_file


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for rigid body sampling."""

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        cns_script = RECIPE_PATH / "cns" / "rigidbody.cns"
        super().__init__(order, path, initial_params, cns_script)

    @classmethod
    def confirm_installation(cls):
        """Confirm module is installed."""
        return

    def run(self, **params):
        """Execute module."""
        log.info("Running [rigidbody] module")

        super().run(params)

        # Pool of jobs to be executed by the CNS engine
        jobs = []

        # Get the models generated in previous step
        #  The models from topology come in a dictionary
        input_dic = {}
        for i, model_dic in enumerate(self.previous_io.output):
            input_dic[i] = []
            for key in model_dic:
                input_dic[i].append(model_dic[key])

        if len(input_dic) == 1:
            _msg = ("Only one molecule found. To perform docking you must"
                    " input each chain as a separate molecule.")
            self.finish_with_error(_msg)

        if not self.params['crossdock']:
            # docking should be paired
            # A1-B1, A2-B2, A3-B3, etc

            # check if all ensembles contain the same number of models
            sub_lists = iter(input_dic.values())
            _len = len(next(sub_lists))
            if not all(len(sub) == _len for sub in sub_lists):
                _msg = ('With crossdock=false, the number of models inside each'
                        ' ensemble must be the same')
                self.finish_with_error(_msg)

            # prepare pairwise combinations
            models_to_dock = [values for values in zip(*input_dic.values())]

        elif self.params['crossdock']:
            # All combinations should be sampled
            # A1-B1, A1-B2, A3-B3, A2-B1, A2-B2, etc
            # prepare combinations as cartesian product
            models_to_dock = [values for values in product(*input_dic.values())]

        # How many times each combination should be sampled,
        #  cannot be smaller than 1
        sampling_factor = int(params['sampling'] / len(models_to_dock))
        if sampling_factor < 1:
            self.finish_with_error('Sampling is smaller than the number'
                                   ' of model combinations '
                                   f'#model_combinations={len(models_to_dock)},'
                                   f' sampling={params["sampling"]}.')

        # Prepare the jobs
        idx = 1
        structure_list = []
        for combination in models_to_dock:

            for _i in range(sampling_factor):
                inp_file = generate_docking(
                    idx,
                    combination,
                    self.path,
                    self.recipe_str,
                    self.params,
                    ambig_fname=self.params['ambig_fname'],
                    )

                log_fname = Path(self.path, f"rigidbody_{idx}.out")
                output_pdb_fname = Path(self.path, f"rigidbody_{idx}.pdb")

                # Create a model for the expected output
                model = PDBFile(output_pdb_fname, path=self.path)
                # FIXME: Investigate this topology structure and simplify
                topology_list = [e.topology for e in combination]
                model.topology = list(chain(*topology_list))
                structure_list.append(model)

                job = CNSJob(
                    inp_file,
                    log_fname,
                    cns_folder=self.cns_folder_path,
                    modpath=self.path,
                    config_path=self.params['config_path'],
                    cns_exec=self.params['cns_exec'],
                    )
                jobs.append(job)

                idx += 1

        # Run CNS engine
        log.info(f"Running CNS engine with {len(jobs)} jobs")
        engine = Scheduler(jobs, ncores=self.params['ncores'])
        engine.run()
        log.info("CNS engine has finished")

        # Get the weights according to CNS parameters
        _weight_keys = \
            ('w_vdw', 'w_elec', 'w_desolv', 'w_air', 'w_bsa')
        weights = {e: self.params[e] for e in _weight_keys}

        not_present = []
        for model in structure_list:
            if not model.is_present():
                not_present.append(model.full_name)

            # Score the model
            haddock_score = \
                HaddockModel(model.full_name).calc_haddock_score(**weights)

            model.score = haddock_score

        # Check for generated output
        if len(not_present) == len(structure_list):
            # fail if not all expected files are found
            self.finish_with_error("No models were generated.")

        if not_present:
            # also fail is some are not found
            # Note: we can add some fault tolerancy here,
            #  and only finish if a given % of models were not generated
            self.finish_with_error(f"Several models were not generated"
                                   f" {not_present}")

        # Save module information
        io = ModuleIO()
        io.add(structure_list, "o")
        io.save(self.path)
