"""HADDOCK3 module for flexible refinement."""
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
from haddock.libs.libontology import Format, ModuleIO, PDBFile
from haddock.libs.libparallel import Scheduler
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


def generate_flexref(
        identifier,
        input_file,
        step_path,
        recipe_str,
        defaults,
        ambig=None,
        ):
    """Generate the .inp file that will run the docking."""
    # prepare the CNS header that will read the input

    # read the default parameters
    default_params = load_workflow_params(defaults)
    param, top, link, topology_protonation, \
        trans_vec, tensor, scatter, \
        axis, water_box = generate_default_header()

    # for element in input_files:
    pdb = Path(input_file.path, input_file.file_name)
    psf_list = []
    for psf in input_file.topology:
        psf_list.append(Path(psf.path, psf.file_name))

    # pdb_list.append(str(pdb))
    # psf_list.append(str(psf))

    input_str = prepare_multiple_input([pdb], psf_list)

    if ambig:
        ambig_str = load_ambig(ambig)
    else:
        ambig_str = ""

    output_pdb_filename = step_path / f'flexref_{identifier}.pdb'
    output = f"{linesep}! Output structure{linesep}"
    output += (f"eval ($output_pdb_filename="
               f" \"{output_pdb_filename}\"){linesep}")
    output += (f"eval ($count="
               f" {identifier}){linesep}")
    inp = default_params + param + top + input_str + output \
        + topology_protonation + ambig_str + recipe_str

    inp_file = step_path / f'flexref_{identifier}.inp'
    with open(inp_file, 'w') as fh:
        fh.write(inp)

    return inp_file


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for flexible refinement."""

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        cns_script = RECIPE_PATH / "cns" / "flexref.cns"
        super().__init__(order, path, initial_params, cns_script)

    @classmethod
    def confirm_installation(cls):
        """Confirm module is installed."""
        return

    def run(self, **params):
        """Execute module."""
        log.info("Running [flexref] module")

        super().run(params)

        # Pool of jobs to be executed by the CNS engine
        jobs = []

        # Get the models generated in previous step
        models_to_refine = [
            p
            for p in self.previous_io.output
            if p.file_type == Format.PDB
            ]

        first_model = models_to_refine[0]
        topologies = first_model.topology

        refined_structure_list = []
        for idx, model in enumerate(models_to_refine, start=1):
            inp_file = generate_flexref(
                idx,
                model,
                self.path,
                self.recipe_str,
                self.params,
                ambig=self.params.get('ambig', None),
                )

            out_file = self.path / f"flexref_{idx}.out"
            structure_file = self.path / f"flexref_{idx}.pdb"
            refined_structure_list.append(structure_file)

            job = CNSJob(
                inp_file,
                out_file,
                cns_folder=self.cns_folder_path,
                modpath=self.path,
                config_path=self.params['config_path'],
                cns_exec=self.params['cns_exec'],
                )

            jobs.append(job)

        # Run CNS engine
        log.info(f"Running CNS engine with {len(jobs)} jobs")
        engine = Scheduler(jobs, ncores=self.params['ncores'])
        engine.run()
        log.info("CNS engine has finished")

        # Get the weights from the defaults
        _weight_keys = \
            ('w_vdw', 'w_elec', 'w_desolv', 'w_air', 'w_bsa')
        weights = {e: self.params[e] for e in _weight_keys}

        excepted = self.evaluate_expected(
            refined_structure_list,
            weights,
            topologies,
            )

        # Save module information
        io = ModuleIO()
        io.add(refined_structure_list)
        io.add(expected, "o")
        io.save(self.path)
