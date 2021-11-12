"""HADDOCK3 rigid-body docking module."""
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


def generate_docking(
        identifier,
        input_files,
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

    # input_files is the ontology, unwrap it
    pdb_list = []
    psf_list = []
    for element in input_files:
        pdb = Path(element.path, element.file_name)
        psf = Path(element.path, element.topology.file_name)

        pdb_list.append(str(pdb))
        psf_list.append(str(psf))

    input_str = prepare_multiple_input(pdb_list, psf_list)

    if ambig:
        ambig_str = load_ambig(ambig)
    else:
        ambig_str = ""

    output_pdb_filename = step_path / f'rigidbody_{identifier}.pdb'
    output = f"{linesep}! Output structure{linesep}"
    output += (f"eval ($output_pdb_filename="
               f" \"{output_pdb_filename}\"){linesep}")
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
        models_to_dock = [
            p
            for p in self.previous_io.output
            if p.file_type == Format.PDB
            ]

        # TODO: Make the topology aquisition generic,
        # here its expecting this module
        # to be preceeded by topology
        topologies = [
            p
            for p in self.previous_io.output
            if p.file_type == Format.TOPOLOGY
            ]

        # Sampling
        structure_list = []
        for idx in range(params['sampling']):
            inp_file = generate_docking(
                idx,
                models_to_dock,
                self.path,
                self.recipe_str,
                self.params,
                ambig=self.params.get('ambig', None),
                )

            out_file = self.path / f"rigidbody_{idx}.out"
            structure_file = self.path / f"rigidbody_{idx}.pdb"
            structure_list.append(structure_file)

            job = CNSJob(
                inp_file,
                out_file,
                cns_folder=self.cns_folder_path,
                cns_exec=self.params['cns_exec'],
                )

            jobs.append(job)

        # Run CNS engine
        log.info(f"Running CNS engine with {len(jobs)} jobs")
        engine = Scheduler(jobs, ncores=self.params['ncores'])
        engine.run()
        log.info("CNS engine has finished")

        # Get the weights according to CNS parameters
        _weight_keys = \
            ('w_vdw_0', 'w_elec_0', 'w_desolv_0', 'w_air_0', 'w_bsa_0')
        weights = {e: self.params[e] for e in _weight_keys}

        expected = []
        not_found = []
        for model in structure_list:
            if not model.exists():
                not_found.append(model.name)

            haddock_score = HaddockModel(model).calc_haddock_score(**weights)

            pdb = PDBFile(model, path=self.path)
            pdb.score = haddock_score
            pdb.topology = topologies
            expected.append(pdb)

        if not_found:
            # Check for generated output,
            # fail if not all expected files are found
            self.finish_with_error("Several files were not generated:"
                                   f" {not_found}")

        # Save module information
        io = ModuleIO()
        io.add(structure_list)
        io.add(expected, "o")
        io.save(self.path)
