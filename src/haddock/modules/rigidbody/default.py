"""HADDOCK3 rigid-body docking module"""
import logging
from os import linesep
from pathlib import Path
from haddock.modules import BaseHaddockModule
from haddock.cns.engine import CNSJob, CNSEngine
from haddock.cns.util import (load_axis, load_ff_parameters, load_ff_topology,
                              load_link, load_protonation_state, load_scatter,
                              load_tensor, load_trans_vectors, load_waterbox,
                              prepare_output, generate_default_header)
from haddock.cns.util import load_workflow_params, prepare_multiple_input
# from haddock.cns.topology import get_topology_header
from haddock.ontology import Format, ModuleIO, PDBFile

logger = logging.getLogger(__name__)


def generate_docking(identifier, input_files, step_path, recipe_str, defaults):
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

    output_pdb_filename = step_path / f'rigidbody_{identifier}.pdb'
    output = f"{linesep}! Output structure{linesep}"
    # output += (f"eval ($input_psf_filename="
    #            f" \"{input_psf_filename}\"){linesep}")
    output += (f"eval ($output_pdb_filename="
               f" \"{output_pdb_filename}\"){linesep}")
    inp = default_params + param + top + input_str + output \
        + topology_protonation + recipe_str

    inp_file = step_path / f'rigidbody_{identifier}.inp'
    with open(inp_file, 'w') as fh:
        fh.write(inp)

    return inp_file


class HaddockModule(BaseHaddockModule):

    def __init__(self, stream, order, path):
        self.stream = stream
        recipe_path = Path(__file__).resolve().parent.absolute()
        cns_script = recipe_path / "cns" / "rigidbody.cns"
        defaults = recipe_path / "cns" / "rigidbody.toml"
        super().__init__(order, path, cns_script, defaults)

    def run(self, module_information):
        logger.info("Running [rigidbody] module")

        # Pool of jobs to be executed by the CNS engine
        jobs = []

        # Get the models generated in previous step
        models_to_dock = [p for p in self.previous_io.output if p.file_type == Format.PDB]

        # xSampling
        structure_list = []
        for idx in range(module_information['sampling']):
            inp_file = generate_docking(idx,
                                        models_to_dock,
                                        self.path,
                                        self.recipe_str,
                                        self.defaults)

            out_file = self.path / f"rigidbody_{idx}.out"
            structure_file = self.path / f"rigidbody_{idx}.pdb"
            structure_list.append(structure_file)

            job = CNSJob(inp_file, out_file, cns_folder=self.cns_folder_path)

            jobs.append(job)

        # Run CNS engine
        logger.info(f"Running CNS engine with {len(jobs)} jobs")
        engine = CNSEngine(jobs)
        engine.run()
        logger.info("CNS engine has finished")

        # Check for generated output, fail it not all expected files are found
        expected = []
        not_found = []
        for model in structure_list:
            if not model.exists():
                not_found.append(model.name)
            pdb = PDBFile(model, path=self.path)
            expected.append(pdb)
        if not_found:
            self.finish_with_error("Several files were not generated:"
                                   f" {not_found}")

        # Save module information
        io = ModuleIO()
        io.add(structure_list)
        io.add(expected, "o")
        io.save(self.path)
