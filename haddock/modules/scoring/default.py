import logging
from os import linesep
from pathlib import Path
from haddock.modules import BaseHaddockModule
from haddock.cns.engine import CNSJob, CNSEngine
from haddock.cns.util import load_recipe_params, prepare_input
from haddock.cns.topology import get_topology_header
from haddock.ontology import Format, ModuleIO


logger = logging.getLogger(__name__)


def generate_scoring(input_pdb, course_path, recipe_str, defaults):
    general_param = load_recipe_params(defaults)

    param, top, link, topology_protonation, trans_vec, tensor, scatter, axis, water_box = get_topology_header()

    output_pdb_filename = course_path / input_pdb.name
    input_abs_path = input_pdb.resolve().parent.absolute()
    output_psf_filename = course_path / f'{input_pdb.stem}.{Format.TOPOLOGY}'
    output = f'{linesep}! Output structure{linesep}'
    output += f"eval ($output_psf_filename= \"{output_psf_filename}\"){linesep}"
    output += f"eval ($output_pdb_filename= \"{output_pdb_filename}\"){linesep}"

    input_str = prepare_input(str(input_pdb.resolve().absolute()), input_abs_path)

    inp = general_param + param + top + input_str + output + topology_protonation + recipe_str

    output_inp_filename = course_path / f'{input_pdb.stem}.{Format.CNS_INPUT}'
    with open(output_inp_filename, 'w') as output_handler:
        output_handler.write(inp)

    return output_inp_filename


class HaddockModule(BaseHaddockModule):

    def __init__(self, order, path):
        recipe_path = Path(__file__).resolve().parent.absolute()
        cns_script = recipe_path / "cns" / "scoring.cns"
        cns_defaults = recipe_path / "cns" / "scoring.toml"
        super().__init__(order, path, cns_script, cns_defaults)

    def run(self, module_information):
        logger.info("Running [scoring] module")

        # Pool of jobs to be executed by the CNS engine
        jobs = []

        # Get the models generated in previous step
        models_to_score = [o[0] for o in self.previous_io.output if o[1] == Format.PDB]
        for input_pdb_filename in models_to_score:
            input_pdb = self.previous_path() / input_pdb_filename
            scoring_filename = generate_scoring(input_pdb, self.path, self.recipe_str, self.defaults)
            output_filename = self.path / f'{input_pdb.stem}_scoring.out'
            jobs.append(CNSJob(scoring_filename, output_filename, cns_folder=self.cns_folder_path))

        # Run CNS engine
        logger.info(f'Running CNS engine with {len(jobs)} jobs')
        engine = CNSEngine(jobs)
        engine.run()
        logger.info('CNS engine has finished')

        # Check for generated output, fail it not all expected files are found
        expected = []
        not_found = []
        for model in models_to_score:
            model = Path(model)
            if not (self.path / model).is_file():
                not_found.append(model)
            expected.append((model, Format.PDB))
        if not_found:
            self.finish_with_error(f'Several files were not generated: {not_found}')

        # Save module information
        io = ModuleIO()
        for model in models_to_score:
            io.add(model, Format.PDB)
        for expected_file, output_format in expected:
            io.add(expected_file.name, output_format, 'o')
        io.save(self.path)
