"""HADDOCK3 scoring module"""
import logging
from os import linesep
from pathlib import Path

from haddock.libs.libcns import (
    generate_default_header,
    load_workflow_params,
    prepare_single_input,
    )
from haddock.libs.libparallel import Scheduler
from haddock.libs.libsubprocess import CNSJob
from haddock.modules import BaseHaddockModule
from haddock.libs.libontology import Format, ModuleIO, PDBFile


logger = logging.getLogger(__name__)


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


def generate_scoring(model, course_path, recipe_str, defaults):
    general_param = load_workflow_params(defaults)

    param, top, _, topology_protonation, _, _, _, _, _ = \
        generate_default_header()

    output_pdb_filename = course_path / Path(model.file_name)
    input_abs_path = Path(model.path).resolve()
    input_pdb_filename = input_abs_path / model.file_name
    input_psf_filename = (Path(model.topology.path) /
                          Path(model.topology.file_name))
    output = f"{linesep}! Output structure{linesep}"
    output += (f"eval ($input_psf_filename="
               f" \"{input_psf_filename}\"){linesep}")
    output += (f"eval ($output_pdb_filename="
               f" \"{output_pdb_filename}\"){linesep}")

    input_str = prepare_single_input(str(input_pdb_filename))

    inp = general_param + param + top + input_str + output \
        + topology_protonation + recipe_str

    output_inp_filename = (course_path /
                           f"{Path(model.file_name).stem}.{Format.CNS_INPUT}")
    with open(output_inp_filename, "w") as output_handler:
        output_handler.write(inp)

    return output_inp_filename


class HaddockModule(BaseHaddockModule):

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        cns_script = RECIPE_PATH / "cns" / "scoring.cns"
        super().__init__(order, path, initial_params, cns_script)

    @classmethod
    def confirm_installation(cls):
        return

    def run(self, **params):
        logger.info("Running [scoring] module")

        super().run(params)

        # Pool of jobs to be executed by the CNS engine
        jobs = []

        # Get the models generated in previous step
        models_to_score = [p for p in self.previous_io.output if p.file_type == Format.PDB]

        for model in models_to_score:
            scoring_filename = generate_scoring(model,
                                                self.path,
                                                self.recipe_str,
                                                self.params,
                                                )
            output_filename = (self.path /
                               f"{Path(model.file_name).stem}_scoring.out")

            job = CNSJob(
                scoring_filename,
                output_filename,
                cns_folder=self.cns_folder_path,
                cns_exec=self.params['cns_exec'],
                )

            jobs.append(job)

        # Run CNS engine
        logger.info(f"Running CNS engine with {len(jobs)} jobs")
        engine = Scheduler(jobs, ncores=self.params['ncores'])
        engine.run()
        logger.info("CNS engine has finished")

        # Check for generated output, fail it not all expected files are found
        expected = []
        not_found = []
        for model in models_to_score:
            model_path = Path(self.path, model.file_name)
            if model_path.is_file():
                expected.append(
                    PDBFile(
                        model.file_name,
                        topology=model.topology,
                        path=self.path
                        ))
            else:
                not_found.append(model.file_name)

        if not_found:
            self.finish_with_error("Several files were not generated:"
                                   f" {not_found}")

        # Save module information
        io = ModuleIO()
        io.add(models_to_score)
        io.add(expected, "o")
        io.save(self.path)
