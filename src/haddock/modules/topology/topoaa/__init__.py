"""CNS Topology creation and management module"""
import logging
import shutil
from pathlib import Path
from haddock.modules import BaseHaddockModule
from haddock.structure import Molecule, make_molecules
from haddock.pdbutil import PDBFactory
from haddock.libs.libsubprocess import CNSJob
from haddock.cns.util import (
    load_workflow_params,
    generate_default_header,
    prepare_output,
    prepare_single_input,
    )
from haddock.defaults import TOPOLOGY_PATH
from haddock.error import StepError
from haddock.libs.libparallel import Scheduler
from haddock.ontology import ModuleIO, Format, PDBFile, TopologyFile


logger = logging.getLogger(__name__)


def generate_topology(input_pdb, step_path, recipe_str, defaults,
                      protonation=None):
    """Generate a HADDOCK topology file from input_pdb"""
    general_param = load_workflow_params(defaults)

    param, top, link, topology_protonation, \
        trans_vec, tensor, scatter, \
        axis, water_box = generate_default_header(protonation)

    abs_path = input_pdb.resolve().parent.absolute()
    output_pdb_filename = abs_path / (f'{input_pdb.stem}_'
                                      f'haddock{input_pdb.suffix}')
    output_psf_filename = abs_path / (f'{input_pdb.stem}_'
                                      f'haddock.{Format.TOPOLOGY}')
    output = prepare_output(output_psf_filename, output_pdb_filename)

    input_str = prepare_single_input(str(input_pdb.resolve().absolute()))

    inp = general_param + param + top + input_str + output + link \
        + topology_protonation + trans_vec + tensor + scatter + axis \
        + water_box + recipe_str

    output_inp_filename = abs_path / f'{input_pdb.stem}.{Format.CNS_INPUT}'
    with open(output_inp_filename, 'w') as output_handler:
        output_handler.write(inp)

    return output_inp_filename


class HaddockModule(BaseHaddockModule):

    def __init__(self, order, path):
        recipe_path = Path(__file__).resolve().parent
        cns_script = recipe_path / "cns" / "generate-topology.cns"
        defaults = recipe_path / "cns" / "generate-topology.toml"
        super().__init__(order, path, cns_script, defaults)

    def run(self, **params):
        logger.info("Running [allatom] module")
        logger.info("Generating topologies")

        molecules = make_molecules(params['molecules'])

        # Pool of jobs to be executed by the CNS engine
        jobs = []

        models = []
        for i, molecule in enumerate(molecules):
            logger.info(f"{i + 1} - {molecule.file_name}")

            # Copy the molecule to the step folder
            step_molecule_path = Path(self.path, molecule.file_name.name)
            shutil.copyfile(molecule.file_name, step_molecule_path)

            # Split models
            logger.info(f"Split models if needed for {step_molecule_path}")
            ens = PDBFactory.split_ensemble(step_molecule_path)
            splited_models = sorted(ens)

            # Sanitize the different PDB files
            for model in splited_models:
                logger.info(f"Sanitizing molecule {model.name}")
                models.append(model)
                PDBFactory.sanitize(model, overwrite=True)

                # Prepare generation of topologies jobs
                topology_filename = generate_topology(model,
                                                      self.path,
                                                      self.recipe_str,
                                                      self.defaults)
                logger.info("Topology CNS input created in"
                            f" {topology_filename}")

                # Add new job to the pool
                output_filename = Path(
                    model.resolve().parent,
                    f"{model.stem}.{Format.CNS_OUTPUT}",
                    )

                job = CNSJob(topology_filename,
                             output_filename,
                             cns_folder=self.cns_folder_path)

                jobs.append(job)

        # Run CNS engine
        logger.info(f"Running CNS engine with {len(jobs)} jobs")
        engine = Scheduler(jobs)
        engine.run()
        logger.info("CNS engine has finished")

        # Check for generated output, fail it not all expected files
        #  are found
        expected = []
        not_found = []
        for model in models:
            model_name = model.stem
            processed_pdb = (self.path / f"{model_name}_haddock.{Format.PDB}")
            if not processed_pdb.is_file():
                not_found.append(processed_pdb.name)
            processed_topology = (self.path /
                                  f"{model_name}_haddock"
                                  f".{Format.TOPOLOGY}")
            if not processed_topology.is_file():
                not_found.append(processed_topology.name)
            topology = TopologyFile(processed_topology,
                                    path=self.path)
            expected.append(topology)
            expected.append(PDBFile(processed_pdb,
                                    topology,
                                    path=self.path))
        if not_found:
            self.finish_with_error("Several files were not generated:"
                                   f" {not_found}")

        # Save module information
        io = ModuleIO()
        for model in models:
            io.add(PDBFile(model))
        io.add(expected, "o")
        io.save(self.path)
