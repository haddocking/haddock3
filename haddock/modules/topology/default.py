"""CNS Topology creation and management module"""
import logging
import shutil
from pathlib import Path
from haddock.modules import BaseHaddockModule
from haddock.structure import Molecule
from haddock.pdbutil import PDBFactory
from haddock.cns.topology import generate_topology
from haddock.cns.engine import CNSJob, CNSEngine
from haddock.ontology import ModuleIO, Format, PDBFile, TopologyFile
from haddock.error import RecipeError
from haddock.defaults import TOPOLOGY_PATH


logger = logging.getLogger(__name__)


class HaddockModule(BaseHaddockModule):

    def __init__(self, order, path):
        recipe_path = Path(__file__).resolve().parent.absolute()
        cns_script = recipe_path / "cns" / "generate-topology.cns"
        defaults = recipe_path / "cns" / "generate-topology.toml"
        super().__init__(order, path, cns_script, defaults)

    def run(self, module_information):
        logger.info("Running [topology] module")

        logger.info("Generating topologies")

        # Pool of jobs to be executed by the CNS engine
        jobs = []

        try:
            molecules = HaddockModule._map(module_information["molecules"][0])
        except KeyError:
            self.finish_with_error("No molecules found in recipe")
        except RecipeError as re:
            self.finish_with_error(re)

        for order, molecule in enumerate(molecules):
            logger.info(f"{order+1} - {molecule.file_name}")

            # Get the molecule path and copy it to the course path
            molecule_orig_file_name = self.path.parent / molecule.file_name
            molecule_file_name = self.path / molecule.file_name
            shutil.copyfile(molecule_orig_file_name, molecule_file_name)

            # Split models
            logger.info(f"Split models if needed for {molecule.file_name}")
            models = sorted(PDBFactory.split_ensemble(molecule_file_name))

            # Sanitize the different PDB files
            for model in models:
                logger.info(f"Sanitizing molecule {model.name}")
                PDBFactory.sanitize(model, overwrite=True)

                # Prepare generation of topologies jobs
                topology_filename = generate_topology(model,
                                                      self.path,
                                                      self.recipe_str,
                                                      self.defaults)
                logger.info("Topology CNS input created in"
                            f" {topology_filename}")

                # Add new job to the pool
                output_filename = (model.resolve().parent.absolute()
                                   / f"{model.stem}.{Format.CNS_OUTPUT}")
                jobs.append(CNSJob(topology_filename,
                                   output_filename,
                                   cns_folder=self.cns_folder_path))

            # Run CNS engine
            logger.info(f"Running CNS engine with {len(jobs)} jobs")
            engine = CNSEngine(jobs)
            engine.run()
            logger.info("CNS engine has finished")

            # Check for generated output, fail it not all expected files
            #  are found
            expected = []
            not_found = []
            for model in models:
                model_name = model.stem
                processed_pdb = (self.path /
                                 f"{model_name}_haddock.{Format.PDB}")
                if not processed_pdb.is_file():
                    not_found.append(processed_pdb.name)
                processed_topology = (self.path /
                                      f"{model_name}_haddock"
                                      f".{Format.TOPOLOGY}")
                if not processed_topology.is_file():
                    not_found.append(processed_topology.name)
                topology = TopologyFile(processed_topology,
                                        path=(self.path / TOPOLOGY_PATH))
                expected.append(topology)
                expected.append(PDBFile(processed_pdb,
                                        topology,
                                        path=(self.path / TOPOLOGY_PATH)))
            if not_found:
                self.finish_with_error("Several files were not generated:"
                                       f" {not_found}")

            # Save module information
            io = ModuleIO()
            for model in models:
                io.add(PDBFile(model))
            io.add(expected, "o")
            io.save(self.path)

    @staticmethod
    def _map(raw_data):
        molecules = []
        for _, data in raw_data.items():
            file_name = data[0]["file"]
            try:
                segid = data[0]["segid"]
            except KeyError:
                segid = None
            molecule = Molecule(file_name, segid)
            molecules.append(molecule)
        return molecules
