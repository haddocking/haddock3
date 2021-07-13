"""Running Lightdock as a module"""
import logging
import shutil
import subprocess
from pathlib import Path
from haddock.modules import BaseHaddockModule
from haddock.ontology import Format, ModuleIO, PDBFile
from haddock.pdbutil import PDBFactory
from haddock.modules import working_directory
from haddock.defaults import NUM_CORES


logger = logging.getLogger(__name__)


class HaddockModule(BaseHaddockModule):

    def __init__(self, order, path):
        recipe_path = Path(__file__).resolve().parent.absolute()
        defaults = recipe_path / "sampling.toml"
        super().__init__(order, path, defaults=defaults)

    def run(self, module_information):
        logger.info("Running [sampling-lightdock] module")

        # Apply module information to defaults
        self.patch_defaults(module_information)

        # Get the models generated in previous step
        models_to_score = [p for p in self.previous_io.output if p.file_type == Format.PDB]

        # Check if multiple models are provided
        if len(models_to_score) > 1:
            _msg = "Only one model allowed in LightDock sampling module"
            self.finish_with_error(_msg)

        model = models_to_score[0]
        # Check if chain IDs are present
        _path = Path(model.path, model.file_name)
        segids, chains = PDBFactory.identify_chainseg(_path)
        if set(segids) != set(chains):
            logger.info("No chain IDs found, using segid information")
            PDBFactory.swap_segid_chain(Path(model.path) / model.file_name,
                                        self.path / model.file_name)
        else:
            # Copy original model to this working path
            shutil.copyfile(Path(model.path) / model.file_name, self.path /
                            model.file_name)

        model_with_chains = self.path / model.file_name
        # Split by chain
        new_models = PDBFactory.split_by_chain(model_with_chains)
        if model_with_chains in new_models:
            self.finish_with_error(f"Input {model_with_chains} cannot be"
                                   " split by chain")

        # Receptor and ligand PDB structures
        rec_chain = self.defaults["params"]["receptor_chains"][0]
        lig_chain = self.defaults["params"]["ligand_chains"][0]
        receptor_pdb_file = (f"{Path(model.file_name).stem}_"
                             f"{rec_chain}.{Format.PDB}")
        ligand_pdb_file = (f"{Path(model.file_name).stem}_"
                           f"{lig_chain}.{Format.PDB}")

        # Setup
        logger.info("Running LightDock setup")
        with working_directory(self.path):
            swarms = self.defaults["params"]["swarms"]
            glowworms = self.defaults["params"]["glowworms"]
            noxt = self.defaults["params"]["noxt"]
            noh = self.defaults["params"]["noh"]
            cmd = (f"lightdock3_setup.py {receptor_pdb_file}"
                   f" {ligand_pdb_file} -s {swarms} -g {glowworms}")
            if noxt:
                cmd += " --noxt"
            if noh:
                cmd += " --noh"
            subprocess.call(cmd, shell=True)

        # Simulation
        logger.info("Running LightDock simulation")
        with working_directory(self.path):
            steps = self.defaults["params"]["steps"]
            scoring = self.defaults["params"]["scoring"]
            cores = NUM_CORES
            cmd = f"lightdock3.py setup.json {steps} -c {cores} -s {scoring}"
            subprocess.call(cmd, shell=True)

        # Clustering

        # Ranking
        logger.info("Generating ranking")
        with working_directory(self.path):
            steps = self.defaults["params"]["steps"]
            swarms = self.defaults["params"]["swarms"]
            cmd = f"lgd_rank.py {swarms} {steps}"
            subprocess.call(cmd, shell=True)

        # Generate top, requires a hack to use original structures (H, OXT,
        #  etc.)
        logger.info("Generating top structures")
        with working_directory(self.path):
            # Save structures, needs error control
            shutil.copyfile(self.path / receptor_pdb_file, self.path /
                            f"tmp_{receptor_pdb_file}")
            shutil.copyfile(self.path / ligand_pdb_file, self.path /
                            f"tmp_{ligand_pdb_file}")
            shutil.copy(self.path / receptor_pdb_file, self.path /
                        f"lightdock_{receptor_pdb_file}")
            shutil.copy(self.path / ligand_pdb_file, self.path /
                        f"lightdock_{ligand_pdb_file}")
            # Create top
            steps = self.defaults["params"]["steps"]
            top = self.defaults["params"]["top"]
            cmd = (f"lgd_top.py {receptor_pdb_file} {ligand_pdb_file}"
                   f" rank_by_scoring.list {top}")
            subprocess.call(cmd, shell=True)

        # Tidy top files
        expected = []
        top = self.defaults["params"]["top"]
        for i in range(top):
            file_name = f"top_{i+1}.{Format.PDB}"
            tidy_file_name = f"haddock_top_{i+1}.{Format.PDB}"
            PDBFactory.tidy(self.path / file_name, self.path / tidy_file_name)
            expected.append(PDBFile(tidy_file_name,
                                    topology=model.topology,
                                    path=(self.path / tidy_file_name)))

        # Save module information
        io = ModuleIO()
        io.add(models_to_score)
        io.add(expected, "o")
        io.save(self.path)
