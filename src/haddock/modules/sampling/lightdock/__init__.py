"""Running Lightdock as a module"""
import logging
import shutil
import subprocess

from pathlib import Path

from haddock.libs import libpdb
from haddock.libs.libontology import Format, ModuleIO, PDBFile
from haddock.libs.libutil import check_subprocess
from haddock.modules import BaseHaddockModule, working_directory


logger = logging.getLogger(__name__)

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


class HaddockModule(BaseHaddockModule):

    def __init__(
            self,
            order,
            path,
            *ignore,
            initial_params=DEFAULT_CONFIG,
            **everything,
            ):
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm this module is installed."""
        check_subprocess('lightdock3.py -h')

    def run(self, **params):
        logger.info("Running [sampling-lightdock] module")

        super().run(params)

        # Get the models generated in previous step
        models_to_score = [p for p in self.previous_io.output if p.file_type == Format.PDB]

        # Check if multiple models are provided
        if len(models_to_score) > 1:
            _msg = "Only one model allowed in LightDock sampling module"
            self.finish_with_error(_msg)

        model = models_to_score[0]
        # Check if chain IDs are present
        _path = Path(model.path, model.file_name)
        segids, chains = libpdb.identify_chainseg(_path)
        if set(segids) != set(chains):
            logger.info("No chain IDs found, using segid information")
            libpdb.swap_segid_chain(Path(model.path) / model.file_name,
                                        self.path / model.file_name)
        else:
            # Copy original model to this working path
            shutil.copyfile(Path(model.path) / model.file_name, self.path /
                            model.file_name)

        model_with_chains = self.path / model.file_name
        # Split by chain
        new_models = libpdb.split_by_chain(model_with_chains)
        if model_with_chains in new_models:
            self.finish_with_error(f"Input {model_with_chains} cannot be"
                                   " split by chain")

        # Receptor and ligand PDB structures
        rec_chain = self.params["receptor_chains"][0]
        lig_chain = self.params["ligand_chains"][0]
        receptor_pdb_file = (f"{Path(model.file_name).stem}_"
                             f"{rec_chain}.{Format.PDB}")
        ligand_pdb_file = (f"{Path(model.file_name).stem}_"
                           f"{lig_chain}.{Format.PDB}")

        # Setup
        logger.info("Running LightDock setup")
        with working_directory(self.path):
            swarms = self.params["swarms"]
            glowworms = self.params["glowworms"]
            noxt = self.params["noxt"]
            noh = self.params["noh"]
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
            steps = self.params["steps"]
            scoring = self.params["scoring"]
            cores = self.params['ncores'] or 1
            cmd = f"lightdock3.py setup.json {steps} -c {cores} -s {scoring}"
            subprocess.call(cmd, shell=True)

        # Clustering

        # Ranking
        logger.info("Generating ranking")
        with working_directory(self.path):
            steps = self.params["steps"]
            swarms = self.params["swarms"]
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
            steps = self.params["steps"]
            top = self.params["top"]
            cmd = (f"lgd_top.py {receptor_pdb_file} {ligand_pdb_file}"
                   f" rank_by_scoring.list {top}")
            subprocess.call(cmd, shell=True)

        # Tidy top files
        expected = []
        top = self.params["top"]
        for i in range(top):
            file_name = f"top_{i+1}.{Format.PDB}"
            tidy_file_name = f"haddock_top_{i+1}.{Format.PDB}"
            libpdb.tidy(self.path / file_name, self.path / tidy_file_name)
            expected.append(PDBFile(tidy_file_name,
                                    topology=model.topology,
                                    path=self.path))

        # Save module information
        io = ModuleIO()
        io.add(models_to_score)
        io.add(expected, "o")
        io.save(self.path)
