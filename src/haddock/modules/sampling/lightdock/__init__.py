"""LightDock integration sampling module."""
import shutil
import subprocess
from pathlib import Path

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, FilePath
from haddock.libs import libpdb
from haddock.libs.libio import working_directory
from haddock.libs.libontology import Format, PDBFile
from haddock.libs.libutil import check_subprocess
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 Lightdock module."""

    name = RECIPE_PATH.name

    def __init__(
            self,
            order: int,
            path: Path,
            *ignore: Any,
            initial_params: FilePath = DEFAULT_CONFIG,
            **everything: Any,
            ) -> None:
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm this module is installed."""
        check_subprocess('lightdock3.py -h')

    def _run(self) -> None:
        """Execute module."""
        # Get the models generated in previous step
        models: list[PDBFile] = [
            p
            for p in self.previous_io.output
            if p.file_type == Format.PDB
            ]

        # Check if multiple models are provided
        if len(models) > 1:
            _msg = "Only one model allowed in LightDock sampling module"
            self.finish_with_error(_msg)

        model = models[0]
        # Check if chain IDs are present
        _path = Path(model.path, model.file_name)
        segids, chains = libpdb.identify_chainseg(_path)
        if set(segids) != set(chains):
            log.info("No chain IDs found, using segid information")
            libpdb.swap_segid_chain(
                Path(model.path, model.file_name),
                Path(self.path, model.file_name),
                )
        else:
            # Copy original model to this working path
            shutil.copyfile(
                Path(model.path, model.file_name),
                Path(self.path, model.file_name),
                )

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
        log.info("Running LightDock setup")
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
        log.info("Running LightDock simulation")
        with working_directory(self.path):
            steps = self.params["steps"]
            scoring = self.params["scoring"]
            cores = self.params['ncores'] or 1
            cmd = f"lightdock3.py setup.json {steps} -c {cores} -s {scoring}"
            subprocess.call(cmd, shell=True)

        # Clustering

        # Ranking
        log.info("Generating ranking")
        with working_directory(self.path):
            steps = self.params["steps"]
            swarms = self.params["swarms"]
            cmd = f"lgd_rank.py {swarms} {steps}"
            subprocess.call(cmd, shell=True)

        # Generate top, requires a hack to use original structures (H, OXT,
        #  etc.)
        log.info("Generating top structures")
        with working_directory(self.path):
            # Save structures, needs error control
            shutil.copyfile(
                Path(self.path, receptor_pdb_file),
                Path(self.path, f"tmp_{receptor_pdb_file}"),
                )
            shutil.copyfile(
                Path(self.path, ligand_pdb_file),
                Path(self.path, f"tmp_{ligand_pdb_file}")
                )
            shutil.copy(
                Path(self.path, receptor_pdb_file),
                Path(self.path, f"lightdock_{receptor_pdb_file}"),
                )
            shutil.copy(
                Path(self.path, ligand_pdb_file),
                Path(self.path, f"lightdock_{ligand_pdb_file}"),
                )
            # Create top
            steps = self.params["steps"]
            top = self.params["top"]
            cmd = (f"lgd_top.py {receptor_pdb_file} {ligand_pdb_file}"
                   f" rank_by_scoring.list {top}")
            subprocess.call(cmd, shell=True)

        # Tidy top files
        expected: list[PDBFile] = []
        top = self.params["top"]
        for i in range(top):
            file_name = f"top_{i+1}.{Format.PDB}"
            tidy_file_name = f"haddock_top_{i+1}.{Format.PDB}"
            libpdb.tidy(self.path / file_name, self.path / tidy_file_name)
            expected.append(PDBFile(tidy_file_name,
                                    topology=model.topology,
                                    path=self.path))

        self.output_models = models
        self.export_io_models()
