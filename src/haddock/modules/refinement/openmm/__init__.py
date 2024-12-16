"""OpenMM Molecular Dynamics refinement module for HADDOCK3.

The potential of OpenMM can be exploited to perform potentially different
tasks, such as:

* Run MD simulation for each model from previous step;
* Refine the models in the middle of a docking run. For example, it can be used
  to refine the models coming from a `[rigidbody]` module before `[flexref]` is
  executed, or to replace the `[mdref]` step.
* Generate conformers prior to their use in a thorough docking run.

To get a list of all possible parameters, run:

>>> haddock3-cfg -m openmm

Module workflow:

* Generate openmm topology and fix atoms
* Build solvation box
* Equilibration solvation box restraining the protein
* Run MD simulation: increase temperature, run MD, reduce temperature.
* Either generate an ensemble of multiple frames or return the last frame.

This module will refine all models coming from the previous workflow step and
send them to the next step in the workflow. If you want to use other modules
such as `flexref` or `emref` after the OpenMM module, you need to recreate the
topologies by simply adding a `[topoaa]` step in the workflow.
See examples in `examples/thirdparty/openmm` folder.
"""
import os
import shutil

from contextlib import suppress
from pathlib import Path
from subprocess import run as subprocrun

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.exceptions import ThirdPartyIntallationError
from haddock.libs.libontology import PDBFile
from haddock.modules import BaseHaddockModule, get_engine

# allow general testing when OpenMM is not installed
with suppress(ImportError):
    from haddock.modules.refinement.openmm.openmm import OPENMM


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 OpenMM module."""

    name = RECIPE_PATH.name

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        super().__init__(order, path, initial_params)

    def create_directories(self) -> dict[str, str]:
        """Create the necessary directories and provides the paths."""
        directory_list = (
            "pdbfixer",
            "solvation_boxes",
            "intermediates",
            "md_raw_output",
            "openmm_output",
            "simulation_stats",
            )

        directory_dict: dict[str, str] = {}
        for dir in directory_list:
            self.log(f"Creating directory {dir}")
            os.mkdir(dir)
            directory_dict[dir] = dir
        
        return directory_dict
    
    def remove_directories(self) -> None:
        """Remove unnecessary directories full of heavy files."""
        directory_list = (
            "pdbfixer",
            "md_raw_output",
            "solvation_boxes",
            )
        for dir in directory_list:
            self.log(f"Removing temporary directory {dir}")
            shutil.rmtree(dir)
        return None

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm installation of openmm and pdfixer.

        Raises
        ------
        ThirdPartyIntallationError
            When OpenMM is not installed
        ThirdPartyIntallationError
            When OpenMM pdbfixer is not installed
        """
        try:
            import openmm
        except ModuleNotFoundError:
            raise ThirdPartyIntallationError(
                "OpenMM is not installed."
                )
        try:
            import pdbfixer
        except ModuleNotFoundError:
            raise ThirdPartyIntallationError(
                "OpenMM pdbfixer is not installed."
                )
        return None
    
    @staticmethod
    def set_max_cpu(nbcpu: int) -> None:
        from openmm import Platform
        cpu_platform = Platform.getPlatformByName('CPU')
        cpu_platform.setPropertyDefaultValue('Threads', str(nbcpu))

    def _run(self) -> None:
        """Execute module."""
        # Retrieve previous models
        previous_models = self.previous_io.retrieve_models(
            individualize=True,
            )

        # create directories
        directory_dict = self.create_directories()

        # Limit cpu usage
        self.set_max_cpu(self.params["ncores"])

        # Build list of OPENMM jobs
        openmm_jobs: list[OPENMM] = []
        for i, model_to_be_simulated in enumerate(previous_models, start=1):
            # Create a openmm job
            openmm_job_i = OPENMM(
                i,
                model_to_be_simulated,
                Path("."),
                directory_dict,
                self.params,
                )
            # Hold it
            openmm_jobs.append(openmm_job_i)

        # running jobs
        scheduling_engine = get_engine(self.params["mode"], self.params)
        scheduler = scheduling_engine(openmm_jobs)
        scheduler.run()

        # Retrieve generated models
        output_pdbs = list(Path(directory_dict["openmm_output"]).glob("*.pdb"))

        # Check if at least one output file has been generated
        if len(output_pdbs) == 0:
            self.finish_with_error(
                "No output models generated. "
                "Check Openmm Execution and logfile."
                )

        # deleting unnecessary directories
        self.log("Removing unnecessary directories...")
        # self.remove_directories()

        # Add small message to the user
        self.log("Completed OpenMM module run.")
        self.log(
            "If you want to continue the haddock3 workflow after "
            "the OpenMM module, the next module should be `[topoaa]`, "
            "to rebuild the CNS molecular topologies."
            )
        
        # Setting the output variable
        self.output_models = [
            PDBFile(openmmout)
            for openmmout in sorted(output_pdbs)
            ]
        # Generating standardized haddock3 outputs
        self.export_io_models()
