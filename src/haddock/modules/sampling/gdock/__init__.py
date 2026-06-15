"""gdock integration sampling module."""
import shutil
from pathlib import Path

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, FilePath
from haddock.libs import libpdb
from haddock.libs.libontology import PDBFile
from haddock.modules import BaseHaddockModule
from haddock.modules.sampling.gdock.gdock import (
    GdockWrapper,
    gdock_is_available,
    parse_restraints,
)


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 gdock module."""

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
        """Confirm if the module is ready to use."""
        if not gdock_is_available():
            raise Exception(
                "You are trying to use the `gdock` module but it is not"
                " available, please check the installation instructions"
            )

    def _run(self) -> None:
        """Execute module."""
        models = self.previous_io.retrieve_models()

        if len(models) > 1:
            self.finish_with_error("Only one model allowed in gdock sampling module")

        model = models[0]

        # Make sure chain IDs are present, copy a local working copy
        _path = Path(model.path, model.file_name)
        segids, chains = libpdb.identify_chainseg(_path)
        if set(segids) != set(chains):
            log.info("No chain IDs found, using segid information")
            libpdb.swap_segid_chain(
                Path(model.path, model.file_name),
                Path(self.path, model.file_name),
            )
        elif Path(model.path, model.file_name).resolve() != Path(
            self.path, model.file_name
        ).resolve():
            shutil.copyfile(
                Path(model.path, model.file_name),
                Path(self.path, model.file_name),
            )

        model_with_chains = self.path / model.file_name

        # Split the complex into receptor and ligand structures
        new_models = libpdb.split_by_chain(model_with_chains)
        if model_with_chains in new_models:
            self.finish_with_error(
                f"Input {model_with_chains} cannot be split by chain"
            )

        rec_chain = self.params["receptor_chains"][0]
        lig_chain = self.params["ligand_chains"][0]
        receptor_pdb_file = self.path / f"{Path(model.file_name).stem}_{rec_chain}.pdb"
        ligand_pdb_file = self.path / f"{Path(model.file_name).stem}_{lig_chain}.pdb"

        # Convert restraints, if provided, to gdock's residue pair format
        restraints = None
        ambig_fname = self.params["ambig_fname"]
        if ambig_fname:
            restraints = parse_restraints(ambig_fname, rec_chain, lig_chain)

        log.info("Running gdock")
        gdock_wrapper = GdockWrapper(
            receptor_pdb_file=receptor_pdb_file,
            ligand_pdb_file=ligand_pdb_file,
            restraints=restraints,
            max_generations=self.params["max_generations"],
            seed=self.params["seed"],
        )
        gdock_wrapper.run()
        poses = gdock_wrapper.save_poses(self.path, self.params["sampling"])

        expected: list[PDBFile] = []
        for pose in poses:
            expected.append(
                PDBFile(
                    pose["file_name"],
                    topology=model.topology,
                    path=self.path,
                    score=pose["fitness"],
                    unw_energies={
                        "vdw": pose["vdw"],
                        "elec": pose["elec"],
                        "desolv": pose["desolv"],
                        "air": pose["air"],
                    },
                )
            )

        self.output_models = expected
        self.export_io_models()
