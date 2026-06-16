"""gdock integration sampling module."""

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

    def _identify_receptor_ligand(
        self, combination: tuple[PDBFile, ...], rec_chain: str, lig_chain: str
    ) -> tuple[PDBFile, PDBFile]:
        """Identify the receptor and ligand PDBFile in a model combination."""
        by_chain: dict[str, PDBFile] = {}
        for pdb in combination:
            _, chains = libpdb.identify_chainseg(Path(pdb.path, pdb.file_name))
            for chain in chains:
                by_chain.setdefault(chain, pdb)

        receptor = by_chain.get(rec_chain)
        ligand = by_chain.get(lig_chain)
        if receptor is None or ligand is None:
            self.finish_with_error(
                f"Could not find receptor chain '{rec_chain}' and ligand"
                f" chain '{lig_chain}' among the input models"
                f" {[c.file_name for c in combination]}"
            )

        return receptor, ligand  # type: ignore

    def _run(self) -> None:
        """Execute module."""
        # Each combination is a tuple with one PDBFile per input molecule,
        #  e.g. (receptor, ligand), as produced by `topoaa`.
        models_to_dock = self.previous_io.retrieve_models()

        rec_chain = self.params["receptor_chains"][0]
        lig_chain = self.params["ligand_chains"][0]

        # Convert restraints, if provided, to gdock's residue pair format
        restraints = None
        ambig_fname = self.params["ambig_fname"]
        if ambig_fname:
            restraints = parse_restraints(ambig_fname, rec_chain, lig_chain)

        sampling_factor = max(1, self.params["sampling"] // len(models_to_dock))

        log.info(f"ncores={self.params['ncores']} sampling={self.params['sampling']}")

        expected: list[PDBFile] = []
        for idx, combination in enumerate(models_to_dock, start=1):
            receptor, ligand = self._identify_receptor_ligand(
                combination, rec_chain, lig_chain
            )

            gdock_wrapper = GdockWrapper(
                receptor_pdb_file=Path(receptor.path, receptor.file_name),
                ligand_pdb_file=Path(ligand.path, ligand.file_name),
                restraints=restraints,
                max_generations=self.params["max_generations"],
                number_of_individuals=self.params["number_of_individuals"],
                seed=self.params["seed"],
                sampling=sampling_factor,
            )
            gdock_wrapper.run()
            models = gdock_wrapper.save_models(".", prefix=f"gdock_{idx}")

            for m in models:
                expected.append(
                    PDBFile(
                        m["file_name"],
                        topology=[receptor.topology, ligand.topology],
                        path=self.path,
                        score=m["fitness"],
                        unw_energies={
                            "vdw": m["vdw"],
                            "elec": m["elec"],
                            "desolv": m["desolv"],
                            "air": m["air"],
                        },
                    )
                )

        self.output_models = expected
        self.export_io_models()
