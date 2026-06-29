"""gdock integration sampling module."""

import time
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

    def _split_combination(
        self, combination: tuple[PDBFile, ...]
    ) -> tuple[PDBFile, str, PDBFile, str]:
        """Return (receptor, rec_chain, ligand, lig_chain) from a combination.

        The molecule order follows the `molecules` list in the run config, so
        index 0 is always the receptor and index 1 is always the ligand.
        Chain IDs are read directly from the files.
        """
        if len(combination) != 2:
            self.finish_with_error(
                f"gdock requires exactly 2 input molecules per combination, "
                f"got {len(combination)}: {[c.file_name for c in combination]}"
            )

        receptor, ligand = combination[0], combination[1]

        _, rec_chains = libpdb.identify_chainseg(
            Path(receptor.path, receptor.file_name)
        )
        _, lig_chains = libpdb.identify_chainseg(Path(ligand.path, ligand.file_name))

        if not rec_chains or not lig_chains:
            self.finish_with_error(
                f"Could not detect chain IDs in receptor '{receptor.file_name}' "
                f"or ligand '{ligand.file_name}'"
            )

        return receptor, rec_chains[0], ligand, lig_chains[0]  # type: ignore

    def _run(self) -> None:
        """Execute module."""
        # Each combination is a tuple with one PDBFile per input molecule,
        #  e.g. (receptor, ligand), as produced by `topoaa`.
        models_to_dock = self.previous_io.retrieve_models(
            crossdock=self.params["crossdock"]
        )

        sampling_factor = max(1, self.params["sampling"] // len(models_to_dock))

        log.info(f"ncores={self.params['ncores']} sampling={self.params['sampling']}")

        expected: list[PDBFile] = []
        for idx, combination in enumerate(models_to_dock, start=1):
            receptor, rec_chain, ligand, lig_chain = self._split_combination(
                combination
            )

            # Convert restraints using chain IDs read from the actual files.
            restraints = None
            ambig_fname = self.params["ambig_fname"]
            if ambig_fname:
                restraints = parse_restraints(ambig_fname, rec_chain, lig_chain)

            gdock_wrapper = GdockWrapper(
                receptor_pdb_file=Path(receptor.path, receptor.file_name),
                ligand_pdb_file=Path(ligand.path, ligand.file_name),
                restraints=restraints,
                max_generations=self.params["max_generations"],
                number_of_individuals=self.params["number_of_individuals"],
                ncores=self.params["ncores"],
                seed=self.params["seed"],
                sampling=sampling_factor,
            )
            t0 = time.monotonic()
            gdock_wrapper.run()
            elapsed = time.monotonic() - t0
            log.info(
                f"gdock run {idx}: {receptor.file_name} + {ligand.file_name}"
                f" generations={gdock_wrapper.result['generationsRun']}"
                f" models={len(gdock_wrapper.result['models'])}"
                f" elapsed={elapsed:.1f}s"
            )
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
