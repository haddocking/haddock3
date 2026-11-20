"""DisVis integration sampling module.

Generates the complexes of two molecules by a full 6D search of rotations
and translations with DisVis, keeping those consistent with most of the
ambiguous restraints, optionally reduced to the centers of their clusters.

When the input molecules are ensembles, a separate search is run for each
combination of their models (all combinations with `crossdock = true`, else
models paired by index), and up to `sampling` models are passed on for each
combination: e.g. 1 receptor and a 10-model ligand ensemble give 10 searches
and up to 10 x `sampling` models.
"""

import time
from pathlib import Path

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, FilePath
from haddock.libs.libontology import PDBFile
from haddock.modules import BaseHaddockModule
from haddock.modules.sampling.disvis.disvis import (
    cluster_centers,
    disvis_is_available,
    read_heavy_atoms,
    search,
    write_model,
)


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 disvis module."""

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
        if not disvis_is_available():
            raise Exception(
                "You are trying to use the `disvis` module but disvis>=3.2.0"
                " is not available, install it with"
                " `pip install 'haddock3[disvis]'` (Python >= 3.11)"
            )

    def _run(self) -> None:
        """Execute module."""
        from disvis.helpers import parse_restraints

        ambig_fname = self.params["ambig_fname"]
        if not ambig_fname:
            self.finish_with_error("disvis requires an `ambig_fname` restraints file")

        # Each combination is a tuple with one PDBFile per input molecule,
        #  e.g. (receptor, ligand), as produced by `topoaa`.
        models_to_dock = self.previous_io.retrieve_models(
            crossdock=self.params["crossdock"]
        )

        expected: list[PDBFile] = []
        for idx, combination in enumerate(models_to_dock, start=1):
            if len(combination) != 2:
                self.finish_with_error(
                    "disvis requires exactly 2 input molecules per combination, "
                    f"got {len(combination)}: {[c.file_name for c in combination]}"
                )
            receptor, ligand = combination
            receptor_file = Path(receptor.path, receptor.file_name)
            ligand_file = Path(ligand.path, ligand.file_name)
            fixed = read_heavy_atoms(receptor_file)
            scanning = read_heavy_atoms(ligand_file)
            restraints = parse_restraints(
                ambig_fname,
                fixed,
                scanning,
                log.info,
                self.params["min_maxdis"],
            )

            t0 = time.monotonic()
            complexes, rotmat, threshold, nfound = search(
                fixed,
                scanning,
                restraints,
                directory=".",
                angle=self.params["angle"],
                voxelspacing=self.params["voxelspacing"],
                interaction_radius=self.params["interaction_radius"],
                max_clash=self.params["max_clash"],
                min_interaction=self.params["min_interaction"],
                percentage=self.params["percentage"],
                max_complexes=self.params["max_complexes"],
                ncores=self.params["ncores"],
                gpu=self.params["gpu"],
            )
            log.info(
                f"disvis run {idx}: {receptor.file_name} + {ligand.file_name}"
                f" {nfound} complexes consistent with at least {threshold}"
                f" of {len(restraints)} restraints, {len(complexes.rotation)} kept"
                f" elapsed={time.monotonic() - t0:.1f}s"
            )

            selected = list(range(len(complexes.rotation)))
            if self.params["cluster"]:
                selected, cutoff = cluster_centers(
                    complexes,
                    rotmat,
                    scanning,
                    self.params["cluster_target"],
                    self.params["min_population"],
                )
                log.info(
                    f"disvis run {idx}: {len(selected)} clusters of at least"
                    f" {self.params['min_population']} complexes within"
                    f" {cutoff:.2f} A RMSD"
                )

            for k, i in enumerate(selected[: self.params["sampling"]], start=1):
                file_name = f"disvis_{idx}_{k}.pdb"
                write_model(
                    receptor_file,
                    ligand_file,
                    scanning.center,
                    rotmat[complexes.rotation[i]],
                    complexes.translation[i],
                    file_name,
                )
                expected.append(
                    PDBFile(
                        file_name,
                        topology=[receptor.topology, ligand.topology],
                        path=self.path,
                        # Lower is better: the most consistent restraints
                        #  first, then the largest clusters (stable sorts).
                        score=-int(complexes.consistent[i]),
                    )
                )

        if not expected:
            self.finish_with_error(
                "disvis found no complexes (or clusters) to pass on: try a"
                " finer search, a lower `percentage` or `min_population`"
            )
        self.output_models = expected
        self.export_io_models()
