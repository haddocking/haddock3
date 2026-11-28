"""DisVis integration sampling module.

Generates the complexes of two molecules by a full 6D search of rotations
and translations with DisVis, keeping those consistent with most of the
ambiguous restraints, optionally reduced to the centers of their clusters.

When the input molecules are ensembles, a separate search is run for each
combination of their models (all combinations with `crossdock = true`, else
models paired by index). Their complexes are pooled: the threshold of
consistent restraints is taken over all searches, and they are clustered
together, so that up to `sampling` models are passed on in total.

The models are written largest cluster first, and their score is minus the
number of restraints their complex (the cluster center) is consistent with.
When there are fewer clusters than `sampling`, `fill_clusters = true` adds
other members of the clusters after the centers, up to `sampling`: the
second member of each cluster, then the third, and so on. On the HADDOCK3
examples with fewer clusters than `sampling`, this improved the top-ranked
final cluster in four of five cases (but protein-peptide); it needs more
benchmarking.
A `seletop` step directly after `disvis` therefore selects the models
consistent with the most restraints, and among those (often many with the
same number) the largest clusters first, since its sort keeps the order of
equal scores. This number of restraints says little about the quality of a
model: to rank the models on the HADDOCK score before `seletop`, a short
`emref` can be added, with the restraints and few minimisation steps, e.g.::

    [disvis]
    ambig_fname = "ambig.tbl"
    [emref]
    ambig_fname = "ambig.tbl"
    nemsteps = 50
    [seletop]
    select = 200

On one test case (E2A-HPR, with a 10-model HPR ensemble) this ranked a few
more near-native models at the top than the disvis score did (4 instead of
2 of the top 10), while `emscoring`, without the restraints, ranked fewer.
This needs more benchmarking.

For a small molecule, such as a ligand of a few tens of atoms, lower
`min_interaction` (e.g. to 100): the default minimum volume of interaction
(300 A^3) may not be reached, and then no complex is found.
"""

import os
import time
from pathlib import Path

import numpy as np

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import Any, FilePath
from haddock.libs.libontology import PDBFile
from haddock.modules import BaseHaddockModule
from haddock.modules.sampling.disvis.disvis import (
    cluster_complexes,
    disvis_is_available,
    fill_from_clusters,
    read_heavy_atoms,
    pool,
    search,
    superposition,
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
        # The heavy atoms of the input models, by file name
        self._read: dict[Path, Any] = {}

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm if the module is ready to use."""
        if not disvis_is_available():
            raise Exception(
                "You are trying to use the `disvis` module but disvis>=3.4.2"
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
        combinations = self.previous_io.retrieve_models(
            crossdock=self.params["crossdock"]
        )
        for combination in combinations:
            if len(combination) != 2:
                self.finish_with_error(
                    "disvis requires exactly 2 input molecules per combination, "
                    f"got {len(combination)}: {[c.file_name for c in combination]}"
                )

        # One search per combination, pooled before the selection
        parts = []
        for idx, (receptor, ligand) in enumerate(combinations, start=1):
            fixed, scanning = self._heavy_atoms(receptor), self._heavy_atoms(ligand)
            restraints = parse_restraints(
                ambig_fname, fixed, scanning, log.info, self.params["min_maxdis"]
            )
            # The result files of DisVis, one folder per search
            directory = Path(f"disvis_search_{idx}")
            directory.mkdir(exist_ok=True)
            t0 = time.monotonic()
            complexes, rotmat = search(
                fixed,
                scanning,
                restraints,
                directory=directory,
                angle=self.params["angle"],
                voxelspacing=self.params["voxelspacing"],
                interaction_radius=self.params["interaction_radius"],
                max_clash=self.params["max_clash"],
                min_interaction=self.params["min_interaction"],
                percentage=self.params["percentage"],
                max_complexes=self.params["max_complexes"],
                ncores=self.params["ncores"],
                gpu=self.params["gpu"],
                inputs={
                    "receptor": os.path.relpath(
                        Path(receptor.path, receptor.file_name)
                    ),
                    "ligand": os.path.relpath(Path(ligand.path, ligand.file_name)),
                    "restraints": str(ambig_fname),
                    "min_maxdis": self.params["min_maxdis"],
                },
            )
            log.info(
                f"disvis search {idx}: {receptor.file_name} + {ligand.file_name}:"
                f" at most {complexes.most()} of {len(restraints)} restraints,"
                f" {len(complexes.rotation)} complexes collected"
                f" elapsed={time.monotonic() - t0:.1f}s"
            )
            parts.append(complexes)

        # Without clustering, the `sampling` complexes consistent with the
        #  most restraints
        max_complexes = self.params["max_complexes"]
        if not self.params["cluster"]:
            max_complexes = min(max_complexes, self.params["sampling"])
        complexes, part, threshold, nfound = pool(
            parts, len(rotmat), self.params["percentage"], max_complexes
        )
        log.info(
            f"disvis: {nfound} complexes of {len(parts)} searches consistent with"
            f" at least {threshold} restraints, {len(complexes.rotation)} kept"
        )

        selected = list(range(len(complexes.rotation)))
        if self.params["cluster"]:
            clusters, cutoff = self._cluster(complexes, part, rotmat, combinations)
            log.info(
                f"disvis: {len(clusters)} clusters of at least"
                f" {self.params['min_population']} complexes within"
                f" {cutoff:.2f} A RMSD"
            )
            selected = [int(members[0]) for members in clusters]
            if self.params["fill_clusters"] and len(clusters) < self.params["sampling"]:
                selected = fill_from_clusters(
                    clusters,
                    complexes.consistent,
                    complexes.key,
                    self.params["sampling"],
                )
                log.info(
                    f"disvis: {len(selected) - len(clusters)} other cluster"
                    " members added"
                )

        expected: list[PDBFile] = []
        for k, i in enumerate(selected[: self.params["sampling"]], start=1):
            receptor, ligand = combinations[part[i]]
            file_name = f"disvis_{k}.pdb"
            write_model(
                Path(receptor.path, receptor.file_name),
                Path(ligand.path, ligand.file_name),
                self._heavy_atoms(ligand).center,
                rotmat[complexes.rotation[i]],
                complexes.translation[i],
                file_name,
            )
            expected.append(
                PDBFile(
                    file_name,
                    topology=[receptor.topology, ligand.topology],
                    path=".",
                    # Lower is better: the most consistent restraints
                    #  first, then the largest clusters (stable sorts).
                    score=-int(complexes.consistent[i]),
                )
            )

        if not expected:
            self.finish_with_error(
                "disvis found no complexes (or clusters) to pass on: try a"
                " finer search, a lower `percentage` or `min_population`, or"
                " for a small molecule a lower `min_interaction`"
            )
        self.output_models = expected
        self.export_io_models()

    def _heavy_atoms(self, model: PDBFile) -> Any:
        """Return the heavy atoms of an input model, read once."""
        fname = Path(model.path, model.file_name)
        if fname not in self._read:
            self._read[fname] = read_heavy_atoms(fname)
        return self._read[fname]

    def _cluster(
        self,
        complexes: Any,
        part: np.ndarray,
        rotmat: np.ndarray,
        combinations: list,
    ) -> tuple[list[np.ndarray], float]:
        """Cluster the pooled complexes on the RMSD of the ligand.

        The poses of the receptor models after the first are superposed on
        it, and the ligand models can differ (an ensemble), with the same
        RMSD atoms (CA and P, else heavy atoms).
        """
        from disvis.complexes import rmsd_atoms

        receptors, receptor = np.unique(
            [str(Path(r.path, r.file_name)) for r, _ in combinations],
            return_inverse=True,
        )
        ligands, ligand = np.unique(
            [str(Path(m.path, m.file_name)) for _, m in combinations],
            return_inverse=True,
        )
        receptor_coor = [rmsd_atoms(self._read[Path(r)]) for r in receptors]
        ligand_coor = [rmsd_atoms(self._read[Path(m)]) for m in ligands]
        if any(
            len({c.shape for c in coor}) > 1 for coor in (receptor_coor, ligand_coor)
        ):
            self.finish_with_error(
                "disvis clusters the complexes of ensembles on their CA and P"
                " (else heavy) atoms, which differ between the models"
            )
        frames = [superposition(c, receptor_coor[0]) for c in receptor_coor]
        rotation = np.asarray([frames[r][0] for r in receptor[part]])
        shift = np.asarray([frames[r][1] for r in receptor[part]])
        return cluster_complexes(
            np.einsum("nij,nj->ni", rotation, complexes.translation) + shift,
            rotation @ rotmat[complexes.rotation],
            np.stack(
                [c - self._read[Path(m)].center for c, m in zip(ligand_coor, ligands)]
            ),
            ligand[part],
            self.params["cluster_target"],
            self.params["min_population"],
        )
