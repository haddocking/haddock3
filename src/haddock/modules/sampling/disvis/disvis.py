"""Wrapper around the `disvis` package."""

import logging
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from pdbtools.pdb_reatom import run as renumber_atoms

from haddock import log
from haddock.core.typing import Any, FilePath, Optional
from haddock.libs.libalign import kabsch


class _HaddockLog(logging.Handler):
    """Pass the records of the disvis logger on to the haddock3 log."""

    def emit(self, record: logging.LogRecord) -> None:
        log.handle(record)


# The messages of disvis (on its `disvis` logger) go to the haddock3 log only
_disvis_log = logging.getLogger("disvis")
_disvis_log.addHandler(_HaddockLog())
_disvis_log.setLevel(logging.INFO)
_disvis_log.propagate = False


def disvis_is_available() -> bool:
    """Check whether a `disvis` with `write_results` is importable."""
    try:
        import disvis.main
    except ImportError:
        return False
    # disvis >= 3.4.0, also with poses of several ligand models (3.3.0)
    return hasattr(disvis.main, "write_results")


def read_heavy_atoms(pdb_file: FilePath) -> Any:
    """Read a structure with DisVis, without its hydrogens.

    Hydrogens enlarge the clash volumes and the number of restraint atom
    pairs, so the search is done on heavy atoms only.
    """
    from disvis import PDB

    return PDB.fromfile(str(pdb_file)).select("element", "H", "!=")


def search(
    receptor: Any,
    ligand: Any,
    restraints: list,
    directory: FilePath,
    angle: float = 9.72,
    voxelspacing: float = 2.0,
    interaction_radius: float = 3.0,
    max_clash: float = 200,
    min_interaction: float = 300,
    percentage: float = 75,
    max_complexes: int = 50000,
    ncores: int = 1,
    gpu: bool = False,
    inputs: Optional[dict] = None,
) -> tuple[Any, np.ndarray]:
    """Run the DisVis search and write its result files to `directory`.

    The search is run as by the `disvis` command line, on `ncores`
    processes or on the GPU, and the files are those of the command
    (accessible_interaction_space.mrc, accessible_complexes.out,
    violations.out and z-score.out) but complexes.out, and args.json with
    the options and `inputs` (e.g. the input files). Returns the complexes
    collected, those consistent with at least `percentage` percent of the
    most restraints found (at most `max_complexes`, see
    `disvis.disvis.Complexes`), and the rotation matrices.
    """
    import disvis.main
    from disvis.helpers import dump_args
    from disvis.rotations import proportional_orientations, quat_to_rotmat

    quaternions, weights, _ = proportional_orientations(angle)
    rotmat = quat_to_rotmat(quaternions)

    # The options of the `disvis` command line read by its search functions
    args = SimpleNamespace(
        voxelspacing=voxelspacing,
        interaction_radius=interaction_radius,
        max_clash=max_clash,
        min_interaction=min_interaction,
        interaction_restraints_cutoff=max(len(restraints) - 3, 1),
        interaction_selection=None,
        occupancy_analysis=False,
        write_complexes=percentage,
        max_complexes=max_complexes,
        nproc=ncores,
        progress=False,
        directory=str(directory),
    )
    search_args = (rotmat, weights, restraints, None, None, args)

    if gpu:
        from disvis.disvis import PYOPENCL

        if not PYOPENCL:
            raise RuntimeError(
                "GPU acceleration of disvis requires pyopencl and pyvkfft"
            )
        from disvis.helpers import get_queue

        queue = get_queue() or get_queue("all")
        if queue is None:
            raise RuntimeError("No OpenCL device was found for disvis")
        results = disvis.main.run_disvis(queue, receptor, ligand, *search_args)
    elif ncores > 1:
        results = disvis.main.mp_cpu_disvis(receptor, ligand, *search_args)
    else:
        results = disvis.main.run_disvis(None, receptor, ligand, *search_args)

    # complexes.out is not written: the complexes are pooled and clustered
    disvis.main.write_results(
        results,
        rotmat,
        ligand,
        len(restraints),
        None,
        None,
        SimpleNamespace(**{**vars(args), "write_complexes": None, "delimiter": None}),
    )
    args.angle = angle
    dump_args(
        SimpleNamespace(**vars(args), **(inputs or {})), Path(directory, "args.json")
    )
    return results.complexes, rotmat


def pool(
    parts: list, nrot: int, percentage: float, max_complexes: int
) -> tuple[Any, np.ndarray, int, int]:
    """Pool the complexes of several searches and select them.

    The complexes kept are those consistent with at least `percentage`
    percent of the most restraints found in all searches (see
    `disvis.disvis.Complexes.select`): each search collected all those
    above its own, lower or equal, threshold. At most `max_complexes` are
    kept, those consistent with the most restraints first, and a sample of
    the same from run to run (the smallest sampling keys) of those
    consistent with the last number of restraints kept. Returns them, with
    `rotation` the index in the `nrot` rotations, the search of each
    complex, the threshold and the number of complexes found above it.
    """
    from disvis.disvis import Complexes

    pooled = Complexes.concatenate(parts, [k * nrot for k in range(len(parts))])
    selected, threshold, nfound = pooled.select(percentage)
    keep = np.lexsort((selected.key, -selected.consistent))[:max_complexes]
    keep.sort()
    kept = Complexes(
        selected.rotation[keep],
        selected.translation[keep],
        selected.consistent[keep],
        selected.key[keep],
        selected.counts,
    )
    part, kept.rotation = np.divmod(kept.rotation, nrot)
    return kept, part, threshold, nfound


def superposition(
    coor: np.ndarray, reference: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return the rotation S and translation s best superposing coor on reference.

    The superposed coordinates are `coor @ S.T + s`.
    """
    center, reference_center = coor.mean(axis=0), reference.mean(axis=0)
    rotation = kabsch(coor - center, reference - reference_center).T
    return rotation, reference_center - rotation @ center


def cluster_complexes(
    translations: np.ndarray,
    rotations: np.ndarray,
    coor: np.ndarray,
    model: np.ndarray,
    cluster_target: float = 6.0,
    min_population: int = 2,
) -> tuple[list[np.ndarray], float]:
    """Cluster the complexes on the RMSD of the ligand.

    `coor` gives the RMSD atoms of each ligand model, relative to its
    center, and `model` the ligand model of each complex (see
    `disvis.complexes.Poses`). Returns the clusters (indices of their
    complexes, center first), largest first, and the RMSD cutoff, tuned for clusters of at
    most `cluster_target` A, as `disvis-complexes --cluster`.
    """
    from disvis.complexes import Poses, cluster, tune_cutoff

    poses = Poses(translations, rotations, coor, model)
    cutoff = tune_cutoff(poses, cluster_target, min_population)
    return cluster(poses, cutoff, min_population), cutoff


def fill_from_clusters(
    clusters: list[np.ndarray], consistent: np.ndarray, key: np.ndarray, total: int
) -> list[int]:
    """Return the cluster centers, then other members, at most `total`.

    The other members follow in rounds over the clusters, in their order
    (largest first): the second member of each cluster, then the third, and
    so on; within a cluster, those consistent with the most restraints
    first (and the smallest sampling keys).
    """
    members = [m[1:][np.lexsort((key[m[1:]], -consistent[m[1:]]))] for m in clusters]
    selected = [int(m[0]) for m in clusters]
    for k in range(max(map(len, members), default=0)):
        if len(selected) >= total:
            break
        selected += [int(m[k]) for m in members if len(m) > k]
    return selected[:total]


def write_model(
    receptor_file: FilePath,
    ligand_file: FilePath,
    center: np.ndarray,
    rotation: np.ndarray,
    translation: np.ndarray,
    output_file: FilePath,
) -> None:
    """Write the receptor and the ligand placed as in a complex.

    The atoms of `ligand_file`, hydrogens included, are placed at
    `rotation @ (x - center) + translation`, with `center` that of the
    ligand of the search. The other columns are kept as in the input files
    (segids included, which CNS uses), and the atom serials renumbered.
    """
    receptor_lines = _atom_lines(receptor_file)
    ligand_lines = _atom_lines(ligand_file)
    coor = np.array(
        [[float(line[i : i + 8]) for i in (30, 38, 46)] for line in ligand_lines]
    )
    coor = (coor - center) @ rotation.T + translation
    placed = [
        f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
        for line, (x, y, z) in zip(ligand_lines, coor, strict=True)
    ]
    lines = receptor_lines + placed + ["END\n"]
    Path(output_file).write_text("".join(renumber_atoms(lines, 1)))


def _atom_lines(pdb_file: FilePath) -> list[str]:
    """Return the ATOM and HETATM lines of a PDB file."""
    with open(pdb_file) as fin:
        return [line for line in fin if line.startswith(("ATOM", "HETATM"))]
