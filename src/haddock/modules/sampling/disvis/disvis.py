"""Wrapper around the `disvis` package."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
from pdbtools.pdb_reatom import run as renumber_atoms

from haddock.core.typing import Any, FilePath


def disvis_is_available() -> bool:
    """Check whether a `disvis` with complexes collection is importable."""
    try:
        from disvis.complexes import Poses  # noqa: F401
    except ImportError:
        return False
    return True


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
) -> tuple[Any, np.ndarray, int, int]:
    """Run the DisVis search and select the complexes.

    The search is run as by the `disvis` command line, on `ncores`
    processes or on the GPU. Returns the complexes consistent with at least
    `percentage` percent of the most restraints found (at most
    `max_complexes`, see `disvis.disvis.Complexes.select`), the rotation
    matrices, the threshold of consistent restraints and the number of
    complexes found above it.
    """
    from disvis.main import mp_cpu_disvis, run_disvis
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
        results = run_disvis(queue, receptor, ligand, *search_args)
    elif ncores > 1:
        results = mp_cpu_disvis(receptor, ligand, *search_args)
    else:
        results = run_disvis(None, receptor, ligand, *search_args)

    selected, threshold, nfound = results.complexes.select(percentage, max_complexes)
    return selected, rotmat, threshold, nfound


def cluster_centers(
    complexes: Any,
    rotmat: np.ndarray,
    ligand: Any,
    cluster_target: float = 6.0,
    min_population: int = 2,
) -> tuple[list[int], float]:
    """Cluster the complexes on the RMSD of the ligand.

    Returns the indices of the cluster centers, largest cluster first, and
    the RMSD cutoff, tuned for clusters of at most `cluster_target` A, as
    `disvis-complexes --cluster`.
    """
    from disvis.complexes import Poses, cluster, rmsd_atoms, tune_cutoff

    poses = Poses(
        complexes.translation,
        rotmat[complexes.rotation],
        rmsd_atoms(ligand) - ligand.center,
    )
    cutoff = tune_cutoff(poses, cluster_target, min_population)
    clusters = cluster(poses, cutoff, min_population)
    return [int(members[0]) for members in clusters], cutoff


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
