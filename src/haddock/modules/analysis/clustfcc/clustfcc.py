from pathlib import Path
from typing import Any

from fccpy import read_pdb
from fccpy.clustering import disjoint_taylor_butina
from fccpy.contacts import (
    BY_RESIDUE,
    get_intermolecular_contacts,
    hash_many,
    read_contacts,
    write_contacts,
)
from fccpy.similarity import build_matrix, fcc, read_matrix, write_matrix

from haddock.libs.libontology import PDBFile


def calculate_contacts(pdb_list: list[PDBFile], cutoff: float, path: str) -> list[Path]:
    """Calculate contacts."""
    contact_list = []
    for pdb in pdb_list:
        pdb_fname = Path(pdb.rel_path)
        contact_fname = Path(pdb.file_name.replace(".pdb", ".con"))
        s = read_pdb(pdb_fname)
        clist = get_intermolecular_contacts(structure=s, max_distance=cutoff)
        write_contacts(atom_pairs=clist, filepath=contact_fname)
        contact_list.append(contact_fname)
    return contact_list


def calculate_fcc_matrix(contact_list: list[Path], matrix_outfile: Path) -> Path:
    """Calculate fcc matrix."""
    unique = True
    selector1 = selector2 = BY_RESIDUE
    contact_list.sort(key=lambda x: int(str(x).split("_")[-1].split(".")[0]))
    clist = [list(read_contacts(f)) for f in contact_list]
    clist_hashed = [
        hash_many(c, unique=unique, selector1=selector1, selector2=selector2)
        for c in clist
    ]

    idxs, sims = build_matrix(clist_hashed, metric=fcc)

    matrix_outfile.touch()

    write_matrix(idxs, sims, matrix_outfile)

    return matrix_outfile


def cluster_fcc(matrix_fname: Path, similarity: float, minsize: int) -> dict[Any, int]:
    """Cluster fcc."""
    idxs, sims = read_matrix(filepath=matrix_fname)
    labels = disjoint_taylor_butina(idxs, sims, eps=similarity, minsize=minsize)

    # write_clusters(labels, cluster_outfile)

    # return cluster_outfile
    return labels
