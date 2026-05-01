"""RMSD calculation and output for the rmsdfilter module."""

import copy
from math import isnan
from pathlib import Path

import numpy as np

from haddock import log
from haddock.libs.libalign import (
    ALIGNError,
    calc_rmsd,
    centroid,
    get_atoms,
    kabsch,
    load_coords,
)
from haddock.libs.libontology import PDBFile


class RMSDFilter:
    """Compute RMSD between model and reference."""

    def __init__(
        self,
        identificator: int,
        model: PDBFile,
        reference: Path,
        path: Path,
        params: dict,
        align_func,
        ref_id: int = 1,
    ) -> None:
        self.identificator = identificator
        self.model = model
        self.reference = reference
        self.path = path
        self.params = params
        self.align_func = align_func
        self.ref_id = ref_id
        self.rmsd = float("nan")

    def run(self) -> "RMSDFilter":
        """Compute RMSD and return self."""
        allatoms = self.params["allatoms"]
        keep_hetatm = self.params["keep_hetatm"]

        # map model residue numbers onto reference residue numbers
        try:
            model2ref_numbering, model2ref_chain_dict = self.align_func(
                self.reference, self.model, self.path
            )
        except ALIGNError:
            log.warning(
                f"Alignment failed between {self.reference} "
                f"and {self.model}, skipping..."
            )
            # deepcopy so the scheduler gets an independent object with rmsd=nan
            return copy.deepcopy(self)

        atoms = get_atoms(self.model, full=allatoms)
        atoms.update(get_atoms(self.reference, full=allatoms))

        ref_coord_dic, _ = load_coords(
            self.reference,
            atoms,
            keep_hetatm=keep_hetatm,
        )
        try:
            # numbering_dic remaps model residues to reference numbering so that
            # coordinate keys are directly comparable between the two dicts
            mod_coord_dic, _ = load_coords(
                self.model,
                atoms,
                numbering_dic=model2ref_numbering,
                model2ref_chain_dict=model2ref_chain_dict,
                keep_hetatm=keep_hetatm,
            )
        except ALIGNError as e:
            log.warning(e)
            return copy.deepcopy(self)

        # Only superpose on shared atoms between model and reference
        common_keys = ref_coord_dic.keys() & mod_coord_dic.keys()
        if not common_keys:
            log.warning(
                f"No common atoms found between {self.reference} and {self.model}"
            )
            return copy.deepcopy(self)

        Q = np.asarray([ref_coord_dic[k] for k in common_keys])
        P = np.asarray([mod_coord_dic[k] for k in common_keys])

        # Kabsch superposition: centre both structures Q (ref) and P (model),
        # find optimal rotation U, apply it to P, then compute RMSD on the superposed coordinates.
        Q = Q - centroid(Q)
        P = P - centroid(P)
        U = kabsch(P, Q)
        P = np.dot(P, U)

        self.rmsd = calc_rmsd(P, Q)
        return copy.deepcopy(self)

# helper functions to not clutter init.py
def _sort_key(row: tuple, sortby: str) -> float:
    """Return the sort value for a (model, rmsd) row."""
    model, rmsd = row
    if sortby == "score":
        val = model.score
        if val is None or (isinstance(val, float) and isnan(val)):
            return float("inf")
        return val
    return float("inf") if isnan(rmsd) else rmsd


def collect_rmsd_map(results: list) -> dict[int, float]:
    """Filter per-job results to minimum RMSD per model across all references."""
    rmsd_map: dict[int, float] = {}
    for job in results:
        current = rmsd_map.get(job.identificator, float("nan"))
        job_rmsd = job.rmsd
        if isnan(current):
            rmsd_map[job.identificator] = job_rmsd
        elif not isnan(job_rmsd):
            rmsd_map[job.identificator] = min(current, job_rmsd)
    return rmsd_map


def build_sorted_rows(
    rmsd_map: dict[int, float],
    id_to_model: dict[int, PDBFile],
    sortby: str,
    sort_ascending: bool,
) -> list[tuple]:
    """Build a list of (model, rmsd) pairs sorted by sortby."""
    rows = [(id_to_model[i], rmsd_map[i]) for i in sorted(rmsd_map.keys())]
    rows.sort(key=lambda row: _sort_key(row, sortby), reverse=not sort_ascending)
    return rows


def write_rmsdfilter_ss(
    rows: list[tuple],
    filtered: list[PDBFile],
    percent_filtered: float,
    threshold: float,
    fname: str = "rmsdfilter_ss.tsv",
) -> None:
    """Write _ss_.tsv."""
    with open(fname, "w") as fh:
        fh.write(
            f"# RMSD filtering threshold is set to {threshold:.3f} Å; "
            f"{len(filtered)} model(s) were kept; {percent_filtered:.2f}% were filtered out. "
            "This file contains all models for user information.\n"
        )
        fh.write("model\tscore\trmsd\n")
        for model, rmsd in rows:
            score_str = (
                f"{model.score:.3f}"
                if model.score is not None
                and not (isinstance(model.score, float) and isnan(model.score))
                else "nan"
            )
            rmsd_str = f"{rmsd:.3f}" if not isnan(rmsd) else "nan"
            fh.write(f"{model.rel_path}\t{score_str}\t{rmsd_str}\n")


def write_rmsdfilter_multiref(
    results: list,
    id_to_model: dict[int, PDBFile],
    fname: str = "rmsdfilter_ss_multiref.tsv",
) -> None:
    """Write multiref file: one row per (model, reference) pair, sorted by model name then ref_id."""
    multiref_rows = sorted(
        results,
        key=lambda job: (str(id_to_model[job.identificator].rel_path), job.ref_id),
    )
    with open(fname, "w") as fh:
        fh.write("model\tref_id\trmsd\n")
        for job in multiref_rows:
            model = id_to_model[job.identificator]
            rmsd_str = f"{job.rmsd:.3f}" if not isnan(job.rmsd) else "nan"
            fh.write(f"{model.rel_path}\t{job.ref_id}\t{rmsd_str}\n")
