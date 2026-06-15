"""Wrapper around the `gdock` package."""

import re
from pathlib import Path

from pdbtools.pdb_chainxseg import run as chain_to_seg
from pdbtools.pdb_reatom import run as renumber_atoms

from haddock.core.typing import FilePath, Optional
from haddock.libs.librestraints import extract_restraint_entries


_SELECTION_RE = re.compile(
    r"\(\s*(?:resid?\s+(\d+)\s+and\s+segid\s+(\w+)"
    r"|segid\s+(\w+)\s+and\s+resid?\s+(\d+))\s*\)",
    re.IGNORECASE,
)


def gdock_is_available() -> bool:
    """Check whether the `gdock` package is importable."""
    try:
        import gdock  # noqa: F401
    except ImportError:
        raise
    return True


def _extract_selections(entry: str) -> list[tuple[str, int]]:
    """Extract `(segid, resi)` selections from a TBL `assign` entry, in order."""
    selections: list[tuple[str, int]] = []
    for match in _SELECTION_RE.finditer(entry):
        resi_a, segid_a, segid_b, resi_b = match.groups()
        if resi_a is not None:
            selections.append((segid_a, int(resi_a)))
        else:
            selections.append((segid_b, int(resi_b)))
    return selections


def parse_restraints(
    tbl_filename: FilePath,
    receptor_chain: str,
    ligand_chain: str,
) -> list[tuple[int, int]]:
    """Convert a TBL ambiguous restraints file into gdock restraint pairs.

    Each `assign` entry is treated as one anchor selection followed by one
    or more (possibly OR-combined) partner selections, which is the standard
    HADDOCK AIR format. A pair `(receptor_resseq, ligand_resseq)` is emitted
    for every anchor/partner combination spanning `receptor_chain` and
    `ligand_chain`.

    Parameters
    ----------
    tbl_filename : FilePath
        Path to the AIR.tbl restraints file.
    receptor_chain : str
        Chain ID of the receptor.
    ligand_chain : str
        Chain ID of the ligand.

    Returns
    -------
    list[tuple[int, int]]
        Sorted list of `(receptor_resseq, ligand_resseq)` pairs.
    """
    pairs: set[tuple[int, int]] = set()
    for entry in extract_restraint_entries(tbl_filename):
        selections = _extract_selections(entry)
        if len(selections) < 2:
            continue

        anchor_chain, anchor_resi = selections[0]
        for chain, resi in selections[1:]:
            if anchor_chain == receptor_chain and chain == ligand_chain:
                pairs.add((anchor_resi, resi))
            elif anchor_chain == ligand_chain and chain == receptor_chain:
                pairs.add((resi, anchor_resi))

    return sorted(pairs)


class GdockWrapper:
    """Wrapper to run `gdock`'s genetic algorithm docking and save its poses."""

    def __init__(
        self,
        receptor_pdb_file: FilePath,
        ligand_pdb_file: FilePath,
        restraints: Optional[list[tuple[int, int]]] = None,
        max_generations: int = 250,
        seed: int = 42,
    ) -> None:
        self.receptor_pdb_file = Path(receptor_pdb_file)
        self.ligand_pdb_file = Path(ligand_pdb_file)
        self.restraints = restraints
        self.max_generations = max_generations
        self.seed = seed
        self.result: Optional[dict] = None

    def run(self) -> None:
        """Run gdock's docking pipeline."""
        import gdock

        receptor_pdb = self.receptor_pdb_file.read_text()
        ligand_pdb = self.ligand_pdb_file.read_text()

        self.result = gdock.dock(
            receptor_pdb,
            ligand_pdb,
            restraints=self.restraints,
            max_generations=self.max_generations,
            seed=self.seed,
        )

    def save_poses(
        self, output_dir: FilePath, top: int, prefix: str = "gdock"
    ) -> list[dict]:
        """Write the top-ranked poses as receptor+ligand complex PDB files.

        gdock only returns the transformed ligand structure for each pose,
        so the (unchanged) receptor structure is prepended to reconstruct
        the full complex.

        Returns
        -------
        list[dict]
            One entry per saved pose with keys `file_name`, `fitness`,
            `vdw`, `elec`, `desolv` and `air`.
        """
        if self.result is None:
            raise RuntimeError("`run` must be called before `save_poses`")

        # Drop the receptor's trailing END record: it would otherwise
        # terminate CNS's coordinate reading before the ligand atoms below it.
        receptor_lines = [
            line
            for line in self.receptor_pdb_file.read_text().splitlines(keepends=True)
            if line.strip() != "END"
        ]
        receptor_pdb = "".join(receptor_lines)
        poses = sorted(self.result["poses"], key=lambda p: p["rank"])[:top]

        saved = []
        for pose in poses:
            file_name = f"{prefix}_{pose['rank']}.pdb"
            # gdock's output PDB for the ligand pose lacks the segID column
            # (cols 73-76), which CNS relies on to map atoms to the PSF. Fill
            # it in from the chain ID column before assembling the complex.
            ligand_pdb = "".join(chain_to_seg(pose["pdb"].splitlines(keepends=True)))
            complex_pdb = receptor_pdb + ligand_pdb
            # gdock numbers the ligand pose atoms starting from 1, which
            # collides with the receptor's atom serials and confuses CNS's
            # PDB reader (it expects strictly increasing serial numbers).
            complex_pdb = "".join(
                renumber_atoms(complex_pdb.splitlines(keepends=True), 1)
            )
            Path(output_dir, file_name).write_text(complex_pdb)
            saved.append(
                {
                    "file_name": file_name,
                    "fitness": pose["fitness"],
                    "vdw": pose["vdw"],
                    "elec": pose["elec"],
                    "desolv": pose["desolv"],
                    "air": pose["air"],
                }
            )

        return saved
