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
    """Wrapper to run `gdock`'s genetic algorithm docking and save its models."""

    def __init__(
        self,
        receptor_pdb_file: FilePath,
        ligand_pdb_file: FilePath,
        restraints: Optional[list[tuple[int, int]]] = None,
        max_generations: int = 250,
        number_of_individuals: int = 150,
        ncores: int = 1,
        seed: int = 42,
        sampling: Optional[int] = None,
    ) -> None:
        self.receptor_pdb_file = Path(receptor_pdb_file)
        self.ligand_pdb_file = Path(ligand_pdb_file)
        self.restraints = restraints
        self.max_generations = max_generations
        self.number_of_individuals = number_of_individuals
        self.ncores = ncores
        self.seed = seed
        self.sampling = sampling
        self.result: Optional[dict] = None
        self.converged_early: bool = False

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
            population_size=self.number_of_individuals,
            ncores=self.ncores,
            seed=self.seed,
            sampling=self.sampling,
        )
        self.converged_early = self.result.get("convergedEarly", False)

    def save_models(self, output_dir: FilePath, prefix: str = "gdock") -> list[dict]:
        """Write all ranked models as PDB files to `output_dir`.

        Each model PDB already contains the full receptor+ligand complex as
        returned by gdock. The segID column (cols 73-76) is filled from the
        chain ID before writing so that CNS can map atoms to the PSF, and atom
        serials are renumbered to ensure they are strictly increasing.

        Returns
        -------
        list[dict]
            One entry per saved model with keys `file_name`, `fitness`,
            `vdw`, `elec`, `desolv` and `air`.
        """
        if self.result is None:
            raise RuntimeError("`run` must be called before `save_models`")

        saved = []
        for model in self.result["models"]:
            file_name = f"{prefix}_{model['rank']}.pdb"
            # Fill in the segID column (cols 73-76) from chain ID: CNS relies
            # on segIDs to map atoms to the PSF topology.
            complex_pdb = "".join(chain_to_seg(model["pdb"].splitlines(keepends=True)))
            # Renumber atoms so serials are strictly increasing across the
            # combined receptor+ligand complex, as CNS's PDB reader requires.
            complex_pdb = "".join(
                renumber_atoms(complex_pdb.splitlines(keepends=True), 1)
            )
            Path(output_dir, file_name).write_text(complex_pdb)
            saved.append(
                {
                    "file_name": file_name,
                    "fitness": model["fitness"],
                    "vdw": model["vdw"],
                    "elec": model["elec"],
                    "desolv": model["desolv"],
                    "air": model["air"],
                }
            )

        return saved
