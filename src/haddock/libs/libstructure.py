"""Molecular data structures."""
import os
import string
from pathlib import Path

from pdbtools import pdb_chain, pdb_chainxseg, pdb_segxchain, pdb_tidy

from haddock.core.exceptions import MoleculeError
from haddock.libs.libpdb import slc_chainid, slc_segid


class Molecule:
    """Input molecule, usually a PDB file."""

    def __init__(self, file_path):
        # the rest of the code is too dependent on the Path API
        assert isinstance(file_path, Path), "`file_name` must be pathlib.Path"

        self.file_path = file_path

    def read_lines(self):
        """Read lines of Molecule file to attribute `lines`."""
        self.lines = [
            _l + os.linesep  # adds os.linesep for compatiblity with pdb-tools
            for _l in self.file_path.read_text().split(os.linesep)
            ]

    def read_chain(self):
        """
        Read the chain ID of the molecule.

        Requires `.read_lines()` before.

        Saves chainID to attribute `chaindid`.
        """
        # read the chain ID
        chainids = set(
            line[slc_chainid].strip()
            for line in self.lines
            if line.startswith(('ATOM', 'HETATM', 'ANISOU'))
            )

        if len(chainids) > 1:
            raise ValueError(f'PDB files can have only one chain. {chainids}')
        self.chainid = chainids.pop()

    def read_segid(self):
        """
        Read the seg ID of the molecule.

        Requires `.read_lines()` before.

        Saves segID to attribute `segid`.
        """
        # read the segid
        segids = set(
            line[slc_segid].strip()
            for line in self.lines
            if line.startswith(('ATOM', 'HETATM', 'ANISOU'))
            )

        if len(segids) > 1:
            raise ValueError('PDB files can have only one segid.')
        self.segid = segids.pop()


def make_molecules(paths):
    """Get input molecules from the data stream."""
    return list(map(Molecule, paths))


def clean_chainID_segID(molecules):
    """
    Clean molecules chainID and segID.

    Reads Molecules chainID and segID. If one of the two is missing,
    fill the missing gap with the information of the other.

    If both chainID and segID are missing, fills both with an uppercase
    letter what as not been used yet by any of the `molecules`.

    Modified lines are saved in molecules.lines attribute.

    New lines are corrected by `pdbtools.pdb_tidy`.

    Parameters
    ----------
    molecules : list
        A list of :class:`Molecules` objects.

    Raises
    ------
    MolecureError
        If chainID and segID differ within the same molecule.
    """
    for mol in molecules:
        mol.read_lines()
        mol.read_chain()
        mol.read_segid()

    for mol in molecules:
        if mol.chainid != mol.segid:
            raise MoleculeError(
                'chainID and segID differ for molecule: '
                f'{str(mol.file_path)!r}')

    # i don't use sets beause i want to keep order
    chainids = list(mol.chainid for mol in molecules)
    uc = list(string.ascii_uppercase)
    remaining_chains = [c for c in uc if c not in chainids][::-1]

    for mol in molecules:
        if not mol.chainid and mol.segid:
            lines = pdb_segxchain(mol.lines)

        elif mol.chainid and not mol.segid:
            lines = pdb_chainxseg(mol.lines)

        elif not mol.chainid and not mol.segid:
            _chains = pdb_chain.run(mol.lines, remaining_chains.pop())
            lines = pdb_chainxseg.run(_chains)

        elif mol.chainid != mol.segid:
            # repeats the error...just in case
            # this will likely drop when we write tests
            raise MoleculeError(
                'chainID and segID differ for molecule: '
                f'{str(mol.file_path)!r}')

        else:
            # nothing to change
            mol.lines = [_l for _l in pdb_tidy.run(mol.lines) if _l.strip()]
            continue

        mol.lines = [_l for _l in pdb_tidy.run(lines) if _l.strip()]
