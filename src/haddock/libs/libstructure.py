"""Molecular data structures."""
import os
import string
from pathlib import Path

from pdbtools import pdb_chainxseg, pdb_segxchain, pdb_chain


slc_record = slice(0, 6)
slc_serial = slice(6, 11)
slc_name = slice(12, 16)
slc_altLoc = slice(16, 17)
slc_resName = slice(17, 20)
slc_chainID = slice(21, 22)
slc_resSeq = slice(22, 26)
slc_iCode = slice(26, 27)
slc_x = slice(30, 38)
slc_y = slice(38, 46)
slc_z = slice(46, 54)
slc_occ = slice(54, 60)
slc_temp = slice(60, 66)
slc_segid = slice(72, 76)
slc_element = slice(76, 78)
slc_model = slice(78, 80)


class Molecule:
    """Input molecule, usually a PDB file."""

    def __init__(self, file_path):
        # the rest of the code is too dependent on the Path API
        assert isinstance(file_path, Path), "`file_name` must be pathlib.Path"

        self.file_path = file_path

    def read_lines(self):
        self.lines = self.file_path.read_text().split(os.linesep)

    def read_chain(self):
        # read the chain ID
        chainids = set(
            line[slc_chainID].strip()
            for line in self.lines
            if line.startswith(('ATOM', 'HETATM', 'ANISOU'))
            )

        if len(chainids) > 1:
            raise ValueError(f'PDB files can have only one chain. {chainids}')
        self.chainid = chainids.pop()

    def read_segid(self):
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
    """Clean molecules chainID and segID."""

    for mol in molecules:
        mol.read_lines()
        mol.read_chain()
        mol.read_segid()

    chainids = list(mol.chainid for mol in molecules)
    segids = list(mol.segid for mol in molecules)

    if chainids != segids:
        raise ValueError

    # i don't use sets beause i want to keep order
    uc = list(string.ascii_uppercase)
    remaining_chains = [c for c in uc if c not in chainids][::-1]

    for mol in molecules:
        if not mol.chainid and mol.segid:
            mol.lines = list(pdb_segxchain(mol.lines))

        elif mol.chainid and not mol.segid:
            mol.lines = list(pdb_chainxseg(mol.lines))

        elif not mol.chainid and not mol.segid:
            _chains = pdb_chain.run(mol.lines, remaining_chains.pop())
            mol.lines = list(pdb_chainxseg.run(_chains))

        elif mol.chainid != mol.segid:
            raise ValueError('ChainID differs from segID')
