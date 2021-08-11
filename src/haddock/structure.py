"""Molecular data structures"""
from pathlib import Path


class Molecule:
    """Input molecule, usually a PDB file."""
    def __init__(self, file_name, segid=None):
        # the rest of the code is too dependent on the Path API
        assert isinstance(file_name, Path), "`file_name` must be pathlib.Path"

        self.file_name = file_name
        self.segid = segid


def make_molecules(paths):
    """Get input molecules from the data stream."""
    return list(map(Molecule, paths))
