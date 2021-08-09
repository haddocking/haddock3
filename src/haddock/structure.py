"""Molecular data structures"""


class Molecule:
    """Input molecule, usually a PDB file."""
    def __init__(self, file_name, segid=None):
        self.file_name = file_name
        self.segid = segid
