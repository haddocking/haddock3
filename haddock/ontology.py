"""This module describes the Haddock3 ontology used for communicating between modules"""

import json
from os import linesep
from enum import Enum
from pathlib import Path
from haddock.defaults import MODULE_IO_FILE


class Format(Enum):
    """Input and Output possible formats"""
    PDB = "pdb"
    PDB_ENSEMBLE = "pdb"
    CNS_INPUT = "inp"
    CNS_OUTPUT = "out"
    TOPOLOGY = "psf"

    def __str__(self):
        return str(self.value)


class ModuleIO:
    """This class represents an object for intercommunicating modules and exchange
    input/output information"""
    def __init__(self):
        self.input = []
        self.output = []

    def add(self, filename, format=Format.PDB, mode='i'):
        """Add a given filename as input or output"""
        f = Path(filename).resolve().name
        if mode == 'i':
            self.input.append((f, format.name))
        else:
            self.output.append((f, format.name))

    def save(self, path, filename=MODULE_IO_FILE):
        """Save Input/Output needed files by this module to disk"""
        with open(path / filename, 'w') as output_handler:
            to_save = {'input': self.input, 'output': self.output}
            json.dump(to_save, output_handler, indent=4)
        return path / filename

    def load(self, filename):
        """Load the content of a given IO filename"""
        with open(filename) as json_file:
            content = json.load(json_file)
            self.input = content['input']
            for pair in self.input:
                pair[1] = Format[pair[1]]
            self.output = content['output']
            for pair in self.output:
                pair[1] = Format[pair[1]]

    def __repr__(self):
        return f'Input: {self.input}{linesep}Output: {self.output}'
