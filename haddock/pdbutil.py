import os
from pathlib import Path
from pdbtools.pdb_splitmodel import split_model
from haddock.data.topology import Topology
from haddock.modules import working_directory


class PDBFactory:
    """A factory class to deal with PDB logic"""

    __to_remove = ['REMAR', 'CTERB', 'CTERA', 'NTERA', 'NTERB', 'CONECT']
    __to_rename = {'HSD': 'HIS', 'HSE': 'HIS', 'HID': 'HIS', 'HIE': 'HIS',
                   'WAT ': 'TIP3', ' 0.00969': ' 0.00   '}
    __to_keep = Topology.get_supported()

    @staticmethod
    def split_ensemble(pdb_file_path):
        """"Split a PDB file in multiple structures if different models are found"""
        new_models = []
        abs_path = Path(pdb_file_path).resolve().parent.absolute()
        with open(pdb_file_path) as input_handler:
            with working_directory(abs_path):
                split_model(input_handler)

        basename = Path(pdb_file_path)
        new_models = list(abs_path.glob(f'{basename.stem}_*{basename.suffix}'))
        return new_models

    @staticmethod
    def sanitize(pdb_file_path, overwrite=True):
        """Sanitize a PDB file"""
        good_lines = []
        with open(pdb_file_path) as input_handler:
            for line in input_handler:
                line = line.rstrip(os.linesep)
                # Ignoring lines containing any tag from __to_remove
                if not any([tag in line for tag in PDBFactory.__to_remove]):
                    for tag, new_tag in PDBFactory.__to_rename.items():
                        line = line.replace(tag, new_tag)
                    # check if this residue is known
                    res = line[17:20].strip()
                    if res and res in PDBFactory.__to_keep:
                        good_lines.append(line)
            if len(good_lines) and good_lines[-1] != 'END':
                good_lines.append('END')

        if overwrite:
            with open(pdb_file_path, 'w') as output_handler:
                for line in good_lines:
                    output_handler.write(line + os.linesep)
            return pdb_file_path
        else:
            basename = Path(pdb_file_path)
            abs_path = Path(pdb_file_path).resolve().parent.absolute()
            new_pdb_file = abs_path / f'{basename.stem}_cleaned{basename.suffix}'
            with open(new_pdb_file, 'w') as output_handler:
                for line in good_lines:
                    output_handler.write(line + os.linesep)
            return new_pdb_file

    @staticmethod
    def identify_chainseg(pdb_file_path):
        """"Return segID OR chainID"""
        segid_l = []
        with open(pdb_file_path) as input_handler:
            for line in input_handler:
                if line.startswith('ATOM  '):
                    try:
                        segid = line[72:76].strip()[:1]
                    except IndexError:
                        segid = ''
                    try:
                        chainid = line[21].strip()[:1]
                    except IndexError:
                        chainid = ''

                    if segid:
                        segid_l.append(segid)
                    elif chainid:
                        segid_l.append(segid)

        segid_l = list(set(segid_l))
        segid_l.sort()
        return segid_l
