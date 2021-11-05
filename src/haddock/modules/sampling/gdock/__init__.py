"""HADDOCK3 gdock integration module"""
import logging
import os
import re
import subprocess
import sys

from pathlib import Path

from haddock.libs import libpdb
from haddock.libs.libontology import Format, ModuleIO, PDBFile
from haddock.libs.libutil import check_subprocess
from haddock.modules import BaseHaddockModule, working_directory


logger = logging.getLogger(__name__)


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")




def ambig2dic(ambig_f):
    """Read an ambig.tbl file and convert it to a dictionary"""
    ambig_regex = r"resid\s*(\d*)\s*and\s*segid\s*(\w)"
    ambig_dic = {}
    with open(ambig_f) as fh:
        for line in fh.readlines():
            matches = re.finditer(ambig_regex, line)
            for m in matches:
                resid = int(m.group(1))
                chain = m.group(2)
                if chain not in ambig_dic:
                    ambig_dic[chain] = []

                ambig_dic[chain].append(resid)
    return ambig_dic


class HaddockModule(BaseHaddockModule):

    def __init__(self, order, path, initial_params=DEFAULT_CONFIG):
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm this module is installed."""
        gdock_path = os.environ['GDOCK_PATH']
        gdock_exec = Path(gdock_path, 'gdock.py')
        check_subprocess(f'{sys.executable} {gdock_exec}')

    def run(self, **params):
        logger.info("Running [gdock] module")

        super().run(params)

        try:
            gdock_path = os.environ['GDOCK_PATH']
        except KeyError:
            self.finish_with_error('GDOCK_PATH not defined')

        gdock_exec = Path(gdock_path, 'gdock.py')
        if not gdock_exec.exists():
            self.finish_with_error(f'{gdock_exec} not found')

        # Get the models generated in previous step
        models_to_dock = [p for p in self.previous_io.output if p.file_type == Format.PDB]

        if '00_topoaa' not in Path(models_to_dock[0].path).stem:
            _msg = 'This module must come after Topology generation'
            self.finish_with_error(_msg)

        topologies = [p for p in self.previous_io.output if p.file_type == Format.TOPOLOGY]

        input_a = Path(models_to_dock[0].path, models_to_dock[0].file_name)
        input_b = Path(models_to_dock[1].path, models_to_dock[1].file_name)

        input = {'A': input_a, 'B': input_b}
        # Check if chain IDs are present
        for chain in input:
            pdb = input[chain]
            chain_pdb = Path(self.path, pdb.name)
            segids, chains = libpdb.identify_chainseg(pdb)
            if set(segids) != set(chains):
                logger.info("No chain IDs found, using segid information")
                libpdb.swap_segid_chain(pdb, chain_pdb)
            input[chain] = chain_pdb

        # convert ambig to list
        ambig_dic = ambig2dic(self.params.get('ambig', None))

        input_toml = '' + os.linesep
        input_toml += '[main]' + os.linesep
        input_toml += 'identifier = "gdock-integration"' + os.linesep

        # this is needed because 'ncores' is defined in BaseHaddockModule
        # by default as None
        ncores = self.params['ncores'] or 1
        input_toml += f'number_of_processors = {ncores}' + os.linesep

        input_toml += '[restraints]' + os.linesep

        for chain in ambig_dic:
            reslist = list(set(ambig_dic[chain]))
            input_toml += f'{chain} = {reslist}' + os.linesep

        input_toml += '[molecules]' + os.linesep
        input_toml += f'A = \"{input["A"]}\"' + os.linesep
        input_toml += f'B = \"{input["B"]}\"' + os.linesep
        input_toml += os.linesep

        with working_directory(self.path):
            with open('input.toml', 'w') as inp_fh:
                inp_fh.write(input_toml)

            cmd = f'{sys.executable} {gdock_exec} --dry input.toml'

            subprocess.call(cmd, shell=True)

        # retrieve the structures
        output_structures = []
        structure_folder = Path(self.path, 'gdock-integration/structures')
        for model in structure_folder.glob('*pdb'):
            pdb = PDBFile(model, path=model.parent)
            pdb.score = .0
            pdb.topology = topologies
            output_structures.append(pdb)

        io = ModuleIO()
        io.add(output_structures, "o")
        io.save(self.path)
