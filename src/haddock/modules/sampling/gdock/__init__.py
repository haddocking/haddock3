"""gdock integration sampling module."""
import os
import re
import shlex
import subprocess
import sys
from pathlib import Path

from pdbtools import pdb_tidy

from haddock import log
from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.typing import FilePath
from haddock.libs import libpdb
from haddock.libs.libontology import PDBFile
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, MODULE_DEFAULT_YAML)


def ambig2dic(ambig_f: FilePath) -> dict[str, list[int]]:
    """Read an ambig.tbl file and convert it to a dictionary."""
    ambig_regex = r"resid\s*(\d*)\s*and\s*segid\s*(\w)"
    ambig_dic: dict[str, list[int]] = {}
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
    """HADDOCK3 gdock module."""

    name = RECIPE_PATH.name

    def __init__(self,
                 order: int,
                 path: Path,
                 initial_params: FilePath = DEFAULT_CONFIG) -> None:
        super().__init__(order, path, initial_params)

    @classmethod
    def confirm_installation(cls) -> None:
        """Confirm this module is installed."""
        gdock_path = os.environ['GDOCK_PATH']
        gdock_exec = Path(gdock_path, 'gdock.py')
        cmd = f'{sys.executable} {gdock_exec}'
        p = subprocess.run(shlex.split(cmd), capture_output=True)
        # out = p.stdout.decode('utf-8')
        err = p.stderr.decode('utf-8')
        if "error: the following arguments are required: input_file" in err:
            # all good :)
            pass
        else:
            raise Exception('gdock is not installed properly')

    def _run(self) -> None:
        """Execute module."""
        try:
            gdock_path = os.environ['GDOCK_PATH']
        except KeyError:
            self.finish_with_error('GDOCK_PATH not defined')

        gdock_exec = Path(gdock_path, 'gdock.py')
        if not gdock_exec.exists():
            self.finish_with_error(f'{gdock_exec} not found')

        # Get the models generated in previous step
        models_to_dock: list[PDBFile] = self.previous_io.retrieve_models()[0]
        topologies = [e.topology for e in models_to_dock]

        input_a = models_to_dock[0].rel_path
        input_b = models_to_dock[1].rel_path

        input = {'A': input_a, 'B': input_b}
        # Check if chain IDs are present
        for chain in input:
            pdb = input[chain]
            chain_pdb = Path(self.path, pdb.name)
            segids, chains = libpdb.identify_chainseg(pdb)
            if set(segids) != set(chains):
                log.info("No chain IDs found, using segid information")
                libpdb.swap_segid_chain(pdb, chain_pdb)
            if chain_pdb.exists():
                input[chain] = chain_pdb

        # convert ambig to list
        ambig_dic = ambig2dic(self.params['ambig_fname'])

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

        # ===============
        # Placeholder, this is not yet implemented in gdock
        # input_toml += '[parameters]' + os.linesep
        # input_toml += 'population_size = 100' + os.linesep
        # input_toml += 'max_number_of_generations = 50' + os.linesep
        # ===============

        with open('input.toml', 'w') as inp_fh:
            inp_fh.write(input_toml)

        cmd = f'{sys.executable} {gdock_exec} --dry input.toml'

        subprocess.call(cmd, shell=True)

        # retrieve the structures
        self.output_models: list[PDBFile] = []
        structure_folder = Path('gdock-integration/structures')
        for model in structure_folder.glob('*pdb'):
            # Make sure the output is tidy, this should be handled
            #  by gdock, but check it here to be sure
            with open(model, 'r') as fin:
                lines = list(
                    pdb_tidy.run(fin, strict=False)
                    )  # be explicit in the `strict`

            with open(model, 'w') as fout:
                fout.write(''.join(lines))

            pdb = PDBFile(model)
            pdb.score = .0
            pdb.topology = topologies
            self.output_models.append(pdb)

        self.export_io_models()
