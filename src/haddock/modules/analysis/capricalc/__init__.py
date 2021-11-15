"""Calculate CAPRI metrics."""
import os
import shutil
import tempfile
from pathlib import Path

from Bio import PDB
from pdbtools import pdb_segxchain

from haddock import bin_path, log
from haddock.gear.config_reader import read_config
from haddock.libs.libontology import Format, ModuleIO
from haddock.libs.libparallel import Scheduler
from haddock.libs.libsubprocess import Job
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")


def calc_irmsd(reference, pdb_list, contact_exec, cutoff=10.0):
    """Calculate the Interface RMSD."""
    irmsd_dic = {}
    interface_list = identify_interface(reference,
                                        cutoff,
                                        contact_exec)

    bio_parser = PDB.PDBParser(QUIET=True)
    ref_structure = bio_parser.get_structure("reference", reference)

    ref_atoms = []
    for model in ref_structure:
        for chain in model:
            for res in chain:
                res_id = f'{chain.id}.{res.get_id()[1]}'
                if res_id in interface_list:
                    for atom in res:
                        if atom.get_id() == 'CA':
                            ref_atoms.append(atom)

    for pdb in pdb_list:
        pdb = Path(pdb.path, pdb.file_name)
        # make sure the model has chainID
        pdb = add_chain_from_segid(pdb)

        structure = bio_parser.get_structure("model", pdb)

        model_atoms = []
        for model in structure:
            for chain in model:
                for res in chain:
                    res_id = f'{chain.id}.{res.get_id()[1]}'
                    if res_id in interface_list:
                        for atom in res:
                            if atom.get_id() == 'CA':
                                model_atoms.append(atom)

        super_imposer = PDB.Superimposer()
        super_imposer.set_atoms(ref_atoms, model_atoms)
        super_imposer.apply(structure.get_atoms())

        irmsd = super_imposer.rms
        irmsd_dic[pdb] = irmsd
    return irmsd_dic


# def calc_fnat():
#     pass


# def calc_lrmsd():
#     pass


# def calc_ilrmsd():
#     pass


def identify_interface(pdb_f, cutoff, contact_exec, np=1):
    """Calculate the interface of a complex."""
    log.info('Calculating Interface')
    temp_f = tempfile.NamedTemporaryFile(delete=False, mode='w+t')
    interface_job = Job(pdb_f, temp_f.name, contact_exec, cutoff)
    interface_engine = Scheduler([interface_job],
                                 ncores=np)
    interface_engine.run()
    interface_reslist = []

    with open(interface_job.output, 'r') as fh:
        for line in fh.readlines():
            data = line.split()
            res_i = int(data[0])
            chain_i = data[1]
            # atom_i = data[2]
            res_j = int(data[3])
            chain_j = data[4]
            # atom_j = data[5]
            # distance = float(data[6])
            int_i = f'{chain_i}.{res_i}'
            int_j = f'{chain_j}.{res_j}'
            if int_i not in interface_reslist:
                interface_reslist.append(f'{chain_i}.{res_i}')
            if int_j not in interface_reslist:
                interface_reslist.append(f'{chain_j}.{res_j}')
    os.unlink(temp_f.name)
    return interface_reslist


def add_chain_from_segid(pdb_path):
    """Replace the chainID with the segID."""
    temp_f = tempfile.NamedTemporaryFile(delete=False,
                                         mode='w+t')
    with open(pdb_path) as fh:
        for line in list(pdb_segxchain.run(fh)):
            temp_f.writelines(line)
    temp_f.close()
    # REPLACE!
    new_pdb_path = shutil.move(temp_f.name, pdb_path)
    return new_pdb_path


class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module to calculate the CAPRI metrics."""

    def __init__(
            self,
            order,
            path,
            *ignore,
            init_params=DEFAULT_CONFIG,
            **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm if contact executable is compiled."""
        dcfg = read_config(DEFAULT_CONFIG)
        exec_path = Path(bin_path, dcfg['contact_exec'])

        if not os.access(exec_path, mode=os.F_OK):
            raise Exception(f'Required {str(exec_path)} file does not exist.')

        if not os.access(exec_path, mode=os.X_OK):
            raise Exception(f'Required {str(exec_path)} file is not executable')

        return

    def run(self, **params):
        """Execute module."""
        log.info("Running [capricalc] module")

        super().run(params)

        # This is not the same as FCC uses (:
        contact_exec = Path(bin_path, self.params['contact_exec'])

        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            self.finish_with_error('This module cannot come after one'
                                   ' that produced an iterable')

        models_to_calc = [
            p
            for p in self.previous_io.output
            if p.file_type == Format.PDB
            ]

        reference = self.params['reference']
        if not reference:
            # This was not defined, we cannot move on
            _msg = 'reference structure not defined'
            self.finish_with_error(_msg)

        irmsd_dic = calc_irmsd(reference,
                               models_to_calc,
                               contact_exec,
                               cutoff=self.params['irmsd_cutoff'])

        # Temporary
        for k in irmsd_dic:
            v = irmsd_dic[k]
            print(k, v)

        selected_models = models_to_calc
        io = ModuleIO()
        io.add(selected_models, "o")
        io.save(self.path)
