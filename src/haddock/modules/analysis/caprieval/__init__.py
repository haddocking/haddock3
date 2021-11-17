"""Calculate CAPRI metrics."""
import copy
import os
import shutil
import tempfile
from pathlib import Path

import numpy as np
from fccpy import read_pdb
from fccpy.contacts import get_intermolecular_contacts
from pdbtools import pdb_segxchain

from haddock import bin_path, log
from haddock.gear.config_reader import read_config
from haddock.libs.libontology import Format, ModuleIO
from haddock.modules import BaseHaddockModule


RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.cfg")

# FIXME: Add this as a parameter
BACKBONE_ATOMS = ['C', 'CA', 'N', 'O']


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


def centroid(X):
    """Get the centroid."""
    return X.mean(axis=0)


def kabsch(P, Q):
    """Find the rotation matrix using Kabsch algorithm."""
    # Covariance matrix
    C = np.dot(np.transpose(P), Q)
    # use SVD
    V, S, W = np.linalg.svd(C)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]
    # Create Rotation matrix U
    U = np.dot(V, W)
    return U


def calc_rmsd(V, W):
    """Calculate the RMSD from two vectors."""
    diff = np.array(V) - np.array(W)
    N = len(V)
    return np.sqrt((diff * diff).sum() / N)


def read_res(pdb_f):
    """Read residue numbers in a PDB file."""
    res_dic = {}
    with open(pdb_f, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('ATOM'):
                chain = line[21]
                resnum = int(line[22:26])
                if chain not in res_dic:
                    res_dic[chain] = []
                if resnum not in res_dic[chain]:
                    res_dic[chain].append(resnum)
    return res_dic


# Debug only
def write_coords(output_name, coor_list):
    """Add a dummy atom to a PDB file according to a list of coordinates."""
    with open(output_name, 'w') as fh:
        for i, dummy_coord in enumerate(coor_list):
            atom_num = f'{i}'.rjust(4, ' ')
            resnum = f'{i}'.rjust(3, ' ')
            dum_x = f'{dummy_coord[0]:.3f}'.rjust(7, ' ')
            dum_y = f'{dummy_coord[1]:.3f}'.rjust(7, ' ')
            dum_z = f'{dummy_coord[2]:.3f}'.rjust(7, ' ')
            dummy_line = (f'ATOM   {atom_num}  H   DUM X {resnum}   '
                          f'  {dum_x} {dum_y} {dum_z}  1.00  1.00   '
                          '        H  ' + os.linesep)
            fh.write(dummy_line)


def load_contacts(pdb_f, cutoff):
    """Load residue-based contacts."""
    con_list = []
    structure = read_pdb(pdb_f)
    for atom_i, atom_j in get_intermolecular_contacts(structure, cutoff):
        con = (atom_i.chain, atom_i.resid, atom_j.chain, atom_j.resid)
        con_list.append(con)
    return set(con_list)


def load_coords(pdb_f, filter_resdic=None, atoms=None):
    """Load coordinates from PDB."""
    C = []
    chain_dic = {}
    idx = 0
    with open(pdb_f, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('ATOM'):

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                resnum = int(line[22:26])
                chain = line[21]

                if chain not in chain_dic:
                    chain_dic[chain] = []

                atom_name = line[12:16].strip()
                if atoms:
                    if atom_name not in atoms:
                        continue

                if filter_resdic:

                    if chain in filter_resdic:
                        if resnum in filter_resdic[chain]:
                            # OK
                            C.append(np.asarray([x, y, z], dtype=float))
                            chain_dic[chain].append(idx)
                            idx += 1
                else:
                    # append anyway
                    C.append(np.asarray([x, y, z], dtype=float))
                    chain_dic[chain].append(idx)
                    idx += 1

    chain_ranges = {}
    for chain in chain_dic:
        min_idx = min(chain_dic[chain])
        max_idx = max(chain_dic[chain])
        chain_ranges[chain] = (min_idx, max_idx)

    return np.asarray(C), chain_ranges


def identify_interface(pdb_f, cutoff):
    """Identify the interface."""
    pdb = read_pdb(pdb_f)
    interface_resdic = {}
    for atom_i, atom_j in get_intermolecular_contacts(pdb, cutoff):
        if atom_i.chain not in interface_resdic:
            interface_resdic[atom_i.chain] = []
        if atom_j.chain not in interface_resdic:
            interface_resdic[atom_j.chain] = []

        if atom_i.resid not in interface_resdic[atom_i.chain]:
            interface_resdic[atom_i.chain].append(atom_i.resid)
        if atom_j.resid not in interface_resdic[atom_j.chain]:
            interface_resdic[atom_j.chain].append(atom_j.resid)
    return interface_resdic


class CAPRI:
    """CAPRI class."""

    def __init__(self, reference, model_list, atoms):
        self.reference = reference
        self.model_list = []
        self.irmsd_dic = {}
        self.lrmsd_dic = {}
        self.ilrmsd_dic = {}
        self.fnat_dic = {}
        self.atoms = atoms

        for struct in model_list:
            pdb_f = Path(struct.path, struct.file_name)
            pdb_w_chain = add_chain_from_segid(pdb_f)
            self.model_list.append(pdb_w_chain)

    def irmsd(self, cutoff=10.):
        """Calculate the I-RMSD."""
        # Identify interface
        ref_interface_resdic = identify_interface(self.reference, cutoff)

        # Load interface coordinates
        Q, _ = load_coords(self.reference,
                           filter_resdic=ref_interface_resdic,
                           atoms=self.atoms)
        # Move to centroids
        Q = Q - centroid(Q)

        for model in self.model_list:

            P, _ = load_coords(model,
                               filter_resdic=ref_interface_resdic,
                               atoms=self.atoms)

            P = P - centroid(P)
            U = kabsch(P, Q)
            P = np.dot(P, U)
            i_rmsd = calc_rmsd(P, Q)
            # write_coords('ref.pdb', P)
            # write_coords('model.pdb', Q)
            self.irmsd_dic[model] = i_rmsd

        return self.irmsd_dic

    def lrmsd(self, receptor_chain, ligand_chain):
        """Calculate the L-RMSD."""
        ref_resdic = read_res(self.reference)

        ref_receptor_resdic = copy.deepcopy(ref_resdic)
        del ref_receptor_resdic[ligand_chain]

        # Get reference coordinates
        Q_all, chain_ranges = load_coords(self.reference,
                                          filter_resdic=ref_resdic,
                                          atoms=self.atoms)

        receptor_start = chain_ranges[receptor_chain][0]
        receptor_end = chain_ranges[receptor_chain][1]
        Q_receptor = Q_all[receptor_start:receptor_end]

        # loop goes here
        model = self.model_list[0]
        for model in self.model_list:

            P_all, _ = load_coords(model,
                                   filter_resdic=ref_resdic,
                                   atoms=self.atoms)

            receptor_start = chain_ranges[receptor_chain][0]
            receptor_end = chain_ranges[receptor_chain][1]
            P_receptor = P_all[receptor_start:receptor_end]

            # write_coords('ref_ori.pdb', Q_all)
            # write_coords('model_ori.pdb', P_all)

            # Center receptors and get rotation matrix
            Q_receptor_centroid = centroid(Q_receptor)
            Q_receptor -= Q_receptor_centroid
            P_receptor_centroid = centroid(P_receptor)
            P_receptor -= P_receptor_centroid
            U_receptor = kabsch(P_receptor, Q_receptor)

            # Center complexes in the receptor centroids
            P_all -= P_receptor_centroid
            Q_all -= Q_receptor_centroid

            # Apply rotation to complex
            #  - complex are aligned on the receptor
            P_all = np.dot(P_all, U_receptor)

            # write_coords('ref.pdb', Q_all)
            # write_coords('model.pdb', P_all)

            # Identify the ligand coordinates
            ligand_start = chain_ranges[ligand_chain][0]
            ligand_end = chain_ranges[ligand_chain][1]

            Q_ligand = Q_all[ligand_start:ligand_end]
            P_ligand = P_all[ligand_start:ligand_end]

            # write_coords('ref_ligand.pdb', Q_ligand)
            # write_coords('model_ligand.pdb', P_ligand)

            # Calculate the RMSD of the ligands
            l_rmsd = calc_rmsd(P_ligand, Q_ligand)

            self.lrmsd_dic[model] = l_rmsd

        return self.lrmsd_dic

    def ilrmsd(self, ligand_chain, cutoff):
        """Calculate the Interface Ligand RMSD."""
        ref_resdic = read_res(self.reference)
        # Identify interface
        ref_interface_resdic = identify_interface(self.reference, cutoff)

        # Load interface coordinates
        Q, chain_ranges = load_coords(self.reference,
                                      filter_resdic=ref_resdic,
                                      atoms=self.atoms)
        Q_int, _ = load_coords(self.reference,
                               filter_resdic=ref_interface_resdic,
                               atoms=self.atoms)
        # Move to centroids
        Q_int_centroid = centroid(Q_int)
        Q_int = Q_int - Q_int_centroid

        for model in self.model_list:

            Q_all = copy.deepcopy(Q)

            P_all, _ = load_coords(model,
                                   filter_resdic=ref_resdic,
                                   atoms=self.atoms)

            P_int, _ = load_coords(model,
                                   filter_resdic=ref_interface_resdic,
                                   atoms=self.atoms)

            P_int_centroid = centroid(P_int)
            P_int = P_int - P_int_centroid

            # find the rotation that minimizes the interface rmsd
            U_int = kabsch(P_int, Q_int)

            P_all -= P_int_centroid
            Q_all -= Q_int_centroid

            # apply this rotation to the model
            P_all = np.dot(P_all, U_int)

            # Calculate the rmsd of the ligand
            ligand_start = chain_ranges[ligand_chain][0]
            ligand_end = chain_ranges[ligand_chain][1]

            Q_ligand = Q_all[ligand_start:ligand_end]
            P_ligand = P_all[ligand_start:ligand_end]

            # write_coords('ref.pdb', P_ligand)
            # write_coords('model.pdb', Q_ligand)

            # this will be the interface-ligand-rmsd
            i_l_rmsd = calc_rmsd(P_ligand, Q_ligand)
            self.ilrmsd_dic[model] = i_l_rmsd

        return self.ilrmsd_dic

    def fnat(self, cutoff=5.0):
        """Calculate the frequency of native contacts."""
        ref_contacts = load_contacts(self.reference, cutoff)
        for model in self.model_list:
            model_contacts = load_contacts(model, cutoff)
            intersection = ref_contacts & model_contacts
            fnat = len(intersection) / float(len(ref_contacts))
            self.fnat_dic[model] = fnat
        return self.fnat_dic

    def output(self, output_f):
        """Output the CAPRI results to a .tsv file."""
        sep = '\t'
        with open(output_f, 'w') as fh:
            header = 'model' + sep
            if self.fnat_dic:
                header += 'fnat' + sep
            if self.irmsd_dic:
                header += 'irmsd' + sep
            if self.lrmsd_dic:
                header += 'lrmsd' + sep
            if self.ilrmsd_dic:
                header += 'ilrmsd' + sep

            header += os.linesep
            fh.write(header)

            for model in self.model_list:
                row = f'{model.name}' + sep
                if model in self.fnat_dic:
                    row += f'{self.fnat_dic[model]:.3f}' + sep
                if model in self.irmsd_dic:
                    row += f'{self.irmsd_dic[model]:.2f}' + sep
                if model in self.lrmsd_dic:
                    row += f'{self.lrmsd_dic[model]:.2f}' + sep
                if model in self.ilrmsd_dic:
                    row += f'{self.ilrmsd_dic[model]:.2f}' + sep
                row += os.linesep
                fh.write(row)


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
        log.info("Running [caprieval] module")
        log.info('Calculating CAPRI metrics...')

        super().run(params)

        # Get the models generated in previous step
        if type(self.previous_io) == iter:
            self.finish_with_error('This module cannot come after one'
                                   ' that produced an iterable')

        models_to_calc = [
            p
            for p in self.previous_io.output
            if p.file_type == Format.PDB
            ]

        reference = Path(self.params['reference'])
        if not reference:
            # This was not defined, we cannot move on
            _msg = 'reference structure not defined'
            self.finish_with_error(_msg)

        log.info(f'Using {reference} as reference structure')

        capri = CAPRI(reference,
                      models_to_calc,
                      atoms=self.params['atoms'])

        if self.params['fnat']:
            log.info(' Calculating FNAT')
            capri.fnat(cutoff=self.params['fnat_cutoff'])

        if self.params['irmsd']:
            log.info(' Calculating I-RMSD')
            capri.irmsd(cutoff=self.params['irmsd_cutoff'])

        if self.params['lrmsd']:
            log.info(' Calculating L-RMSD')
            capri.lrmsd(receptor_chain=self.params['receptor_chain'],
                        ligand_chain=self.params['ligand_chain'])

        if self.params['ilrmsd']:
            log.info(' Calculating I-L-RMSD')
            capri.ilrmsd(ligand_chain=self.params['ligand_chain'],
                         cutoff=self.params['irmsd_cutoff'])

        output_fname = Path(self.path, 'capri.tsv')
        log.info(f'Saving output to {output_fname}')
        capri.output(output_fname)

        selected_models = models_to_calc
        io = ModuleIO()
        io.add(selected_models, "o")
        io.save(self.path)
