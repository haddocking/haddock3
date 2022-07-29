"""PDB related functions."""
from fccpy import read_pdb
from fccpy.contacts import get_intermolecular_contacts

from haddock.modules.design.optint.modules.residues import AA_TO_CODE


def pdb_to_dict(pdb_f):
    """Read the PDB into a dictionary."""
    coord_dic = {}
    with open(pdb_f) as fh:
        for line in fh:
            if line.startswith("ATOM"):
                chain = line[21]
                if chain not in coord_dic:
                    coord_dic[chain] = {}
                resnum = int(line[22:26])
                resname_three_letter = line[17:20]
                if resnum not in coord_dic[chain]:
                    coord_dic[chain][resnum] = AA_TO_CODE[resname_three_letter]
    return coord_dic


def get_interface_dict(pdb_f, cutoff=5.0):
    """Load the PDB file and return a dictionary with interface information."""
    coord_dic = {}
    interface_dict = {}
    with open(pdb_f) as fh:
        for line in fh:
            if line.startswith("ATOM"):
                chain = line[21]
                if chain not in coord_dic:
                    coord_dic[chain] = {}
                x = float(line[31:38])
                y = float(line[39:46])
                z = float(line[47:54])
                resnum = int(line[22:26])
                resname_three_letter = line[17:20]
                # atom_name = line[12:16].strip()
                if resnum not in coord_dic[chain]:
                    coord_dic[chain][resnum] = []
                coord_dic[chain][resnum].append((x, y, z))

                if chain not in interface_dict:
                    interface_dict[chain] = {}

                interface_dict[chain][resnum] = AA_TO_CODE[resname_three_letter]

    pdb = read_pdb(pdb_f)
    _resdic = {}
    for atom_i, atom_j in get_intermolecular_contacts(pdb, cutoff):

        if atom_i.chain not in _resdic:
            _resdic[atom_i.chain] = []
        if atom_j.chain not in _resdic:
            _resdic[atom_j.chain] = []

        if atom_i.resid not in _resdic[atom_i.chain]:
            _resdic[atom_i.chain].append(atom_i.resid)
        if atom_j.resid not in _resdic[atom_j.chain]:
            _resdic[atom_j.chain].append(atom_j.resid)

    for chain in interface_dict:
        for resnum in interface_dict[chain]:
            if resnum not in _resdic[chain]:
                # This residue is not in the interface
                #  set it to None
                interface_dict[chain][resnum] = None

    return interface_dict
