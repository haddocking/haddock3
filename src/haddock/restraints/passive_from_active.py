"""haddock3-restraints passive_from_active subcommand.

Given a list of active_residues and a PDB structure, it will return a list of
surface exposed passive residues within a 6.5A radius from the active residues.

When provided with a list of surface residues, it will filter the list for those
that are within 6.5A from the active residues.

Usage:
    haddock3-restraints passive_from_active <pdb_file> <active_list> [-c <chain_id>] [-s <surface_list>]
"""
from pathlib import Path

import numpy as np
import sys
from Bio.PDB import PDBParser, NeighborSearch


def add_pass_from_act_arguments(pass_from_act_subcommand):
    """Add arguments to the pass_from_act subcommand."""
    pass_from_act_subcommand.add_argument(
        "structure",
        type=str,
        help="input PDB structure.",
        )

    pass_from_act_subcommand.add_argument(
        "active_list",
        help="List of active residues IDs (int) separated by commas",
        type=str,
        )

    pass_from_act_subcommand.add_argument(
        "-c",
        "--chain-id",
        help="Chain id to be used in the PDB file (default: All)",
        required=False,
        )

    pass_from_act_subcommand.add_argument(
        "-s",
        "--surface-list",
        help="List of surface residues IDs (int) separated by commas",
        required=False,
        type=str,
        )

    return pass_from_act_subcommand

# Scaling factors for relative ASA
# Calculated using extended ALA-X-ALA peptides
# Taken from NACCESS
rel_asa = {
    'total':
        {
            'ALA': 107.95,
            'CYS': 134.28,
            'ASP': 140.39,
            'GLU': 172.25,
            'PHE': 199.48,
            'GLY': 80.10,
            'HIS': 182.88,
            'ILE': 175.12,
            'LYS': 200.81,
            'LEU': 178.63,
            'MET': 194.15,
            'ASN': 143.94,
            'PRO': 136.13,
            'GLN': 178.50,
            'ARG': 238.76,
            'SER': 116.50,
            'THR': 139.27,
            'VAL': 151.44,
            'TRP': 249.36,
            'TYR': 212.76,
        },
    'bb':
        {
            'ALA': 38.54,
            'CYS': 37.53,
            'ASP': 37.70,
            'GLU': 37.51,
            'PHE': 35.37,
            'GLY': 47.77,
            'HIS': 35.80,
            'ILE': 37.16,
            'LYS': 37.51,
            'LEU': 37.51,
            'MET': 37.51,
            'ASN': 37.70,
            'PRO': 16.23,
            'GLN': 37.51,
            'ARG': 37.51,
            'SER': 38.40,
            'THR': 37.57,
            'VAL': 37.16,
            'TRP': 38.10,
            'TYR': 35.38,
        },
    'sc':
        {
            'ALA': 69.41,
            'CYS': 96.75,
            'ASP': 102.69,
            'GLU': 134.74,
            'PHE': 164.11,
            'GLY': 32.33,
            'HIS': 147.08,
            'ILE': 137.96,
            'LYS': 163.30,
            'LEU': 141.12,
            'MET': 156.64,
            'ASN': 106.24,
            'PRO': 119.90,
            'GLN': 140.99,
            'ARG': 201.25,
            'SER': 78.11,
            'THR': 101.70,
            'VAL': 114.28,
            'TRP': 211.26,
            'TYR': 177.38,
        }
}

def get_surface_resids(structure, cutoff=15):
    """
    Calls freesasa using its Python API and returns
    per-residue accessibilities.
    """
    try:
        from freesasa import Classifier, structureFromBioPDB, calc
    except ImportError as err:
        print('[!] The binding affinity prediction tools require the \'freesasa\' Python API', file=sys.stderr)
        raise ImportError(err)

    asa_data, rsa_data, rel_main_chain, rel_side_chain = {}, {}, {}, {}
    _rsa = rel_asa['total']
    _rsa_bb = rel_asa['bb']
    _rsa_sc = rel_asa['sc']

    #classifier = Classifier(config_path)
    classifier = Classifier()

    struct = structureFromBioPDB(structure, classifier, )
    result = calc(struct)

    # iterate over all atoms to get SASA and residue name
    for idx in range(struct.nAtoms()):
        atname = struct.atomName(idx).strip()
        resname = struct.residueName(idx)
        resid = int(struct.residueNumber(idx))
        chain = struct.chainLabel(idx)
        at_uid = (chain, resname, resid, atname)
        res_uid = (chain, resname, resid)

        asa = result.atomArea(idx)
        asa_data[at_uid] = asa
        # add asa to residue
        rsa_data[res_uid] = rsa_data.get(res_uid, 0) + asa

        if atname in ('C', 'N', 'O'):
            rel_main_chain[res_uid] = rel_main_chain.get(res_uid, 0) + asa
        else:
            rel_side_chain[res_uid] = rel_side_chain.get(res_uid, 0) + asa

    # convert total asa ro relative asa
    rsa_data.update((res_uid, asa / _rsa[res_uid[1]]) for res_uid, asa in rsa_data.items())
    rel_main_chain.update((res_uid, asa / _rsa_bb[res_uid[1]] * 100) for res_uid, asa in rel_main_chain.items())
    rel_side_chain.update((res_uid, asa / _rsa_sc[res_uid[1]] * 100) for res_uid, asa in rel_side_chain.items())

    # We format to fit the pipeline
    resid_access = {}
    for res_uid, access in rel_main_chain.items():
        resid_access[res_uid[2]] = {'side_chain_rel': rel_side_chain.get(res_uid), 'main_chain_rel': access}
    surface_resids = [r for r, v in resid_access.items() if v['side_chain_rel'] >= cutoff or
                      v['main_chain_rel'] >= cutoff]
    return surface_resids


def passive_from_active(structure, active_list, chain_id=None, surface_list=""):
    """Get the passive residues."""

    # Parse the PDB file
    if Path(structure).exists():
        try:
            p = PDBParser(QUIET=True)
            s = p.get_structure('pdb', structure)
        except Exception as e:
            print('Error while parsing the PDB file: {0}'.format(e))
            sys.exit(1)
    else:
        print('File not found: {0}'.format(structure))
        sys.exit(1)

    try:
        if chain_id:
            atom_list = [a for a in s[0][chain_id].get_atoms()]
        else:
            atom_list = [a for a in s[0].get_atoms()]
    except KeyError as e:
        print('Chain {0} does not exist in the PDB file {1}, please enter a proper chain id'.
              format(chain_id, structure))
        sys.exit(1)

    try:
        active_list = [int(res) for res in active_list.split(',')]
        act_atoms = [a.get_coord() for a in atom_list if a.parent.id[1] in active_list]
    except:
        print('The list of active residues must be provided as a comma-separated list of integers')
        sys.exit(1)

    try:
        if surface_list:
            surface_list = [int(res) for res in surface_list.split(',')]
        else:
            surface_list = get_surface_resids(s)
    except Exception as e:
        print("There was an error while calculating surface residues: {}".format(e))
        sys.exit(1)

    ns = NeighborSearch(atom_list)
    neighbors = []
    for a in act_atoms:
        neighbors.append(ns.search(a, 6.5, "R"))  # HADDOCK used 6.5A as default

    passive_list = set()
    for n in neighbors:
        for r in n:
            passive_list.add(r.id[1])
    tmp = passive_list & set(surface_list)
    passive_list = tmp - set(active_list)
    print(' '.join([str(r) for r in sorted(passive_list)]))
    return 