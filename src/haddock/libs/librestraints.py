"""restraints-related functions."""

import numpy as np
from haddock import log
from haddock.libs.libalign import get_atoms, load_coords


def get_unique_resids(models):
    """
    Parameters
    ----------
    models : list of models

    Returns
    -------
    chains : dict
        dictionary of all the unique resids for each chain
    """
    chains = {}
    for mod in models:
        # getting atoms and coords
        atoms = get_atoms(mod)
        coords, coords_range = load_coords(mod, atoms)
        # getting unique chains
        chain_keys = [atom[0] for atom in coords.keys()]
        unique_chain_keys = set(chain_keys)
        # getting unique residues
        for ch in unique_chain_keys:
            res_keys = [atom[1] for atom in coords.keys() if atom[0] == ch]
            unique_res_keys = set(res_keys)
            if ch not in chains.keys():
                chains[ch] = list(unique_res_keys)
            else:
                chains[ch].extend(list(unique_res_keys))
    return chains

def validate_ambig_fname(ambig_fname, models):
    """
    validate ambig_fname against input models.

    At first the function extracts the restraints from ambig_fname and the
    unique residue ids for each chain from the models. Then it checks for the
    existence of each restraint in such dictionary.
    
    Parameters
    ----------
    ambig_fname : str or Path
        restraints filename
    models : list
        list of models
    """
    log.debug(f"validating {ambig_fname}")
    restraints = parse_tbl(ambig_fname)
    resid_dict = get_unique_resids(models)
    # checking if restraints exist
    found, not_found = 0, 0
    for r_one, r_two in restraints:
        if r_one[0] in resid_dict.keys() and r_two[0] in resid_dict.keys():
            if r_one[1] in resid_dict[r_one[0]] and r_two[1] in resid_dict[r_two[0]]:
                found += 1
                continue
        not_found += 1
    if not_found != 0:
        log.warning(f"{not_found} restraint were not valid for models {models} in {ambig_fname}. {found} are valid.")
    else:
        log.debug(f"{not_found} restraint were not valid for models {models} in {ambig_fname}. {found} are valid.")
    if found == 0:
        raise Exception(f"No valid restraints are available for models {models} in {ambig_fname}. Aborting")
    return


def parse_tbl(ambig_fname):
    """
    parse .tbl file
    
    Parameters
    ----------
    ambig_fname : str or Path
        restraints filename
    """
    file_content = open(ambig_fname, "r").read().split("\n")
    restr_pairs = []
    for ln in file_content:
        if ln != "":
            splt_ln = ln.split()
            if splt_ln[0] == "assign":
                first_chain = splt_ln[-1].rstrip(')')
                first_resid = int(splt_ln[3])
                first_partner = (first_chain, first_resid)
            elif splt_ln[0] == "(" and len(splt_ln) > 1:
                if splt_ln[1] == "resid":
                    sec_chain = splt_ln[-1].rstrip(')')
                    sec_resid = int(splt_ln[2])
                    second_partner = (sec_chain, sec_resid)
                    restr_pairs.append((first_partner, second_partner))
    return restr_pairs