"""restraints-related functions."""

from haddock import log
from haddock.libs.libalign import get_atoms, load_coords


def parse_tbl(ambig_fname):
    """
    Parse .tbl file.
    
    Parameters
    ----------
    ambig_fname : str or Path
        restraints filename
    
    Returns
    -------
    restr_pairs : list of tuples
        list of restraints
    """
    file_content = open(ambig_fname, "r").read().split("\n")
    restr_pairs = []
    for ln in file_content:
        if ln != "":
            splt_ln = ln.split()
            if splt_ln[0] == "assign":  # active residue
                first_chain = splt_ln[6].rstrip(')')
                first_resid = int(splt_ln[3])
                first_partner = (first_chain, first_resid)
            elif splt_ln[0] == "(" and len(splt_ln) > 1:
                if splt_ln[1] == "resid":  # partner residue
                    sec_chain = splt_ln[5].rstrip(')')
                    sec_resid = int(splt_ln[2])
                    second_partner = (sec_chain, sec_resid)
                    restr_pairs.append((first_partner, second_partner))
    return restr_pairs


def get_unique_resids(models):
    """
    Get unique residues for each chain in input models.

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
    Validate ambig_fname against input models.

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
    res_dict = get_unique_resids(models)
    # checking if restraints exist
    found, not_found = 0, 0
    for rone, rtwo in restraints:
        if rone[0] in res_dict.keys() and rtwo[0] in res_dict.keys():
            if rone[1] in res_dict[rone[0]] and rtwo[1] in res_dict[rtwo[0]]:
                found += 1
                continue
        not_found += 1
    msg = (
        f"{not_found} restraints not valid for {models} in {ambig_fname}."
        f" {found} are valid."
        )
    if not_found != 0:
        log.warning(msg)
    else:
        log.debug(msg)
    if found == 0:
        raise Exception(f"No valid restraints for {models} in {ambig_fname}.")
    return True
