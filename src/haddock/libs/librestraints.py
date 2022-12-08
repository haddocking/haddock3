"""restraints-related functions."""

import numpy as np
import time
from haddock.libs.libalign import get_atoms, load_coords


def validate_ambig_fname(ambig_fname, models):
    """validates ambig_fname"""
    st_time = time.time()
    print(f"validating {ambig_fname}")
    restraints = parse_tbl(ambig_fname)
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
    # checking if restraints exist
    found, not_found = 0, 0
    for res_one, res_two in restraints:
        if res_one[0] in chains.keys() and res_two[0] in chains.keys():
            if res_one[1] in chains[res_one[0]] and res_two[1] in chains[res_two[0]]:
                found += 1
                continue
        not_found += 1
    print(f"found {found} vs not_found {not_found} restraints")
    if found == 0:
        raise Exception("No valid restraints are available. Aborting")
    print(f"time for validating{time.time() - st_time}")
    return


def parse_tbl(restr_file):
    """parsing tbl file"""
    file_content = open(restr_file, "r").read().split("\n")
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