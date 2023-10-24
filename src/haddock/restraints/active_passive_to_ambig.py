"""haddock3-restraints active_passive_to_ambig subcommand.

Given two files containing active (in the first line) and passive (second line)
to be used by HADDOCK, this command gives in output the corresponding ambig.tbl
 file.

Usage:
    haddock3-restraints active_passive_to_ambig file_actpass_one 

An example content for file_actpass_one is
    72 73 74 75 81 83 84 89 90 92 94 96 97 98 115 116 117
    3 24 46 47 48 50 66 76 77 79 80 82 86 87 88 91 93 95 118 119 120
"""
from pathlib import Path


def add_actpass_to_ambig_arguments(actpass_to_ambig_subcommand):
    """Add arguments to the score subcommand."""
    actpass_to_ambig_subcommand.add_argument(
        "actpass_one",
        type=str,
        help="First actpass file",
        )

    actpass_to_ambig_subcommand.add_argument(
        "actpass_two",
        type=str,
        help="Second actpass file",
        )

    actpass_to_ambig_subcommand.add_argument(
        "--segid-one",
        type=str,
        help="Segid to use for the first model",
        default="A"
        )
    
    actpass_to_ambig_subcommand.add_argument(
        "--segid-two",
        type=str,
        help="Segid to use for the second model",
        default="B"
        )

    return actpass_to_ambig_subcommand


def parse_actpass_file(actpass_file):
    """Parse actpass file
    
    Parameters
    ----------
    actpass_file: str or Path
        path to actpass_file

    Returns
    -------
    active: list
        list of active residues
    passive: list
        list of passive residues
    """
    
    if Path(actpass_file).exists() is False:
        raise Exception(f"actpass file {actpass_file} does not exist.")

    lines = open(actpass_file, "r").readlines()
    nlines = len(lines)
    if nlines != 2:
        raise Exception(f"actpass file {actpass_file} does not have two lines (counted {nlines})")
    active, passive = [[int(x) for x in line.split()] for line in lines]
    return active, passive


def active_passive_to_ambig(active1, passive1, active2, passive2, segid1='A', segid2='B'):
    """Convert active and passive residues to Ambiguous Interaction Restraints

    Parameters
    ----------
    active1 : list
        List of active residue numbers of the first segid

    passive1 : list
        List of passive residue numbers of the first segid

    passive2 : list
        List of passive residue numbers of the second segid

    active2 : list
        List of active residue numbers of the second segid
    
    active2 : list
        List of passive residue numbers of the second segid

    segid1 : string
        Segid to use for the first model

    segid2 : string
        Segid to use for the second model

    """

    all1 = active1 + passive1
    all2 = active2 + passive2

    for resi1 in active1:
        print('assign (resi {:d} and segid {:s})'.format(resi1, segid1))
        print('(')
        c = 0
        for resi2 in all2:
            print('       (resi {:d} and segid {:s})'.format(resi2, segid2))
            c += 1
            if c != len(all2):
                print('        or')

        print(') 2.0 2.0 0.0\n')
            
    for resi2 in active2:
        print('assign (resi {:d} and segid {:s})'.format(resi2, segid2))
        print('(')
        c = 0
        for resi1 in all1:
            print('       (resi {:d} and segid {:s})'.format(resi1, segid1))
            c += 1
            if c != len(all1):
                print('        or')

        print(') 2.0 2.0 0.0\n')

def actpass_to_ambig(actpass_one, actpass_two, segid_one, segid_two):
    """generate ambig from two actpass files."""

    active1, passive1 = parse_actpass_file(actpass_one)
    active2, passive2 = parse_actpass_file(actpass_two)
    active_passive_to_ambig(active1, passive1, active2, passive2, segid_one, segid_two)
    #
    return 