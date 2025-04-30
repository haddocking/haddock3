"""haddock3-restraints active_passive_to_ambig subcommand.

Given two files containing active (in the first line) and passive (second line)
 residues to be used by HADDOCK, this command gives in output the corresponding
 ambig.tbl file.

Usage:
    haddock3-restraints active_passive_to_ambig file_actpass_one file_actpass_two [--segid-one] [--segid-two]

An example content for file_actpass_one is
    72 73 74 75 81 83 84 89 90 92 94 96 97 98 115 116 117
    3 24 46 47 48 50 66 76 77 79 80 82 86 87 88 91 93 95 118 119 120
"""
from haddock.libs.librestraints import active_passive_to_ambig, parse_actpass_file


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


def actpass_to_ambig(actpass_one, actpass_two, segid_one, segid_two):
    """Generate ambig from two actpass files.
    
    Parameters
    ----------
    actpass_one : str
        path to first actpass file
    
    actpass_two : str
        path to second actpass file
    
    segid_one : str
        segid to use for the first model
    
    segid_two : str
        segid to use for the second model
    """
    active1, passive1 = parse_actpass_file(actpass_one)
    active2, passive2 = parse_actpass_file(actpass_two)
    # Check if there are active residues in at least one of the two files
    if active1 == [] and active2 == []:
        raise ValueError(
            f"No active residues found in {actpass_one} and {actpass_two}. "
             "No restraints will be generated."
            )
    active_passive_to_ambig(active1, passive1, active2, passive2, segid_one, segid_two)
    
    return
