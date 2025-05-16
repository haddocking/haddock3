"""haddock3-restraints random_removal subcommand.

Given an input ambiguous file (.tbl), this subcommand will generate an archive
containing multiple tbl files, each containing a subset of the initial ones.

The subset is tuned by the optional argument --ratio (-r)
The number of generated files in the archive is tuned by the argument --nb-tbl (-n)
The random seed can be tuned using --seed (-s), for reproducibility issues

Usage:
    haddock3-restraints random_removal <tbl_file> [-r <ratio>] [-n <nb-tbl>] [-s <seed>]

positional arguments:
  tblfile               input tbl restraint file.

options:
  -h, --help            show this help message and exit
  -r RATIO, --ratio RATIO
                        Ratio of restraints to be randomly removed.
  -s SEED, --seed SEED  Pseudo-random seed.
  -n NB_TBL, --nb-tbl NB_TBL
                        Number of ambiguous files to generate in the archive.
"""

import os
import sys
import tarfile
from io import BytesIO
from pathlib import Path

from haddock.core.typing import Union
from haddock.libs.librestraints import get_restraint_subset


def add_rand_removal_arguments(rand_removal_subcommand):
    """Add arguments to the random removal subcommand."""
    rand_removal_subcommand.add_argument(
        "tblfile",
        type=str,
        help="input tbl restraint file.",
        )
    
    rand_removal_subcommand.add_argument(
        "-r",
        "--ratio",
        help="Ratio of restraints to be randomly removed. (default: %(default)s)",
        required=False,
        default=0.5,
        type=float,
        )

    rand_removal_subcommand.add_argument(
        "-s",
        "--seed",
        help="Pseudo-random seed. (default: %(default)s)",
        required=False,
        default=917,  # Same as the one in various modules (iniseed)
        )

    rand_removal_subcommand.add_argument(
        "-n",
        "--nb-tbl",
        help="Number of ambiguous files to generate in the archive. (default: %(default)s)",
        required=False,
        type=int,
        default=10,
        )

    return rand_removal_subcommand


def random_removal(
        tblfile: Union[str, Path],
        ratio: float,
        nb_tbl: int = 10,
        seed: int = 917,
        ) -> Path:
    """Generate an archive containing the randomly removed restraints.

    Parameters
    ----------
    tblfile : Union[str, Path]
        Path the the input ambiguous file
    ratio : float
        Ration of restraints to be removed
    nb_tbl : int, optional
        Number of ambig files to generate in the archive, by default 10
    seed : int, optional
        Initial random seed, by default 917

    Returns
    -------
    Union[str, Path]
        Path to the generated archive containing `nb_tbl` restraints in it.
    """
    if nb_tbl < 1:
        sys.exit(
            "Number of restraints files to generate must be "
            f">= 1 (now set to {nb_tbl})"
        )
    # Initiate the tbl archive holding all the restraints
    tbl_archive_fpath = Path(
        Path(tblfile).parent.resolve(),
        Path(tblfile).name + ".tgz",
        )
    tbl_archive = tarfile.open(tbl_archive_fpath, "w:gz")

    # Initiate restraints subset generator
    try:
        subsets_restraints = get_restraint_subset(
            tblfile,
            ratio,
            seed=seed,
            )
    except ValueError as e:
        sys.exit(e)

    # Loop over number of tbl file to generate
    for i in range(1, nb_tbl + 1):
        # Obtain a subset of restraints
        subset_restraints = next(subsets_restraints)
        # Combine them to a string containing the restraints
        subset_restraints_str = os.linesep.join(subset_restraints)
        # Create tarinfo object
        tarinfo = tarfile.TarInfo(f"{Path(tblfile).stem}_rr{i}.tbl")
        tarinfo.size = len(subset_restraints_str)
        # Add this file to the archive
        tbl_archive.addfile(
            tarinfo,
            BytesIO(subset_restraints_str.encode("utf-8")),
            )
    # Closing things
    tbl_archive.close()
    subsets_restraints.close()
    # Return the generated archive path
    return tbl_archive_fpath


def main(
        tblfile: Union[str, Path],
        ratio: float,
        nb_tbl: int = 10,
        seed: int = 917,
        ) -> None:
    """Simple wrapper of the random_removal function.

    Parameters
    ----------
    tblfile : Union[str, Path]
        Path the the input ambiguous file
    ratio : float
        Ration of restraints to be removed
    nb_tbl : int, optional
        Number of ambig files to generate in the archive, by default 10
    seed : int, optional
        Initial random seed, by default 917

    Returns
    -------
    Union[str, Path]
        Path to the generated archive containing `nb_tbl` restraints in it.
    """
    archive_fpath = random_removal(tblfile, ratio, nb_tbl=nb_tbl, seed=seed)
    print(
        f"Generated archive path: {archive_fpath}!{os.linesep}"
        f"Note that resulting restraints might be redundant{os.linesep}"
        f"By using this file, we suggest to:{os.linesep}"
        f"- turn off random removal (random_removal = false){os.linesep}"
        "- turn on the previous ambig in later stages (previous_ambig = true)"
        )
