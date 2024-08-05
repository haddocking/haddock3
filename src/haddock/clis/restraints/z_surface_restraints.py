"""haddock3-restraints z-surface-restraints subcommand.

Generate both z-restraints and corresponding z-surfaces based on
input pdb structure and residue selection.

Usage:
    haddock3-restraints z-surface-restraints
        --pdb       <path/to/the/structure.pdb>
        --residues  <list of coma separated residue index>
        --output    <base/path/where/to/output/data>
        --spacing   <spacing_between_two_beads>
        --x-size    <size_in_x_dim>
        --y-size    <size_in_y_dim>
        --z-padding <additional_z_padding_between_two_extrema_plans>

e.g:
haddock3-restraints z-surface-restraints
  --pdb mystructure.pdb
  --residues 1,2,3 7,8,9
  --spacing 20
  --x-size 200
  --y-size 200
  --z-padding 5
  --output myZrestaints
"""

import logging
import os
from pathlib import Path

from haddock.core.typing import Generator, Optional, Union
from haddock.libs.librestraints import read_structure, calc_euclidean


# As this script is a subcommand of the `haddock3-restraints` cli,
# it requires its own options and arguments that are managed here.
def add_z_surf_restraints_arguments(z_surf_restraints_subcommand):
    """Add arguments to the z_plan subcommand."""
    z_surf_restraints_subcommand.add_argument(
        "--pdb", "-p",
        help="Path to a pdb file.",
        required=True,
        default='',
        type=str,
        )

    z_surf_restraints_subcommand.add_argument(
        "--residues", "-r",
        help=(
            "List of comma separated residues (can be multiple selections). "
            "Example 1,2,3 7,8,9 for two selections."
            ),
        required=False,
        default=[],
        nargs='+',
        type=str,
        )
    
    z_surf_restraints_subcommand.add_argument(
        "--output", "-o",
        help=(
            "Base output path. This script will generate two files, "
            "therefore no extention needed here"
            ),
        required=False,
        default=None,
        type=str,
        )

    z_surf_restraints_subcommand.add_argument(
        "--spacing", "-s",
        type=float,
        help="Spacing between two beads (A)",
        required=False,
        default=20,
        )

    z_surf_restraints_subcommand.add_argument(
        "--x-size",
        "-x",
        help="Size of the plan in X dimension (A)",
        required=False,
        default=100,
        type=float,
        )
    
    z_surf_restraints_subcommand.add_argument(
        "--y-size",
        "-y",
        help="Size of the plan in Y dimension",
        required=False,
        default=100,
        type=float,
        )
    
    z_surf_restraints_subcommand.add_argument(
        "--z-padding",
        "-z",
        help="Additional padding between two external plans.",
        required=False,
        default=5.0,
        type=float,
        )
    
    z_surf_restraints_subcommand.add_argument(
        "--log_level",
        default='INFO',
        choices=('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'),
        help="Logging level",
        required=False,
        )

    return z_surf_restraints_subcommand


def setup_logging(log_level: str = "INFO") -> None:
    """Set log level and format."""
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s L%(lineno)d %(levelname)s - %(message)s',
        datefmt='%d/%m/%Y %H:%M:%S',
        )


############################
# SET OF USEFULL FUNCTIONS #
############################
def load_selected_resiudes_coords(
        pdb_fpath: Union[str, Path],
        selections: dict[str, list[int]],
        ) -> tuple[dict[str, list[tuple[float, float, float]]], list[str], list[str]]:
    """Load coordinates of selected residues.

    Parameters
    ----------
    pdb_fpath : Union[str, Path]
        Path to the PDB file to be parsed.
    selections : dict[str, list[int]]
        Dictionary holding the various residues indices for each selection.

    Returns
    -------
    selection_coords: dict[str, list[tuple[float, float, float]]]
        Dictionary holding the various Calpha coordinates for each selection.
    """
    # Load structures atom coordinates
    pdb_atoms = read_structure(pdb_fpath)
    # Set paring variables
    selection_coords: dict[str, list[tuple[float, float, float]]] = {}
    selected_chains: list[str] = []
    selected_atoms: list[str] = []
    # Loop over atoms
    for chain, resi, atname, coords in pdb_atoms:
        # Simplify the probleme to Calpha/Phosphates/BackBone atoms only
        # FIXME: maybe change P to C1 or C9 ?
        if not atname in ("CA", "P", "BB", ):
            continue
        # Loop over selections
        for selection, sele_resis in selections.items():
            # Check if residue of interest
            if resi in sele_resis:
                # Make sure the key is created
                selection_coords.setdefault(selection, [])
                # Hold this coordinates
                selection_coords[selection].append(coords)
                # Also add chain and atom type
                selected_chains.append(chain)
                selected_atoms.append(atname)
    set_selected_chains = list(set(selected_chains))
    set_selected_atoms = list(set(selected_atoms))
    return selection_coords, set_selected_chains, set_selected_atoms


def compute_barycenter(
        resi_coords: list[tuple[float, float, float]],
        ) -> tuple[float, float, float]:
    """Compute center of mass of multiple resiudes coordinates.

    Parameters
    ----------
    resi_coords : list[tuple[float, float, float]]
        List of Calpha coordinates.

    Returns
    -------
    barycenter : tuple[float, float, float]
        Corrdinates of the center of mass.
    """
    # Combine coordinates
    xs, ys, zs = [], [], []
    for coords in resi_coords:
        xs.append(coords[0])
        ys.append(coords[1])
        zs.append(coords[2])
    # Compute average
    x_avg = sum(xs) / len(xs)
    y_avg = sum(ys) / len(ys)
    z_avg = sum(zs) / len(zs)
    # Return barycenter
    barycenter = (x_avg, y_avg, z_avg)
    return barycenter


def load_selections(residues_lists: list[str]) -> dict[str, list[int]]:
    """Split and cast residues from an initial string to list.

    Parameters
    ----------
    residues_lists : list[str]
        List of strings containing coma separated resiudes indices.

    Returns
    -------
    selections: dict[str, list[int]]
        Dictionary of resiudes indices.
    """
    selections: dict[str, list[int]] = {}
    for listid, str_resiudes in enumerate(residues_lists, start=1):
        resid_indices: list[int] = []
        selection_key = f"selection_{listid}"
        for strresid in str_resiudes.split(','):
            try:
                resid = int(strresid)
            except Exception as _e:
                msg = f"Could not cast residue {strresid} from {selection_key}"
                logging.warning(msg)
            else:
                resid_indices.append(resid)
        if resid_indices == []:
            err_msg = f"Not considering {selection_key} as it is empty !"
            logging.error(err_msg)
        else:
            selections[selection_key] = resid_indices
    return selections


def get_z_coords(
        select_coords: dict[str, list[tuple[float, float, float]]],
        padding: float = 5.0,
        ) -> dict[str, float]:
    """Generate z-coordinates from selection of residues.

    Here the idea is to find the most distant points between selections,
    and project it on a Z axis to be able to later orient the protein.

    Parameters
    ----------
    selection_coords: dict[str, list[tuple[float, float, float]]]
        Dictionary holding the various Calpha coordinates for each selection.
    padding : float, optional
        Extra padding (in Angstrom) of z-coordinate, by default 10

    Returns
    -------
    selection_z : dict[str, float]
        Z coodrinate for each selection.
    """
    # Default when less than 1 selection was made
    if select_coords == {}:
        return {"z_1": 0}
    elif len(select_coords.keys()) == 1:
        return {s: 0 for s in select_coords.keys()}

    # Compute geometrical center
    select_centers = {
        select: compute_barycenter(select_resids_coords)
        for select, select_resids_coords in select_coords.items()
        }
    # Compute distances
    dists: dict[str, dict[str, float]] = {s: {} for s in select_centers.keys()}
    max_dist: float = -1
    max_dist_keys: list[str] = []
    for select, center in select_centers.items():
        for select2, center2 in select_centers.items():
            if select == select2:
                continue
            # Compute dist
            dist = calc_euclidean(center, center2)
            # Hold data
            dists[select][select2] = dist
            dists[select2][select] = dist
            # Define max
            if max_dist < dist:
                max_dist = dist
                max_dist_keys = [select, select2]
    # Initiate z-coords boundaries
    max_z = (padding + max_dist) / 2
    # Compute Z coords for external surfaces
    selection_z = {
        max_dist_keys[0]: max_z,
        max_dist_keys[1]: -max_z,
        }
    # Compute Z coords for internal surfaces (when nb. selections >= 3)
    for select, dist in dists.items():
        # If the selection is part of the external surfaces
        if select in max_dist_keys:
            continue
        # Compute location of z coordinate
        dist_to_upper = max_z - dist[max_dist_keys[0]]
        dist_to_lower = dist[max_dist_keys[1]] - max_z
        select_z_coord = (dist_to_upper + dist_to_lower) / 2
        # Hold data
        selection_z[select] = select_z_coord
    
    return selection_z


def gen_z_restraints(
        res_select: dict[str, list[int]],
        selection_z: dict[str, float],
        rest_dist: float = 7.5,
        segids: list[str] = ["A"],
        atome_types: list[str] = ["CA"],
        ) -> str:
    """Generate set of z ambiguous restraints according to residue selections.

    Parameters
    ----------
    res_select : dict[str, list[int]]
        Dictionary holding the various residues indices for each selection.
    selection_z : dict[str, float]
        Z coodrinate for each selection.
    rest_dist : float, optional
        Upper boundary (in Angstrom) of satisfied restraints, by default 7.5

    Returns
    -------
    all_restraints : str
        A string containing the AIR restraints.
    """
    # Gather all coordinates
    all_z_coords = [zcoord for zcoord in selection_z.values()]
    # Point min and max values
    minz = min(all_z_coords)
    maxz = max(all_z_coords)
    # Initiate restraints holder
    restraints: list[str] = []
    # z_padding variable is used as a padding for the Z beads selection.
    # This is meant for the CNS selection method using lt(lower than)
    # and gt (geater than) rather than equal to.
    z_padding: float = 0.1
    # Compile chain selection
    chain_selection_string = _compile_multiple_cns_selections(
        "segid", segids,
        )
    # Compile atome selection
    atom_selection_string = _compile_multiple_cns_selections(
        "name", atome_types,
        )
    # Loop over selections
    for select in res_select.keys():
        # Point data
        residues = res_select[select]
        z_coord = selection_z[select]
        # Add comment
        list_residues = ",".join([str(r) for r in residues])
        restraints.append(f"! z restraints for {select}: {list_residues}")
        # Loop over residues selection
        for resid in residues:
            # Compute lower/greater than z-coord
            lt_coord = z_coord + z_padding
            gt_coord = z_coord - z_padding
            # If the lower Z coordinate plan
            if z_coord == minz:
                rest = (
                    f"assign (resid {resid:>7d} and "
                    f"{atom_selection_string} and {chain_selection_string}) "
                    f"(name SHA and attr z lt {lt_coord:>-8.2f}) "
                    f"{rest_dist:.1f} {rest_dist:.1f} 0.0"
                    )
            # If the upper Z coordinate plan
            elif z_coord == maxz:
                rest = (
                    f"assign (resid {resid:>7d} and "
                    f"{atom_selection_string} and {chain_selection_string}) "
                    f"(name SHA and attr z gt {gt_coord:>-8.2f}) "
                    f"{rest_dist:.1f} {rest_dist:.1f} 0.0"
                    )
            # If in between lower and upper plans (when nb. plans >= 3)
            else:
                rest = (
                    f"assign (resid {resid:>7d} and "
                    f"{atom_selection_string} and {chain_selection_string}) "
                    f"(name SHA and attr z lt {lt_coord:>-8.2f} "
                    f"and attr z gt {gt_coord:>-8.2f}) "
                    f"{rest_dist:.1f} {rest_dist:.1f} 0.0"
                    )
            restraints.append(rest)
    all_restraints = os.linesep.join(restraints)
    return all_restraints


def _compile_multiple_cns_selections(
        selection_key: str,
        selections: list[str],
        logical_operator: str = 'OR',
        ) -> str:
    """Generate a selection from multiple ones using logical operator.

    Parameters
    ----------
    selection_key : str
        Name of the selection key. (e.g.: segid, name, resid, ...)
    selections : list[str]
        List of selections.

    Returns
    -------
    combined_selection : str
        The combined selection.
    """
    assert len(selections) >= 1
    if len(selections) == 1:
        return f"{selection_key} {selections[0]}"
    compiled_selections = [
        f"{selection_key} {select}"
        for select in selections
        ]
    joined_selections = f" {logical_operator} ".join(compiled_selections)
    combined_selection = f"({joined_selections})"
    return combined_selection


def output_data(
        restraints: str,
        plans: str,
        output: Optional[Union[str, Path]] = None,
        ) -> tuple[str, str]:
    """Write output files.

    Parameters
    ----------
    restraints : str
        String containing the ambiguous restraints.
    plans : str
        String containing shape beads coordinates as PDB file.
    output : Optional[Union[str, Path]], optional
        Base output path, by default None

    Returns
    -------
    restraints_fpath: str
        Path to the generated AIRs.
    beadplans_fpath : str
        Path to the generated PDB file containing beads.
    """
    # Define base output path if not given
    if not output:
        output = 'Zrestraints'
    # Write restraints
    restraints_fpath = f"{output}.tbl"
    with open(restraints_fpath, 'w') as filout:
        filout.write(restraints)
    # Write restraints
    beadplans_fpath = f"{output}_beads.pdb"
    with open(beadplans_fpath, 'w') as filout:
        filout.write(plans)
    # Return filepaths
    return restraints_fpath, beadplans_fpath


def gen_bead_plans(
        spacing: float = 40,
        x_size: float = 200,
        y_size: float = 200,
        z_coords: Optional[list[float]] = None,
        ) -> str:
    """Generate multiple bead plans.

    Parameters
    ----------
    spacing : float, optional
        Spacing (in Angstrom) between beads in same dimension, by default 40
    x_size : float, optional
        Width (in Angstrom) of the plan, by default 200
    y_size : float, optional
        Height (in Angstrom) of the plan, by default 200
    z_coords : Optional[list[float]], optional
        List of z-coordinates where to generate plans, by default None

    Returns
    -------
    bead_plans : str
        A PDB file containing multiple plans.
    """
    # Presets
    bead_plans: str = ''
    resindex: int = 0
    if not z_coords:
        z_coords = [0]
    # Loop over z-coords
    for z in z_coords:
        plan, resindex = bead_plan(
            spacing=spacing,
            x_size=x_size,
            y_size=y_size,
            z_coord=z,
            resindex=resindex,
            )
        bead_plans += plan
    return bead_plans
    
    
def bead_plan(
        spacing: float = 40,
        x_size: float = 200,
        y_size: float = 200,
        z_coord: float = 0,
        resindex: int = 0,
        ) -> tuple[str, int]:
    """Generate a PDB plan made of beads.

    Parameters
    ----------
    spacing : float, optional
        Spacing (in Angstrom) between beads in same dimension, by default 40
    x_size : float, optional
        Width (in Angstrom) of the plan, by default 200
    y_size : float, optional
        Height (in Angstrom) of the plan, by default 200
    z_coord : float, optional
        Z-coordinate where to generate the plan, by default 0
    resindex : int, optional
        From which resiude to start the , by default 0

    Returns
    -------
    plan : str
        The PDB plan made of beads.
    resindex : int
        Index of the last residue index added.
    """
    plan_beads: list[str] = []
    # Loop over x coords
    for x_coord in step_coords(x_size, spacing):
        # Loop over y coords
        for y_coord in step_coords(y_size, spacing):
            resindex += 1
            # Generate new bead
            bead = shape_bead(x_coord, y_coord, z_coord, resindex)
            plan_beads.append(bead)
    # Finalize plan
    plan = ''.join(plan_beads)
    return plan, resindex


def step_coords(_size: float, _spacing: float) -> Generator[float, None, None]:
    """Generate set of evenly spaced coordinates between of defined size.

    Parameters
    ----------
    size : float
        Size (in Angstrom) to be sampled
    spacing : float
        Spacing between each coordinate

    Return
    ------
    Generator[float, None, None]
        1D coodinate of current position.
    """
    # Convert to absolute value
    size = abs(_size)
    spacing = abs(_spacing)
    # Check if size is indeed greater than spacing
    if spacing > size:
        logging.error("Size must be greater than spacing!")
        raise ValueError
    return _step_coords(size, spacing)


def _step_coords(size: float, spacing: float) -> Generator[float, None, None]:
    """Generate set of evenly spaced coordinates between of defined size.

    Parameters
    ----------
    size : float
        Size (in Angstrom) to be sampled
    spacing : float
        Spacing between each coordinate

    Yields
    ------
    float
        1D coodinate of current position.
    """
    # Define initial position
    coord = - (size / 2)
    # Define oversized position
    oversized = (size + spacing) / 2
    # Loop until oversized
    while coord <= oversized:
        yield coord
        coord += spacing


def shape_bead(
        x: float,
        y: float,
        z: float,
        resindex: int,
        chain: str = "S",
        atindex: int = 1,
        bfactor: float = 1.00,
        ) -> str:
    """Generate a PDB shape bead.

    Parameters
    ----------
    x : float
        x coordinate of the bead
    y : float
        y coordinate of the bead
    z : float
        z coordinate of the bead
    resindex : int
        Residue index
    chain : str, optional
        Chain id, by default "S"
    atindex : int, optional
        Atome index, by default 1
    bfactor : float, optional
        B-factor of the bead, by default 1.00

    Returns
    -------
    bead : str
        A valid PDB shape bead.
    """
    bead = f"ATOM  {atindex:>5d}  SHA SHA {chain}{resindex:>4d}    {x:-8.3f}{y:-8.3f}{z:-8.3f}  1.00{bfactor:6.2f}         SHA  {os.linesep}"  # noqa : E501
    return bead


def _get_ideal_restraint_dist(spacing: float) -> float:
    """Computes ideal restraint distance based on spacing.

    Basically, want to return (spacing / 2) - 2

    Parameters
    ----------
    spacing : int
        Spacing between beads.

    Returns
    -------
    float
        Effective distance restraint.
    """
    # Compute de distance
    effective_distance_restraint: float = (spacing / 2) - 2
    # Make sure it is not out of allowed boundaries
    bounded_distance_restraint = max(effective_distance_restraint, 2.0)
    return bounded_distance_restraint


def main(
        pdb: Union[str, Path],
        residues: Optional[list[str]] = None,
        output: Optional[str] = None,
        spacing: float = 40,
        x_size: float = 200,
        y_size: float = 200,
        z_padding: float = 5.0,
        log_level: str = "INFO",
        ) -> tuple[str, str]:
    """Generate both z-restraints and z-surface beads from residue selection.

    Parameters
    ----------
    pdb : Union[str, Path]
        Path to the PDB file to be parsed.
    residues : list[str]
        List of strings containing coma separated resids, by default None
    spacing : float, optional
        Spacing (in Angstrom) between beads in same dimension, by default 40
    x_size : float, optional
        Width (in Angstrom) of the plan, by default 200
    y_size : float, optional
        Height (in Angstrom) of the plan, by default 200

    Returns
    -------
    tuple[Union[str, Path], Union[str, Path]]
        Paths to the Z_restraints.tbl and Z_surface.pdb
    """
    # Load residue selection
    res_select = load_selections(residues)
    # Load corresponding coordinates
    select_coords, chainids, atomtypes = load_selected_resiudes_coords(
        pdb, res_select
        )
    # Compute z-coordinates for each selection
    selection_z = get_z_coords(select_coords, padding=z_padding)
    # Generate corresponding z-surfaces
    plans = gen_bead_plans(
        spacing=spacing,
        x_size=x_size,
        y_size=y_size,
        z_coords=[zcoord for zcoord in selection_z.values()],
        )
    # Compute ideal restraint distance based on spacing
    restraint_distance = _get_ideal_restraint_dist(spacing)
    # Generate corresponding z-surface restraints
    restraints = gen_z_restraints(
        res_select,
        selection_z,
        rest_dist=restraint_distance,
        segids=chainids,
        atome_types=atomtypes,
        )
    # Output data
    restraints_tbl, plans_pdb = output_data(restraints, plans, output=output)
    return restraints_tbl, plans_pdb


gen_z_surface_restraints = main


############################
# COMMAND LINE ENTRY POINT #
############################
if __name__ == "__main__":
    import argparse
    # Command line interface parser
    ap = argparse.ArgumentParser(
        prog="haddock3-restraints",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    add_z_surf_restraints_arguments(ap)
    args = vars(ap.parse_args())
    setup_logging(log_level=args['log_level'])
    # Launch main
    restraints_fpath, plan_s_fpath = gen_z_surface_restraints(
        args['residues'],
        residues=args['residues'],
        output=args['output'],
        spacing=args['spacing'],
        x_size=args['x_size'],
        y_size=args['y_size'],
        z_padding=args['z_padding'],
        log_level=args["log_level"],
        )
    logging.info(restraints_fpath)
    logging.info(plan_s_fpath)
