"""haddock3-restraints z-surface-restraints subcommand.

Generate both z-restraints and corresponding z-surfaces based on
input pdb structure and residue selection.

Usage:
    haddock3-restraints z-surface-restraints
        --pdb      <path/to/the/structure.pdb>
        --residues <list of coma separated residue index>
        --output   <base/path/where/to/output/data>
        --spacing  <spacing_between_two_beads>
        --x-size   <size_in_x_dim>
        --y-size   <size_in_y_dim>

e.g:
haddock3-restraints z-surface-restraints
  --pdb mystructure.pdb
  --residues 1,2,3 7,8,9
  --spacing 40
  --x-size 200
  --y-size 200
  --output myZrestaints
"""


import logging
import os
from pathlib import Path
from haddock.core.typing import Union, Optional
from haddock.clis.restraints.z_surfaces import gen_bead_plans
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
        help="List of coma separated residues (can be multiple selections).",
        required=False,
        default=[],
        nargs='+',
        type=str,
        )
    
    z_surf_restraints_subcommand.add_argument(
        "--output", "-o",
        help="Base output path.",
        required=False,
        default=None,
        type=str,
        )

    z_surf_restraints_subcommand.add_argument(
        "--spacing", "-s",
        type=float,
        help="Spacing between two beads",
        required=False,
        default=40,
        )

    z_surf_restraints_subcommand.add_argument(
        "--x-size",
        "-x",
        help="Size of the plan in X dimension",
        required=False,
        default=200,
        type=float,
        )
    
    z_surf_restraints_subcommand.add_argument(
        "--y-size",
        "-y",
        help="Size of the plan in Y dimension",
        required=False,
        default=200,
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


############################
# SET OF USEFULL FUNCTIONS #
############################
def load_selected_resiudes_coords(
        pdb_fpath: Union[str, Path],
        selections: dict[str, list[int]],
        ) -> dict[str, list[tuple[float, float, float]]]:
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
    # Loop over atoms
    selection_coords: dict[str, list[tuple[float, float, float]]] = {}
    for _chain, resi, aname, coords in pdb_atoms:
        # Simplify the probleme to Calpha only
        if aname != 'CA':
            continue
        # Loop over selections
        for selection, sele_resis in selections.items():
            # Check if residue of interest
            if resi in sele_resis:
                # Cehck if key was created
                if selection not in selection_coords.keys():
                    selection_coords[selection] = []
                # Hold this coordinates
                selection_coords[selection].append(coords)
    return selection_coords


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
    selections = {
        f"selection_{listid}": [int(resid) for resid in resiudes.split(',')]
        for listid, resiudes in enumerate(residues_lists, start=1)
        }
    return selections


def get_z_coords(
        select_coords: dict[str, list[tuple[float, float, float]]],
        padding: float = 10,
        ) -> dict[str, float]:
    """Generate z-coordinates from selection of residues.

    Here the idea is to find the most distant points between selection,
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
        return {"z": 0}
    elif len(select_coords.keys()) == 1:
        return {s: 0 for s in select_coords.keys()}

    # Compute center of mass
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
    # Compute Z coords
    selection_z = {
        max_dist_keys[0]: max_z,
        max_dist_keys[1]: -max_z,
        }
    for select, dist in dists.items():
        if select in max_dist_keys:
            continue
        # Compute z coord
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
        z_padding: float = 1.0,
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
    z_padding : float, optional
        Padding (in Angstrom) for selection of z coordinates, by default 1.0

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
            if z_coord == minz:
                rest = (
                    "assign "
                    f"(resid {resid:>7d} and name CA and segid A) "
                    f"(name SHA and attr z lt {lt_coord:>-8.2f}) "
                    f"{rest_dist:.1f} {rest_dist:.1f} 0.0"
                    )
            elif z_coord == maxz:
                rest = (
                    "assign "
                    f"(resid {resid:>7d} and name CA and segid A) "
                    f"(name SHA and attr z gt {gt_coord:>-8.2f}) "
                    "{rest_dist:.1f} {rest_dist:.1f} 0.0"
                    )
            else:
                rest = (
                    "assign "
                    f"(resid {resid:>7d} and name CA and segid A) "
                    f"(name SHA and attr z lt {lt_coord:>-8.2f} "
                    f"and attr z gt {gt_coord:>-8.2f}) "
                    "{rest_dist:.1f} {rest_dist:.1f} 0.0"
                    )
            restraints.append(rest)
    all_restraints = os.linesep.join(restraints)
    return all_restraints


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


def main(
        pdb: Union[str, Path],
        residues: Optional[list[str]] = None,
        output: Optional[str] = None,
        spacing: float = 40,
        x_size: float = 200,
        y_size: float = 200,
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
        _description_
    """
    # Load residue selection
    res_select = load_selections(residues)
    # Load corresponding coordinates
    select_coords = load_selected_resiudes_coords(pdb, res_select)
    # Compute z-coordinates for each selection
    selection_z = get_z_coords(select_coords)
    # Generate corresponding z-surfaces
    plans = gen_bead_plans(
        spacing=spacing,
        x_size=x_size,
        y_size=y_size,
        z_coords=[zcoord for zcoord in selection_z.values()],
        )
    # Generate corresponding z-surface restraints
    restraints = gen_z_restraints(res_select, selection_z)
    # Output data
    restraints_tbl, plans_pdb = output_data(restraints, plans, output=output)
    return restraints_tbl, plans_pdb


gen_z_surfrace_restraints = main


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
    logging.basicConfig(
        level=args['log_level'],
        format='%(asctime)s L%(lineno)d %(levelname)s - %(message)s',
        datefmt='%d/%m/%Y %H:%M:%S',
        )
    # Launch main
    restraints, plan_s = gen_z_surfrace_restraints(
        args['residues'],
        residues=args['residues'],
        output=args['output'],
        spacing=args['spacing'],
        x_size=args['x_size'],
        y_size=args['y_size'],
        )
