"""haddock3-restraints z-surface subcommand.

Generate a set of shape beads on the x, y axis

Usage:
    haddock3-restraints z-surface
        --spacing  <spacing_between_two_beads>
        --x-size   <size_in_x_dim>
        --y-size   <size_in_y_dim>
        --z-coords <list of z positions where to generate surfaces>
"""


import logging
import os

from haddock.core.typing import Optional, Generator


# As this script is a subcommand of the `haddock3-restraints` cli,
# it requires its own options and arguments that are managed here.
def add_z_plan_arguments(z_plan_subcommand):
    """Add arguments to the z_plan subcommand."""
    z_plan_subcommand.add_argument(
        "--spacing", "-s",
        type=float,
        help="Spacing between two beads",
        required=False,
        default=40,
        )

    z_plan_subcommand.add_argument(
        "--x-size",
        "-x",
        help="Size of the plan in X dimension",
        required=False,
        default=200,
        type=float,
        )
    
    z_plan_subcommand.add_argument(
        "--y-size",
        "-y",
        help="Size of the plan in Y dimension",
        required=False,
        default=200,
        type=float,
        )
    
    z_plan_subcommand.add_argument(
        "--z-coords",
        "-z",
        help="List of z levels where to generate plans.",
        required=False,
        default=[0],
        nargs='+',
        type=float,
        )
    
    z_plan_subcommand.add_argument(
        "--log_level",
        default='INFO',
        choices=('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'),
        help="Logging level",
        required=False,
        )

    return z_plan_subcommand


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


def step_coords(size: float, spacing: float) -> Generator[float]:
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
    oversized = (size / 2) + spacing
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


if __name__ == "__main__":
    import argparse
    # Command line interface parser
    ap = argparse.ArgumentParser(
        prog="haddock3-restraints",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        )
    add_z_plan_arguments(ap)
    args = vars(ap.parse_args())
    logging.basicConfig(
        level=args['log_level'],
        format='%(asctime)s L%(lineno)d %(levelname)s - %(message)s',
        datefmt='%d/%m/%Y %H:%M:%S',
        )
    # Run main function
    plan_s = gen_bead_plans(
        spacing=args['spacing'],
        x_size=args['x_size'],
        y_size=args['y_size'],
        z_coords=args['z_coords'],
        )
    print(plan_s)  # noqa : T201
