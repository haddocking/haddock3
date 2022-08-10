"""
A simple tool to calculate the HADDOCK-score of a complex.

Usage:
    haddock3-score complex.pdb
"""
import argparse
import sys


ap = argparse.ArgumentParser(
    prog="haddock3-score",
    description=__doc__,
    )

ap.add_argument("pdb_file", help="Input PDB file")

ap.add_argument(
    "--full",
    action="store_true",
    help="Print all energy components",
    )

ap.add_argument(
    "--outputpdb",
    action="store_true",
    help="Save the output PDB file (minimized structure)",
    )

ap.add_argument(
    "--outputpsf",
    action="store_true",
    help="Save the output PSF file (topology)",
    )

ap.add_argument(
    "-k"
    "--keep-all",
    dest="keep_all",
    action="store_true",
    help="Keep the whole run folder."
    )


def _ap():
    return ap


def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = load_args(ap)
    main(**vars(cmd))


def maincli():
    """Execute main client."""
    cli(ap, main)


def main(
        pdb_file,
        full=False,
        outputpdb=False,
        outputpsf=False,
        keep_all=False,
        ):
    """Calculate the score of a complex."""
    import logging
    import shutil
    from pathlib import Path

    from haddock import log
    from haddock.gear.haddockmodel import HaddockModel
    from haddock.gear.zerofill import zero_fill
    from haddock.libs.libio import working_directory
    from haddock.libs.libworkflow import WorkflowManager

    log.setLevel(logging.ERROR)

    input_pdb = Path(pdb_file).resolve()

    run_dir = Path("haddock-score-client")
    run_dir.mkdir()
    zero_fill.set_zerofill_number(2)

    params = {
        "topoaa": {"molecules": [input_pdb]},
        "emscoring": {},
        }

    print("> starting calculations...")  # noqa: T201
    with working_directory(run_dir):
        workflow = WorkflowManager(
            workflow_params=params,
            start=0,
            run_dir=run_dir,
            )

        workflow.run()

    minimized_mol = Path(run_dir, "1_emscoring", "emscoring_1.pdb")
    haddock_score_component_dic = HaddockModel(minimized_mol).energies

    vdw = haddock_score_component_dic["vdw"]
    elec = haddock_score_component_dic["elec"]
    desolv = haddock_score_component_dic["desolv"]
    air = haddock_score_component_dic["air"]
    # bsa = haddock_score_component_dic["bsa"]

    # emscoring is equivalent to itw
    haddock_score_itw = \
        1.0 * vdw \
        + 0.2 * elec \
        + 1.0 * desolv \
        + 0.1 * air

    print("> HADDOCK-score = (1.0 * vdw) + (0.2 * elec) + (1.0 * desolv) + (0.1 * air)")  # noqa: T201, E501
    print(f"> HADDOCK-score (emscoring) = {haddock_score_itw:.4f}")  # noqa: T201, E501

    if outputpdb:
        shutil.copy(
            Path("run_dir", "1_emscoring", "emscoring_1.pdb"),
            Path(f"{input_pdb.name}_hs.pdb"),
            )

    if outputpsf:
        shutil.copy(
            Path(run_dir, "0_topoaa", f'{input_pdb.name}_haddock.psf'),
            Path(f"{input_pdb.name}_hs.pdb"),
            )

    if not keep_all:
        shutil.rmtree(run_dir)


if __name__ == "__main__":
    sys.exit(maincli())
