"""
A simple tool to calculate the HADDOCK-score of a complex.

You can pass to the command-line any parameter accepted by the `emscoring`
module. For this, use the ``-p`` option writing the name of the parameters
followed by the desired value. Write booleans with capital letter.

Use the ``haddock3-cfg`` command-line to obtain the list of parameters for
the ``emscoring`` module.

Usage::

    haddock3-score complex.pdb
    haddock3-score complex.pdb -p nemsteps 50
    haddock3-score complex.pdb -p nemsteps 50 w_air 1
    haddock3-score complex.pdb -p nemsteps 50 w_air 1 electflag True

"""
import argparse
import sys
import tempfile

from haddock.core.typing import (
    Any,
    ArgumentParser,
    Callable,
    FilePath,
    Namespace,
    )
from haddock.libs.libcli import _ParamsToDict


ap = argparse.ArgumentParser(
    prog="haddock3-score",
    description=__doc__,
    )

ap.add_argument("pdb_file", help="Input PDB file")

ap.add_argument(
    "--run_dir",
    default="haddock-score-client",
    type=str,
    required=False,
    help="Run directory name.",
    )

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
    "-k" "--keep-all",
    dest="keep_all",
    action="store_true",
    help="Keep the whole run folder.",
    )

ap.add_argument(
    "-p",
    "--other-params",
    dest="other_params",
    help=(
        "Any other parameter of the `emscoring` module. "
        "For example: -p nemsteps 1000. You can give any number of "
        "parameters."
    ),
    action=_ParamsToDict,
    default={},
    nargs="*",
    )


def _ap() -> ArgumentParser:
    return ap


def load_args(ap: ArgumentParser) -> Namespace:
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap: ArgumentParser, main: Callable[..., None]) -> None:
    """Command-line interface entry point."""
    cmd = vars(load_args(ap))
    kwargs = cmd.pop("other_params")
    main(**cmd, **kwargs)


def maincli() -> None:
    """Execute main client."""
    cli(ap, main)


def main(
        pdb_file: FilePath,
        run_dir: FilePath,
        full: bool = False,
        outputpdb: bool = False,
        outputpsf: bool = False,
        keep_all: bool = False,
        **kwargs: Any,
        ) -> None:
    """
    Calculate the score of a complex using the ``emscoring`` module.

    Parameters
    ----------
    pdb_file : str or pathlib.Path
        The path to the PDB containing the complex.

    full : bool
        Print all energy components.

    outputpdb : bool
        Save the PDB file resulting from the scoring calculation.

    outputpsf : bool
        Save  the PSF file resulting from the scoring calculation; this
        is the CNS topology created before running the calculation.

    keep_all : bool
        Keep the whole temporary run folder. ``haddock3-score`` creates
        a temporary run folder where the calculations are performed. If
        ``keep_all`` is True, this folder is **not** deleted after when
        the calculation finishes.

    kwargs : any
        Any additional arguments that will be passed to the ``emscoring``
        module.
    """
    import os
    import logging
    import shutil
    from contextlib import suppress
    from pathlib import Path

    from haddock import log
    from haddock.gear.haddockmodel import HaddockModel
    from haddock.gear.yaml2cfg import read_from_yaml_config
    from haddock.gear.zerofill import zero_fill
    from haddock.libs.libio import working_directory
    from haddock.libs.libworkflow import WorkflowManager
    from haddock.modules.scoring.emscoring import DEFAULT_CONFIG

    log.setLevel(logging.ERROR)

    input_pdb = Path(pdb_file).resolve()
    if not input_pdb.exists():
        sys.exit(f"* ERROR * Input PDB file {str(input_pdb)!r} does not exist")

    # config all parameters are correctly spelled.
    default_emscoring = read_from_yaml_config(DEFAULT_CONFIG)
    ems_dict = default_emscoring.copy()
    n_warnings = 0
    for param, value in kwargs.items():
        if param not in default_emscoring:
            sys.exit(
                f"* ERROR * Parameter {param!r} is not a "
                f"valid `emscoring` parameter.{os.linesep}"
                f"Valid emscoring parameters are: {', '.join(sorted(default_emscoring))}"
                )
        if value != default_emscoring[param]:
            print(
                f"* ATTENTION * Value ({value}) of parameter {param} different from default ({default_emscoring[param]})"
            )  # noqa:E501
            # get the type of default value
            default_type = type(default_emscoring[param])
            # convert the value to the same type
            if default_type == bool:
                if value.lower() not in ["true", "false"]:
                    sys.exit(f"* ERROR * Boolean parameter {param} should be True or False")
                value = value.lower() == "true"
            else:
                value = default_type(value)
            ems_dict[param] = value
            n_warnings += 1
    if n_warnings != 0:
        print(
            "* ATTENTION * Non-default parameter values were used. "
            "They should be properly reported if the output "
            "data are used for publication."
            )
        print(f"used emscoring parameters: {ems_dict}")

    # create run directory
    run_dir = Path(run_dir)
    with suppress(FileNotFoundError):
        shutil.rmtree(run_dir)
    run_dir.mkdir()
    zero_fill.set_zerofill_number(2)

    # create temporary file
    with tempfile.NamedTemporaryFile(prefix=input_pdb.stem, suffix=".pdb") as tmp:

        # create a copy of the input pdb
        input_pdb_copy = Path(tmp.name)
        shutil.copy(input_pdb, input_pdb_copy)
    
        params = {
            "topoaa": {"molecules": [input_pdb_copy]},
            "emscoring": ems_dict,
        }

        print("> starting calculations...")

        # run workflow
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
    bsa = haddock_score_component_dic["bsa"]

    # emscoring is equivalent to itw
    haddock_score_itw = (
        ems_dict["w_vdw"] * vdw
        + ems_dict["w_elec"] * elec
        + ems_dict["w_desolv"] * desolv
        + ems_dict["w_air"] * air
        + ems_dict["w_bsa"] * bsa
        )

    print(
        "> HADDOCK-score ="
        f" ({ems_dict['w_vdw']} * vdw)"
        f" + ({ems_dict['w_elec']} * elec)"
        f" + ({ems_dict['w_desolv']} * desolv)"
        f" + ({ems_dict['w_air']} * air)"
        f" + ({ems_dict['w_bsa']} * bsa)"
        )
    print(f"> HADDOCK-score (emscoring) = {haddock_score_itw:.4f}")

    if full:
        print(f"> vdw={vdw},elec={elec},desolv={desolv},air={air},bsa={bsa}")

    if outputpdb:
        outputpdb_name = Path(f"{input_pdb.stem}_hs.pdb")
        print(f"> writing {outputpdb_name}")
        shutil.copy(
            Path(run_dir, "1_emscoring", "emscoring_1.pdb"),
            outputpdb_name,
            )

    if outputpsf:
        outputpsf_name = Path(f"{input_pdb.stem}_hs.psf")
        print(f"> writing {outputpsf_name}")
        shutil.copy(
            Path(run_dir, "0_topoaa", f"{input_pdb_copy.stem}_haddock.psf"),
            outputpsf_name,
            )

    if not keep_all:
        shutil.rmtree(run_dir)
    else:
        print(
            'The folder where the calculations were performed was kept.'
            f' See folder: {run_dir}'
            )


if __name__ == "__main__":
    sys.exit(maincli())  # type: ignore
