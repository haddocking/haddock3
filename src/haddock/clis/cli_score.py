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
import ast
import sys


class _ParamsToDict(argparse.Action):
    """
    Convert command-line parameters in an argument to a dictionary.

    Example
    -------

    Where ``-x`` is an optional argument of the command-line client
    interface.

        >>> par1 1 par2 'my name' par3 [1,2,3] par4 True
        >>> {'par1': 1, 'par2': 'my name', 'par3': [1, 2, 3]}

    """

    def __call__(self, parser, namespace, ivalues, option_string=None):
        """Execute."""
        params = ivalues[::2]
        values = ivalues[1::2]

        if len(params) != len(values):
            raise parser.error(
                "The parameters and value pairs "
                "do not match for argument `-p`"
                )

        param_dict = {}
        for k, v in zip(params, values):
            try:
                param_dict[k] = ast.literal_eval(v)
            except (ValueError, TypeError, SyntaxError):
                raise parser.error(f"Parameter {k} with invalid value {v}")

        setattr(namespace, self.dest, param_dict)


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

ap.add_argument(
    "-p",
    "--other-params",
    dest="other_params",
    help=(
        "Any other parameter of the `emscoring` module."
        "For example: -p nemsteps 1000. You can give any number of "
        "parameters."
        ),
    action=_ParamsToDict,
    default={},
    nargs="*",
    )


def _ap():
    return ap


def load_args(ap):
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap, main):
    """Command-line interface entry point."""
    cmd = vars(load_args(ap))
    kwargs = cmd.pop("other_params")
    main(**cmd, **kwargs)


def maincli():
    """Execute main client."""
    cli(ap, main)


def main(
        pdb_file,
        full=False,
        outputpdb=False,
        outputpsf=False,
        keep_all=False,
        **kwargs,
        ):
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
        sys.exit(f'* ERROR * Input PDB file {str(input_pdb)!r} does not exist')

    # config all parameters are correctly spelled.
    default_emscoring = read_from_yaml_config(DEFAULT_CONFIG)
    ems_dict = default_emscoring.copy()
    n_warnings = 0
    for param in kwargs:
        if param not in default_emscoring:
            sys.exit(f'* ERROR * Parameter {param!r} is not a valid `emscoring` parameter')  # noqa:E501
        if kwargs[param] != default_emscoring[param]:
            print(f"* ATTENTION * Value ({kwargs[param]}) of parameter {param} different from default ({default_emscoring[param]})")  # noqa:E501
            ems_dict[param] = kwargs[param]
            n_warnings += 1
    
    if n_warnings != 0:
        print("* ATTENTION * Non-default parameter values were used. They should be properly reported if the output data are used for publication.")  # noqa:E501
    
    params = {
        "topoaa": {"molecules": [input_pdb]},
        "emscoring": ems_dict,
        }

    print("> starting calculations...")

    run_dir = Path("haddock-score-client")
    with suppress(FileNotFoundError):
        shutil.rmtree(run_dir)
    run_dir.mkdir()
    zero_fill.set_zerofill_number(2)

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
    haddock_score_itw = \
        ems_dict["w_vdw"] * vdw \
        + ems_dict["w_elec"] * elec \
        + ems_dict["w_desolv"] * desolv \
        + ems_dict["w_air"] * air \
        + ems_dict["w_bsa"] * bsa

    print(f"""> HADDOCK-score = ({ems_dict['w_vdw']} * vdw) + ({ems_dict['w_elec']} * elec) + ({ems_dict['w_desolv']} * desolv) + ({ems_dict['w_air']} * air) + ({ems_dict['w_bsa']} * bsa)""")  # noqa: E501
    print(f"> HADDOCK-score (emscoring) = {haddock_score_itw:.4f}")

    if full:
        print(f"> vdw={vdw},elec={elec},desolv={desolv},air={air},bsa={bsa}")

    if outputpdb:
        outputpdb_name = Path(f"{input_pdb.name}_hs.pdb")
        print(f"> writing {outputpdb_name}")
        shutil.copy(
            Path(run_dir, "1_emscoring", "emscoring_1.pdb"),
            outputpdb_name,
            )

    if outputpsf:
        outputpsf_name = Path(f"{input_pdb.name}_hs.psf")
        print(f"> writing {outputpsf_name}")
        shutil.copy(
            Path(run_dir, "0_topoaa", f'{input_pdb.name}_haddock.psf'),
            outputpsf_name,
            )

    if not keep_all:
        shutil.rmtree(run_dir)
    else:
        print(
            'The folder where the calculations where performed was kept. See '
            'folder: haddock-scoring-client'
            )


if __name__ == "__main__":
    sys.exit(maincli())
