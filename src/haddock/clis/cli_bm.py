"""
Prepare HADDOCK3 benchmark configuration files and job scripts.

Creates HADDOCK3 configuration files and job files. Details on each
parameter are explained in the `-h` menu.

USAGE:
    haddock3-bm -h
    haddock3-bm <BM folder> <results folder> --job-sys [OPTION]

A `BM folder` is a folder with the characteristics of:
    https://github.com/haddocking/BM5-clean
"""
import argparse
import shutil
import string
import sys
from functools import partial
from pathlib import Path

from haddock import log


capital_and_digits = tuple(string.digits + string.ascii_uppercase)


def _dir_path(path):
    path = Path(path)
    if not path.is_dir():
        _p = str(path.resolve())
        raise argparse.ArgumentTypeError(f'{_p!r} is not a directory.')
    return path


def _is_valid(f, cap_and_dig=capital_and_digits):
    """Assert if directory is a valid model directory."""
    _is_valid = \
        f.name.startswith(cap_and_dig) \
        and f.is_dir()
    return _is_valid


def create_cfg_test_daemon(
        run_dir,
        receptor_f,
        ligand_f,
        *ignore,
        **everythinelse):
    """
    Create HADDOCK3 configuration file just for topology.

    This function is usefull to use test the benchmark daemon.
    """
    cfg_str = \
f"""
run_dir = {str(run_dir)!r}
ncores = 2

molecules = [
    {str(receptor_f)!r},
    {str(ligand_f)!r}
    ]

[topoaa]

"""
    return cfg_str


def create_cfg_scn_1(run_dir, receptor_f, ligand_f, ambig_f):
    """
    Create HADDOCK3 configuration file.

    Parameters
    ----------
    run_dir : path or str
        Path to the run directory; where run results will be saved.

    receptor_f : Path or str
        Absolute path pointing to the receptor PDB file.

    ligand_f : Path or str
        Absolute path pointing to the ligand PDB file.

    ambig_f : Path or str
        Absolute path pointing to the `ambig.tbl` file.

    Return
    ------
    str
        The HADDOCK3 configuration file for benchmarking.
    """
    cfg_str = \
f"""
run_dir = {str(run_dir)!r}
ncores = 48

molecules = [
    {str(receptor_f)!r},
    {str(ligand_f)!r}
    ]

[topoaa]

[rigidbody]
ambig = {str(ambig_f)!r}
sampling = 1000
noecv = false

[seletop]
select = 200

[flexref]
ambig = {str(ambig_f)!r}

[mdref]
ambig = {str(ambig_f)!r}
"""
    return cfg_str


def create_cfg_scn_2(run_dir, receptor_f, ligand_f, ambig_f):
    """
    Create HADDOCK3 configuration file, scenario #2.

    Parameters
    ----------
    run_dir : path or str
        Path to the run directory; where run results will be saved.

    receptor_f : Path or str
        Absolute path pointing to the receptor PDB file.

    ligand_f : Path or str
        Absolute path pointing to the ligand PDB file.

    ambig_f : Path or str
        Absolute path pointing to the `ambig.tbl` file.

    Return
    ------
    str
        The HADDOCK3 configuration file for benchmarking.
    """
    cfg_str = \
f"""
run_dir = {str(run_dir)!r}
ncores = 48

molecules = [
    {str(receptor_f)!r},
    {str(ligand_f)!r}
    ]

[topoaa]
autohis = true

[rigidbody]
ambig = {str(ambig_f)!r}
sampling = 10000
noecv = false
cmrest = true

[flexref]
ambig = {str(ambig_f)!r}
cool1_steps = 500
#sampling = 400

[mdref]
ambig = {str(ambig_f)!r}
#sampling = 400
cool1_steps = 10
noecv = false
"""
    return cfg_str


def create_torque_job(job_name, scenario, path, **job_params):
    """
    Create HADDOCK3 Alcazer job file.

    Parameters
    ----------
    job_name : str
        The name of the job.

    **job_params
        According to `job_setup`.

    Return
    ------
    str
        Torque-based job file for HADDOCK3 benchmarking.
    """
    header = \
f"""#!/usr/bin/env tcsh
#PBS -N {job_name}-{scenario}-BM5
#PBS -q medium
#PBS -l nodes=1:ppn=48
#PBS -S /bin/tcsh
#PBS -o {str(path)}/haddock.out
#PBS -e {str(path)}/haddock.err
"""
    return job_setup(header, path, **job_params)


def create_slurm_job(job_name, scenario, root_path, run_path, **job_params):
    """
    Create HADDOCK3 Slurm Batch job file.

    Parameters
    ----------
    job_name : str
        The name of the job.

    **job_params
        According to `job_setup`.

    Return
    ------
    str
        Slurm-based job file for HADDOCK3 benchmarking.
    """
    header = \
f"""#!/usr/bin/env bash
#SBATCH -J {job_name}-{scenario}-BM5
#SBATCH -p medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH --output=/dev/null
#SBATCH --error=logs/{scenario}-haddock.err
#SBATCH --workdir={root_path}
"""
    return job_setup(header, run_path, **job_params)


def job_setup(header, run_path, conf_f):
    """
    Write body for the job script.

    Parameters
    ----------
    header : str
        Configures queue manager system.

    run_path : Path or str
        Path where results will be stored. Must be synchronized with the
        `run_dir` path in the configuration file.

    conf_f : Path or str
        Path to the configuration file.
    """
    conda_sh = Path(
        Path(sys.executable).parents[2],
        'etc', 'profile.d', 'conda.sh'
        )

    job = \
f"""{header}

source {str(conda_sh)}
conda activate haddock3

haddock3 {str(conf_f)}
"""
    return job


job_systems = {
    'torque': create_torque_job,
    'slurm': create_slurm_job,
    }

benchmark_scenarios = {
    'true-interface': create_cfg_scn_1,
    'center-of-mass': create_cfg_scn_2,
    }

test_scenarios = {
    'test-daemon': create_cfg_test_daemon,
    }

scenarios_acronym = {
    'true-interface': 'ti',
    'center-of-mass': 'com',
    'test-daemon': 'tdmn',
    }

ap = argparse.ArgumentParser(description='Setup HADDOCK3 benchmark.')
ap.add_argument(
    "benchmark_path",
    help=(
        'Location of BM5 folder. '
        'This folder should have a subfolder for each model. '
        'We expected each subfolder to have the name of a PDBID. '
        'That is, folder starting with capital letters or numbers '
        'will be considered.'
        ),
    type=_dir_path,
    )

ap.add_argument(
    "output_path",
    help="Where the executed benchmark will be stored.",
    type=_dir_path,
    )

ap.add_argument(
    '--job-sys',
    dest='job_sys',
    help='The system where the jobs will be run. Default `slurm`.',
    choices=list(job_systems.keys()),
    default='slurm',
    )

ap.add_argument(
    '-td',
    '--test-daemon',
    dest='test_daemon',
    help=(
        'Creates configuration files just with the `topology` module. '
        'Useful to test the benchmark daemon.'
        ),
    action='store_true',
    )


def process_target(source_path, result_path, create_job_func, scenarios):
    """
    Process each model example for benchmarking.

    Parameters
    ----------
    source_path : Path
        The folder path to the target.

    result_path : Path
        The path where the results for the different scenarios will be
        saved.

    create_job_func : callable
        A function to create the job script. This is an argument because
        there are several queue systems available. See `job_systems`.
    """
    pdb_id = source_path.name

    # Define the root directory
    root_p = Path(result_path, pdb_id).resolve()
    root_p.mkdir(exist_ok=True, parents=True)

    # Define the input directory
    #  all files required for the different scenarios
    #  should be inside this folder
    input_p = Path(root_p, 'input')
    input_p.mkdir(exist_ok=True, parents=True)

    ligand = shutil.copy(Path(source_path, f'{pdb_id}_l_u.pdb'), input_p)
    receptor = shutil.copy(Path(source_path, f'{pdb_id}_r_u.pdb'), input_p)
    ambig_tbl = shutil.copy(Path(source_path, 'ambig.tbl'), input_p)
    target = shutil.copy(Path(source_path, 'ana_scripts', 'target.pdb'), input_p)  # noqa: E501

    ligand = Path(ligand).relative_to(root_p)
    receptor = Path(receptor).relative_to(root_p)
    ambig_tbl = Path(ambig_tbl).relative_to(root_p)
    target = Path(target).relative_to(root_p)

    # Define the job folder
    #  all job files should be here
    job_p = Path(root_p, 'jobs')
    job_p.mkdir(exist_ok=True, parents=True)

    # Define the logs folder
    log_p = Path(root_p, 'logs')
    log_p.mkdir(exist_ok=True, parents=True)

    for scn_name, scn_func in scenarios.items():

        run_folder = Path(result_path, pdb_id, f'run-{scn_name}').resolve()
        run_folder = run_folder.relative_to(root_p)

        cfg_file = Path(input_p, f'{scn_name}.cfg')
        job_file = Path(job_p, f'{scn_name}.job')

        cfg_str = scn_func(
            run_folder,
            receptor,
            ligand,
            ambig_tbl,
            )

        job_str = create_job_func(
            pdb_id,
            scenario=scenarios_acronym[scn_name],
            root_path=root_p,
            run_path=run_folder,
            conf_f=cfg_file.relative_to(root_p),
            )

        cfg_file.write_text(cfg_str)
        job_file.write_text(job_str)

    return


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


def main(benchmark_path, output_path, job_sys='slurm', test_daemon=False):
    """
    Create configuration and job scripts for HADDOCK3 benchmarking.

    Developed for https://github.com/haddocking/BM5-clean

    Parameters
    ----------
    benchmark_path : Path
        The path to the benchmark models folder. In BM5-clean would be
        the 'HADDOCK-ready' folder.

    output_path : Path
        Where the results will be saved. A subfolder for each model in
        `benchmark_path` will be created.

    job_sys : str
        A key for `job_systems` dictionary.
    """
    log.info('*** Preparing Benchnark scripts')

    _ = (f for f in benchmark_path.glob('*') if _is_valid(f))
    source_folders = sorted(_)
    log.info(f'* Creating benchmark jobs for {len(source_folders)} targets')

    _scenarios = test_scenarios if test_daemon else benchmark_scenarios

    pe = partial(
        process_target,
        result_path=output_path,
        create_job_func=job_systems[job_sys],
        scenarios=_scenarios
        )

    for source_path in source_folders:
        pe(source_path)

    log.info('* done')


if __name__ == '__main__':
    sys.exit(maincli())
