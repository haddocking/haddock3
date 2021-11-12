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
autohis = true

[rigidbody]
ambig = {str(ambig_f)!r}
sampling = 400
noecv = false

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


def create_torque_job(job_name, path, **job_params):
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
#PBS -N {job_name}-BM5
#PBS -q medium
#PBS -l nodes=1:ppn=48
#PBS -S /bin/tcsh
#PBS -o {str(path)}/haddock.out
#PBS -e {str(path)}/haddock.err
"""
    return job_setup(header, path, **job_params)


def create_slurm_job(job_name, path, **job_params):
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
        Slurm-based job file for HADDOCK3 benchmarking.
    """
    header = \
f"""#!/usr/bin/env bash
#SBATCH -J {job_name}-BM5
#SBATCH -p medium
#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH -o {str(path)}/haddock.out
#SBATCH -e {str(path)}/haddock.err
"""
    return job_setup(header, path, **job_params)


def job_setup(header, path, conf_f):
    """
    Write body for the job script.

    Parameters
    ----------
    header : str
        Configures queue manager system.

    path : Path or str
        Path where results will be stored. Must be synchronized with the
        `run_dir` path in the configuration file.

    conf_f : Path or str
        Path to the configuration file.
    """
    job = \
f"""{header}

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate haddock3

cd {str(path)}
touch RUNNING
haddock3 {str(conf_f)}
rm RUNNING

#if ( -f {str(path)}/run-ti/ss.stats ) then
#    touch DONE
#else
#    touch FAIL
#endif"""
    return job


job_systems = {
    'torque': create_torque_job,
    'slurm': create_slurm_job,
    }

scenarios = {
    'true-interface': create_cfg_scn_1,
    'center-of-mass': create_cfg_scn_2,
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
    "results_path",
    help="Where results will be stored.",
    type=_dir_path,
    )

ap.add_argument(
    '--job-sys',
    dest='job_sys',
    help='The system where the jobs will be run. Default `torque`.',
    choices=list(job_systems.keys()),
    default='torque',
    )


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


def process_example(source_path, results_path, create_job_func):
    """
    Process each model example for benchmarking.

    Parameters
    ----------
    source_path : Path
        The folder path to the model example.

    results_path : Path
        The path where the results for the different scenarios will be
        saved.

    create_job_func : callable
        A function to create the job script. This is an argument because
        there are several queue systems available. See `job_systems`.


    """
    pdb_id = source_path.name

    ligand = Path(source_path, f'{pdb_id}_l_u.pdb').resolve()
    receptor = Path(source_path, f'{pdb_id}_r_u.pdb').resolve()

    results_p = Path(results_path, pdb_id).resolve()
    results_p.mkdir(exist_ok=True, parents=True)

    for scn_name, scn_func in scenarios.items():

        model_folder = Path(results_path, pdb_id, scn_name).resolve()
        model_folder.mkdir(exist_ok=True, parents=True)

        cfg_file = Path(model_folder, 'run.cfg').resolve()
        job_file = Path(model_folder, 'run.job').resolve()

        run_folder = Path(results_p, scn_name, 'run').resolve()

        cfg_str = scn_func(
            run_folder,
            receptor,
            ligand,
            Path(source_path, 'ambig.tbl'),
            )

        job_str = create_job_func(
            pdb_id,
            path=model_folder,
            conf_f=cfg_file,
            )

        cfg_file.write_text(cfg_str)
        job_file.write_text(job_str)

    return


def main(benchmark_path, results_path, job_sys='torque'):
    """
    Create configuration and job scripts for HADDOCK3 benchmarking.

    Developed for https://github.com/haddocking/BM5-clean

    Parameters
    ----------
    benchmark_path : Path
        The path to the benchmark models folder. In BM5-clean would be
        the 'HADDOCK-ready' folder.

    results_path : Path
        Where the results will be saved. A subfolder for each model in
        `benchmark_path` will be created.

    job_sys : str
        A key for `job_systems` dictionary.
    """
    log.info('*** Creating benchmark scripts')

    _ = (f for f in benchmark_path.glob('*') if _is_valid(f))
    source_folders = sorted(_)
    log.info(f'* creating {len(source_folders)} benchmark jobs')

    pe = partial(
        process_example,
        results_path=results_path,
        create_job_func=job_systems[job_sys],
        )

    for source_path in source_folders:
        pe(source_path)

    log.info('* done')


if __name__ == '__main__':
    sys.exit(maincli())
