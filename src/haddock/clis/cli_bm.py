"""
Prepare HADDOCK3 benchmark configuration files and job scripts.

Creates HADDOCK3 configuration files and job files.

USAGE:
    haddock3-bm <BM folder> <results folder> --job-sys [OPTION]

A `BM folder` is a folder with the characteristics of:
    https://github.com/haddocking/BM5-clean
"""
import argparse
import string
import sys
from pathlib import Path


capital_and_digits = tuple(string.digits + string.ascii_uppercase)


def _dir_path(path):
    path = Path(path)
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f'{path!r} is not a directory.')
    return path


def _is_valid(f, cap_and_dig=capital_and_digits):
    """Assert if directory is a valid model directory."""
    _is_valid = \
        f.name.startswith(cap_and_dig) \
        and f.is_dir()
    return _is_valid



def create_cfg(run_dir, receptor_f, ligand_f, ambig_f):
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
sampling = 1000
noecv = false

[flexref]
ambig = {str(ambig_f)!r}
cool1_steps = 500
#sampling = 400

[mdref]
ambig = {(ambig_f)!r}
#sampling = 400
cool1_steps = 10
noecv = false
"""
    return cfg_str


def create_torque_job(job_name, **job_params):
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
#PBS -o haddock.out
#PBS -e haddock.err
"""
    return job_setup(header, **job_params)


def create_slurm_job(job_name, **job_params):
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
#SBATCH -o haddock.out
#SBATCH -e haddock.err
"""
    return job_setup(header, **job_params)


def job_setup(header, path, conf_f):
    job = \
"""{header}

source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate haddock3

touch RUNNING
haddock3 {str(cfg_f)}
rm RUNNING

if ( -f {path}run-ti/ss.stats ) then
    touch DONE
else
    touch FAIL
endif"""
    return job


job_system = {
    'torque': create_torque_job,
    'slurm': create_slurm_job,
    }


ap = argparse.ArgumentParser(description='Setup HADDOCK3 benchmark.')
ap.add_argument(
    "benchmark_path",
    help=(
        'Location of BM5 folder. '
        'This folder should have a subfolder for each model.'
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
    choices=list(job_system.keys()),
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


def main():
    """Run main logic."""
    cmd = ap.parse_args()

    _ = (f for f in cmd.benchmark_path.glob('*') if _is_valid(f))
    source_folders = sorted(_)

    for source_path in source_folders:
        pdb_id = source_path.name

        ligand = Path(source_path, f'{pdb_id}_l_u.pdb').resolve()
        receptor = Path(source_path, f'{pdb_id}_r_u.pdb').resolve()

        results_p = Path(cmd.results_path, pdb_id).resolve()
        results_p.mkdir(exist_ok=True, parents=True)

        cfg_file = Path(cmd.results_path, pdb_id, 'run.cfg').resolve()
        job_file = Path(cmd.results_path, pdb_id, 'run.job').resolve()

        #
        cfg_str = create_cfg(
            results_p.resolve(),
            receptor,
            ligand,
            Path(source_path, 'ambig.tbl'),
            )

        job_str = job_system[cmd.job_sys](results_p, cfg_file, pdb_id)
        #

        cfg_file.write_text(cfg_str)
        job_file.write_text(job_str)


if __name__ == '__main__':
    sys.exit(maincli())
