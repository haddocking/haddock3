"""
Prepare HADDOCK3 benchmark configuration files and job scripts.

Creates HADDOCK3 configuration files and job files. Details on each
parameter are explained in the `-h` menu.

There's also a test flag that generates jobs only with `topology`
creation. This feature helps testing the `haddock3-dmn` client.

The state of the jobs is identified by a file tag:
    - AVAILABLE
    - RUNNING
    - DONE
    - FAIL

At start, all jobs have the AVAILABLE tag, and this tag is upgraded as
the job completes. To know which jobs are in each state, navigate to the
<output dir> and grep the tags, for example:

    find . -name "AVAILABLE"
    find . -name "RUNNING"
    find . -name "DONE"
    find . -name "FAIL"

Jobs are identified as FAIL if there are messages in the stderr file.

USAGE:
    haddock3-bm -h
    haddock3-bm <BM dir> <output dir> --job-sys <option>
    haddock3-bm <BM dir> <output dir> --job-sys <option> -n <num cores>
    haddock3-bm <BM dir> <output dir> --job-sys <option> -n <num cores> -td

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


def create_job(
        create_job_header,
        create_job_body,
        create_job_tail,
        job_name_prefix,
        scenario_name,
        job_name_suffix,
        queue_name,
        ncores,
        work_dir,
        run_dir,
        #
        config_file,
        ):
    """."""
    # create job header
    job_name = f'{job_name_prefix}-{scenario_name}-{job_name_suffix}'
    std_out = str(Path('logs', 'haddock.out'))
    std_err = str(Path('logs', 'haddock.err'))

    job_header = create_job_header(
        job_name,
        work_dir=work_dir,
        stdout_path=std_out,
        stderr_path=std_err,
        queue=queue_name,
        ncores=ncores,
        )

    available_flag = str(Path(run_dir, 'AVAILABLE'))
    running_flag = str(Path(run_dir, 'RUNNING'))
    done_flag = str(Path(run_dir, 'DONE'))
    fail_flag = str(Path(run_dir, 'FAIL'))

    job_body = create_job_body(available_flag, running_flag, config_file)

    job_tail = create_job_tail(std_err, done_flag, fail_flag)

    return job_header + job_body + job_tail


def create_torque_header(
        job_name,
        work_dir,
        stdout_path,
        stderr_path,
        queue='medium',
        ncores=48,
        ):
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
#PBS -N {job_name}
#PBS -q {queue}
#PBS -l nodes=1:ppn={str(ncores)}
#PBS -S /bin/tcsh
#PBS -o {stdout_path}
#PBS -e {stderr_path}
#PBS -wd {work_dir}
"""
    return header


def create_slurm_header(
        job_name,
        work_dir,
        stdout_path,
        stderr_path,
        queue='medium',
        ncores=48,
        ):
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
#SBATCH -J {job_name}
#SBATCH -p {queue}
#SBATCH --nodes=1
#SBATCH --tasks-per-node={str(ncores)}
#SBATCH --output={stdout_path}
#SBATCH --error={stderr_path}
#SBATCH --workdir={work_dir}
"""
    return header


def setup_haddock3_job(available_flag, running_flag, conf_f):
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
        Path(sys.executable).parents[3],
        'etc', 'profile.d', 'conda.sh'
        )

    job_body = \
f"""
source {str(conda_sh)}
conda activate haddock3

rm {str(available_flag)}
touch {str(running_flag)}
haddock3 {str(conf_f)}
rm {str(running_flag)}
"""
    return job_body


def process_job_execution_status(stderr, done, fail):
    """Add execution status tail to job."""
    tail = \
f"""
if [ -s {stderr} ]; then
    # the file is not empty
    touch {str(fail)}
else
    # the file is empty
    touch {str(done)}
fi
"""
    return tail


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
        there are several queue systems available. See
        `create_job_header_funcs`.
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
        run_folder.mkdir(parents=True, exist_ok=True)

        available_file = Path(run_folder, 'AVAILABLE')
        available_file.touch()

        run_folder_rel = run_folder.relative_to(root_p)

        run_job_folder = Path(run_folder_rel, 'run')

        cfg_file = Path(input_p, f'{scn_name}.cfg')
        job_file = Path(job_p, f'{scn_name}.job')

        cfg_str = scn_func(
            run_job_folder,
            receptor,
            ligand,
            ambig_tbl,
            )

        job_str = create_job_func(
            job_name_prefix=pdb_id,
            scenario_name=scenarios_acronym[scn_name],
            work_dir=root_p,
            run_dir=run_folder_rel,
            config_file=cfg_file.relative_to(root_p),
            )

        cfg_file.write_text(cfg_str)
        job_file.write_text(job_str)

    return


create_job_header_funcs = {
    'torque': create_torque_header,
    'slurm': create_slurm_header,
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

ap = argparse.ArgumentParser(
    prog='HADDOCK3 benchmark setup.',
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
    )

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
    help="Where the prepared jobs and the executed benchmark will be stored.",
    type=_dir_path,
    )

ap.add_argument(
    '--job-sys',
    dest='job_sys',
    help='The system where the jobs will be run. Default `slurm`.',
    choices=list(create_job_header_funcs.keys()),
    default='slurm',
    )

ap.add_argument(
    '-td',
    '--test-daemon',
    dest='test_daemon',
    help=(
        'Creates configuration files just with the `topology` module. '
        'Useful to test the benchmark daemon. Default `False`.'
        ),
    action='store_true',
    )

ap.add_argument(
    '-qu',
    '--queue-name',
    dest='queue_name',
    help=(
        'The name of the queue where to send the jobs. '
        'Default \'medium\'.'
        ),
    default='medium',
    )

ap.add_argument(
    '-n',
    '--ncores',
    help='Maximum number of processors to use. Default: 48',
    default=48,
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


def main(
        benchmark_path,
        output_path,
        job_sys='slurm',
        ncores=48,
        queue_name='medium',
        test_daemon=False,
        ):
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
        A key for `create_job_header_funcs` dictionary. These relate to
        the queue managing software installed in your system. Examples
        are 'slurm' and 'torque'.

    ncores : int
        The number of CPUs to use in the created jobs.

    queue_name : str
        The name of the queue where the jobs will be sent. This depends
        on your system configurations.

    test-daemon : bool
        If `True`, generates short jobs where only the `topology` will
        be created. This facilitates testing the `haddoc3 benchmark
        daemon`.
    """
    log.info('*** Preparing Benchnark scripts')

    _ = (f for f in benchmark_path.glob('*') if _is_valid(f))
    source_folders = sorted(_)
    log.info(f'* Creating benchmark jobs for {len(source_folders)} targets')

    _scenarios = test_scenarios if test_daemon else benchmark_scenarios

    _create_job_func = partial(
        create_job,
        create_job_header=create_job_header_funcs[job_sys],
        create_job_body=setup_haddock3_job,
        create_job_tail=process_job_execution_status,
        queue_name=queue_name,
        ncores=ncores,
        job_name_suffix='BM5',
        )

    pe = partial(
        process_target,
        result_path=output_path,
        create_job_func=_create_job_func,
        scenarios=_scenarios
        )

    for source_path in source_folders:
        pe(source_path)

    log.info('* done')


if __name__ == '__main__':
    sys.exit(maincli())
