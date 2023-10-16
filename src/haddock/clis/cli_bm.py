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
<output dir> and search for the tags, for example::

    find . -name AVAILABLE
    find . -name RUNNING
    find . -name DONE
    find . -name FAIL

Jobs are identified as FAIL if there are messages in the stderr file.

Finally, a daemon job file is created to facilitate the usage of the
daemon without directly using the `haddock3-dmn` client.

Usage::

    haddock3-bm -h
    haddock3-bm <BM dir> <output dir> [OPTIONS]
    haddock3-bm <BM dir> <output dir> --workload-manager <option>
    haddock3-bm <BM dir> <output dir> --workload-manager <option> -n <num cores>
    haddock3-bm <BM dir> <output dir> --workload-manager <option> -n <num cores> -td

A `BM folder` is a folder with the characteristics of:
    https://github.com/haddocking/BM5-clean

For more information read our benchmark tutorial at `docs/benchmark.tut`
in HADDOCK3 repository site: https://github.com/haddocking/haddock3
"""  # noqa: E501
import argparse
import shutil
import string
import sys
from functools import partial
from pathlib import Path

from haddock import log
from haddock.core.typing import (
    Any,
    ArgumentParser,
    Callable,
    FilePath,
    Namespace,
    Union,
)
from haddock.libs.libhpc import create_job_header_funcs


# first character allowed for benchmark test cases, we use digits and
# upper cases because we consider BM test cases folders are named after
# their PDBID code.
capital_and_digits = tuple(string.digits + string.ascii_uppercase)


# client helper functions
def _dir_path(path: FilePath) -> Path:
    path = Path(path)
    if not path.is_dir():
        _p = str(path.resolve())
        raise argparse.ArgumentTypeError(f"{_p!r} is not a directory.")
    return path


def _is_valid(
    f: Path, cap_and_dig: Union[str, tuple[str, ...]] = capital_and_digits
) -> bool:
    """Assert if directory is a valid model directory."""
    _is_valid = f.name.startswith(cap_and_dig) and f.is_dir()
    return _is_valid


def get_conda_path() -> Path:
    """Get conda source path."""
    return Path(Path(sys.executable).parents[3], "etc", "profile.d", "conda.sh")


def create_cfg_test_daemon(
    run_dir: FilePath,
    receptor_f: FilePath,
    ligand_f: FilePath,
    *ignore: Any,
    **everythinelse: Any,
) -> str:
    """
    Create HADDOCK3 configuration file that only generates the topology.

    This function is usefull to use test the benchmark daemon.
    """
    cfg_str = f"""
run_dir = {str(run_dir)!r}
ncores = 2

molecules = [
    {str(receptor_f)!r},
    {str(ligand_f)!r}
    ]

[topoaa]

"""
    return cfg_str


# FIXME: This should not be hardcoded here
def create_cfg_ti(
    run_dir: FilePath,
    receptor_f: FilePath,
    ligand_f: FilePath,
    ambig_f: FilePath,
    target_f: FilePath,
) -> str:
    """
    Create HADDOCK3 configuration file for the first scenario.

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
    cfg_str = f"""
run_dir = {str(run_dir)!r}
ncores = 48

molecules = [
    {str(receptor_f)!r},
    {str(ligand_f)!r}
    ]

[topoaa]

[rigidbody]
ambig_fname = {str(ambig_f)!r}
sampling = 1000
noecv = false

[caprieval]
reference = {str(target_f)!r}

[seletop]
select = 200

[flexref]
ambig_fname = {str(ambig_f)!r}
noecv = true

[caprieval]
reference = {str(target_f)!r}

[mdref]
ambig_fname = {str(ambig_f)!r}
noecv = true

[caprieval]
reference = {str(target_f)!r}
"""
    return cfg_str


def create_job(
    create_job_header: Callable[..., str],
    create_job_body: Callable[[FilePath, FilePath, FilePath], str],
    create_job_tail: Callable[[FilePath, FilePath, FilePath], str],
    job_name_prefix: str,
    scenario_name: str,
    job_name_suffix: str,
    queue_name: str,
    ncores: int,
    work_dir: Path,
    run_dir: Path,
    #
    config_file: FilePath,
) -> str:
    """
    Create the job file.

    The jobs is created by assembling three parts: the job header,
    the body, and the final tail (post execution process).

    The different parameters will be injected in the respective job
    creation functions.

    Parameters
    ----------
    create_job_header : callable
        The function that will create the header.

    create_job_body : callable
        The function that will create the job body.

    create_job_tail: callable
        The function that will create the job tail.

    job_name_prefix : str
        A prefix for the job name. Normally this is the name of the job
        test case, for example the PDB ID.
        Injected in `create_job_header`.

    scenario_name : str
        The name of the benchmark scenario.
        Injected in `create_job_header`.

    job_name_suffix : str
        An additional suffix for the job name. Normally, `BM5`.
        Injected in `create_job_header`.

    queue_name : str
        The name of the queue. Injected in `create_job_header`.

    ncores : int
        The number of cpu cores to use in the jobs. Injected in
        `create_job_header`.

    work_dir : pathlib.Path
        The working dir of the example. That is, the directory where
        `input`, `jobs`, and `logs` reside. Injected in `create_job_header`.

    run_dir : pathlib.Path
        The running directory of the scenario.

    config_file : pathlib.Path
        Path to the scenario configuration file.
        Injected in `create_job_body`.

    Returns
    -------
    str
        The job file in the form of string.
    """
    # create job header
    job_name = f"{job_name_prefix}-{scenario_name}-{job_name_suffix}"
    std_out = str(Path("logs", "haddock.out"))
    std_err = str(Path("logs", "haddock.err"))

    job_header = create_job_header(
        job_name,
        work_dir=work_dir,
        stdout_path=std_out,
        stderr_path=std_err,
        queue=queue_name,
        ncores=ncores,
    )

    available_flag = str(Path(run_dir, "AVAILABLE"))
    running_flag = str(Path(run_dir, "RUNNING"))
    done_flag = str(Path(run_dir, "DONE"))
    fail_flag = str(Path(run_dir, "FAIL"))

    job_body = create_job_body(available_flag, running_flag, config_file)

    job_tail = create_job_tail(std_err, done_flag, fail_flag)

    return job_header + job_body + job_tail


def setup_haddock3_job(
    available_flag: FilePath, running_flag: FilePath, conf_f: FilePath
) -> str:
    """
    Write body for the job script.

    Parameters
    ----------
    available_flag : pathlib.Path
        Path where to generate the `AVAILABLE` file tag. Relative to the
        job header `workdir`.

    running_flag : pathlib.Path
        Path where to generate the `RUNNING` file tag. Relative to the
        job header `workdir`.

    conf_f : Path or str
        Path to the configuration file. Relative to the job header
        `workdir`.

    Returns
    -------
    str
        The body of the job file.
    """
    conda_sh = get_conda_path()

    job_body = f"""
source {str(conda_sh)}
conda activate haddock3

rm {str(available_flag)}
touch {str(running_flag)}
haddock3 {str(conf_f)}
rm {str(running_flag)}
"""
    return job_body


def process_job_execution_status(
    stderr: FilePath, done_tag: FilePath, fail_tag: FilePath
) -> str:
    """
    Add execution status tail to job.

    If job completes properly adds `DONE`, else adds `FAIL`.

    Parameters
    ----------
    done_tag : pathlib.Path
        Path where to generate the `DONE` file tag. Relative to the
        job header `workdir`.

    fail_flag : pathlib.Path
        Path where to generate the `FAIL` file tag. Relative to the
        job header `workdir`.
    """
    tail = f"""
if [ -s {stderr} ]; then
    # the file is not empty
    touch {str(fail_tag)}
else
    # the file is empty
    touch {str(done_tag)}
fi
"""
    return tail


def process_target(
    source_path: Path,
    result_path: FilePath,
    create_job_func: Callable[..., str],
    scenarios: dict[str, Callable[..., str]],
) -> None:
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
    input_p = Path(root_p, "input")
    input_p.mkdir(exist_ok=True, parents=True)

    # copy the required files to the `input` path
    ligand = shutil.copy(Path(source_path, f"{pdb_id}_l_u.pdb"), input_p)
    receptor = shutil.copy(Path(source_path, f"{pdb_id}_r_u.pdb"), input_p)
    ambig_tbl = shutil.copy(Path(source_path, "ambig.tbl"), input_p)
    target = shutil.copy(
        Path(source_path, "ana_scripts", "target.pdb"), input_p
    )  # noqa: E501

    # make all paths relative
    ligand = Path(ligand).relative_to(root_p)
    receptor = Path(receptor).relative_to(root_p)
    ambig_tbl = Path(ambig_tbl).relative_to(root_p)
    target = Path(target).relative_to(root_p)

    # Define the job folder
    #  all job files should be here
    job_p = Path(root_p, "jobs")
    job_p.mkdir(exist_ok=True, parents=True)

    # Define the logs folder
    log_p = Path(root_p, "logs")
    log_p.mkdir(exist_ok=True, parents=True)

    # for each scenario...
    for scn_name, scn_func in scenarios.items():
        # creates a scenario run folder
        run_folder = Path(result_path, pdb_id, f"run-{scn_name}").resolve()
        run_folder.mkdir(parents=True, exist_ok=True)

        # creates the AVAILABLE tag
        available_file = Path(run_folder, "AVAILABLE")
        available_file.touch()

        run_folder_rel = run_folder.relative_to(root_p)

        # the actual scenario results will be inside the `run` folder
        # which is inside the `run_folder` itself.
        run_job_folder = Path(run_folder_rel, "run")

        cfg_file = Path(input_p, f"{scn_name}.cfg")
        job_file = Path(job_p, f"{scn_name}.job")

        # creates the HADDOCK3 configuration file for the scenario
        cfg_str = scn_func(
            run_job_folder,
            receptor,
            ligand,
            ambig_tbl,
            target,
        )

        # creates the job for the scenario
        job_str = create_job_func(
            job_name_prefix=pdb_id,
            scenario_name=scenarios_acronym[scn_name],
            work_dir=root_p,
            run_dir=run_folder_rel,
            config_file=cfg_file.relative_to(root_p),
        )

        # saves files to disk
        cfg_file.write_text(cfg_str)
        job_file.write_text(job_str)

    return


def make_daemon_job(
    create_job_func: Callable[..., str],
    workdir: FilePath,
    target_dir: FilePath,
    job_name: str = "HD3-dmn",
    stdout_path: FilePath = "daemon.out",
    stderr_path: FilePath = "daemon.err",
    queue: str = "short",
) -> str:
    """Make a daemon-ready job."""
    job_header = create_job_func(
        job_name=job_name,
        work_dir=workdir,
        stdout_path=stdout_path,
        stderr_path=stderr_path,
        queue=queue,
        ncores=1,
    )

    conda_sh = get_conda_path()

    job = f"""{job_header}

source {str(conda_sh)}
conda activate haddock3

haddock3-dmn {str(target_dir)}
"""
    return job


# helper dictionaries

# the different scenarios covered
benchmark_scenarios = {
    "true-interface": create_cfg_ti,
    # 'center-of-mass': create_cfg_scn_2,
}

test_scenarios = {
    "test-daemon": create_cfg_test_daemon,
}

scenarios_acronym = {
    "true-interface": "ti",
    "center-of-mass": "com",
    "test-daemon": "tdmn",
}

# prepares the command-line client arguments
ap = argparse.ArgumentParser(
    prog="HADDOCK3 benchmark setup.",
    description=__doc__,
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

ap.add_argument(
    "benchmark_path",
    help=(
        "Location of BM5 folder. "
        "This folder should have a subfolder for each model. "
        "We expected each subfolder to have the name of a PDBID. "
        "That is, folder starting with capital letters or numbers "
        "will be considered."
    ),
    type=_dir_path,
)

ap.add_argument(
    "output_path",
    help="Where the prepared jobs and the executed benchmark will be stored.",
    type=_dir_path,
)

ap.add_argument(
    "-wm",
    "--workload-manager",
    dest="workload_manager",
    help="The system where the jobs will be run. Default `slurm`.",
    choices=list(create_job_header_funcs.keys()),
    default="slurm",
)

ap.add_argument(
    "-td",
    "--test-daemon",
    dest="test_daemon",
    help=(
        "Creates configuration files just with the `topology` module. "
        "Useful to test the benchmark daemon. Default `False`."
    ),
    action="store_true",
)

ap.add_argument(
    "-qu",
    "--queue-name",
    dest="queue_name",
    help=("The name of the queue where to send the jobs. " "Default 'medium'."),
    default="medium",
)

ap.add_argument(
    "-n",
    "--ncores",
    help="Maximum number of processors to use in the jobs. Default: 48",
    default=48,
)

ap.add_argument(
    "-s",
    "--suffix",
    help=(
        "A common suffix for all jobs. Defaults to 'BM5'. Avoid using "
        "long names because the job-name has a limited amount of characters."
    ),
    default="BM5",
    type=str,
)


def _ap() -> ArgumentParser:
    return ap


# client helper functions
def load_args(ap: ArgumentParser) -> Namespace:
    """Load argument parser args."""
    return ap.parse_args()


def cli(ap: ArgumentParser, main: Callable[..., None]) -> None:
    """Command-line interface entry point."""
    cmd = load_args(ap)
    main(**vars(cmd))


def maincli() -> None:
    """Execute main client."""
    cli(ap, main)


def main(
    benchmark_path: Path,
    output_path: FilePath,
    workload_manager: str = "slurm",
    ncores: int = 48,
    queue_name: str = "medium",
    test_daemon: bool = False,
    suffix: str = "BM5",
) -> None:
    """
    Create configuration and job scripts for HADDOCK3 benchmarking.

    Developed for https://github.com/haddocking/BM5-clean

    The parameters defined here are the same as defined in the client
    arguments.

    This is the main function of the client. If you want to run the benchmark
    creation routine withOUT using the command line and instead importing its
    functionalities and setting it up from another pythong script, you should
    import this function.

    >>> from haddock.clis.cli_bm import main

    Parameters
    ----------
    benchmark_path : Path
        The path to the benchmark models folder. In BM5-clean would be
        the 'HADDOCK-ready' folder.

    output_path : Path
        Where the results will be saved. A subfolder for each model in
        `benchmark_path` will be created.

    workload_manager : str
        A key for `create_job_header_funcs` dictionary. These relate to
        the queue managing software installed in your system. Examples
        are 'slurm' and 'torque'.

    ncores : int
        The number of CPUs to use in the created jobs.

    queue_name : str
        The name of the queue where the jobs will be sent. This depends
        on your system configurations.

    test_daemon : bool
        If `True`, generates short jobs where only the `topology` will
        be created. This facilitates testing the `haddoc3 benchmark
        daemon`.

    suffix : str
        A common suffix for all jobs. Avoid using more than three chars.
    """
    log.info("*** Preparing Benchnark scripts")

    _ = (f for f in benchmark_path.glob("*") if _is_valid(f))
    source_folders = sorted(_)
    log.info(f"* Creating benchmark jobs for {len(source_folders)} targets")

    # which scenarios to use?
    _scenarios = test_scenarios if test_daemon else benchmark_scenarios

    # prepares a `create_job_func` with predefined parameters.
    _create_job_func = partial(
        create_job,
        create_job_header=create_job_header_funcs[workload_manager],
        create_job_body=setup_haddock3_job,
        create_job_tail=process_job_execution_status,
        queue_name=queue_name,
        ncores=ncores,
        job_name_suffix=suffix,
    )

    # prepares a `process_target` function with predefined parameters
    # before entering the for loop
    pe = partial(
        process_target,
        result_path=output_path,
        create_job_func=_create_job_func,
        scenarios=_scenarios,
    )

    # for each benchmark test case...
    for source_path in source_folders:
        pe(source_path)

    # makes a job to run the daemon as a job.
    dmn_job = make_daemon_job(
        create_job_header_funcs[workload_manager],
        Path.cwd(),
        output_path,
        queue=queue_name,
    )

    Path(output_path, "hd3-daemon.job").write_text(dmn_job)

    log.info("* done")

    return


if __name__ == "__main__":
    sys.exit(maincli())  # type: ignore
