from pathlib import Path
import argparse
import toml
import os
import subprocess
import logging
import shutil
import sys

bm_log = logging.getLogger('bm_log')
bm_log.setLevel(logging.INFO)
ch = logging.StreamHandler()
formatter = logging.Formatter(' %(asctime)s L%(lineno)d'
                              ' %(levelname)s %(message)s')
ch.setFormatter(formatter)
bm_log.addHandler(ch)


def load_workflows(workflow_path):
    """Load the workflow toml and return a dictionary"""
    workflow_dic = {}
    for workflow in Path(workflow_path).glob('*.toml'):
        with open(workflow) as fh:
            workflow_dic[workflow.stem] = ''.join(fh.readlines())
    return workflow_dic


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--setup",
        help="Setup only",
        action="store_true",
        default=False)
    parser.add_argument(
        "--force",
        help="Force re-running if a previous execution is found",
        action="store_true",
        default=False)
    parser.add_argument(
        "--skip",
        help="Skip the execution if a previous execution is found",
        action="store_true",
        default=False)
    parser.add_argument(
        "input",
        help="Input file containing instruction to run the benchmark (.toml)")

    args = parser.parse_args()

    if args.force and args.skip:
        bm_log.error('--force and --skip are mutually exclusive.')
        sys.exit()

    params = toml.load(args.input)

    bm5 = Path(params['bm5_path'], 'HADDOCK-ready')
    output_path = Path(params['output_path'])
    to_be_executed = params

    workflow_dic = load_workflows(params['workflow_path'])

    # Setup the benchmark
    bm_log.info(f'Setting up the benchmark at {output_path}')
    folders_to_skip = ['ana_scripts', 'data', 'scripts']
    for folder in bm5.glob('*'):
        if folder.is_dir() and folder.name not in folders_to_skip:
            bm_log.debug(f' Working on {folder.name}')
            # create the .toml for the different workflows
            receptor_fname = f'{folder.name}_r_u.pdb'
            ligand_fname = f'{folder.name}_l_u.pdb'
            ambig_fname = 'ambig.tbl'
            reference_folder = Path(folder, 'ana_scripts')

            target_path = Path(output_path, folder.name)
            target_path.mkdir(exist_ok=True, parents=True)

            shutil.copy(Path(folder, receptor_fname), target_path)
            shutil.copy(Path(folder, ligand_fname), target_path)
            shutil.copy(Path(folder, ambig_fname), target_path)
            shutil.copytree(reference_folder,
                            Path(target_path, 'ana_scripts'),
                            dirs_exist_ok=True)

            for workflow in workflow_dic:
                if workflow not in params['workflows']:
                    continue
                input_toml = ''
                input_toml += f'run_dir = "run-{workflow}"' + os.linesep
                input_toml += (f'molecules = [ "{receptor_fname}",'
                               f' "{ligand_fname}" ]' + os.linesep)
                input_toml += workflow_dic[workflow] + os.linesep

                input_f = Path(target_path, f'{workflow}.toml')
                with open(input_f, 'w') as inp_fh:
                    inp_fh.write(input_toml)

    bm_log.info('Benchmark setup complete')
    if args.setup:
        bm_log.info('Running with --setup, we are done.')
        sys.exit()

    # Run the benchmark
    bm_log.info('Running benchmark')
    for folder in output_path.glob('*'):

        os.chdir(folder)

        for workflow in folder.glob('*.toml'):
            output_f = Path(f'{workflow.stem}.out')
            error_f = Path(f'{workflow.stem}.err')

            bm_log.info(f'Running {folder.name} - {workflow.stem} workflow')

            if output_f.exists() or error_f.exists():
                _msg = (f' Previous execution detected at'
                        f' {folder}/run-{workflow.stem}')
                bm_log.warning(_msg)

                if args.force:
                    bm_log.warning(f'  Removing run-{workflow.stem}')
                    shutil.rmtree(f'run-{workflow.stem}')
                    bm_log.warning(f'  Re-running run-{workflow.stem}')

                elif args.skip:
                    bm_log.info('  Skipping...')
                    continue

                else:
                    _msg = '  Clean the directory or run with --force/--skip'
                    bm_log.error(_msg)
                    sys.exit()

            cmd = f'haddock3 {workflow}'
            p = subprocess.Popen(cmd,
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 shell=True)
            out, err = p.communicate()

            with open(output_f, 'w') as out_fh:
                out_fh.write(out.decode('utf8'))

            with open(error_f, 'w') as err_fh:
                err_fh.write(err.decode('utf8'))

            p.kill()

            if 'error' in err.decode('utf-8').lower():
                bm_log.warning(f'  Failed {folder.name} - {workflow.stem}')
                bm_log.error(f'Benchmarking failed, check {error_f.resolve()}')
                sys.exit()
            else:
                bm_log.info(f'  Completed {folder.name} - {workflow.stem}')
