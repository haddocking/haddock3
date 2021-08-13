from pathlib import Path
import argparse
import toml
import os
import subprocess


def load_workflows(workflow_path):
    """Load the workflow toml and return a dictionary"""
    workflow_dic = {}
    for workflow in Path(workflow_path).glob('*.toml'):
        with open(workflow) as fh:
            workflow_dic[workflow.stem] = ''.join(fh.readlines())
    return workflow_dic


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")

    args = parser.parse_args()
    params = toml.load(args.input)

    bm5 = Path(params['bm5_path'], 'HADDOCK-ready')
    ouput_path = Path(params['output_path'])
    to_be_executed = params

    workflow_dic = load_workflows(params['workflow_path'])

    # Setup the benchmark
    folders_to_skip = ['ana_scripts', 'data', 'scripts']
    for folder in bm5.glob('*'):
        if folder.is_dir() and folder.name not in folders_to_skip:
            # create the .toml for the different workflows
            receptor_fname = f'{folder.name}_r_u.pdb'
            ligand_fname = f'{folder.name}_l_u.pdb'

            for workflow in workflow_dic:
                if workflow not in params['workflows']:
                    continue
                input_toml = ''
                input_toml += f'run_dir = "run-{workflow}"' + os.linesep
                input_toml += (f'molecules = [ "{receptor_fname}",'
                               f' "{ligand_fname}" ]' + os.linesep)
                input_toml += workflow_dic[workflow] + os.linesep

                input_f = Path(folder, f'{workflow}.toml')
                with open(input_f, 'w') as inp_fh:
                    inp_fh.write(input_toml)

    # Run the benchmark
    for folder in bm5.glob('*'):
        if folder.is_dir() and folder.name not in folders_to_skip:
            os.chdir(folder)
            for workflow in workflow_dic:
                cmd = f'haddock3 {workflow}.toml'
                subprocess.call(cmd, shell=True)
