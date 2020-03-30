# Setup the Benchmark
import argparse
import os
import glob

def create_toml(receptor_f, ligand_f, pdb_name):
    toml_str = f'''title = "HADDOCK3 Setup file"
[molecules]
mol1 = '{receptor_f}'
mol2 = '{ligand_f}'
[restraints]
ambig = 'ambig.tbl'
[identifier]
run = 'ti'
[execution_parameters]
scheme = 'parallel'
nproc = 48
reference.analysis  = '/home/rodrigo/projects/BM5-clean/HADDOCK-ready/{pdb_name}/ana_scripts/target.pdb'
# Stage specific parameters
[stage]
[stage.topology]
recipe='default'
[stage.rigid_body]
recipe='rigid-body-minimal.cns'
sampling = 1000
params.noecv = false
[stage.semi_flexible]
recipe='default'
sampling = 400
params.noecv = false
[stage.water_refinement]
recipe='default'
sampling = 400
params.noecv = false'''
    return toml_str


def create_job(path, toml_f, job_name):
    job = f'''#!/usr/bin/env tcsh
#PBS -N {job_name}-BM5
#PBS -q medium
#PBS -l nodes=1:ppn=48
#PBS -S /bin/tcsh
setenv HADDOCK3 /home/rodrigo/software/haddock3
setenv PYTHONPATH $PYTHONPATH\:$HADDOCK3
setenv CNS_EXE /home/software/science/cns/cns_solve_1.31-UU/intel-x86_64bit-linux/bin/cns
echo $PYTHONPATH
echo $HADDOCK3
echo $CNS_EXE
cd {path}
touch RUNNING
/home/rodrigo/miniconda3/bin/python /home/rodrigo/software/haddock3/bin/haddock3.py {toml_f} > haddock.out
rm RUNNING
if ( -f {path}run-ti/ss.stats ) then
    touch DONE
else
    touch FAIL
endif'''
    return job


if __name__ == '__main__':
    ''' Setup BM5 for Haddock3 benchmarking '''

    parser = argparse.ArgumentParser(description='Setup your HADDOCK run')
    parser.add_argument("benchmark_path", help="Location of BM5")
    args = parser.parse_args()

    benchmark_path = args.benchmark_path

    if not os.path.isdir(benchmark_path):
        print('Could not find benchmark patch')
        exit()

    blacklisted_folders = ['ana_scripts','data', 'scripts']
    target_folder = [f for f in glob.glob(f'{benchmark_path}/HADDOCK-ready/*/') if f not in blacklisted_folders]

    for target_path in target_folder:
        pdb_id = target_path.split('/')[-2]
        ligand = f'{pdb_id}_l_u.pdb'
        receptor = f'{pdb_id}_r_u.pdb'
        toml_file = f'{target_path}{pdb_id}.toml'
        job_file = f'{target_path}{pdb_id}.job'
        toml_str = create_toml(receptor, ligand, pdb_id)
        with open(toml_file, 'w') as fh:
            fh.write(toml_str)

        job_str = create_job(target_path, toml_file, pdb_id)
        with open(job_file, 'w') as fh:
            fh.write(job_str)
