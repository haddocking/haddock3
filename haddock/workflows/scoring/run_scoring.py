# convert the analysis scripts to python
# import os
import sys
import toml
import json
from datetime import datetime
from haddock.workflows.scoring.analysis.ana import Ana
from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import InputGenerator
from haddock.modules.worker.distribution import JobCreator
from haddock.modules.functions import *
from haddock.run_haddock import generate_topology

etc_folder = get_full_path('haddock', 'etc')
with open(f'{etc_folder}/default.json', 'r') as fh:
    default_recipes = json.load(fh)
fh.close()


def run_analysis(pdb_l):
    ana = Ana(pdb_l)

    # ana.retrieve_structures()

    ana.extract_energies()

    ana.calculate_haddock_score()

    # ana.cluster(cutoff=0.75, threshold=1)
    # ana.cluster(cutoff=0.75, threshold=2)
    ana.cluster(cutoff=0.75, threshold=4)

    # ana.cluster(cutoff=0.65, threshold=1)
    # ana.cluster(cutoff=0.65, threshold=2)
    ana.cluster(cutoff=0.65, threshold=4)

    # ana.run_fastcontact()

    # ana.run_dfire()

    # ana.run_dockq()

    ana.output()


def greeting():
    now = datetime.now().replace(second=0, microsecond=0)
    python_version = sys.version
    return (f'''##############################################
#                                            #
#         Starting HADDOCK v3.0beta1         #
#                                            #
#             EXPERIMENTAL BUILD             #
#                                            #
#              Scoring Workflow              #
#                                            #
##############################################

Starting HADDOCK on {now}

HADDOCK version: 3.0 beta 1
Python {python_version}
''')


def score_models(run_param):
    print('++ Running scoring workflow')

    recipe_name = run_param['stage']['scoring']['recipe']
    if recipe_name == 'default':
        recipe_name = default_recipes['scoring']

    recipe = f'scoring/template/{recipe_name}'
    if not os.path.isfile(recipe):
        print('+ ERROR: Template recipe for scoring not found')

    scoring_gen = InputGenerator(recipe_file=recipe,
                                 input_folder='scoring')

    jobs = JobCreator(job_id='scoring',
                      job_folder='scoring')

    input_list = []
    # match pdbs and psfs
    for pdbf in glob.glob('topology/*pdb'):
        psff = pdbf.replace('.pdb', '.psf')
        if os.path.isfile(psff):
            input_list.append((pdbf, psff))

    job_counter = 1
    for mol in input_list:
        pdb_struct = mol[0]
        psf_struct = mol[1]
        output_struct = pdb_struct.split('/')[1].split('.')[0]
        input_f = scoring_gen.generate(protonation_dic={},
                                       output_pdb=True,
                                       output_psf=True,
                                       input_pdb=pdb_struct,
                                       input_psf=psf_struct,
                                       output_fname=output_struct)

        jobs.delegate(job_num=job_counter,
                      input_file_str=input_f)

        job_counter += 1

    cns = CNS()
    cns.run(jobs)
    output = retrieve_output(jobs)
    pdb_list = [e[0] for e in output.values()]
    return pdb_list


if __name__ == '__main__':

    print(greeting())

    run_f = 'data/run.toml'
    if not os.path.isfile(run_f):
        print('+ ERROR: data/run.toml not found, make sure you are in the correct folder.')
        exit()

    run_parameters = toml.load(run_f)
    molecules = get_begin_molecules('data/')

    # 1 Generate topologies
    begin_models = generate_topology(molecules, run_parameters)

    rescored = score_models(run_parameters)

    run_analysis(rescored)