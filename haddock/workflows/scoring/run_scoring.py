import sys
import toml
import json
import argparse
from datetime import datetime

from bin.haddock3 import generate_topology
from haddock.modules.analysis.ana import Ana
from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import InputGenerator
from haddock.modules.setup import Setup
from haddock.modules.worker.distribution import JobCreator
from haddock.modules.functions import *
from haddock.version import CURRENT_VERSION

etc_folder = get_full_path('haddock', 'etc')
with open(f'{etc_folder}/default.json', 'r') as fh:
    default_recipes = json.load(fh)
fh.close()

start = None


def greeting():
    global start
    start = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
    python_version = sys.version
    print(f'''##############################################
#                                            #
#              Starting HADDOCK              #
#             EXPERIMENTAL BUILD             #
#                                            #
#              Scoring Workflow              #
#                                            #
##############################################

Starting HADDOCK {CURRENT_VERSION} on {start}

Python {python_version}
''')


def score_models(run_param):
    print('\n++ Running scoring workflow')

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
    # FIXME: Find a better way to retrieve molecules from topology
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


def preprocess_scoring(raw_molecule_dic):

    print('\n+ Pre-processing input molecules')

    p = PDB()

    # Split the ensamble in multiple models
    print('++ Treating ensambles')
    ensemble_dic = p.treat_ensemble(raw_molecule_dic)

    # Deal with chains (WIP)
    print('++ Analyzing chains')
    organized_dic = p.organize_chains(ensemble_dic)

    # Clean problematic parts
    print('++ Cleaning problematic parts')
    clean_molecule_dic = p.sanitize(organized_dic)

    # Deal with chainIDs and segIDs
    print('++ Adding segID (from chainID)')
    clean_pdb_list = []
    for pdb in clean_molecule_dic:
        clean_pdb = p.fix_id(pdb, priority='chain')
        clean_pdb_list.append(clean_pdb)

    return clean_pdb_list


def adieu():
    end = datetime.now().strftime('%d/%m/%Y %H:%M:%S')
    duration = datetime.strptime(end, '%d/%m/%Y %H:%M:%S') - datetime.strptime(start, '%d/%m/%Y %H:%M:%S')
    salut = bye()
    print(f'''
 Finished at {end} ({duration})

{salut}
''')


if __name__ == '__main__':

    # Initialize ======================================================================================================#

    greeting()

    parser = argparse.ArgumentParser(description='Setup your HADDOCK Scoring run')
    parser.add_argument("run_file", help="The run file containing the parameters of your scoring run (.toml)")
    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_usage()
        exit()

    # Setup & pre-process =============================================================================================#

    setup_dictionary = toml.load(args.run_file)
    run_id = setup_dictionary['identifier']['run']
    if type(run_id) == int:
        run_folder = f'run{run_id}'
    else:
        run_folder = f'run-{run_id}'

    if not os.path.isdir(run_folder):
        print(f'+ Setting up {run_folder}')
        s = Setup(setup_dictionary)
        s.prepare_folders()
        s.configure_recipes()
        # this one is not being used anywhere
        _ = preprocess_scoring(setup_dictionary['molecules'])

    print(f'\n+ Executing {run_folder}')
    os.chdir(run_folder)

    run_f = 'data/run.toml'
    if not os.path.isfile(run_f):
        print('+ ERROR: data/run.toml not found, make sure you are in the correct folder.')
        exit()

    run_parameters = toml.load(run_f)
    molecules = get_begin_molecules('data/')

    # Generate Topologies =============================================================================================#

    begin_models = generate_topology(molecules, run_parameters)

    # Run HADDOCK's CNS scoring recipe ================================================================================#

    rescored = score_models(run_parameters)

    # Analysis ========================================================================================================#

    ana = Ana(rescored)

    reference_structure = setup_dictionary['execution_parameters']['reference']['analysis']
    chain_num_reference = setup_dictionary['execution_parameters']['reference']['numbering']

    method = setup_dictionary['clustering']['method']
    detail_flag = setup_dictionary['clustering']['details']
    clustering_params = setup_dictionary['clustering']['params']

    # HADDOCK-Score
    ana.extract_energies()
    ana.calculate_haddock_score()

    # ana.match_renumber(chain_num_reference)

    # Clustering
    if method != 'fcc':
        print(f'+ ERROR: Clustering method {method} not supported')
        exit()

    for cutoff, threshold in clustering_params:
        ana.cluster(cutoff=cutoff, size=threshold)

    # Third-party
    ana.run_fastcontact()
    ana.run_dfire()
    ana.run_dockq(reference_structure)

    # Output ==========================================================================================================#
    ana.output()

    # Done! ===========================================================================================================#

    adieu()
