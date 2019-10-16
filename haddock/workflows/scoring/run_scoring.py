# convert the analysis scripts to python
import os
import sys
from datetime import datetime
import toml

from haddock.modules.functions import get_begin_molecules
from haddock.workflows.scoring.analysis.ana import Ana
from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import RecipeGenerator
from haddock.modules.structure.utils import PDB
from haddock.modules.worker.distribution import JobCreator
# import haddock.workflows.scoring.config as config
import config
import glob


def run_scoring(ensamble_f, suffix):

	pdb = PDB()
	recipe_gen = RecipeGenerator()
	jobs = JobCreator()
	cns = CNS()

	# Prepare PDBs
	pdb.prepare(ensamble_f)

	# Generate the Recipe
	#  - Protonation dictionary
	#
	recipe = recipe_gen.generate(protonation_dic=pdb.protonation_dic, out_suffix=suffix)

	jobs.delegate(recipe, pdb.model_list)

	cns.run(jobs)

	return glob.glob(f'structures/*{suffix}.pdb')


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

	ana.run_dockq()

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


def score_models(mol_dic, run_params):
	pass


if __name__ == '__main__':

	print(greeting())

	run_f = 'data/scoring.toml'
	if not os.path.isfile(run_f):
		print('+ ERROR: data/run.toml not found, make sure you are in the correct folder.')
		exit()

	run_parameters = toml.load(run_f)
	models = get_begin_molecules('data/')

	# TODO: Refactor scoring CNS script

	#
	# input_f = config.param_dic['input']['ensemble']
	#
	# converted_pdb_list = run_scoring(input_f, '_conv')
	# run_analysis(converted_pdb_list)
