import random
import sys
from datetime import datetime
from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import RecipeGenerator
from haddock.modules.docking.it0 import RigidBody
from haddock.modules.structure.utils import PDB
from haddock.modules.worker.distribution import JobCreator
from haddock.modules.functions import *
import config


def logo():
	now = datetime.now().replace(second=0, microsecond=0)
	version = sys.version
	print(f'''##############################################
#                                            #
#         Starting HADDOCK v3.0beta1         #
#                                            #
#             EXPERIMENTAL BUILD             #
#                                            #                               
##############################################

Starting HADDOCK on {now}

HADDOCK version: 3.0 beta 1
Python {version}
''')


def pre_process():
	molecule_dic = config.param_dic['input']['molecules']

	p = PDB()

	# TODO: Edit chainID/segID
	treated_dic = treat_ensemble(molecule_dic)

	# update data structure ?
	config.param_dic['input']['molecules'] = treated_dic

	for mol in treated_dic:
		pdb_list = treated_dic[mol]
		p.sanitize(pdb_list)

	return treated_dic


def generate_topology(molecule_dictionary):
	print('++ Generating topologies')
	recipe_gen = RecipeGenerator()
	jobs = JobCreator()
	cns = CNS()

	job_counter = 1
	for mol in molecule_dictionary:
		for input_strc in molecule_dictionary[mol]:
			output_strct = input_strc.split('/')[1].split('.')[0]
			recipe = recipe_gen.generate(recipe_file=config.param_dic['topology_recipe'],
			                             molecule_id=mol,
			                             protonation_dic={},
			                             prefix_folder='begin/',
			                             pdb=input_strc,
			                             psf=None,
			                             output_name=output_strct)

			jobs.delegate(job_num=job_counter,
			              job_id='generate',
			              recipe_str=recipe,
			              job_folder='begin')
			job_counter += 1

	cns.run(jobs)
	output = retrieve_output(jobs)
	topology_dic = molecularize(output)
	return topology_dic


def run_it0(model_dic):
	print('\n++ Running it0')
	supported_modules = []
	recipe = config.param_dic['it0']['recipe']
	print(f'+ Recipe: {recipe}')
	complex_list = None

	if '.cns' not in recipe:
		# Its a module, look for it
		if recipe not in supported_modules:
			print(f'+ ERROR: {recipe} not supported.')
			exit()

	else:
		# Its a HADDOCK recipe
		it0 = RigidBody()
		job_dic = it0.init(model_dic)
		complex_list = it0.run(job_dic)
		file_list = it0.output(complex_list)

	# recipe_gen = RecipeGenerator()
	# jobs = JobCreator()
	# cns = CNS()

	return complex_list


def bye():
	salutations = ['Tot ziens!', 'Good bye!', 'Até logo!', 'Ciao!', 'Au revoir!', 'Adéu-siau!', 'Agur!']
	return '\n' + random.choice(salutations)


if __name__ == '__main__':
	logo()

	# 0 Initialize
	init()

	# 1 Pre-process initial structures
	mol_dic = pre_process()

	# 2 Generate topologies
	begin_models = generate_topology(mol_dic)

	# 3 Docking

	# Input:
	#  molecule dictionary, keys=molecule, values=list of tuples containing .psf and .pdb
	#   example:
	# {
	#  'mol1': [('begin/mol1_1.psf', 'begin/mol1_1.pdb')],
	#  'mol2': [('begin/mol2_1.psf', 'begin/mol2_1.pdb'),
	#           ('begin/mol2_2.psf', 'begin/mol2_2.pdb')]
	#  'mol3': [('begin/mol3_1.psf', 'begin/mol3_1.pdb')]
	#   }

	# Output:
	#  List of complexes, sorted according to ranking
	#   example:
	# ['complex_42.pdb', 'complex_10.pdb', 'complex_23.pdb']

	rigid_complexes = run_it0(begin_models)
	# semi_flexible_complexes = run_it1(begin_models)
	# water_refined_complexes = run_itw(begin_models)
	# md_models = run_md(models)

	print(bye())
