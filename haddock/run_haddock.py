import glob
import random
import sys
import toml
import os
from datetime import datetime
from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import InputGenerator
from haddock.modules.docking.it0 import RigidBody
# from haddock.modules.structure.utils import PDB
from haddock.modules.worker.distribution import JobCreator
from haddock.modules.functions import *


# import config


def greeting():
	now = datetime.now().replace(second=0, microsecond=0)
	python_version = sys.version
	return (f'''##############################################
#                                            #
#         Starting HADDOCK v3.0beta1         #
#                                            #
#             EXPERIMENTAL BUILD             #
#                                            #
##############################################

Starting HADDOCK on {now}

HADDOCK version: 3.0 beta 1
Python {python_version}
''')


def generate_topology(mol_dic, run_param):
	print('++ Generating topologies')

	recipe = 'topology/template/' + run_param['stage']['topology']['recipe']
	if not os.path.isfile(recipe):
		print('+ ERROR: Template recipe for topology not found')

	topo_gen = InputGenerator(recipe_file=recipe,
	                          input_folder='topology')

	jobs = JobCreator(job_id='generate',
	                  job_folder='topology')

	job_counter = 1
	for mol in mol_dic:
		for input_strc in mol_dic[mol]:
			output_strct = input_strc.split('/')[1].split('.')[0]
			input_f = topo_gen.generate(protonation_dic={},
			                            pdb=True,
			                            psf=True,
			                            input_data=input_strc,
			                            output_fname=output_strct)

			jobs.delegate(job_num=job_counter,
			              recipe_str=input_f)

			job_counter += 1

	cns = CNS()
	cns.run(jobs)
	output = retrieve_output(jobs)
	topology_dic = molecularize(output)
	return topology_dic


def run_it0(model_dic, run_param):
	print('\n++ Running it0')
	supported_modules = []
	# recipe = config.param_dic['it0']['recipe']
	recipe = 'topology/template/' + run_param['stage']['rigid_body']['recipe']
	# print(f'+ Recipe: {recipe}')
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


if __name__ == '__main__':
	print(greeting())

	# 0 Identify input
	#

	# TODO: Add USAGE
	run_parameters = toml.load('data/run.toml')

	molecules = get_begin_molecules('data/')


	# 1 Generate topologies
	begin_models = generate_topology(molecules, run_parameters)

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

	rigid_complexes = run_it0(begin_models, run_parameters)
	# rigid_complexes = run_lightdock(begin_models)
	# semi_flexible_complexes = run_it1(begin_models)
	# water_refined_complexes = run_itw(begin_models)
	# md_models = run_md(models)

	print(bye())
