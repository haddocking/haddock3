import sys
import toml
import json
from datetime import datetime
from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import InputGenerator
from haddock.modules.docking.it0 import RigidBody
from haddock.modules.worker.distribution import JobCreator
from haddock.modules.functions import *
import argparse

etc_folder = get_full_path('haddock', 'etc')
with open(f'{etc_folder}/default.json', 'r') as fh:
	default_recipes = json.load(fh)
fh.close()


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

	recipe_name = run_param['stage']['topology']['recipe']
	if recipe_name == 'default':
		recipe_name = default_recipes['topology']

	recipe = f'topology/template/{recipe_name}'
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
			                            output_pdb=True,
			                            output_psf=True,
			                            input_pdb=input_strc,
			                            input_psf=None,
			                            output_fname=output_strct)

			jobs.delegate(job_num=job_counter,
			              input_file_str=input_f)

			job_counter += 1

	cns = CNS()
	cns.run(jobs)
	output = retrieve_output(jobs)
	topology_dic = molecularize(output)
	return topology_dic


def run_it0(model_dic, run_param):
	print('\n++ Running it0')
	file_list = []
	supported_modules = []

	recipe_name = run_param['stage']['rigid_body']['recipe']
	if recipe_name == 'default':
		recipe_name = default_recipes['rigid_body']

	recipe = f'rigid_body/template/{recipe_name}'
	if not os.path.isfile(recipe):
		print('+ ERROR: Template recipe for topology not found')

	if '.cns' not in recipe:
		# Its a module, look for it
		if recipe not in supported_modules:
			print(f'+ ERROR: {recipe} not supported.')
			exit()

	else:
		# Its a HADDOCK recipe
		it0 = RigidBody()
		job_dic = it0.init(recipe, model_dic, run_param)
		complex_list = it0.run(job_dic)
		file_list = it0.output(complex_list)

	return file_list


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

	# 2 Dock

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
