from haddock.modules.cns.engine import CNS
from haddock.modules.cns.input import RecipeGenerator
from haddock.modules.docking.it0 import RigidBody
from haddock.modules.structure.utils import PDB
from haddock.modules.worker.distribution import JobCreator
from haddock.modules.functions import *
import config


def pre_process():
	molecule_dic = config.param_dic['input']['molecules']

	p = PDB()
	treated_dic = treat_ensemble(molecule_dic)

	# update data structure ?
	config.param_dic['input']['molecules'] = treated_dic

	for mol in treated_dic:
		pdb_list = treated_dic[mol]
		p.sanitize(pdb_list)

	return treated_dic


def generate_topology(molecule_dictionary):
	recipe_gen = RecipeGenerator()
	jobs = JobCreator()
	cns = CNS()

	job_counter = 1
	for mol in molecule_dictionary:
		for input_strc in molecule_dictionary[mol]:
			recipe = recipe_gen.generate(recipe_file=config.param_dic['topology_recipe'],
			                             molecule_id=mol,
			                             protonation_dic={},
			                             prefix_folder='begin/',
			                             output_name='')

			jobs.delegate(job_num=job_counter,
			              job_id='generate',
			              recipe_str=recipe,
			              input_model=input_strc)
			job_counter += 1

	cns.run(jobs)
	output = retrieve_output(jobs)
	topology_dic = molecularize(output)
	return topology_dic


def run_it0(model_dic):
	supported_modules = []
	recipe = config.param_dic['it0']['recipe']
	complex_list = None

	if '.cns' not in recipe:
		# Its a module, look for it
		if recipe not in supported_modules:
			print(f'+ ERROR: {recipe} not supported.')
			exit()
	# else:
	# 	exec =
	else:
		# Its a HADDOCK recipe

		it0 = RigidBody(model_dic)
		it0.init()
		it0.run()
		complex_list = it0.retrieve()
	#
	# recipe_gen = RecipeGenerator()
	# jobs = JobCreator()
	# cns = CNS()

	return complex_list


if __name__ == '__main__':
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
