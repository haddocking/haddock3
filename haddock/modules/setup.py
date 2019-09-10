import os
import sys
import json
from utils.files import get_full_path
from haddock.modules.worker.recipe import RecipeComposer


class Setup:

	def __init__(self, setup_dic):
		self.setup_dic = setup_dic
		self.protocol_path = get_full_path('haddock', 'protocols')

		run_id = self.setup_dic['identifier']['run']
		self.run_dir = None
		if type(run_id) == int:
			self.run_dir = f"run{run_id}"
		elif type(run_id) == str:
			run_id = run_id.replace(' ', '-')
			self.run_dir = f"run-{run_id}"
		else:
			print('+ ERROR: Run identifier can only be string or integer')
			exit()

	def prepare_folders(self):
		""" Create folder structure and copy significant input files """

		if len(self.setup_dic['molecules']) >= 20:
			print('+ ERROR: Too many molecules')
			exit()

		# move input molecules to correct places
		if not os.path.isdir(self.run_dir):
			os.system(f'mkdir {self.run_dir}')
		else:
			print(f'+ WARNING: {self.run_dir} already present')
		# Exit?

		data_dir = f'{self.run_dir}/data'
		if not os.path.isdir(data_dir):
			os.system(f'mkdir {data_dir}')

		os.system(f'cp {sys.argv[1]} {self.run_dir}/data/')

		for mol_id in self.setup_dic['molecules']:
			molecule = self.setup_dic['molecules'][mol_id]

			if molecule == mol_id + '.pdb':
				print("+ ERROR: Name mol_X.pdb not supported, please rename your molecule.")
				exit()

			# ensembles will be treated later
			os.system(f'cp {molecule} {self.run_dir}/data/{mol_id}_1.pdb')

			self.setup_dic['molecules'][mol_id] = f'{self.run_dir}/data/{mol_id}_1.pdb'

		try:
			ambig_fname = self.setup_dic['restraints']['ambig']
			os.system(f'cp {ambig_fname} {self.run_dir}/data/ambig.tbl')
		except KeyError:
			pass

		try:
			unambig_fname = self.setup_dic['restraints']['unambig']
			os.system(f'cp {unambig_fname} {self.run_dir}/data/unambig.tbl')
		except KeyError:
			pass

		try:
			hbond_fname = self.setup_dic['restraints']['hbond']
			os.system(f'cp {hbond_fname} {self.run_dir}/data/hbond.tbl')
		except KeyError:
			pass
		try:
			dihe_fname = self.setup_dic['restraints']['dihedrals']
			os.system(f'cp {dihe_fname} {self.run_dir}/data/dihe.tbl')
		except KeyError:
			pass

		# FIXME: What should this function return?
		return True

	def configure_recipes(self):
		""" Prepare recipes with parameters from setup dictionary """

		# default_recipes = {'topology': 'generate-topology.cns',
		#                    'rigid_body': 'it0.cns'}

		stage_dic = self.setup_dic['stage']
		for stage in stage_dic:
			recipe_id = stage_dic[stage]['recipe']
			# check if recipe is valid
			# if recipe_id == 'default':
			# 	recipe_id = default_recipes[stage]

			recipe_params_file = self.protocol_path + '/' + recipe_id.split('.')[0] + '.json'
			if os.path.isfile(recipe_params_file):
				with open(recipe_params_file, 'r') as f:
					default_parameter_dic = json.load(f)
				f.close()
			else:
				print(f'+ ERROR: Default parameters not found for {recipe_id}')
				exit()

			custom_parameter_dic = dict([(a, stage_dic[stage][a]) for a in stage_dic[stage] if a != 'recipe'])

			# TODO: Save the custom json file alongside the template
			# parameter_dic = self.merge_parameters(default_parameter_dic, custom_parameter_dic)

			# Generate the recipe and place it in the appropriate place
			r = RecipeComposer(recipe_id, default_parameter_dic, custom_parameter_dic)
			complete_recipe = r.compose()

			if not os.path.isdir(f'{self.run_dir}/{stage}'):
				os.mkdir(f'{self.run_dir}/{stage}')
				if not os.path.isdir(f'{self.run_dir}/{stage}/template'):
					os.mkdir(f'{self.run_dir}/{stage}/template')

			with open(f'{self.run_dir}/{stage}/template/{recipe_id}', 'w') as fh:
				fh.write(complete_recipe)
			fh.close()

		# FIXME: What should this function return?
		return True
