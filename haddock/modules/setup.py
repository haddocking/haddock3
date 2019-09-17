import os
import sys
import json
import string

from haddock.modules.structure.utils import PDB
from utils.files import get_full_path
from haddock.modules.worker.recipe import RecipeComposer

etc_folder = get_full_path('haddock', 'etc')


class Setup:

	def __init__(self, setup_dic):
		self.setup_dic = setup_dic
		self.protocol_path = get_full_path('haddock', 'protocols')

		with open(f'{etc_folder}/default.json', 'r') as fh:
			self.default_recipes = json.load(fh)
		fh.close()

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
		""" Create folder structure and copy/edit significant input files """

		if len(self.setup_dic['molecules']) >= 20:
			print('+ ERROR: Too many molecules')
			exit()

		# move input molecules to correct places
		if not os.path.isdir(self.run_dir):
			os.system(f'mkdir {self.run_dir}')
		else:
			print(f'+ ERROR: {self.run_dir} already present')
			# exit()

		data_dir = f'{self.run_dir}/data'
		if not os.path.isdir(data_dir):
			os.system(f'mkdir {data_dir}')

		os.system(f'cp {sys.argv[1]} {self.run_dir}/data/run.toml')

		# Separate by mol and assign segids
		segid_dic = dict([(int(e.split('mol')[1]), {'mol': None, 'segid': None}) for e in self.setup_dic['molecules'] if 'mol' in e])
		for e in self.setup_dic['molecules']:
			if 'mol' in e:
				ident = int(e.split('mol')[1])
				molecule = self.setup_dic['molecules'][e]
				segid_dic[ident]['mol'] = molecule
			if 'segid' in e:
				ident = int(e.split('segid')[1])
				segid = self.setup_dic['molecules'][e]
				segid_dic[ident]['segid'] = segid

		# If segid has not been assigned, check if PDB already has one
		for e in segid_dic:
			structure = segid_dic[e]['mol']
			custom_segid = segid_dic[e]['segid']

			if not custom_segid:
				# Segid not defined in setup, check if it is already present
				molecule_chainseg = PDB.identify_chainseg(structure)

				if molecule_chainseg:
					# keep it
					_ = PDB.add_chainseg(structure, molecule_chainseg)

				if not molecule_chainseg:
					# define sequentially
					molecule_chainseg = string.ascii_uppercase[e-1]
					_ = PDB.add_chainseg(structure, molecule_chainseg)

			if custom_segid:
				_ = PDB.add_chainseg(structure, custom_segid)

		mol_dic = dict([(e, self.setup_dic['molecules'][e]) for e in self.setup_dic['molecules'] if 'mol' in e])
		for mol_id in mol_dic:
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

		stage_dic = self.setup_dic['stage']
		for stage in stage_dic:
			recipe_id = stage_dic[stage]['recipe']
			# TODO: check if recipe is valid
			if recipe_id == 'default':
				recipe_id = self.default_recipes[stage]

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
