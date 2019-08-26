import json
import os
import random
import re
# import haddock.workflows.scoring.config as config
from utils.files import get_full_path
import config


class RecipeGenerator:

	def __init__(self):
		self.recipe = ''

	def generate(self, recipe_file, molecule_id, protonation_dic, prefix_folder, output_name):

		h = HeaderComposer(recipe_file, molecule_id, protonation_dic, prefix_folder, output_name)
		header_string = h.create_header()

		r = RecipeComposer(recipe_file)
		body_string = r.compose()

		self.recipe = header_string + '\n' + body_string

		return self.recipe

	# def generate_it0(self, recipe_file, molecule_id, protonation_dic, prefix_folder, out_suffix):
	#
	# 	h = HeaderComposer(recipe_file, molecule_id, protonation_dic, prefix_folder, out_suffix)
	#
	# 	header_string = h.create_header_it0()
	#
	# 	r = RecipeComposer(recipe_file)
	# 	body_string = r.compose()
	#
	# 	self.recipe = header_string + '\n' + body_string
	#
	# 	return self.recipe


class RecipeComposer:

	def __init__(self, recipe_file):
		self.protocol_path = get_full_path('haddock', 'protocols')
		# self.recipe = self.protocol_path + '/' + config.param_dic['input']['recipe']
		self.recipe = self.protocol_path + '/' + recipe_file

	def compose(self):
		""" Load the CNS recipe identify which modules need to be loaded and compose the recipe body """
		module_regex = r'.*(\@|/)(.*\.cns)'

		_ = self.list_dependencies(self.recipe)

		with open(self.recipe) as f:
			target = f.readlines()
		f.close()

		new = []
		check = True
		module_check = [(i, e) for i, e in enumerate(target) if '.cns' in e if '!' not in e]
		while check:

			for idx, m in module_check:

				# module found, get name and load
				try:
					m = re.findall(module_regex, target[idx])[0][-1]
				except IndexError:
					print(target[idx])
					exit()
				with open(f'{self.protocol_path}/{m}') as f:
					module = f.readlines()
				f.close()

				new = target[:idx] + module + target[idx+1:]

				break

			target = new
			module_check = [(i, e) for i, e in enumerate(target) if '.cns' in e if '!' not in e]
			if module_check:
				check = True
			else:
				check = False

		return ''.join(new)

	@staticmethod
	def identify_modules(cns_file):
		""" Find all modules in this CNS file """
		module_list = []

		if not os.path.isfile(cns_file):
			print(f'+ ERROR: Module {cns_file} not found')
			exit()

		with open(cns_file) as f:
			for pos, l in enumerate(f.readlines()):
				if '!' not in l:
					if '@' in l:
						if '.cns' in l:
							# module detected
							module = l.split('@')[-1].split('/')[-1].split('\n')[0]
							# module_name = f'{self.protocol_path}/{module}'
							module_list.append(module)
							# print(pos, module)
		return module_list

	def list_dependencies(self, target_f):
		""" List the CNS modular dependency of the target file """

		output_list = []

		target_f_name = target_f.split('/')[-1]

		# print(f'+ Creating dependency tree of {target_f_name}')
		# print(f'> {target_f_name}')

		# Lvl 1  Identify modules
		module_list = self.identify_modules(target_f)
		if module_list:
			tbw = ''
			for module in module_list:
				output_list.append(module)
				module_path = f'{self.protocol_path}/{module}'
				tbw += f'  |_ {module} \n'

				# Lvl 2 Dependencies
				# print(module)
				dependency_list = self.identify_modules(module_path)
				if dependency_list:
					for dependency in dependency_list:
						output_list.append(dependency)
						dependency_path = f'{self.protocol_path}/{dependency}'
						tbw += f'     |_ {dependency} \n'

						# Lvl 3 Co-dependency
						co_dependency_list = self.identify_modules(dependency_path)
						if co_dependency_list:
							for co_dependency in co_dependency_list:
								output_list.append(co_dependency)
								co_dependency_path = f'{self.protocol_path}/{co_dependency}'
								tbw += f'        |_ {co_dependency} \n'

								# Lvl 4 Co-co-dependency
								co_co_dependency_list = self.identify_modules(co_dependency_path)
								if co_co_dependency_list:
									for co_co_dependency in co_co_dependency_list:
										output_list.append(co_co_dependency)
										co_co_dependency_path = f'{self.protocol_path}/{co_co_dependency}'
										print(f'          |_ {co_co_dependency}')

										# Lvl 5 Co-co-co-dependency
										# You have gone too far...
										co_co_co_dependency_list = self.identify_modules(co_co_dependency_path)
										if co_co_co_dependency_list:
											for co_co_co_dependency in co_co_co_dependency_list:
												print(f'+ ERROR: Too many dependency levels {target_f_name} > '
												      f'{dependency} > {co_dependency} > {co_co_dependency} > {co_co_co_dependency}')
			# print(tbw)
		return output_list


class HeaderComposer:
	""" Each recipe has a Header with parameters and scoring definitions """

	def __init__(self, recipe_f, molecule_id, prot_dic, prefix_folder, output_name):
		self.protocol_path = get_full_path('haddock', 'protocols')

		self.input_header = ''
		self.ff_param_header = ''
		self.ff_top_header = ''
		self.scoring_header = ''
		self.link_header = ''
		self.protonation_header = ''

		self.header = ''

		self.param_file = config.ini.get('parameters', 'param_file')
		# self.protein_param = config.ini.get('parameters', 'protein_param')
		# self.solvent_param = config.ini.get('parameters', 'solvent_param')
		# self.ion_param = config.ini.get('parameters', 'ion_param')
		# self.ligand_param = config.ini.get('parameters', 'ligand_param')

		self.top_file = config.ini.get('topology', 'top_file')
		# self.protein_top = config.ini.get('topology', 'protein_top')
		# self.solvent_top = config.ini.get('topology', 'solvent_top')
		# self.break_top = config.ini.get('topology', 'break_top')
		# self.ion_top = config.ini.get('topology', 'ion_top')
		# self.ligand_top = config.ini.get('topology', 'ligand_top')

		self.link = config.ini.get('link', 'link')

		self.trans_vectors = {}
		for i in range(51):
			self.trans_vectors[i] = config.ini.get('translation_vectors', f'trans_vector_{i}')

		# self.scoring_params = config.param_dic['scoring-parameters']
		recipe_params = self.protocol_path + '/' + recipe_f.split('.')[0] + '.json'
		self.recipe_params = json.load(open(recipe_params))

		self.mol_id = molecule_id

		self.protonation_dic = prot_dic
		self.output_name = output_name
		self.folder = prefix_folder

	def create_header(self):
		param = self.load_ff_parameters()
		top = self.load_ff_topology()
		recipe_params = self.load_recipe_params()
		link = self.load_link()
		protonation = self.load_protonation_state()
		trans_vec = self.load_trans_vectors()
		output = self.prepare_output()

		self.header = output + param + top + recipe_params + link + protonation + trans_vec

		return self.header

	# def create_header_it0(self):
	# 	param = self.load_ff_parameters()
	# 	top = self.load_ff_topology()
	# 	recipe_params = self.load_recipe_params_it0()
	# 	link = self.load_link()
	# 	protonation = self.load_protonation_state()
	# 	output = self.prepare_output(output_f='dbg')
	#
	# 	self.header = output + param + top + recipe_params + link + protonation
	#
	# 	return self.header

	def prepare_output(self):
		""" Tell the recipe wich should be the output file """
		output_filename = ''
		folder = ''

		if 'output' in self.recipe_params:
			if 'folder' in self.recipe_params['output']:
				folder = self.recipe_params['output']['folder'] + '/'
				if not os.path.isdir(folder):
					os.mkdir(folder)

		output = '\n! Output structure\n'
		if self.output_name:
			output += f"eval ($output_pdb_filename= \"{folder}\" + \"{self.output_name}\" + \".pdb\")\n"
		else:
			if 'output' in self.recipe_params:
				if 'folder' in self.recipe_params['output']:
					folder = self.recipe_params['output']['folder'] + '/'

					if not os.path.isdir(folder):
						os.system(f'mkdir {folder}')

			# 	output_filename += f"{self.recipe_params['output']['folder']}/"
				# if not os.path.isdir(self.recipe_params['output']['folder']):
				# 	os.system(f"mkdir {self.recipe_params['output']['folder']}")

				output_filename += '$file'
				if 'psf' in self.recipe_params['output']:
					output += f"eval ($output_psf_filename= \"{folder}\" + $file_root - \".pdb\" + \".psf\")\n"
				if 'pdb' in self.recipe_params['output']:
					output += f"eval ($output_pdb_filename= \"{folder}\" + $file_root + \".pdb\")\n"
			else:
				print('+ ERROR: No output defined for this recipe')
				exit()

		return output

	def load_protonation_state(self):
		""" Add protonation states to the recipe """
		protonation_header = '\n! Protonation states\n'

		for i, chain in enumerate(self.protonation_dic):
			hise_l = [0] * 10
			hisd_l = [0] * 10
			hisd_counter = 0
			hise_counter = 0
			for res in self.protonation_dic[chain]:
				state = self.protonation_dic[chain][res].lower()
				if state == 'hise':
					hise_l[hise_counter] = res
					hise_counter += 1
				if state == 'hisd':
					hisd_l[hisd_counter] = res
					hisd_counter += 1

			hise_str = ''
			for e in [(i+1, c+1, r) for c, r in enumerate(hise_l)]:
				hise_str += f'eval ($toppar.hise_resid_{e[0]}_{e[1]} = {e[2]})\n'
			hisd_str = ''
			for e in [(i+1, c+1, r) for c, r in enumerate(hisd_l)]:
				hisd_str += f'eval ($toppar.hisd_resid_{e[0]}_{e[1]} = {e[2]})\n'

			protonation_header += hise_str
			protonation_header += hisd_str

		return protonation_header

	def load_ff_parameters(self):
		""" Add force-field specific parameters to its appropriate places in the scoring recipe """
		ff_param_header = '\n! FF parameters\n'
		ff_param_header += 'parameter\n'
		ff_param_header += f'  @@{self.param_file}\n'
		# ff_param_header += f'  @@{self.protein_param}\n'
		# ff_param_header += f'  @@{self.solvent_param}\n'
		# ff_param_header += f'  @@{self.ion_param}\n'
		# ff_param_header += f'  @@{self.ligand_param}\n'
		ff_param_header += 'end\n'

		return ff_param_header

	def load_ff_topology(self):
		""" Add force-field specific topology to its appropriate places in the scoring recipe """
		ff_top_header = '\n! Toplogy\n'
		ff_top_header += 'topology\n'
		ff_top_header += f'  @@{self.top_file}\n'
		# ff_top_header += f'  @@{self.protein_top}\n'
		# ff_top_header += f'  @@{self.solvent_top}\n'
		# ff_top_header += f'  @@{self.break_top}\n'
		# ff_top_header += f'  @@{self.ion_top}\n'
		# ff_top_header += f'  @@{self.ligand_top}\n'
		ff_top_header += 'end\n'

		return ff_top_header

	def load_recipe_params(self):

		recipe_param_header = '\n! Parameters\n'

		for param in self.recipe_params['params']:
			v = self.recipe_params['params'][param]
			if not v:
				# either 0 or empty string
				if isinstance(v, str):
					v = '\"\"'
				if isinstance(v, int):
					v = 0.0
			recipe_param_header += f'eval (${param}={v})\n'

		if 'chain' in self.recipe_params:
			# load molecule specific things
			for mol in self.recipe_params['chain']:
				for param in self.recipe_params['chain'][mol]:
					v = self.recipe_params['chain'][mol][param]
					# this are LOGICAL, which means no quotes
					recipe_param_header += f'eval (${param}_{mol}={v})\n'

			# evaluate($toppar.prot_segid_$nmol = $prot_segid_mol$nmol)
			# evaluate($toppar.fix_origin_$nmol =$fix_origin_mol$nmol)
			# evaluate($toppar.dna_$nmol =$dna_mol$nmol)
			# evaluate($toppar.cyclicpept_$nmol = $cyclicpept_mol$nmol)
			# evaluate($toppar.shape_$nmol = $shape_mol$nmol)
			# evaluate($toppar.cg_$nmol = $cg_mol$nmol)

			pass

		seed = random.randint(100, 999)
		recipe_param_header += f'set seed={seed} end\n'

		return recipe_param_header

	# def load_recipe_params(self):
	# 	# This function must account for ALL different syntaxes
	#
	# 	scoring_header = '\n! Parameters\n'
	#
	# 	if 'flag' in self.recipe_params:
	# 		for flag in self.recipe_params['flags']:
	# 			v = str(self.recipe_params['flags'][flag]).upper()
	# 			scoring_header += f'eval ($Data.flags.{flag} = {v})\n'
	#
	# 	if 'values' in self.recipe_params:
	# 		for value in self.recipe_params['values']:
	# 			v = self.recipe_params['values'][value]
	# 			scoring_header += f'eval (${value}={v})\n'
	#
	# 	if 'ncomp' in self.recipe_params:
	# 		ncomp = self.recipe_params['ncomp']
	# 		scoring_header += f'eval ($data.ncomponents={ncomp})\n'
	#
	# 	if 'segids' in self.recipe_params:
	# 		for i, segid in enumerate(self.recipe_params['segids']):
	# 			scoring_header += f'eval ($Toppar.prot_segid_{i+1}="{segid}")\n'
	#
	# 	if 'seed' in self.recipe_params:
	# 		seed = random.randint(100, 999)
	# 		scoring_header += f'set seed={seed} end\n'
	#
	# 	if 'chain' in self.recipe_params:
	# 		for chain_key in self.recipe_params["chain"][self.mol_id]:
	# 			chain_value = self.recipe_params["chain"][self.mol_id][chain_key]
	# 			scoring_header += f'eval (${chain_key}={chain_value})\n'
	#
	# 	return scoring_header

	def load_link(self):
		self.link_header = '\n! Link file\n'
		self.link_header += f'eval ($link_file = "{self.link}" )\n'

		return self.link_header

	def load_trans_vectors(self):
		trans_header = '\n! Translation vectors\n'
		for vector_id in self.trans_vectors:
			vector_file = self.trans_vectors[vector_id]
			trans_header += f'eval ($trans_vector_{vector_id} = "{vector_file}" )\n'

		return trans_header
