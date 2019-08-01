import os
import re
import haddock.workflows.scoring.config as config
from utils.files import get_full_path


class RecipeGenerator:

	def __init__(self):
		self.recipe = ''

	def generate(self, protonation_dic, out_suffix):
		h = HeaderComposer(protonation_dic, out_suffix)
		header_string = h.create_header()

		r = RecipeComposer()
		body_string = r.compose()

		self.recipe = header_string + '\n' + body_string

		return self.recipe


class RecipeComposer:

	def __init__(self):
		self.protocol_path = get_full_path('haddock', 'protocols')
		self.recipe = self.protocol_path + '/' + config.param_dic['input']['recipe']

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
				m = re.findall(module_regex, target[idx])[0][-1]
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
				if '@' in l:
					if '.cns' in l:
						# module detected
						module = l.split('@')[-1].split('/')[-1].split('\n')[0]
						# module_name = f'{self.protocol_path}/{module}'
						module_list.append(module)
		return module_list

	def list_dependencies(self, target_f):
		""" List the CNS modular dependency of the target file """

		output_list = []

		target_f_name = target_f.split('/')[-1]

		print(f'+ Creating dependency tree of {target_f_name}')
		print(f'> {target_f_name}')

		# Lvl 1  Identify modules
		module_list = self.identify_modules(target_f)
		if module_list:
			tbw = ''
			for module in module_list:
				output_list.append(module)
				module_path = f'{self.protocol_path}/{module}'
				tbw += f'  |_ {module} \n'

				# Lvl 2 Dependencies
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
			print(tbw)
		return output_list


class HeaderComposer:
	""" Each recipe has a Header with parameters and scoring definitions """

	def __init__(self, prot_dic, out_suffix):
		self.input_header = ''
		self.ff_param_header = ''
		self.ff_top_header = ''
		self.scoring_header = ''
		self.link_header = ''
		self.protonation_header = ''

		self.header = ''

		self.protein_param = config.ini.get('parameters', 'protein_param')
		self.solvent_param = config.ini.get('parameters', 'solvent_param')
		self.ion_param = config.ini.get('parameters', 'ion_param')
		self.ligand_param = config.ini.get('parameters', 'ligand_param')

		self.protein_top = config.ini.get('topology', 'protein_top')
		self.solvent_top = config.ini.get('topology', 'solvent_top')
		self.break_top = config.ini.get('topology', 'break_top')
		self.ion_top = config.ini.get('topology', 'ion_top')
		self.ligand_top = config.ini.get('topology', 'ligand_top')

		self.link = config.ini.get('link', 'link')

		self.scoring_params = config.param_dic['scoring-parameters']

		self.protonation_dic = prot_dic
		self.out_suffix = out_suffix

	def create_header(self):
		param = self.load_ff_parameters()
		top = self.load_ff_topology()
		scoring_param = self.load_scoring_parameters()
		link = self.load_link()
		protonation = self.load_protonation_state()
		output = self.prepare_output()

		self.header = output + param + top + scoring_param + link + protonation

		return self.header

	def prepare_output(self):
		""" Tell the recipe wich should be the output file """
		output = '\n! Output structure\n'
		output += f'eval ($filename= $file - ".pdb" + "{self.out_suffix}.pdb")\n'
		output += f'eval ($psfname= $file - ".pdb" + "{self.out_suffix}.psf")\n'

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
		ff_param_header += f'  @@{self.protein_param}\n'
		ff_param_header += f'  @@{self.solvent_param}\n'
		ff_param_header += f'  @@{self.ion_param}\n'
		ff_param_header += f'  @@{self.ligand_param}\n'
		ff_param_header += 'end\n'

		return ff_param_header

	def load_ff_topology(self):
		""" Add force-field specific topology to its appropriate places in the scoring recipe """
		ff_top_header = '\n! Toplogy\n'
		ff_top_header += 'topology\n'
		ff_top_header += f'  @@{self.protein_top}\n'
		ff_top_header += f'  @@{self.solvent_top}\n'
		ff_top_header += f'  @@{self.break_top}\n'
		ff_top_header += f'  @@{self.ion_top}\n'
		ff_top_header += f'  @@{self.ligand_top}\n'
		ff_top_header += 'end\n'

		return ff_top_header

	def load_scoring_parameters(self):
		scoring_header = '\n! Scoring parameters\n'
		for flag in self.scoring_params['flags']:
			v = str(self.scoring_params['flags'][flag]).upper()
			scoring_header += f'eval ($Data.flags.{flag} = {v})\n'

		for value in self.scoring_params['values']:
			v = self.scoring_params['values'][value]
			scoring_header += f'eval (${value}={v})\n'

		ncomp = self.scoring_params['ncomp']
		scoring_header += f'eval ($data.ncomponents={ncomp})\n'

		for i, segid in enumerate(self.scoring_params['segids']):
			scoring_header += f'eval ($Toppar.prot_segid_{i+1}="{segid}")\n'

		return scoring_header

	def load_link(self):
		self.link_header = '\n! Link file\n'
		self.link_header += f'eval ($link_file = "{self.link}" )\n'

		return self.link_header
