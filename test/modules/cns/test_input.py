""" Test the Input generator """
import unittest
from unittest.mock import patch
from utils.files import get_full_path

data_path = get_full_path('test', 'test_data')
param_path = get_full_path('haddock', 'toppar')

with patch("sys.argv", ['', f'{data_path}/scoring.json']):
	from haddock.modules.cns.input import HeaderComposer, RecipeComposer, RecipeGenerator


class TestHeaderComposer(unittest.TestCase):

	def setUp(self):
		prot_dic = {'A': {42: 'HISE'}, 'B': {1: 'HISD', 2: 'HISE'}}
		self.hc = HeaderComposer(prot_dic, '')

	def test_load_ff_parameters(self):
		param = self.hc.load_ff_parameters()
		param_header_str = f'\n! FF parameters\nparameter\n  @@{param_path}/protein-allhdg5-4.param\n  @@{param_path}/water-allhdg5-4.param\n  @@{param_path}/ion.param\n  @@{param_path}/ligand.param\nend\n'

		self.assertEqual(param, param_header_str)

	def test_load_ff_topology(self):
		top = self.hc.load_ff_topology()

		top_header_str = f'\n! Toplogy\ntopology\n  @@{param_path}/protein-allhdg5-4.top\n  @@{param_path}/water-allhdg5-4.top\n  @@{param_path}/protein_break.top\n  @@{param_path}/ion.top\n  @@{param_path}/ligand.top\nend\n'

		self.assertEqual(top, top_header_str)

	def test_load_scoring_parameters(self):
		scoring = self.hc.load_scoring_parameters()

		self.assertEqual(scoring.split('\n')[2], 'eval ($Data.flags.dihed = FALSE)')
		self.assertEqual(scoring.split('\n')[15], 'eval ($disulphide_dist=3)')

	def test_load_link(self):
		link = self.hc.load_link()
		link_header_str = f'\n! Link file\neval ($link_file = "{param_path}/protein-allhdg5-4-noter.link" )\n'

		self.assertEqual(link, link_header_str)

	def test_protonation_dic(self):
		protonation_str = self.hc.load_protonation_state()

		self.assertEqual(protonation_str.split('\n')[2], 'eval ($toppar.hise_resid_1_1 = 42)')
		self.assertEqual(protonation_str.split('\n')[13], 'eval ($toppar.hisd_resid_1_2 = 0)')

	def test_create_header(self):
		header = self.hc.create_header()

		self.assertEqual(header.split('\n')[3], 'eval ($psfname= $file - ".pdb" + ".psf")')
		self.assertEqual(header.split('\n')[42], 'eval ($Toppar.prot_segid_2="B")')
		self.assertEqual(header.split('\n')[50], 'eval ($toppar.hise_resid_1_1 = 42)')


class TestRecipeComposer(unittest.TestCase):

	def setUp(self):
		self.rc = RecipeComposer()

	def test_compose(self):
		body_str = self.rc.compose()

		with open(f'{data_path}/composed-recipe.cns') as f:
			composed = f.readlines()
		f.close()

		composed = ''.join(composed)

		self.assertEqual(body_str, composed)

	def test_identify_modules(self):
		# is it possible to test static methods?
		module_list = self.rc.identify_modules(f'{data_path}/dependency.cns')

		self.assertEqual(module_list, ['def_solv_param.cns'])

	def test_list_dependencies(self):
		dependency_list = self.rc.list_dependencies(f'{data_path}/dependency.cns')

		self.assertEqual(dependency_list, ['def_solv_param.cns'])


class TestRecipeGenerator(unittest.TestCase):

	def setUp(self):
		self.rg = RecipeGenerator()

	def test_generate(self):

		dummy_dic = {}
		generated_recipe_str = self.rg.generate(dummy_dic, '')

		generated_recipe_list = generated_recipe_str.split('\n')

		self.assertEqual(generated_recipe_list[2], 'eval ($filename= $file - ".pdb" + ".pdb")')
		self.assertEqual(generated_recipe_list[6], 'parameter')
		self.assertEqual(generated_recipe_list[31], 'eval ($Data.flags.ncs = FALSE)')


if __name__ == '__main__':
	unittest.main()

