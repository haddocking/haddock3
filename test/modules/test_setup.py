import shutil
import unittest
import os
from haddock.modules.setup import Setup
from haddock.utils.files import get_full_path

data_path = get_full_path('test', 'data')


class TestSetup(unittest.TestCase):
	def setUp(self):

		setup_dic = {'molecules': {'mol1': f'{data_path}/mini.pdb'}, 'identifier': {'run': 1}, 'stage': {'topology': {'recipe': 'default'}}}
		self.Setup = Setup(setup_dic)

	def test_prepare_folders(self):

		check = self.Setup.prepare_folders()
		self.assertTrue(check)
		self.assertTrue(os.path.isdir('run1/'))
		self.assertTrue(os.path.isdir('run1/data'))
		self.assertTrue(os.path.isfile('run1/data/mol1_1.pdb'))
		self.assertTrue(os.path.isfile('run1/data/run.toml'))

	def test_configure_recipes(self):

		check = self.Setup.configure_recipes()

		self.assertTrue(check)
		self.assertTrue(os.path.isfile('run1/topology/template/generate-topology.cns'))

	def tearDown(self):
		shutil.rmtree('run1')


if __name__ == '__main__':
	unittest.main()
