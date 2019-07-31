import filecmp
import math
import unittest
from unittest.mock import patch

from utils.files import get_full_path

data_path = get_full_path('test', 'test_data')
param_path = get_full_path('haddock', 'toppar')

with patch("sys.argv", ['', f'{data_path}/scoring.json']):
	from haddock.workflows.scoring.analysis.ana import Ana, Cluster


class TestAna(unittest.TestCase):

	def setUp(self):
		pdb_list = [f'{data_path}/000001_conv.pdb', f'{data_path}/000002_conv.pdb']
		self.ana = Ana(pdb_list)

	def test_calculate_contacts(self):
		self.ana.calculate_contacts()

		con_list = self.ana.con_list

		self.assertEqual(con_list[0][-19:], 'contacts/000001.con')
		self.assertEqual(con_list[1][-19:], 'contacts/000002.con')

		self.assertTrue(filecmp.cmp(con_list[0], f'{data_path}/000001.con'))
		self.assertTrue(filecmp.cmp(con_list[1], f'{data_path}/000002.con'))

	def test_extract_energies(self):
		self.ana.extract_energies()

		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['haddock-score'], None)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['total'], -809.935)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['bonds'], 123.841)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['angles'], 614.941)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['improper'], 165.975)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['dihe'], 0.0)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['vdw'], -96.4473)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['elec'], -713.488)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['air'], 0.0)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['cdih'], 0.0)
		self.assertTrue(math.isnan(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['sani']))
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['vean'], 0.0)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['dani'], 0.0)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['desolv'], 42.5652)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['bsa'], 4088.2)

	def test_calculate_haddock_score(self):
		self.ana.calculate_haddock_score()

		self.assertEqual(self.ana.structure_dic[f'{data_path}/000001_conv.pdb']['haddock-score'], -196.5797)
		self.assertEqual(self.ana.structure_dic[f'{data_path}/000002_conv.pdb']['haddock-score'], -190.8612)

	def test_run_fastcontact(self):
		pass

	def test_run_dfire(self):
		pass

	def test_run_dockq(self):
		pass

	def test_calc_fcc_matrix(self):
		pass

	def test_cluster(self):
		pass

	def test_output(self):
		pass


class TestCluster(unittest.TestCase):

	def setUp(self):
		pass

	def test_populate(self):
		pass

	def test_add_member(self):
		pass


class TestElement(unittest.TestCase):

	def setUp(self):
		pass

	def test_add_neighbor(self):
		pass

	def test_assign_cluster(self):
		pass


if __name__ == '__main__':
	unittest.main()
