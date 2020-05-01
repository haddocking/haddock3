""" Test the PDB utilities """
import filecmp
import unittest
import os
from shutil import copyfile
from haddock.modules.structure.utils import PDB
from haddock.utils.files import get_full_path

data_path = get_full_path('test', 'data')


class TestPDB(unittest.TestCase):

	def setUp(self):
		self.PDB = PDB()

	def test_treat_ensemble(self):
		copyfile(f'{data_path}/mini_ens.pdb', f'{data_path}/temp_ens.pdb')
		input_pdb_dic = {'mol1': f'{data_path}/temp_ens.pdb'}

		treated_dic = self.PDB.treat_ensemble(input_pdb_dic)
		expected_treated_dic = {'mol1': [f'{data_path}/temp_1.pdb', f'{data_path}/temp_2.pdb']}

		self.assertEqual(treated_dic, expected_treated_dic)
		self.assertTrue(filecmp.cmp(f'{data_path}/temp_1.pdb', f'{data_path}/mini_ens1.pdb'))
		self.assertTrue(filecmp.cmp(f'{data_path}/temp_2.pdb', f'{data_path}/mini_ens2.pdb'))

		os.remove(f'{data_path}/temp_1.pdb')
		os.remove(f'{data_path}/temp_2.pdb')
		os.remove(f'{data_path}/temp_ens.pdb')

	def test_load_structure(self):
		pdb_f = f'{data_path}/miniA.pdb'
		pdb_dic = self.PDB.load_structure(pdb_f)
		expected_pdb_dic = {'A': ['ATOM      2  CA  MET A   1      16.967  12.784   4.338  1.00 10.80      A    C  \n',
							  'ATOM      9  CA  ARG A   2      13.856  11.469   6.066  1.00  8.31      A    C  \n',
							  'ATOM     16  CA  CYS A   3      13.660  10.707   9.787  1.00  5.39      A    C  \n']}
		self.assertEqual(pdb_dic, expected_pdb_dic)

	def test_identify_chains(self):
		pdb_f = f'{data_path}/mini.pdb'
		chain_l = self.PDB.identify_chains(pdb_f)
		expected_chain_l = ['A', 'B','C']
		self.assertEqual(chain_l, expected_chain_l)

	def test_identify_segids(self):
		pdb_f = f'{data_path}/miniA.pdb'
		segid_l = self.PDB.identify_segids(pdb_f)
		expected_segid_l = ['A']
		self.assertEqual(segid_l, expected_segid_l)

	def test_split_models(self):
		ensamble_f = f'{data_path}/mini_ens.pdb'
		model_list = self.PDB.split_models(ensamble_f)
		expected_list = [f'{data_path}/mini_1.pdb', f'{data_path}/mini_2.pdb']

		self.assertEqual(model_list, expected_list, 'Name of list elements differ')
		self.assertTrue(filecmp.cmp(f'{data_path}/mini_1.pdb', f'{data_path}/mini_1.gold'))
		self.assertTrue(filecmp.cmp(f'{data_path}/mini_2.pdb', f'{data_path}/mini_2.gold'))

		for f in model_list:
			os.remove(f)

	def test_fix_id(self):
		nosegid_pdb_f = f'{data_path}/mini.pdb'
		nochain_pdb_f = f'{data_path}/mini_nochain.pdb'

		segid_pdb = self.PDB.fix_id(nosegid_pdb_f, priority='chain', overwrite=False)
		chain_pdb = self.PDB.fix_id(nochain_pdb_f, priority='seg', overwrite=False)

		self.assertTrue(filecmp.cmp(segid_pdb, f'{data_path}/mini_segid.pdb'))
		self.assertTrue(filecmp.cmp(chain_pdb, f'{data_path}/mini_segid.pdb'))

		os.remove(f'{data_path}/mini.pdb_')
		os.remove(f'{data_path}/mini_nochain.pdb_')

	def test_add_chainseg(self):

		copyfile(f'{data_path}/mini.pdb', f'{data_path}/temp.pdb')

		check = self.PDB.add_chainseg(f'{data_path}/temp.pdb', 'A')

		self.assertTrue(check)
		self.assertTrue(filecmp.cmp(f'{data_path}/temp.pdb', f'{data_path}/miniA.pdb'))

		os.remove(f'{data_path}/temp.pdb')

	def test_identify_chainseg(self):

		pdbf = f'{data_path}/miniA.pdb'

		chainseg = self.PDB.identify_chainseg(pdbf)

		self.assertEqual(chainseg, ['A'])

	def test_fix_chainseg(self):
		copyfile(f'{data_path}/mini_1.gold', f'{data_path}/mol1.pdb')
		copyfile(f'{data_path}/mini_2.gold', f'{data_path}/mol2.pdb')

		input_pdb_dic = {'mol1': f'{data_path}/mol1.pdb', 'segid1': 'X', 'mol2': f'{data_path}/mol2.pdb'}

		return_pdb_dic = self.PDB.fix_chainseg(input_pdb_dic)
		expected_return_dic = {'mol1': f'{data_path}/mol1.pdb', 'mol2': f'{data_path}/mol2.pdb'}

		self.assertEqual(return_pdb_dic, expected_return_dic)
		self.assertTrue(filecmp.cmp(f'{data_path}/mol1.pdb', f'{data_path}/miniX.pdb'))

		os.remove(f'{data_path}/mol1.pdb')
		os.remove(f'{data_path}/mol2.pdb')

	def test_sanitize(self):
		copyfile(f'{data_path}/mini.dirty.pdb', f'{data_path}/temp.pdb')

		input_pdb_dic = {'mol1': [f'{data_path}/temp.pdb']}

		model_list = self.PDB.sanitize(input_pdb_dic)

		expected_model_list = [f'{data_path}/temp.pdb']

		self.assertEqual(model_list, expected_model_list)
		self.assertTrue(filecmp.cmp(f'{data_path}/temp.pdb', f'{data_path}/mini.clean.pdb'))

		os.remove(f'{data_path}/temp.pdb')

	def test_count_atoms(self):
		pdb_f = f'{data_path}/mini.pdb'
		atom_count = self.PDB.count_atoms(pdb_f)
		self.assertEqual(atom_count, 3)

	def test_organize_chains(self):
		pass

	def test_replace_chain(self):
		pdb_f = f'{data_path}/mini.pdb'
		newchain_pdb = self.PDB.replace_chain(pdb_f, 'A', 'X', overwrite=False)
		self.assertTrue(filecmp.cmp(newchain_pdb, f'{data_path}/mini_A-X.pdb'))
		os.remove(f'{data_path}/mini.pdb_')

	def test_renumber(self):
		pass

	def test_load_seq(self):
		pass

	def tearDown(self):
		pass
		# os.remove(f'{data_path}/mini.pdb_')

if __name__ == '__main__':
	unittest.main()
