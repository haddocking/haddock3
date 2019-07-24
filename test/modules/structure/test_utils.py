""" Test the PDB utilities """
import filecmp
import unittest
from haddock.modules.structure.utils import PDB
from utils.files import get_full_path

data_path = get_full_path('test', 'test_data')


class TestPDB(unittest.TestCase):

	def setUp(self):
		self.PDB = PDB()

	def test_load_structure(self):
		pdb = f'{data_path}/000042.pdb'
		prot_dic = self.PDB.load_structure(pdb)

		line_a = 'ATOM      1  N   MET A   1      -5.983  45.471  12.068  1.00  0.00A\n'
		line_b = 'ATOM   5673  N   MET B   1      22.877  38.418   1.799  1.00  0.00B\n'
		line_c = 'ATOM   2837  N   MET C   1      -0.321 -40.801 -19.269  1.00  0.00C\n'
		line_d = 'ATOM   5741  N   MET D   1     -19.701 -37.908   5.445  1.00  0.00D\n'

		self.assertEqual(prot_dic['A'][0], line_a)
		self.assertEqual(prot_dic['B'][0], line_b)
		self.assertEqual(prot_dic['C'][0], line_c)
		self.assertEqual(prot_dic['D'][0], line_d)

	def test_identify_chains(self):
		pdb = f'{data_path}/000042.pdb'

		chain_list = ['A', 'B', 'C', 'D']

		self.assertEqual(self.PDB.identify_chains(pdb), chain_list)

	def test_split_models(self):

		ensamble_f = f'{data_path}/T146-ensamble-5.pdb'

		model_list = self.PDB.split_models(ensamble_f)

		self.assertTrue(filecmp.cmp(f'{data_path}/000000.pdb', model_list[0]))
		self.assertTrue(filecmp.cmp(f'{data_path}/000001.pdb', model_list[1]))
		self.assertTrue(filecmp.cmp(f'{data_path}/000002.pdb', model_list[2]))
		self.assertTrue(filecmp.cmp(f'{data_path}/000003.pdb', model_list[3]))
		self.assertTrue(filecmp.cmp(f'{data_path}/000004.pdb', model_list[4]))

	def test_chain2segid(self):

		segless_pdb = f'{data_path}/1crn-chain.pdb'

		fname = self.PDB.chain2segid(segless_pdb)

		self.assertTrue(filecmp.cmp(f'{data_path}/1crn-chain+segid.pdb', fname))

	def test_segid2chain(self):
		chainless_pdb = f'{data_path}/1crn-segid.pdb'

		fname = self.PDB.segid2chain(chainless_pdb)

		self.assertTrue(filecmp.cmp(f'{data_path}/1crn-segid+chain.pdb', fname))

	def test_sanitize(self):

		input_l = [f'{data_path}/1f3g-dirty.pdb']

		clean_list = self.PDB.sanitize(input_l)
		clean_pdb = clean_list[0]
		self.assertTrue(filecmp.cmp(f'{data_path}/1f3g-sanitized.pdb', clean_pdb))


if __name__ == '__main__':
	unittest.main()
