import unittest
from unittest.mock import patch
from utils.files import get_full_path

data_path = get_full_path('test', 'test_data')

with patch("sys.argv", ['', f'{data_path}/scoring.json']):
	from haddock.tools.reduce import run_reduce, analyze_protonation_state


class TestReduce(unittest.TestCase):

	def test_run_reduce(self):
		pdbf = f'{data_path}/test_pdbs/1f3g-sanitized.pdb'
		out, err = run_reduce(pdbf)

		out = out.decode('utf-8').split('\n')

		self.assertEqual(out[0], 'USER  MOD reduce.3.24.130724 H: found=0, std=0, add=1139, rem=0, adj=26')
		self.assertEqual(out[42], 'USER  MOD Single : A 167 LYS NZ  :NH3+    180:sc=       0   (180deg=0)')
		self.assertEqual(out[200], 'ATOM      0  H   GLU A  29      31.237  44.661  46.431  1.00 37.69      A    H   new')

	def test_analyze_protonation_state(self):
		pdbf = f'{data_path}/test_pdbs/1f3g-sanitized.pdb'
		protonation_dic = analyze_protonation_state(pdbf)

		expected_dic = {'A': {75: 'HISE', 90: 'HISD'}}

		self.assertEqual(protonation_dic, expected_dic)


if __name__ == '__main__':
	unittest.main()
