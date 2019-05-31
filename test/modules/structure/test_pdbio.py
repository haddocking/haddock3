"""pdbio library tests"""

import unittest
import os
import shutil
import filecmp
from haddock.modules.structure import pdbio


class TestPDBIO(unittest.TestCase):
    
    def setUp(self):
        self.base_path = os.path.dirname(__file__)
        self.golden_data = os.path.join(self.base_path, 'golden_data', 'pdbio')
        self.scratch = os.path.join(self.base_path, 'scratch_pdbio')
        if not os.path.isdir(self.scratch):
            os.makedirs(self.scratch)

    def tearDown(self):
        try:
            shutil.rmtree(self.scratch)
        except:
            pass

    def test_sanitize(self):
        input_pdb = os.path.join(self.golden_data, 'Target100.pdb')
        output_pdb = os.path.join(self.scratch, 'sanitize_output.pdb')
        reference_pdb = os.path.join(self.golden_data, 'sanitize_Target100.pdb')
        
        pdbio.sanitize(input_pdb, output_pdb)

        self.assertTrue(filecmp.cmp(reference_pdb, output_pdb), 'Files are not identical')


if __name__ == '__main__':
    unittest.main()
