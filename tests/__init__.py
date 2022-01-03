"""Define common test variables."""
from pathlib import Path


test_folder = Path(__file__).resolve().parent
data_folder = Path(test_folder, 'data')

golden_data = Path(test_folder, 'golden_data')

broken_pdb = Path(data_folder, 'broken.pdb')
good_pdb = Path(data_folder, 'good.pdb')
