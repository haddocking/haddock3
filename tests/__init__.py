"""Define common test variables."""
from pathlib import Path


here = Path(__file__).resolve().parent
data_folder = Path(here, 'data')
broken_pdb = Path(data_folder, 'broken.pdb')
good_pdb = Path(data_folder, 'good.pdb')
