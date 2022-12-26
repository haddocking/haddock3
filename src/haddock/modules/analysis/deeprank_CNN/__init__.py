"""deeprank CNN scoring."""
import os
from pathlib import Path

from haddock.modules import BaseHaddockModule

from haddock.modules.analysis.deeprank_CNN.deeprank_CNN import CNN_score

RECIPE_PATH = Path(__file__).resolve().parent
DEFAULT_CONFIG = Path(RECIPE_PATH, "defaults.yaml")

class HaddockModule(BaseHaddockModule):
    """HADDOCK3 module for ranking docking poses with deeprank_CNN"""

    name = RECIPE_PATH.name

    def __init__(self, order, path, *ignore, init_params=DEFAULT_CONFIG,
                 **everything):
        super().__init__(order, path, init_params)

    @classmethod
    def confirm_installation(cls):
        """Confirm if contact executable is compiled."""
        return

    def _run(self):
        '''
        input parameters
        '''
        # Get the models generated in previous step
        models = self.previous_io.retrieve_models(
            individualize=True
            )
        output_dir = os.path.join(models[0].path.split('0_topoaa')[0], '1_deeprank_CNN')
        CNN = CNN_score(pdb_source = models[0].path, chain1 = self.params['chain1'], chain2 = self.params['chain2'], output_dir = os.path.join(output_dir, self.params['output_dir']))
        CNN.create_database()
        CNN.test_CNN()
        result = CNN.analysis_result()
        return result
  