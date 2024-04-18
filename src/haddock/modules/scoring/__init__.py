"""HADDOCK3 modules to score models."""
import pandas as pd

from haddock.core.typing import FilePath
from haddock.modules.base_cns_module import BaseCNSModule
from haddock.modules import BaseHaddockModule


class ScoringModule(BaseHaddockModule):
    """Parent class for Scoring modules."""

    def output(self, output_fname: FilePath, sep: str = "\t") -> None:
        """Save the output in comprehensive tables."""
        # saves scoring data
        sc_data = []
        for pdb in self.output_models:
            sc_data.append([pdb.file_name, pdb.ori_name, pdb.md5, pdb.score])
        
        # converts to pandas dataframe and sorts by score
        df_columns = ["structure", "original_name", "md5", "score"]
        df_sc = pd.DataFrame(sc_data, columns=df_columns)
        df_sc_sorted = df_sc.sort_values(by="score", ascending=True)
        # writes to disk
        df_sc_sorted.to_csv(output_fname, sep=sep, index=False, na_rep="None")

        return

class CNSScoringModule(BaseCNSModule, ScoringModule):
    """Parent class for Scoring modules."""
