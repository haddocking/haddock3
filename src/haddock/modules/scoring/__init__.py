"""HADDOCK3 modules to score models."""
import pandas as pd

from haddock.core.typing import FilePath, Path, Any
from haddock.modules.base_cns_module import BaseCNSModule
from haddock.modules import BaseHaddockModule, PDBFile


class ScoringModule(BaseHaddockModule):
    """Parent class for Scoring modules."""

    def output(
            self,
            output_fname: FilePath,
            sep: str = "\t",
            ascending_sort: bool = True,
            ) -> None:
        r"""Save the output in comprehensive tables.

        Parameters
        ----------
        output_fname : FilePath
            Path to the file where to write scoring data.
        sep : str, optional
            Character used as separator in file, by default "\t"
        ascending_sort : bool, optional
            Should the data be sorted in ascending order, by default True
        """
        # saves scoring data
        sc_data = []
        for pdb in self.output_models:
            sc_data.append([pdb.file_name, pdb.ori_name, pdb.md5, pdb.score])
        
        # converts to pandas dataframe and sorts by score
        df_columns = ["structure", "original_name", "md5", "score"]
        df_sc = pd.DataFrame(sc_data, columns=df_columns)
        df_sc_sorted = df_sc.sort_values(by="score", ascending=ascending_sort)
        # writes to disk
        df_sc_sorted.to_csv(output_fname,
                            sep=sep,
                            index=False,
                            na_rep="None",
                            float_format="%.3f")
        return


class CNSScoringModule(BaseCNSModule, ScoringModule):
    """Parent class for CNS Scoring modules."""

    def per_interface_output(
            self,
            output_fname: FilePath,
            sep: str = "\t",
            ascending_sort: bool = True,
            ) -> None:
        r"""Generate per interface scoring tsv output files.

        Parameters
        ----------
        output_fname : FilePath
            Path to the file where to write scoring data.
        sep : str, optional
            Character used as separator in file, by default "\t"
        ascending_sort : bool, optional
            Should the data be sorted in ascending order, by default True
        """
        # Retrieve all interfaces data for all pdb
        set_interfaces: list[str] = []
        pdb_interfaces_scores: dict[tuple[Any, Any, Any], dict[str, dict[str, float]]] = {}  # noqa : E501
        # Loop over models to recover interfaces
        for pdb in self.output_models:
            # if the pdb does not exist, skip
            if not Path(pdb.file_name).exists():
                continue
            interfaces_scores = self.read_per_interface_scores(pdb)
            reversed_interfaces_scores = {}

            # Hold list of interfaces
            for interface, scores in interfaces_scores.items():
                # Check if reverse chain order present
                split_inter = interface.split('_')
                reverse_interface = f"{split_inter[1]}_{split_inter[0]}"
                if reverse_interface in set_interfaces:
                    reversed_interfaces_scores[reverse_interface] = scores
                # Check if interface present
                if interface not in set_interfaces:
                    set_interfaces.append(interface)

            # Combine with reversed interface scores
            interfaces_scores.update(reversed_interfaces_scores)
            # Hold data
            pdbkey = (pdb.file_name, pdb.ori_name, pdb.md5)
            pdb_interfaces_scores[pdbkey] = interfaces_scores
        
        # Preset output file basename and extension
        output_file = Path(output_fname)
        output_bn = output_file.stem
        ouput_ext = ''.join(output_file.suffixes)
        # Write separated files for all interfaces
        for interface in set_interfaces:
            # Point data
            sc_data = []
            for pdbkey, interfaces_scores in pdb_interfaces_scores.items():
                if interface not in interfaces_scores.keys():
                    continue
                interface_scores = interfaces_scores[interface]
                score = interface_scores['HADDOCKscore']
                sc_data.append([pdbkey[0], pdbkey[1], pdbkey[2], score])
            # Check that the list is not empty
            if len(sc_data) == 0:
                continue
            # converts to pandas dataframe and sorts by score
            df_columns = ["structure", "original_name", "md5", "score"]
            df_sc = pd.DataFrame(sc_data, columns=df_columns)
            df_sc_sorted = df_sc.sort_values(
                by="score",
                ascending=ascending_sort,
                )
            # Generate output filename
            interface_output_fname = f"{output_bn}_{interface}{ouput_ext}"
            # writes to disk
            df_sc_sorted.to_csv(
                interface_output_fname,
                sep=sep,
                index=False,
                na_rep="None",
                float_format="%.3f",
                )
        return

    @staticmethod
    def read_per_interface_scores(pdb: PDBFile) -> dict[str, dict[str, float]]:
        """Read a pdb file and parse per interface scores.

        Parameters
        ----------
        pdb : PDBFile
            A PDBFile object.

        Returns
        -------
        interfaces_scores : dict[str, dict[str, float]]
            Dictionary holding per interfaces scores.
        """
        header = None
        interfaces_scores: dict[str, dict[str, float]] = {}
        with open(pdb.file_name, 'r') as filin:
            for _ in filin:
                if _.startswith('REMARK Interface'):
                    s_ = _.strip().split()[2:]
                    # Extract header
                    if not header:
                        header = s_
                    # Extract data
                    else:
                        chain1 = s_[header.index('Chain1')]
                        chain2 = s_[header.index('Chain2')]
                        haddockscore = float(s_[header.index('HADDOCKscore')])
                        evdw = float(s_[header.index('Evdw')])
                        eelec = float(s_[header.index('Eelec')])
                        edesol = float(s_[header.index('Edesol')])
                        bsa = float(s_[header.index('BSA')])
                        # Combine chains together
                        chains_key = f"{chain1}_{chain2}"
                        # Hold data
                        interfaces_scores[chains_key] = {
                            'HADDOCKscore': haddockscore,
                            'Evdw': evdw,
                            'Eelec': eelec,
                            'Edesol': edesol,
                            'BSA': bsa,
                            }
        return interfaces_scores
