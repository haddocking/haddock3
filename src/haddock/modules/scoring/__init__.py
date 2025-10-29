"""HADDOCK3 modules to score models."""
import pandas as pd

from haddock.core.typing import FilePath, Optional, Path, Any
from haddock.gear.haddockmodel import HaddockModel
from haddock.modules import BaseHaddockModule
from haddock.modules.base_cns_module import BaseCNSModule


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
        sc_data = [
            [pdb.file_name, pdb.ori_name, pdb.md5, pdb.score]
            for pdb in self.output_models
            ]
        
        # converts to pandas dataframe and sorts by score
        df_columns = ["structure", "original_name", "md5", "score"]
        df_sc = pd.DataFrame(sc_data, columns=df_columns)
        df_sc_sorted = df_sc.sort_values(by="score", ascending=ascending_sort)
        # writes to disk
        df_sc_sorted.to_csv(
            output_fname,
            sep=sep,
            index=False,
            na_rep="None",
            float_format="%.3f",
            )
        return


class CNSScoringModule(BaseCNSModule, ScoringModule):
    """Parent class for CNS Scoring modules."""

    def per_interface_output(
            self,
            output_fname: FilePath,
            models: list[HaddockModel],
            sep: str = "\t",
            ascending_sort: bool = True,
            ) -> None:
        r"""Generate per interface scoring tsv output files.

        Parameters
        ----------
        output_fname : FilePath
            Path to the file where to write scoring data.
        models : list[HaddockModel]
            List of HaddockModel object obtained by loading the PDB files.
        sep : str, optional
            Character used as separator in file, by default "\t"
        ascending_sort : bool, optional
            Should the data be sorted in ascending order, by default True
        """
        # Skip the analysis if not desired by the user
        if not self.per_interface_scoring:
            return

        # Retrieve all interfaces data for all pdb
        set_interfaces: list[str] = []
        pdb_interfaces_scores: dict[tuple[Any, Any, Any], dict[str, dict[str, float]]] = {}  # noqa : E501

        # Loop over models to recover interfaces
        for pdb, haddock_model in zip(self.output_models, models):
            # if the pdb does not exist, skip
            if not Path(pdb.file_name).exists():
                continue
            # Make reverse interface checks
            # Score A->B == Score B->A
            reversed_interfaces_scores = {}
            # Hold list of interfaces
            for interface, scores in haddock_model.interface_energies.items():
                # Check if reverse chain order present
                split_inter = interface.split("_")
                reverse_interface = f"{split_inter[1]}_{split_inter[0]}"
                if reverse_interface in set_interfaces:
                    reversed_interfaces_scores[reverse_interface] = scores
                # Check if interface present
                if interface not in set_interfaces:
                    set_interfaces.append(interface)

            # Combine with reversed interface scores
            haddock_model.interface_energies.update(reversed_interfaces_scores)
            # Hold data
            pdbkey = (pdb.file_name, pdb.ori_name, pdb.md5)
            pdb_interfaces_scores[pdbkey] = haddock_model.interface_energies

        # Preset output file basename and extension
        output_file = Path(output_fname)
        output_bn = output_file.stem
        ouput_ext = "".join(output_file.suffixes)
        # Write separated files for all interfaces
        for interface in set_interfaces:
            # Point data
            sc_data = []
            for pdbkey, interfaces_scores in pdb_interfaces_scores.items():
                if interface not in interfaces_scores.keys():
                    continue
                interface_scores = interfaces_scores[interface]
                score = interface_scores["HADDOCKscore"]
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
    
    def extract_interface_combinations(self) -> list[str]:
        """Read interface specific parameters.

        Removes the `interface_combinations` from the parameters as not
        supported by CNS.
        Sets the `per_interface_scoring` to True if `interface_combinations`
        is not empty, as it is required for the interface-scores to be
        present in the PDB file.

        Returns
        -------
        interface_combinations : list[str]
            List of user-defined combinations.
        """
        # Here we pop the parameter as not supported by CNS and only used
        # at the python level for downstream analysis
        interface_combinations = self.params.pop("interface_combinations")
        # Set the per_interface_scoring parameter value as set by the user
        self.per_interface_scoring = self.params["per_interface_scoring"]
        # Check if the parameter is used
        if interface_combinations != []:
            # NOTE: per_interface_scoring must be set to true for the interface
            # scores to be present as REMARK in the header of the PDB file.
            self.params["per_interface_scoring"] = True
        return interface_combinations
    
    def update_pdb_scores(
            self,
            interface_combinations: list[str],
            ) -> tuple[list[HaddockModel], dict[str, list[str]]]:
        """Update the score attributes in the output pdb files.

        Parameters
        ----------
        interface_combinations : list[str]
            Input list of chains to be considered.
            Each list entry must be composed of two chains separated by coma.
            e.g.: 
            []              -> Consider all non-redundant chain pairs
            ["A,H", "A,L"]  -> Consider only the interface scores between A,H and A,L

        Returns
        -------
        output_haddock_models : list[HaddockModel]
            List of HaddockModel for each input pdb, that contain the actual
            scores loaded from the file.
        """
        # Obtain list of user defined interfaces
        desired_interfaces = self.build_interface_sets_combinations(
            interface_combinations
        )
        # Get the weights from the parameters
        _weight_keys = ("w_vdw", "w_elec", "w_desolv", "w_air", "w_bsa")
        weights = {e: self.params[e] for e in _weight_keys}

        interface_errors: dict[str, list[str]] = {}
        # Check for generated output, fail it not all expected files are found
        output_haddock_models: list[HaddockModel] = []
        for pdb in self.output_models:
            if pdb.is_present():
                # Convert pdb file into a HaddockModel to read the scores
                haddock_model = HaddockModel(pdb.file_name)
                # Set the unweighted energy terms
                pdb.unw_energies = haddock_model.energies
                # Compute set of interfaces if some are defined
                if len(desired_interfaces) >= 1:
                    try:
                        score = self.compute_interfaces_score(
                            haddock_model.interface_energies,
                            desired_interfaces,
                        )
                    # In case the output is None, fall back to standard
                    # haddock score that will always be ok to compute
                    except ValueError as interface_error:
                        # Hold that specific error
                        error_msg = str(interface_error)
                        if error_msg not in interface_errors.keys():
                            interface_errors[error_msg] = []
                        interface_errors[error_msg].append(pdb.file_name)
                        # Compute the haddock score instead
                        score = haddock_model.calc_haddock_score(**weights)
                    finally:
                        haddock_score = score
                # Otherwise simply compute the standard haddock score
                else:
                    # Compute the haddock score
                    haddock_score = haddock_model.calc_haddock_score(**weights)
                # Set the score attribute
                pdb.score = haddock_score
                output_haddock_models.append(haddock_model)
        # Log errors
        for error_msg, models in interface_errors.items():
            self.log(
                f"Interface error: '{error_msg}' occured {len(models)} times."
                " Falling back on classic HADDOCKscore."
                )
        return output_haddock_models

    @staticmethod
    def build_interface_sets_combinations(
            interface_combinations: list[str],
            ) -> list[str]:
        """Build desired combinatation of interfaces.

        Parameters
        ----------
        interface_combinations : list[str] | None
            Input list of chains to be considered.
            Each list entry must be composed of two chains separated by coma.
            e.g.: 
            []              -> Consider all non-redundant chain pairs
            ["A,H", "A,L"]  -> Consider only the interface scores between A,H and A,L

        Returns
        -------
        combinations : list[str]
            Unpacked list of interface combinations to consider.
            ["A,H", "A,L"]  -> ["A_H", "A_L"]
        """
        combinations: list[str] = [
            "_".join([c.strip() for c in chains.split(",")])
            for chains in interface_combinations
            if chains.count(",") == 1  # ensures input only contains pairs
        ]
        return combinations

    @staticmethod
    def compute_interfaces_score(
            interface_energies: dict[str, dict[str, float]],
            interface_sets_combinations: list[str],
            ) -> Optional[float]:
        """Compute the sum of selected interface haddock score.

        Parameters
        ----------
        interfaces_scores : dict[str, dict[str, float]]
            Scores of the various interfaces present in the pdb.
        interface_sets_combinations : list[str]
            List of interface combinations to consider
        """
        # Minimum set of interfaces must be >= 1
        if len(interface_sets_combinations) == 0:
            raise ValueError("No input interfaces")

        # Get all desired interfaces scores
        selected_interfaces_scores: list[float] = []
        for interface in interface_sets_combinations:
            if interface in interface_energies.keys():
                interface_score = interface_energies[interface]["HADDOCKscore"]
                selected_interfaces_scores.append(interface_score)
            else:
                # Check if reverse chain order present
                split_inter = interface.split("_")
                reverse_interface = f"{split_inter[1]}_{split_inter[0]}"
                if reverse_interface in interface_energies.keys():
                    interface_score = interface_energies[reverse_interface]["HADDOCKscore"]
                    selected_interfaces_scores.append(interface_score)
        # Sum all desired interfaces scores
        if len(selected_interfaces_scores) >= 1:
            new_score = sum(selected_interfaces_scores)
        else:
            raise ValueError("Selected interface not found")
        return new_score
