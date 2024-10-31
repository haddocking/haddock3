"""Accessibility scoring calculations."""
import copy
from pathlib import Path

from typing import Optional

from haddock.libs.libio import write_nested_dic_to_file
from haddock import log
from haddock.core.typing import ParamDict
from haddock.clis.restraints.calc_accessibility import (
    apply_cutoff,
    get_accessibility,
    )


def calc_acc_score(result_dict, buried_resdic, acc_resdic):
    """Calculate the accessibility score.
    
    The accessibility score is calculated as the sum of the number of
    residues that are not in the correct category (buried or accessible).
    
    Parameters
    ----------
    result_dict : dict
        A dictionary with the results from the accessibility calculation.
        
    buried_resdic : dict
        A dictionary with the buried residues.
    
    acc_resdic : dict
        A dictionary with the accessible residues.
    
    Returns
    -------
    acc_score : int
        The accessibility score.
    """
    # filter the residues
    acc_score: int = 0
    b_viols: dict[str, set] = {}
    a_viols: dict[str, set] = {}
    for ch in result_dict:
        if ch in buried_resdic:
            # for every supposedly buried residue that is not buried,
            # the score should increase by one.
            b_viols[ch] = set(buried_resdic[ch]).intersection(result_dict[ch])
            buried_viol_sc = len(b_viols[ch])
            acc_score += buried_viol_sc
        if ch in acc_resdic:
            # now the opposite logic with the accessible amino acids.
            a_viols[ch] = set(acc_resdic[ch]).difference(result_dict[ch])
            acc_viol_sc = len(a_viols[ch])
            acc_score += acc_viol_sc
    return acc_score, b_viols, a_viols


class AccScore:
    """AccScore class."""

    def __init__(
            self,
            model,
            path,
            buried_resdic,
            acc_resdic,
            cutoff,
            probe_radius,
            identificator,
            ):
        """Initialise AccScore class."""
        self.model = model
        self.path = path
        self.buried_resdic = buried_resdic
        self.acc_resdic = acc_resdic
        self.cutoff = cutoff
        self.data = []
        self.violations = []
        self.probe_radius = probe_radius
        self.violations_data = [self.model.file_name]
        self.identificator = identificator

    def run(self) -> None:
        """Run accessibility calculations."""
        mod_path = str(Path(self.model.path, self.model.file_name))
        try:
            access_data = get_accessibility(mod_path,
                                            probe_radius=self.probe_radius)
            result_dic = apply_cutoff(access_data, self.cutoff)
            acc_sc, b_viols, a_viols = calc_acc_score(result_dic,
                                                      self.buried_resdic,
                                                      self.acc_resdic)
        except AssertionError as e:
            log.warning(f"Error in get_accessibility for {self.model}: {e}.")
            acc_sc, b_viols, a_viols = None, None, None
        self.data = [self.model.file_name, self.model.ori_name, self.model.md5, acc_sc]
        for ch in self.buried_resdic:
            if ch in b_viols and b_viols[ch] != set():
                buried_str = ",".join([str(res) for res in sorted(b_viols[ch])])
                self.violations_data.append(buried_str)
            else:
                self.violations_data.append(None)
        for ch in self.acc_resdic:
            if ch in a_viols and a_viols[ch] != set():
                acc_str = ",".join([str(res) for res in sorted(a_viols[ch])])
                self.violations_data.append(acc_str)
            else:
                self.violations_data.append(None)
        return copy.deepcopy(self)


def extract_data_from_accscore_class(sasascore_objects: list[AccScore],
                                     violations_output_fname: Path
                                     ) -> Optional[dict[int, ParamDict]]:
    """
    Extracts data attributes from a list of AccScore objects into a structured dictionary,
    optionally sorts the data based on a specified key, and writes the sorted data to
    a file.

    Args:
        sasascore_objects (list[AccScore]): List of AccScore objects containing data attributes
                                     to be extracted.
        violations_output_fname (Path): Path to the output file with the sorted violations data

    Returns:
        Optional[dict[int, ParamDict]]: The structured violations_data dictionary
    Raises:
        (Include any specific exceptions the function may raise)
    """
    violations_data: dict[int, ParamDict] = {}
    for i, sasa_obj in enumerate(sasascore_objects):
        nbur = len(sasa_obj.buried_resdic)
        
        # fill in violations
        violations_data[i] = {
            "structure": sasa_obj.data[0],
        }
        for b, bur in enumerate(sasa_obj.buried_resdic):
            bur_key = f"bur_{bur}"
            violations_data[i].update(
                {bur_key: sasa_obj.violations_data[b+1]}
                )
        for a, acc in enumerate(sasa_obj.acc_resdic):
            acc_key = f"acc_{acc}"
            violations_data[i].update(
                {acc_key: sasa_obj.violations_data[nbur+a+1]}
                )

    # violations do not need to be sorted
    write_nested_dic_to_file(violations_data, violations_output_fname)

    return violations_data
