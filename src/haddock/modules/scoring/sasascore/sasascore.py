"""Accessibility scoring calculations."""
from haddock import log
from pathlib import Path
import pandas as pd
from haddock.core.typing import FilePath
from haddock.clis.restraints.calc_accessibility import (
    apply_cutoff,
    get_accessibility,
    )


def prettify_df(output_fname: str, columns: list[str], sortby: Optional[str] = None) -> None:
    """Prettify the output dataframe.

    Parameters
    ----------
    output_fname : str
        The name of the output file.
    
    columns : list
        The columns of the dataframe.
    
    sortby : str
        The column to sort by.
    """
    # dataframe conversion and sorting
    df = pd.read_csv(output_fname, sep="\t", header=None)
    df.columns = columns
    if sortby:
        df = df.sort_values(by=sortby, ascending=True)
    df.to_csv(output_fname, sep="\t", index=False, na_rep="None", float_format="%.0f")
    log.info(f"{output_fname} created.")


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
            # the score should increase by one
            b_viols[ch] = set(buried_resdic[ch]).difference(result_dict[ch])
            buried_viol_sc = len(b_viols[ch])
            acc_score += buried_viol_sc
        if ch in acc_resdic:
            # now the opposite logic with the accessible amino acids.
            a_viols[ch] = set(acc_resdic[ch]).difference(result_dict[ch])
            acc_viol_sc = len(a_viols[ch])
            acc_score += acc_viol_sc
    return acc_score, b_viols, a_viols


class AccScoreJob:
    """A Job dedicated to the parallel running of Accscore jobs."""

    def __init__(
            self,
            accscore_obj):
        """Initialise AccScoreJob."""
        self.accscore_obj = accscore_obj
        self.output = accscore_obj.output_name

    def run(self) -> None:
        """Run this AccScoreJob."""
        log.info(f"core {self.accscore_obj.core}, running AccScore...")
        self.accscore_obj.run()
        self.accscore_obj.output()
        return


class AccScore:
    """AccScore class."""

    def __init__(
            self,
            model_list,
            output_name,
            core,
            path,
            buried_resdic,
            acc_resdic,
            cutoff,
            viol_output_name,
            probe_radius
            ):
        """Initialise AccScore class."""
        self.model_list = model_list
        self.output_name = output_name
        self.core = core
        self.path = path
        self.buried_resdic = buried_resdic
        self.acc_resdic = acc_resdic
        self.cutoff = cutoff
        self.data = []
        self.violations = []
        self.viol_output_name = viol_output_name
        self.probe_radius = probe_radius

    def run(self) -> None:
        """Run accessibility calculations."""
        for mod in self.model_list:
            mod_path = str(Path(mod.path, mod.file_name))
            try:
                access_data = get_accessibility(mod_path,
                                                probe_radius=self.probe_radius)
                result_dic = apply_cutoff(access_data, self.cutoff)
                acc_sc, b_viols, a_viols = calc_acc_score(result_dic,
                                                          self.buried_resdic,
                                                          self.acc_resdic)
            except AssertionError as e:
                log.warning(f"Error in get_accessibility for {mod_path}: {e}.")
                acc_sc, b_viols, a_viols = None, None, None
                continue
            self.data.append([mod.file_name, mod.ori_name, mod.md5, acc_sc])
            # now violations data
            violations_data = [mod.file_name]
            for ch in self.buried_resdic:
                if ch in b_viols:
                    buried_str = ",".join([str(res) for res in sorted(b_viols[ch])])
                    violations_data.append(buried_str)
                else:
                    violations_data.append("None")
            
            for ch in self.acc_resdic:
                if ch in a_viols:
                    acc_str = ",".join([str(res) for res in sorted(a_viols[ch])])
                    violations_data.append(acc_str)
                else:
                    violations_data.append("None")
            self.violations.append(violations_data)
        return

    def output(self) -> None:
        """Write down accessibility scores to file."""
        output_fname = Path(self.path, self.output_name)
        df = pd.DataFrame(self.data)
        df.to_csv(output_fname, sep="\t", index=False, header=False)
        # now we save the violations
        violations_fname = Path(self.path, self.viol_output_name)
        viol_df = pd.DataFrame(self.violations)
        viol_df.to_csv(violations_fname, sep="\t", index=False, header=False)
