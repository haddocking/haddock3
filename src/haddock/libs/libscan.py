"""Shared helpers for the mutagenesis scan modules.

The ``alascan`` (alanine), ``rnascan`` (RNA base) and ``dnascan`` (DNA base
pair) analysis modules share the same scoring, z-score, cluster-output,
b-factor and TSV-writing machinery. The common, module-agnostic pieces live
here; each module keeps only its own mutation logic (which residues/atoms are
mutated, how the interface is scanned, and how mutation jobs are executed).
"""

import io
import os
from contextlib import redirect_stdout
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from haddock import log
from haddock.clis import cli_score
from haddock.libs.libontology import PDBFile


@dataclass
class MutationResult:
    """Result from a single mutation.

    The ``partner_*`` fields are only used by modules that perform paired
    (double) mutations, such as ``dnascan`` (base-pair mutations); they default
    to ``None`` for single-residue mutations (``alascan``, ``rnascan``).
    """

    model_id: str
    chain: str
    resid: int
    ori_resname: str
    target_resname: str
    # components of "mutant_scores": score, vdw, elec, desolv, bsa
    mutant_scores: Tuple[float, float, float, float, float]
    delta_scores: Tuple[float, float, float, float, float]
    success: bool
    error_msg: Optional[str] = None
    partner_chain: Optional[str] = None
    partner_resid: Optional[int] = None
    partner_ori_resname: Optional[str] = None
    partner_target_resname: Optional[str] = None


def get_score_string(
    pdb_f: str,
    run_dir: str,
    outputpdb: bool = False,
    ligand_param_fname: Union[Path, str] = "",
    ligand_top_fname: Union[Path, str] = "",
) -> List[str]:
    """Get score output from cli_score.main.

    Parameters
    ----------
    pdb_f : str
        Path to the pdb file.
    run_dir : str
        Path to the run directory.
    outputpdb : bool, optional
        If True, the output, energy-minimized pdb file will be written.
        Default is False.
    ligand_param_fname : Union[Path, str]
        Path to additional parameter file used by CNS.
    ligand_top_fname : Union[Path, str]
        Path to additional topology file used by CNS.

    Returns
    -------
    out : list[str]
        List of strings with the score output.
    """
    stdout = io.StringIO()
    with redirect_stdout(stdout):
        cli_score.main(
            pdb_f,
            run_dir,
            full=True,
            outputpdb=outputpdb,
            ligand_param_fname=ligand_param_fname,
            ligand_top_fname=ligand_top_fname,
        )
    out = stdout.getvalue().split(os.linesep)
    return out


def calc_score(
    pdb_f: str,
    run_dir: str,
    outputpdb: bool = False,
    ligand_param_fname: Union[Path, str] = "",
    ligand_top_fname: Union[Path, str] = "",
) -> Tuple[float, float, float, float, float]:
    """Calculate the score of a model.

    Parameters
    ----------
    pdb_f : str
        Path to the pdb file.
    run_dir : str
        Path to the run directory.
    outputpdb : bool, optional
        If True, the output, energy-minimized pdb file will be written.
        Default is False.
    ligand_param_fname : Union[Path, str]
        Path to additional parameter file used by CNS
    ligand_top_fname : Union[Path, str]
        Path to additional topology file used by CNS

    Returns
    -------
    Tuple[float, float, float, float, float]
        All the scores
    score : float
        Haddock score.
    vdw : float
        Van der Waals energy.
    elec : float
        Electrostatic energy.
    desolv : float
        Desolvation energy.
    bsa : float
        Buried surface area.

    Raises
    ------
    FileNotFoundError
        Error when the file could not be located.
    """
    scores_strings = get_score_string(
        pdb_f,
        run_dir,
        outputpdb=outputpdb,
        ligand_param_fname=ligand_param_fname,
        ligand_top_fname=ligand_top_fname,
    )
    # Loop over lines to search for the score and components
    for ln in scores_strings:
        if ln.startswith("> HADDOCK-score (emscoring)"):
            score = float(ln.split()[-1])
        if ln.startswith("> vdw"):
            vdw = float(ln.split("vdw=")[1].split(",")[0])
            elec = float(ln.split("elec=")[1].split(",")[0])
            desolv = float(ln.split("desolv=")[1].split(",")[0])
            bsa = float(ln.split("bsa=")[1].split(",")[0])
        if ln.startswith("> writing") and outputpdb:
            outpdb = ln.split()[-1]
            # check if the output pdb file exists
            if not os.path.exists(outpdb):
                raise FileNotFoundError(f"Could not find output pdb file {outpdb}")
    return score, vdw, elec, desolv, bsa


def add_zscores(df_scan_clt, column="delta_score"):
    """Add z-scores to the dataframe.

    Parameters
    ----------
    df_scan : pandas.DataFrame
        Dataframe with the scan results for the model.

    colunm : str
        Column to calculate the z-score.

    Returns
    -------
    df_scan : pandas.DataFrame
        Dataframe with the z-scores added.
    """
    mean_delta = np.mean(df_scan_clt[column])
    std_delta = np.std(df_scan_clt[column])
    if std_delta > 0.0:
        df_scan_clt["z_score"] = (df_scan_clt[column] - mean_delta) / std_delta
    else:
        df_scan_clt["z_score"] = 0.0
    return df_scan_clt


class AddDeltaBFactor:
    """Add the (normalised) delta score into the b-factor column of a PDB.

    Subclasses set :attr:`module_name` (used for the output filename suffix,
    ``<model>_<module_name>.pdb``) and may override :meth:`_residue_keys` to
    control which residues each mutation result is attributed to (e.g. both
    nucleotides of a base pair for ``dnascan``).
    """

    #: filename suffix / module identifier; overridden by subclasses
    module_name: str = "scan"

    def __init__(
        self,
        model: PDBFile,
        path: Path,
        model_results: List[MutationResult],
    ):
        """Initialisation function

        Parameters
        ----------
        model : PDBFile
            PDBfile model to be modified
        path : Path
            Where to write the new file
        """
        self.model = model
        self.path = path
        self.input_results = model_results

    def run(self) -> PDBFile:
        """Perform the addition of delta scores as bfactor in pdb file."""
        self.reorder_results()
        # Define new pdb filename
        model_fname = self.model.file_name.removesuffix(".pdb")
        output_fname = f"{model_fname}_{self.module_name}.pdb"
        # Add delta_score as a bfactor to the model
        self.write_delta_score_to_pdb(output_fname)
        # Update attributes of the model
        self.model.ori_name = self.model.file_name
        self.model.file_name = output_fname
        self.model.full_name = output_fname
        self.model.rel_path = Path("..", Path(self.path).name, output_fname)
        self.model.path = str(Path(".").resolve())
        return self.model

    def _residue_keys(self, mut_result: MutationResult) -> List[Tuple[str, int]]:
        """Return the ``(chain, resid)`` keys a mutation result applies to.

        By default a result applies only to its own residue. Modules that
        perform paired mutations override this to also colour the partner.
        """
        return [(mut_result.chain, mut_result.resid)]

    def reorder_results(self) -> None:
        """Aggregate mutation results into a per-residue delta score map.

        Organise mutation results into a dictionary with ``chain-resid`` keys
        and, as values, the mean delta score over all mutations attributed to
        that residue (e.g. ``{"A-115": 5.2}``), and determine the min and max
        score within the model to normalise the data.
        """
        # First gather every delta score observed for each residue
        per_residue_scores: Dict[str, List[float]] = {}
        for mut_result in self.input_results:
            for chain, resid in self._residue_keys(mut_result):
                chain_res_key = f"{chain}-{resid}"
                per_residue_scores.setdefault(chain_res_key, []).append(
                    mut_result.delta_scores[0]
                )
        # Average the delta scores attributed to each residue
        self.model_results: Dict[str, float] = {}
        all_delta_scores: List[float] = []
        for chain_res_key, delta_scores in per_residue_scores.items():
            mean_delta = float(np.mean(delta_scores))
            self.model_results[chain_res_key] = mean_delta
            all_delta_scores.append(mean_delta)
        # Obtain min and max delta score to be able to normalize later
        if len(all_delta_scores) >= 1:
            self.min_score = min(all_delta_scores)
            self.max_score = max(all_delta_scores)
        else:
            self.min_score = -1
            self.max_score = 1

    def write_delta_score_to_pdb(self, output_path: Path) -> Path:
        """Add delta scores as b-factors in PDB file.

        Parameters
        ----------
        output_path : Path
            Path to the pdb file.

        Returns
        -------
        output_path : Path
            Path to the pdb file with the b-factors added.
        """
        # Input pdb file path
        input_pdbfile = Path(self.model.path, self.model.file_name)
        # Start writting file containing bfactor
        with open(input_pdbfile, "r") as fin, open(output_path, "w") as fout:
            for line in fin:
                if line.startswith("ATOM"):
                    chain = line[21]
                    resid = int(line[22:26])
                    chain_res_key = f"{chain}-{resid}"
                    try:
                        delta_score = self.model_results[chain_res_key]
                        norm_delta = self.normalize_score(delta_score)
                    except KeyError:
                        norm_delta = 0.0
                    # Generate the bfactor string
                    delta_str = f"{norm_delta:.2f}".rjust(6, " ")
                    # Modify the PDB line
                    line = line[:60] + delta_str + line[66:]
                fout.write(line)
        return output_path

    def normalize_score(self, score: float) -> float:
        """Normalise the input score based on observed scores for this model

        In case normalisation cannot be performed, returns 50.0

        Parameters
        ----------
        score : float
            Input score to be normalized

        Returns
        -------
        norm100 : float
            Normalized score between 0 and 100
        """
        try:
            norm = (score - self.min_score) / (self.max_score - self.min_score)
        except ZeroDivisionError:
            norm = 0.5
        norm100 = 100 * norm
        return norm100


#: statistics columns appended after the per-module identity columns in the
#: cluster TSV output (mean/std of each energy term plus the population fraction)
CLUSTER_STAT_COLUMNS = [
    "delta_score",
    "delta_score_std",
    "delta_vdw",
    "delta_vdw_std",
    "delta_elec",
    "delta_elec_std",
    "delta_desolv",
    "delta_desolv_std",
    "delta_bsa",
    "delta_bsa_std",
    "frac_pres",
]


class ClusterOutputer:
    """Base class writing the per-cluster scan TSV (and optional plot).

    Subclasses provide the module-specific identity columns/rows (which depend
    on whether a mutation is a single residue or a base pair) via
    :meth:`_identity_columns` and :meth:`_identity_row`, and set the class
    attributes :attr:`module_name`, :attr:`sort_columns`,
    :attr:`zscore_reference`, :attr:`default_scan_residue` and
    :attr:`plot_func`. The averaging, z-scoring, sorting, file writing and plot
    handling are shared.
    """

    module_name: str = "scan"
    default_scan_residue: str = "residue"
    sort_columns = ["chain", "resid"]
    zscore_reference: str = "residues"
    #: plotting function ``f(df, clt_id, scan_res, offline=...)``; set by subclass
    plot_func = None

    def __init__(
        self,
        cluster_scan_data: Dict[str, Dict[str, Union[float, int]]],
        clt_id: str,
        clt_population: int,
        scan_residue: Optional[str] = None,
        generate_plot: bool = False,
        offline: bool = False,
    ):
        """Initialization function

        Parameters
        ----------
        cluster_scan_data : Dict[str, Dict[str, Union[float, int]]]
            Dictionary containing scan data per mutation identifier.
        clt_id : str
            Cluster identifier
        clt_population : int
            Number of entries in this cluster.
        scan_residue : str, optional
            Label of the scan; defaults to the class' ``default_scan_residue``.
        generate_plot : bool, optional
            Defines if a plot must be generated, by default False
        offline : bool, optional
            Defines if the plot should be functional offline, by default False
        """
        self.cluster_scan_data = cluster_scan_data
        self.clt_id = clt_id
        self.clt_population = clt_population
        self.scanned_residue = (
            scan_residue if scan_residue is not None else self.default_scan_residue
        )
        self.generate_plot = generate_plot
        self.offline = offline

    def _identity_columns(self) -> List[str]:
        """Return the identity column names preceding the statistics columns."""
        raise NotImplementedError

    def _identity_row(self, ident, clt_res_dt) -> list:
        """Return the identity values for one mutation identifier."""
        raise NotImplementedError

    def _extra_header_lines(self) -> List[str]:
        """Return extra comment lines to add to the TSV header (no leading #)."""
        return []

    def _stat_values(self, clt_res_dt) -> list:
        """Return the averaged statistics values for one mutation identifier."""
        return [
            np.mean(clt_res_dt["delta_score"]),
            np.std(clt_res_dt["delta_score"]),
            np.mean(clt_res_dt["delta_vdw"]),
            np.std(clt_res_dt["delta_vdw"]),
            np.mean(clt_res_dt["delta_elec"]),
            np.std(clt_res_dt["delta_elec"]),
            np.mean(clt_res_dt["delta_desolv"]),
            np.std(clt_res_dt["delta_desolv"]),
            np.mean(clt_res_dt["delta_bsa"]),
            np.std(clt_res_dt["delta_bsa"]),
            clt_res_dt["frac_pr"] / self.clt_population,
        ]

    def run(self) -> str:
        """Write cluster scan output to a scan_clt_X.tsv file (and plot).

        Return
        ------
        scan_clt_filename : str
            Name of the tsv file written
        """
        clt_data = []
        for ident, clt_res_dt in self.cluster_scan_data.items():
            clt_data.append(
                self._identity_row(ident, clt_res_dt) + self._stat_values(clt_res_dt)
            )
        df_cols = self._identity_columns() + CLUSTER_STAT_COLUMNS
        df_scan_clt = pd.DataFrame(clt_data, columns=df_cols)
        # adding clt-based Z score
        df_scan_clt = add_zscores(df_scan_clt, "delta_score")
        # Sort rows
        df_scan_clt.sort_values(by=self.sort_columns, inplace=True)
        # Generate output CSV data
        csv_data = io.StringIO()
        df_scan_clt.to_csv(csv_data, index=False, float_format="%.2f", sep="\t")
        csv_data.seek(0)
        # Define output csv filepath
        scan_clt_filename = f"scan_clt_{self.clt_id}.tsv"
        # Write the file
        with open(scan_clt_filename, "w") as fout:
            fout.write(f"{'#' * 80}{os.linesep}")
            fout.write(
                f"# `{self.module_name}` cluster results for cluster "
                f"{self.clt_id}{os.linesep}"
            )
            fout.write(f"# reported values are the average for the cluster{os.linesep}")
            fout.write(f"#{os.linesep}")
            for extra in self._extra_header_lines():
                fout.write(f"# {extra}{os.linesep}")
            fout.write(
                f"# z_score is calculated with respect to the mean values of all "
                f"{self.zscore_reference}{os.linesep}"
            )
            fout.write(f"{'#' * 80}{os.linesep}")
            fout.write(csv_data.read())

        # Generate plot
        self.gen_plot(df_scan_clt)

        return scan_clt_filename

    def gen_plot(self, df_scan_clt: pd.DataFrame) -> None:
        """Generate the scan plot based on provided data.

        Parameters
        ----------
        df_scan_clt : pd.DataFrame
            The data frame containing the data to be plotted.
        """
        if not self.generate_plot or self.plot_func is None:
            return
        try:
            self.plot_func(
                df_scan_clt,
                self.clt_id,
                self.scanned_residue,
                offline=self.offline,
            )
        except Exception as e:
            log.warning(
                f"Could not create interactive plot. The following error occurred {e}"
            )


#: columns of the per-model scan TSV shared by single-residue scans
STANDARD_SCAN_COLUMNS = [
    "chain",
    "res",
    "ori_resname",
    "end_resname",
    "score",
    "vdw",
    "elec",
    "desolv",
    "bsa",
    "delta_score",
    "delta_vdw",
    "delta_elec",
    "delta_desolv",
    "delta_bsa",
]


def standard_scan_row(result: MutationResult) -> list:
    """Build a per-model TSV row for a single-residue mutation result."""
    m_score, m_vdw, m_elec, m_des, m_bsa = result.mutant_scores
    d_score, d_vdw, d_elec, d_des, d_bsa = result.delta_scores
    return [
        result.chain,
        result.resid,
        result.ori_resname,
        result.target_resname,
        m_score,
        m_vdw,
        m_elec,
        m_des,
        m_bsa,
        d_score,
        d_vdw,
        d_elec,
        d_des,
        d_bsa,
    ]


def write_scan_out(
    results: List[MutationResult],
    model_id: str,
    module_name: str,
    sort_columns: List[str],
    row_builder=standard_scan_row,
    columns: Optional[List[str]] = None,
    zscore_reference: str = "residues",
    extra_header_lines=(),
) -> None:
    """Save mutation results for one model to ``scan_{model_id}.tsv``.

    Parameters
    ----------
    results : List[MutationResult]
        List of mutation results from scanning.
    model_id : str
        Identifier for the model used in the filename.
    module_name : str
        Module label used in the file header (``alascan``/``rnascan``/...).
    sort_columns : List[str]
        Column names used to sort the rows.
    row_builder : callable
        Function mapping a successful ``MutationResult`` to a TSV row.
    columns : List[str], optional
        Column names of the rows; defaults to ``STANDARD_SCAN_COLUMNS``.
    zscore_reference : str
        Wording used in the z-score header comment (``residues``/``mutations``).
    extra_header_lines : iterable of str
        Extra comment lines added to the header (without the leading ``#``).

    Notes
    -----
    If ``results`` is empty, no file is created.
    """
    if not results:
        print(f"No scan results for model {model_id}")
        return

    if columns is None:
        columns = STANDARD_SCAN_COLUMNS

    scan_data = []
    native_score = None
    for result in results:
        if result.success:
            scan_data.append(row_builder(result))
            # Get native score from first result (mutant + delta = native)
            if native_score is None:
                native_score = result.mutant_scores[0] + result.delta_scores[0]

    if scan_data:
        df_scan = pd.DataFrame(scan_data, columns=columns)
        # Compute z-scores
        df_scan = add_zscores(df_scan, "delta_score")
        # Sort the rows
        df_scan.sort_values(by=sort_columns, inplace=True)
        # Get csv data as string
        csv_io = io.StringIO()
        df_scan.to_csv(csv_io, index=False, float_format="%.2f", sep="\t")
        csv_io.seek(0)
        # Save to tsv
        output_file = f"scan_{model_id}.tsv"
        with open(output_file, "w") as fout:
            fout.write(f"{'#' * 80}{os.linesep}")
            fout.write(f"# `{module_name}` results for {model_id}{os.linesep}")
            fout.write(f"#{os.linesep}")
            for extra in extra_header_lines:
                fout.write(f"# {extra}{os.linesep}")
            fout.write(f"# native score = {native_score}{os.linesep}")
            fout.write(f"#{os.linesep}")
            fout.write(
                f"# z_score is calculated with respect to the other "
                f"{zscore_reference}{os.linesep}"
            )
            fout.write(f"{'#' * 80}{os.linesep}")
            fout.write(csv_io.read())


def group_scan_by_cluster(
    models: List[PDBFile],
    results_by_model: Dict[str, List[MutationResult]],
    ident_builder,
    metadata_builder=None,
) -> Tuple[Dict[str, Dict[str, Dict[str, Union[float, int]]]], Dict[str, int]]:
    """Group per-model scan data per cluster.

    Parameters
    ----------
    models : List[PDBFile]
        List of input models.
    results_by_model : dict
        Mutation results keyed by model id.
    ident_builder : callable
        Function mapping a ``MutationResult`` to its unique string identifier.
    metadata_builder : callable, optional
        Function mapping a ``MutationResult`` to a dict of static metadata stored
        alongside the accumulated delta lists (used by paired-mutation modules).

    Returns
    -------
    clt_scan : dict
        Scan data for each cluster, grouped by mutation identifier.
    clt_pops : dict
        Number of entries for each cluster.
    """
    clt_scan: Dict[str, Dict[str, Dict[str, Union[float, int]]]] = {}
    clt_pops: Dict[str, int] = {}
    for model in models:
        cl_id = model.clt_id
        # unclustered models have cl_id = None
        if cl_id is None:
            cl_id = "unclustered"
        if cl_id not in clt_scan:
            clt_scan[cl_id] = {}
            clt_pops[cl_id] = 0
        clt_pops[cl_id] += 1
        model_id = model.file_name.removesuffix(".pdb")
        try:
            model_scan_dt = results_by_model[model_id]
        except KeyError:
            continue
        for mut_result in model_scan_dt:
            ident = ident_builder(mut_result)
            if ident not in clt_scan[cl_id]:
                entry: Dict[str, Union[float, int, list]] = {
                    "delta_score": [],
                    "delta_vdw": [],
                    "delta_elec": [],
                    "delta_desolv": [],
                    "delta_bsa": [],
                    "frac_pr": 0,
                }
                if metadata_builder is not None:
                    entry.update(metadata_builder(mut_result))
                clt_scan[cl_id][ident] = entry
            delta_scores = mut_result.delta_scores
            entry = clt_scan[cl_id][ident]
            entry["delta_score"].append(delta_scores[0])
            entry["delta_vdw"].append(delta_scores[1])
            entry["delta_elec"].append(delta_scores[2])
            entry["delta_desolv"].append(delta_scores[3])
            entry["delta_bsa"].append(delta_scores[4])
            entry["frac_pr"] += 1
    return clt_scan, clt_pops
