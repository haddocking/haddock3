"""alascan module."""
import os
import io
import shutil
from pathlib import Path
from contextlib import redirect_stdout
from dataclasses import dataclass
from typing import List, Tuple, Dict

import numpy as np
import pandas as pd


from haddock import log
from haddock.core.typing import Any, Optional, Union
from haddock.libs.libalign import get_atoms, load_coords
from haddock.libs.libontology import PDBFile
from haddock.libs.libplots import make_alascan_plot
from haddock.modules.analysis.caprieval.capri import CAPRI
from haddock.clis import cli_score

ATOMS_TO_BE_MUTATED = ['C', 'N', 'CA', 'O', 'CB']

RES_CODES = dict([
    ("CYS", "C"),
    ("ASP", "D"),
    ("SER", "S"),
    ("GLN", "Q"),
    ("LYS", "K"),
    ("ILE", "I"),
    ("PRO", "P"),
    ("THR", "T"),
    ("PHE", "F"),
    ("ASN", "N"),
    ("GLY", "G"),
    ("HIS", "H"),
    ("LEU", "L"),
    ("ARG", "R"),
    ("TRP", "W"),
    ("ALA", "A"),
    ("VAL", "V"),
    ("GLU", "E"),
    ("TYR", "Y"),
    ("MET", "M"),
    ("ALY", "K"),
    ("ASH", "D"),
    ("CFE", "C"),
    ("CSP", "C"),
    ("CYC", "C"),
    ("CYF", "C"),
    ("CYM", "C"),
    ("DDZ", "A"),
    ("GLH", "E"),
    ("HLY", "P"),
    ("HY3", "P"),
    ("HYP", "P"),
    ("M3L", "K"),
    ("MLY", "K"),
    ("MLZ", "K"),
    ("MSE", "M"),
    ("NEP", "H"),
    ("PNS", "S"),
    ("PTR", "Y"),
    ("SEP", "S"),
    ("TOP", "T"),
    ("TYP", "Y"),
    ("TYS", "Y"),
    ("CIR", "R"),
    ])


@dataclass
class MutationResult:
    """Result from a single mutation."""
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


def mutate(pdb_f, target_chain, target_resid, mut_resname):
    """
    Mutate a residue in a PDB file into a different residue.
    
    Parameters
    ----------
    pdb_f : str
        Path to the pdb file.
    
    target_chain : str
        Chain of the residue to be mutated.
    
    target_resid : int
        Residue number of the residue to be mutated.
    
    mut_resname : str
        Residue name of the residue to be mutated.

    Returns
    -------
    mut_pdb_fname : str
        Path to the mutated pdb file.
    """
    mut_pdb_l = []
    resname = ''
    with open(pdb_f, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('ATOM'):
                chain = line[21]
                resid = int(line[22:26])
                atom_name = line[12:16].strip()
                if target_chain == chain and target_resid == resid:
                    if not resname:
                        resname = line[17:20].strip()
                    if atom_name in ATOMS_TO_BE_MUTATED:
                        # mutate
                        line = line[:17] + mut_resname + line[20:]
                        mut_pdb_l.append(line)
                else:
                    mut_pdb_l.append(line)
    try:
        mut_id = f'{RES_CODES[resname]}{target_resid}{RES_CODES[mut_resname]}'
    except KeyError:
        raise KeyError(f"Could not mutate {resname} into {mut_resname}.")
    mut_pdb_fname = Path(
        pdb_f.name.replace('.pdb', f'-{target_chain}_{mut_id}.pdb'))
    with open(mut_pdb_fname, 'w') as fh:
        fh.write(''.join(mut_pdb_l))
    return mut_pdb_fname


def get_score_string(pdb_f, run_dir, outputpdb=False):
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
    
    Returns
    -------
    out : list
        List of strings with the score output.
    """
    f = io.StringIO()
    with redirect_stdout(f):
        cli_score.main(pdb_f, run_dir, full=True, outputpdb=outputpdb)
    out = f.getvalue().split(os.linesep)
    return out


def calc_score(pdb_f, run_dir, outputpdb=False):
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
    
    Returns
    -------
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
    """
    out_string = get_score_string(pdb_f, run_dir, outputpdb=outputpdb)

    for ln in out_string:
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


def add_zscores(df_scan_clt, column='delta_score'):
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
        df_scan_clt['z_score'] = (df_scan_clt[column] - mean_delta) / std_delta
    else:
        df_scan_clt['z_score'] = 0.0
    return df_scan_clt


class ClusterOutputer():
    """Manage the generation of alascan outputs for cluster-based analysis."""
    def __init__(
            self,
            cluster_scan_data: Dict[str, Dict[str, Union[float, int]]],
            clt_id: str,
            clt_population: int,
            scan_residue: str = "ALA",
            generate_plot: bool = False,
            offline: bool = False,
            ):
        """Initialization function

        Parameters
        ----------
        cluster_scan_data : Dict[str, Dict[str, Union[float, int]]]
            Dictionary containing alascan data per residue identifier.
        clt_id : str
            Cluster identifier
        clt_population : int
            Number of entries in this cluster.
        scan_residue : str, optional
            Residue scanned, by default "ALA"
        generate_plot : bool, optional
            Defines if a plot must be generated, by default False
        offline : bool, optional
            Defines if the plot should be functional offline, by default False
        """
        self.cluster_scan_data = cluster_scan_data
        self.clt_id = clt_id
        self.clt_population = clt_population
        self.scanned_residue = scan_residue
        self.generate_plot = generate_plot
        self.offline = offline

    def run(self) -> str:
        """Wrtie cluster alascan output to scan_clt_X.tsv file,
        including average and stdard deviation data per residue
        and optionally save cluster alascan plots (if generate_plot == True)
        
        Return
        ------
        scan_clt_filename : str
            Name of the tsv file written
        """
        # Gather all data in a list
        clt_data = []
        # Loop over residues
        for ident, clt_res_dt in self.cluster_scan_data.items():
            # Split identifier to retrieve residue data
            chain = ident.split("-")[0]
            resid = int(ident.split("-")[1])
            resname = ident.split("-")[2]
            # Compute averages and stddev and hold data.
            clt_data.append([
                chain,
                resid,
                resname,
                ident,
                np.mean(clt_res_dt['delta_score']),
                np.std(clt_res_dt['delta_score']),
                np.mean(clt_res_dt['delta_vdw']),
                np.std(clt_res_dt['delta_vdw']),
                np.mean(clt_res_dt['delta_elec']),
                np.std(clt_res_dt['delta_elec']),
                np.mean(clt_res_dt['delta_desolv']),
                np.std(clt_res_dt['delta_desolv']),
                np.mean(clt_res_dt['delta_bsa']),
                np.std(clt_res_dt['delta_bsa']),
                clt_res_dt['frac_pr'] / self.clt_population,
                ])
        df_cols = [
            'chain', 'resid', 'resname', 'full_resname',
            'delta_score', 'delta_score_std', 'delta_vdw', 'delta_vdw_std',
            'delta_elec', 'delta_elec_std', 'delta_desolv', 'delta_desolv_std',
            'delta_bsa', 'delta_bsa_std', 'frac_pres',
            ]
        df_scan_clt = pd.DataFrame(clt_data, columns=df_cols)
        # adding clt-based Z score
        df_scan_clt = add_zscores(df_scan_clt, 'delta_score')
        # Sort rows
        df_scan_clt.sort_values(by=['chain', 'resid'], inplace=True)
        # Generate output CSV data
        csv_data = io.StringIO()
        df_scan_clt.to_csv(
            csv_data,
            index=False,
            float_format='%.2f',
            sep="\t")
        csv_data.seek(0)
        # Define output csv filepath
        scan_clt_filename = f"scan_clt_{self.clt_id}.tsv"
        # Write the file
        with open(scan_clt_filename, "w") as fout:
            # add comment to the file
            fout.write(f"{'#' * 80}{os.linesep}")
            fout.write(f"# `alascan` cluster results for cluster {self.clt_id}{os.linesep}")  # noqa E501
            fout.write(f"# reported values are the average for the cluster{os.linesep}")  # noqa E501
            fout.write(f"#{os.linesep}")
            fout.write(f"# z_score is calculated with respect to the mean values of all residues{os.linesep}")  # noqa E501
            fout.write(f"{'#' * 80}{os.linesep}")
            # Write csv data content
            fout.write(csv_data.read())
        
        # Generate plot
        self.gen_alascan_plot(df_scan_clt)

        return scan_clt_filename

    def gen_alascan_plot(self, df_scan_clt: pd.DataFrame) -> None:
        """Generate the alascan plot based on provided data.

        Parameters
        ----------
        df_scan_clt : pd.DataFrame
            The data frame containing the data to be plotted.
        """
        # Check if the plot must be generated
        if not self.generate_plot:
            return

        # Try to plot the data
        try:
            make_alascan_plot(
                df_scan_clt,
                self.clt_id,
                self.scanned_residue,
                offline=self.offline,
                )
        except Exception as e:
            log.warning(
                "Could not create interactive plot. The following error"
                f" occurred {e}"
                )


class AddDeltaBFactor():
    """Add delta score in bfactor column of a PDB."""
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
        output_fname = f"{model_fname}_alascan.pdb"
        # Add delta_score as a bfactor to the model
        self.write_delta_score_to_pdb(output_fname)
        # Update attributes of the model
        self.model.ori_name = self.model.file_name
        self.model.file_name = output_fname
        self.model.full_name = output_fname
        self.model.rel_path = Path("..", Path(self.path).name, output_fname)
        self.model.path = str(Path(".").resolve())
        return self.model
    
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

    def reorder_results(self) -> None:
        """Perform initial data manuputation to simply downstream access.
        
        Organise mutation results into dictionary, with chain-resid keys 
        and values delta score (e.g.: {"A-115": 5}),
        and determine min and max score within the model to normalize data.
        """
        self.model_results: Dict[str, float] = {}
        all_delta_scores: List[float] = []
        # Loop over mutation results
        for mut_result in self.input_results:
            # Create key
            chain_res_key = f"{mut_result.chain}-{mut_result.resid}"
            # Point delta score value
            delta_score = mut_result.delta_scores[0]
            # Hold data
            self.model_results[chain_res_key] = delta_score
            all_delta_scores.append(delta_score)
        # Obtain min and max delta score to be able to normalize later
        if len(all_delta_scores) >= 1:
            self.min_score = min(all_delta_scores)
            self.max_score = max(all_delta_scores)
        else:
            self.min_score = -1
            self.max_score = 1
    
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


def write_scan_out(results: List[MutationResult], model_id: str) -> None:
    """
    Save mutation results per model to tsv file.
    
    Parameters
    ----------
    results : List[MutationResult]
        List of mutation results from scanning (comes from MutationResult)
    model_id : str
        Identifier for the model used in filename.
        
    Returns
    -------
    None
        Function saves results to file `scan_{model_id}.tsv`.
        
    Notes
    -----
    If results list is empty, no file is created.
    """
    if not results:
        print(f"No scan results for model {model_id}")
        return

    # Convert scan output to dataframe
    scan_data = []
    native_score = None
    for result in results:
        if result.success:
            m_score, m_vdw, m_elec, m_des, m_bsa = result.mutant_scores
            d_score, d_vdw, d_elec, d_des, d_bsa = result.delta_scores
            
            scan_data.append([
                result.chain, result.resid, result.ori_resname, 
                result.target_resname, m_score, m_vdw, m_elec, m_des, m_bsa,
                d_score, d_vdw, d_elec, d_des, d_bsa
            ])
            
            # Get native score from first result (mutant + delta = native)
            if native_score is None:
                native_score = result.mutant_scores[0] + result.delta_scores[0]
    
    if scan_data:
        df_columns = ['chain', 'res', 'ori_resname', 'end_resname',
                     'score', 'vdw', 'elec', 'desolv', 'bsa',
                     'delta_score', 'delta_vdw', 'delta_elec',
                     'delta_desolv', 'delta_bsa']
        
        df_scan = pd.DataFrame(scan_data, columns=df_columns)
        # Compute z-scores
        df_scan = add_zscores(df_scan, 'delta_score')
        
        # Sort by chain id, then by residue id
        df_scan.sort_values(by=['chain', 'res'], inplace=True)

        # Get csv data as string
        csv_io = io.StringIO()
        df_scan.to_csv(csv_io, index=False, float_format="%.2f", sep="\t")
        csv_io.seek(0)
        # Save to tsv
        output_file = f"scan_{model_id}.tsv"
        with open(output_file, "w") as fout:
            # Write comments at start of the file
            fout.write(f"{'#' * 80}{os.linesep}")
            fout.write(f"# `alascan` results for {model_id}{os.linesep}")
            fout.write(f"#{os.linesep}")
            fout.write(f"# native score = {native_score}{os.linesep}")
            fout.write(f"#{os.linesep}")
            fout.write(f"# z_score is calculated with respect to the other residues{os.linesep}")
            fout.write(f"{'#' * 80}{os.linesep}")
            # Write csv content
            fout.write(csv_io.read())


def group_scan_by_cluster(
        models: List[PDBFile],
        results_by_model: Dict[str, MutationResult]
        ) -> Tuple[
            Dict[str, Dict[str, Dict[str, Union[float, int]]]],
            Dict[str, int]
            ]:
    """Group models alascan data per cluster.

    Parameters
    ----------
    models : List[PDBFile]
        List of input models

    Returns
    -------
    clt_scan : Dict[str, Dict[str, Dict[str, Union[float, int]]]]
        Dictionary containing alascan data for each cluster, grouped by
        residue identifyer.
    clt_pops: Dict[str, int]
        Dictionary containing number of entries for each cluster.
    """
    # Define holders
    clt_scan: Dict[str, Dict[str, Dict[str, Union[float, int]]]] = {}
    clt_pops: Dict[str, int] = {}
    # Loop over models
    for model in models:
        # Point cluster id
        cl_id = model.clt_id
        # unclustered models have cl_id = None
        if cl_id is None:
            cl_id = "unclustered"
        # Initiate key if this cluster id is encountered for the first time
        if cl_id not in clt_scan:
            clt_scan[cl_id] = {}
            clt_pops[cl_id] = 0
        # Increase the population of that cluster
        clt_pops[cl_id] += 1
        # Point scan results for this model
        model_id = model.file_name.removesuffix(".pdb")
        try:
            model_scan_dt = results_by_model[model_id]
        except KeyError:
            pass
        # Loop over the mutation results
        for mut_result in model_scan_dt:
            # Extract data related to residue information
            chain = mut_result.chain
            res = mut_result.resid
            ori_resname = mut_result.ori_resname
            # Define unique string identifying this residue
            ident = f"{chain}-{res}-{ori_resname}"
            # Create variable with appropriate key
            if ident not in clt_scan[cl_id].keys():
                clt_scan[cl_id][ident] = {
                    "delta_score": [],
                    "delta_vdw": [],
                    "delta_elec": [],
                    "delta_desolv": [],
                    "delta_bsa": [],
                    "frac_pr": 0,
                    }
            # Add data
            delta_scores = mut_result.delta_scores
            clt_scan[cl_id][ident]["delta_score"].append(delta_scores[0])
            clt_scan[cl_id][ident]["delta_vdw"].append(delta_scores[1])
            clt_scan[cl_id][ident]["delta_elec"].append(delta_scores[2])
            clt_scan[cl_id][ident]["delta_desolv"].append(delta_scores[3])
            clt_scan[cl_id][ident]["delta_bsa"].append(delta_scores[4])
            clt_scan[cl_id][ident]["frac_pr"] += 1
    return clt_scan, clt_pops


class InterfaceScanner:
    """Scan interface of a model to get tartget residues and create
    corresponding mutation jobs.
    """

    def __init__(
        self, 
        model: Union[str, Path, Any], 
        mutation_res: str = "ALA", 
        params: Optional[Dict[str, Any]] = None, 
        library_mode: bool = True
    ) -> None:    
        """
        Initialize InterfaceScanner for a single model.
        
        Parameters
        ----------
        model : str, Path, or model object
            HADDOCK model object (if library_mode = False) or 
            Path to PDB file (if library mode = True)
        mutation_res : str
            Target residue for mutation (default: "ALA")
        params : dict, optional
            Additional parameters for interface detection 
            (list on top of alascan/__inint__.py)
        library_mode : bool
            If True, execute mutations sequentially inside InterfaceScanner.run() 
            If False, just prepare jobs - execution will be taken care of in init.py 
            of alascan module by haddock Engine
        """
        self.model = model
        self.mutation_res = mutation_res
        self.library_mode = library_mode
        self.params = params or {}
        self.point_mutations_jobs = []
        self.filter_resdic = {
            key[-1]: value for key, value
            in self.params.items()
            if key.startswith("resdic")
            }
        if isinstance(model, PDBFile):
            self.model_path = model.rel_path
            self.model_id = model.file_name.removesuffix(".pdb")
        else:
            # for library mode
            self.model_path = Path(model)
            self.model_id = self.model_path.stem
      
    def run(self):
        """
        Get interface residues and create mutation jobs.
        If library_mode=True, also execute the mutations sequentially.

        Returns
        -------
        Optional[List[ModelPointMutation]]
            List of mutation jobs if library_mode=False, None if library_mode=True
        """
        try:
            # Calculate native scores
            sc_dir = f"haddock3-score-{self.model_id}-{os.getpid()}"
            try:
                native_scores = calc_score(self.model_path, run_dir=sc_dir, outputpdb=False)
            finally:
                if os.path.exists(sc_dir):
                    shutil.rmtree(sc_dir)
            
            # Load coordinates
            atoms = get_atoms(self.model_path)
            coords, _chain_ranges = load_coords(
                self.model_path,
                atoms,
                add_resname=True,
                )
            
            # Determine target residues: get interface, then apply user filers, if given 
            # Get all interface residues        
            cutoff = self.params.get("int_cutoff", 5.0)
            interface = CAPRI.identify_interface(self.model_path, cutoff=cutoff)

            # get user_chains for the check down the line
            user_chains = self.params.get("chains", [])

            # if user defined target residues, check they are in the interface
            if self.filter_resdic != {'_': []}:
                filtered_interface = {}
                for chain in self.filter_resdic:
                    if chain in interface:
                        # Search for the intersection of user queried residues and interface residues
                        user_res_valid = list(set(self.filter_resdic[chain]).intersection(set(interface[chain])))
                        # If at least one residue must be analyzed, add it to residues to be scanned
                        if user_res_valid:
                            filtered_interface[chain] = user_res_valid
                interface = filtered_interface

            # if (user defined target chains) & (no user target residues) - do use user chains
            elif user_chains:
                interface = {chain: res for chain, res in interface.items() if chain in user_chains}

            # get all atoms of the model to verifiy residue type down the line
            resname_dict = {}
            for chain, resid, _atom, resname in coords.keys():
                key = f"{chain}-{resid}"
                if key not in resname_dict:
                    resname_dict[key] = resname

            # Create mutation 
            output_mutants = self.params.get("output_mutants", False)
            for chain in interface:
                for res in interface[chain]:
                    ori_resname = resname_dict[f"{chain}-{res}"]
                    end_resname = self.mutation_res
                    # Skip if scan_residue is the same as original
                    # (e.g. skip ALA->ALA)
                    if ori_resname != end_resname:
                        job = ModelPointMutation(
                            model_path=self.model_path,
                            model_id=self.model_id,
                            chain=chain,
                            resid=res,
                            ori_resname=ori_resname,
                            target_resname=end_resname,
                            native_scores=native_scores,
                            output_mutants=output_mutants
                        )
                        self.point_mutations_jobs.append(job)
    
            # Execute jobs if in library mode
            if self.library_mode:
                log.info(
                    f"Executing {len(self.point_mutations_jobs)} "
                    f"mutations for {self.model_id}"
                    )
                results = []
                total = len(self.point_mutations_jobs)
                
                for i, job in enumerate(self.point_mutations_jobs, 1):
                    log.info(
                        f"Processing mutation {i}/{total}: "
                        f"{job.chain}:{job.resid} {job.ori_resname}"
                        f"->{job.target_resname}"
                        )
                    result = job.run()
                    results.append(result)
                    if result.success:
                        log.info(f"Delta score: {result.delta_scores[0]:.2f}")
                    else:
                        log.warning(f"Failed: {result.error_msg}")
                
                write_scan_out(results, self.model_id)
            else:
                # return point_mutations_jobs back to alascan/__init__.py
                return self.point_mutations_jobs
                
        except Exception as e:
            log.error(f"Failed to scan model {self.model_id}: {e}")
            raise


class ModelPointMutation:
    """Executes a single point mutation."""
    def __init__(
        self, 
        model_path: Path, 
        model_id: str, 
        chain: str, 
        resid: int, 
        ori_resname: str,
        target_resname: str, 
        native_scores: Tuple[float, float, float, float, float], 
        output_mutants: bool = False
    ) -> None:
        """
        Initialize a single point mutation job.
        
        Parameters
        ----------
        model_path : Path
            Path to the PDB file
        model_id : str
            Identifier for the model
        chain : str
            Chain identifier
        resid : int
            Residue number
        ori_resname : str
            Original residue name
        target_resname : str
            Target residue name for mutation
        native_scores : tuple
            Native model scores (score, vdw, elec, desolv, bsa)
        output_mutants : bool
            Whether to keep mutant PDB files
        """
        self.model_path = Path(model_path)
        self.model_id = model_id
        self.chain = chain
        self.resid = resid
        self.ori_resname = ori_resname
        self.target_resname = target_resname
        self.native_scores = native_scores
        self.output_mutants = output_mutants
    
    def run(self):
        """Execute the point mutation."""
        mutation_id = f"{self.model_id}_{self.chain}{self.resid}{self.target_resname}"
        
        try:
            # Setup working directory
            sc_dir = f"haddock3-score-{mutation_id}"
            os.makedirs(sc_dir, exist_ok=True)
            
            # Perform mutation
            mut_pdb = mutate(self.model_path, self.chain, self.resid, self.target_resname)
            
            # Calculate mutant scores
            mutant_scores = calc_score(mut_pdb, run_dir=sc_dir, outputpdb=self.output_mutants)
            
            # Calculate deltas (native - mutant)
            n_score, n_vdw, n_elec, n_des, n_bsa = self.native_scores
            m_score, m_vdw, m_elec, m_des, m_bsa = mutant_scores
            delta_scores = (n_score - m_score, n_vdw - m_vdw, n_elec - m_elec, 
                          n_des - m_des, n_bsa - m_bsa)

            # Handle output files
            em_mut_pdb = Path(f"{mut_pdb.stem}_hs.pdb")
            if not self.output_mutants:
                # if output_mutants = False, then remove both files
                if os.path.exists(mut_pdb):
                    os.remove(mut_pdb)
                if em_mut_pdb.exists():
                    os.remove(em_mut_pdb)
            else:
                # othervise keep energy-minimized pdb
                if os.path.exists(em_mut_pdb):
                    shutil.move(em_mut_pdb, mut_pdb)
            # clean up scoring dir
            if os.path.exists(sc_dir):
                shutil.rmtree(sc_dir)
            
            return MutationResult(
                model_id=self.model_id, chain=self.chain, resid=self.resid,
                ori_resname=self.ori_resname, target_resname=self.target_resname,
                mutant_scores=mutant_scores, delta_scores=delta_scores, success=True
            )
            
        except Exception as e:
            return MutationResult(
                model_id=self.model_id, chain=self.chain, resid=self.resid,
                ori_resname=self.ori_resname, target_resname=self.target_resname,
                mutant_scores=(0, 0, 0, 0, 0), delta_scores=(0, 0, 0, 0, 0),
                success=False, error_msg=str(e)
            )
