"""alascan module."""
import os
import io
import shutil
from pathlib import Path
from contextlib import redirect_stdout
from dataclasses import dataclass
from haddock.core.typing import Optional, Any
from typing import List, Tuple

import numpy as np
import pandas as pd

from haddock import log
from haddock.libs.libalign import get_atoms, load_coords
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


def mutate(pdb_f, target_chain, target_resnum, mut_resname):
    """
    Mutate a residue in a PDB file into a different residue.
    
    Parameters
    ----------
    pdb_f : str
        Path to the pdb file.
    
    target_chain : str
        Chain of the residue to be mutated.
    
    target_resnum : int
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
                resnum = int(line[22:26])
                atom_name = line[12:16].strip()
                if target_chain == chain and target_resnum == resnum:
                    if not resname:
                        resname = line[17:20].strip()
                    if atom_name in ATOMS_TO_BE_MUTATED:
                        # mutate
                        line = line[:17] + mut_resname + line[20:]
                        mut_pdb_l.append(line)
                else:
                    mut_pdb_l.append(line)
    try:
        mut_id = f'{RES_CODES[resname]}{target_resnum}{RES_CODES[mut_resname]}'
    except KeyError:
        raise KeyError(f"Could not mutate {resname} into {mut_resname}.")
    mut_pdb_fname = Path(
        pdb_f.name.replace('.pdb', f'-{target_chain}_{mut_id}.pdb'))
    with open(mut_pdb_fname, 'w') as fh:
        fh.write(''.join(mut_pdb_l))
    return mut_pdb_fname


def add_delta_to_bfactor(pdb_f, df_scan):
    """Add delta scores as b-factors.
    
    Parameters
    ----------
    pdb_f : str
        Path to the pdb file.

    df_scan : pandas.DataFrame
        Dataframe with the scan results for the model
    
    Returns
    -------
    pdb_f : str
        Path to the pdb file with the b-factors added.
    """
    tmp_pdb_f = pdb_f.replace('.pdb', '_bfactor.pdb')
    max_b, min_b = df_scan["delta_score"].max(), df_scan["delta_score"].min()
    out_pdb_l = []
    with open(pdb_f, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('ATOM'):
                chain = line[21]
                resnum = int(line[22:26])
                norm_delta = 0.0
                # extracting all the elements of df_scan such that
                # chain = chain and res = resnum
                df_scan_subset = df_scan.loc[
                    (df_scan["chain"] == chain) & (df_scan["res"] == resnum)
                    ]
                if df_scan_subset.shape[0] > 0:
                    delta = df_scan_subset["delta_score"].values[0]
                    norm_delta = 100 * (delta - min_b) / (max_b - min_b)

                delta_str = f"{norm_delta:.2f}".rjust(6, " ")
                line = line[:60] + delta_str + line[66:]
            out_pdb_l.append(line)
    with open(tmp_pdb_f, 'w') as out_fh:
        out_fh.write(''.join(out_pdb_l))
    # move tmp_pdb_f to pdb_f
    os.rename(tmp_pdb_f, pdb_f)
    return pdb_f


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


def alascan_cluster_analysis(models):
    """Perform cluster analysis on the alascan data.
    
    Parameters
    ----------
    models : list
        List of models.
    
    path : str
        Path to the run directory.
    """
    clt_scan = {}
    cl_pops = {}
    for native in models:
        cl_id = native.clt_id
        # unclustered models have cl_id = None
        if cl_id is None:
            cl_id = "unclustered"
        # Initiate key if this cluster id is encountered for the first time
        if cl_id not in clt_scan:
            clt_scan[cl_id] = {}
            cl_pops[cl_id] = 0
        # Increase the population of that cluster
        cl_pops[cl_id] += 1
        # read the scan file
        alascan_fname = f"scan_{native.file_name.removesuffix('.pdb')}.tsv"
        df_scan = pd.read_csv(alascan_fname, sep="\t", comment="#")
        # loop over the scan file
        for row_idx in range(df_scan.shape[0]):
            row = df_scan.iloc[row_idx]
            chain = row['chain']
            res = row['res']
            ori_resname = row['ori_resname']
            delta_score = row['delta_score']
            delta_vdw = row['delta_vdw']
            delta_elec = row['delta_elec']
            delta_desolv = row['delta_desolv']
            delta_bsa = row['delta_bsa']
            # add to the cluster data with the ident logic
            ident = f"{chain}-{res}-{ori_resname}"
            # Create variable with appropriate key
            if ident not in clt_scan[cl_id].keys():
                clt_scan[cl_id][ident] = {
                    'delta_score': [],
                    'delta_vdw': [],
                    'delta_elec': [],
                    'delta_desolv': [],
                    'delta_bsa': [],
                    'frac_pr': 0,
                    }
            # Add data
            clt_scan[cl_id][ident]['delta_score'].append(delta_score)
            clt_scan[cl_id][ident]['delta_vdw'].append(delta_vdw)
            clt_scan[cl_id][ident]['delta_elec'].append(delta_elec)
            clt_scan[cl_id][ident]['delta_desolv'].append(delta_desolv)
            clt_scan[cl_id][ident]['delta_bsa'].append(delta_bsa)
            clt_scan[cl_id][ident]['frac_pr'] += 1
    # now average the data for every cluster
    for cl_id in clt_scan:
        scan_clt_filename = f"scan_clt_{cl_id}.tsv"
        log.info(f"Writing {scan_clt_filename}")
        clt_data = []
        # Loop over residues
        for ident in clt_scan[cl_id]:
            # Split identifyer to retrieve residue data
            chain = ident.split("-")[0]
            resnum = int(ident.split("-")[1])
            resname = ident.split("-")[2]
            # Point data for this specific residue
            clt_res_dt = clt_scan[cl_id][ident]
            # Compute averages and stddev and hold data.
            clt_data.append([
                chain,
                resnum,
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
                clt_res_dt['frac_pr'] / cl_pops[cl_id],
                ])
        df_cols = [
            'chain', 'resnum', 'resname', 'full_resname',
            'delta_score', 'delta_score_std', 'delta_vdw', 'delta_vdw_std',
            'delta_elec', 'delta_elec_std', 'delta_desolv', 'delta_desolv_std',
            'delta_bsa', 'delta_bsa_std', 'frac_pres',
            ]
        df_scan_clt = pd.DataFrame(clt_data, columns=df_cols)
        # adding clt-based Z score
        df_scan_clt = add_zscores(df_scan_clt, 'delta_score')

        df_scan_clt.sort_values(by=['chain', 'resnum'], inplace=True)
        df_scan_clt.to_csv(
            scan_clt_filename,
            index=False,
            float_format='%.2f',
            sep="\t")

        # add comment
        fl_content = open(scan_clt_filename, 'r').read()
        with open(scan_clt_filename, 'w') as f:
                f.write(f"{'#' * 80}{os.linesep}")  # noqa E501
                f.write(f"# `alascan` cluster results for cluster {cl_id}{os.linesep}")  # noqa E501
                f.write(f"# reported values are the average for the cluster{os.linesep}")  # noqa E501
                f.write(f"#{os.linesep}")
                f.write(f"# z_score is calculated with respect to the mean values of all residues{os.linesep}")  # noqa E501
                f.write(f"{'#' * 80}{os.linesep}")  # noqa E501
                f.write(fl_content)
    return clt_scan


def generate_alascan_output(models, path):
    """Generate the alascan output files.
    
    Parameters
    ----------
    models : list
        List of models.
    path : str
        Path to the run directory.
    """
    models_to_export = []
    for model in models:
        name = f"{model.file_name.removesuffix('.pdb')}_alascan.pdb"
        # changing attributes
        name_path = Path(name)
        shutil.copy(Path(model.path, model.file_name), name_path) 
        alascan_fname = f"scan_{model.file_name.removesuffix('.pdb')}.tsv"
        # add delta_score as a bfactor to the model
        df_scan = pd.read_csv(alascan_fname, sep="\t", comment="#")
        add_delta_to_bfactor(name, df_scan)
        model.ori_name = model.file_name
        model.file_name = name
        model.full_name = name
        model.rel_path = Path('..', Path(path).name, name)
        model.path = str(Path(".").resolve())
        models_to_export.append(model)
    return models_to_export


def create_alascan_plots(clt_alascan, scan_residue, offline = False):
    """Create the alascan plots."""
    for clt_id in clt_alascan:            
        scan_clt_filename = f"scan_clt_{clt_id}.tsv"
        if not os.path.exists(scan_clt_filename):
            log.warning(f"Could not find {scan_clt_filename}")
            continue
        df_scan_clt = pd.read_csv(
            scan_clt_filename,
            sep="\t",
            comment="#"
            )
        # plot the data
        try:
            make_alascan_plot(
                df_scan_clt,
                clt_id,
                scan_residue,
                offline=offline,
                )
        except Exception as e:
            log.warning(
                "Could not create interactive plot. The following error"
                f" occurred {e}"
                )
    return


def write_scan_out(results: List[Any], model_id: str) -> None:
    """
    Save mutation results per model to tsv file.
    
    Parameters
    ----------
    results : List[Any]
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
        print(f'No scan results for model {model_id}')
        return

    # Convert scan output to dataframe
    scan_data = []
    native_score = None
    for result in results:
        if result.success:
            m_score, m_vdw, m_elec, m_des, m_bsa = result.mutant_scores
            d_score, d_vdw, d_elec, d_des, d_bsa = result.delta_scores
            
            scan_data.append([
                result.chain, result.res_num, result.ori_resname, 
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
        df_scan = add_zscores(df_scan, 'delta_score')
        
        # Sort by chain id, then by residue id
        df_scan.sort_values(by=['chain', 'res'], inplace=True)

        # Save to tsv (per model)
        output_file = f"scan_{model_id}.tsv"
        df_scan.to_csv(output_file, index=False, float_format='%.2f', sep="\t")
        fl_content = open(output_file, 'r').read()
        with open(output_file, 'w') as f:
            f.write(f"{'#' * 70}{os.linesep}")
            f.write(f"# `alascan` results for {model_id}\n")
            f.write(f"#\n")
            f.write(f"# native score = {native_score}\n")
            f.write(f"#\n")
            f.write(f"# z_score is calculated with respect to the other residues\n")
            f.write(f"{'#' * 70}{os.linesep}")
            f.write(fl_content)


@dataclass
class MutationResult:
    """Result from a single mutation."""
    model_id: str
    chain: str
    res_num: int
    ori_resname: str
    target_resname: str
    # components of "mutant_scores": score, vdw, elec, desolv, bsa
    mutant_scores: Tuple[float, float, float, float, float]  
    delta_scores: Tuple[float, float, float, float, float]   
    success: bool
    error_msg: Optional[str] = None


class InterfaceScanner:
    """Scan interface of a model to get tartget residues and create corresponding mutation jobs."""
    
    def __init__(self, model, mutation_res="ALA", params=None, library_mode=True):
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
            if key.startswith("resdic")}
        try:
            self.model_path = model.rel_path
            self.model_id = model.file_name.removesuffix('.pdb')
        except AttributeError:
            # for library mode
            self.model_path = Path(model)
            self.model_id = self.model_path.stem
      
    def run(self):
        """
        Get interface residues and create mutation jobs.
        If library_mode=True, also execute the mutations sequentially.
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
            coords, chain_ranges = load_coords(self.model_path, atoms, add_resname=True)
            
            # Determine target residues - user-given or all interface (same logic as original Scan)
            if self.filter_resdic != {'_': []}:
                # User-specified residues
                interface = {}
                for chain in self.filter_resdic:
                    if chain in chain_ranges:
                        chain_aas = [aa[1] for aa in coords if aa[0] == chain]
                        unique_aas = list(set(chain_aas))
                        # the interface here is the intersection of the
                        # residues in filter_resdic and the residues actually 
                        # present in the model
                        interface[chain] = [
                            res for res in self.filter_resdic[chain]
                            if res in unique_aas]

            else:
                # all interface
                cutoff = self.params.get("int_cutoff", 5.0)
                interface = CAPRI.identify_interface(self.model_path, cutoff=cutoff)
                # in case the user wants to scan only some chains (this is
                # superseded by filter_resdic)
                chains = self.params.get("chains", [])
                if chains != []:
                    interface = {
                        chain: interface[chain] for chain in chains
                        if chain in interface}
            
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
                    # Skip if scan_residue is the same as original (e.g. skip ALA->ALA)
                    if ori_resname != end_resname:
                        job = ModelPointMutation(
                            model_path=self.model_path,
                            model_id=self.model_id,
                            chain=chain,
                            res_num=res,
                            ori_resname=ori_resname,
                            target_resname=end_resname,
                            native_scores=native_scores,
                            output_mutants=output_mutants
                        )
                        self.point_mutations_jobs.append(job)
            
            # Execute jobs if in library mode
            if self.library_mode:
                log.info(f"Executing {len(self.point_mutations_jobs)} mutations for {self.model_id}")
                results = []
                total = len(self.point_mutations_jobs)
                
                for i, job in enumerate(self.point_mutations_jobs, 1):
                    log.info(f"Processing mutation {i}/{total}: {job.chain}:{job.res_num} {job.ori_resname}->{job.target_resname}")
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
        res_num: int, 
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
        res_num : int
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
        self.res_num = res_num
        self.ori_resname = ori_resname
        self.target_resname = target_resname
        self.native_scores = native_scores
        self.output_mutants = output_mutants
    
    def run(self):
        """Execute the point mutation."""
        mutation_id = f"{self.model_id}_{self.chain}{self.res_num}{self.target_resname}"
        
        try:
            # Setup working directory
            sc_dir = f"haddock3-score-{mutation_id}"
            os.makedirs(sc_dir, exist_ok=True)
            
            # Perform mutation
            mut_pdb = mutate(self.model_path, self.chain, self.res_num, self.target_resname)
            
            # Calculate mutant scores
            mutant_scores = calc_score(mut_pdb, run_dir=sc_dir, outputpdb=self.output_mutants)
            
            # Calculate deltas (native - mutant)
            n_score, n_vdw, n_elec, n_des, n_bsa = self.native_scores
            m_score, m_vdw, m_elec, m_des, m_bsa = mutant_scores
            delta_scores = (n_score - m_score, n_vdw - m_vdw, n_elec - m_elec, 
                          n_des - m_des, n_bsa - m_bsa)

            # Handle output files
            em_mut_pdb = Path(mut_pdb.stem + '_hs.pdb')
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
                model_id=self.model_id, chain=self.chain, res_num=self.res_num,
                ori_resname=self.ori_resname, target_resname=self.target_resname,
                mutant_scores=mutant_scores, delta_scores=delta_scores, success=True
            )
            
        except Exception as e:
            return MutationResult(
                model_id=self.model_id, chain=self.chain, res_num=self.res_num,
                ori_resname=self.ori_resname, target_resname=self.target_resname,
                mutant_scores=(0, 0, 0, 0, 0), delta_scores=(0, 0, 0, 0, 0),
                success=False, error_msg=str(e)
            )