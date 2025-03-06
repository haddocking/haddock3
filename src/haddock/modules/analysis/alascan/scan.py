"""alascan module."""
import os
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import io
from contextlib import redirect_stdout

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

def get_score_string(pdb_f, run_dir):
    """Get score output from cli_score.main.

    Parameters
    ----------
    pdb_f : str
        Path to the pdb file.
    
    run_dir : str
        Path to the run directory.
    
    Returns
    -------
    out : list
        List of strings with the score output.
    """
    f = io.StringIO()
    with redirect_stdout(f):
        cli_score.main(pdb_f, run_dir, full=True)
    out = f.getvalue().split(os.linesep)
    return out


def calc_score(pdb_f, run_dir):
    """Calculate the score of a model.

    Parameters
    ----------
    pdb_f : str
        Path to the pdb file.
    run_dir : str
        Path to the run directory.
    
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
    out_string = get_score_string(pdb_f, run_dir)

    for ln in out_string:
        if ln.startswith("> HADDOCK-score (emscoring)"):
            score = float(ln.split()[-1])
        if ln.startswith("> vdw"):
            vdw = float(ln.split("vdw=")[1].split(",")[0])
            elec = float(ln.split("elec=")[1].split(",")[0])
            desolv = float(ln.split("desolv=")[1].split(",")[0])
            bsa = float(ln.split("bsa=")[1].split(",")[0])
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
            cl_id = "-"
        if cl_id not in clt_scan:
            clt_scan[cl_id] = {}
            cl_pops[cl_id] = 1
        else:
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
            if ident not in clt_scan[cl_id]:
                clt_scan[cl_id][ident] = {
                    'delta_score': delta_score,
                    'delta_vdw': delta_vdw,
                    'delta_elec': delta_elec,
                    'delta_desolv': delta_desolv,
                    'delta_bsa': delta_bsa,
                    'frac_pr': 1
                    }
            else:
                clt_scan[cl_id][ident]['delta_score'] += delta_score
                clt_scan[cl_id][ident]['delta_vdw'] += delta_vdw
                clt_scan[cl_id][ident]['delta_elec'] += delta_elec
                clt_scan[cl_id][ident]['delta_desolv'] += delta_desolv
                clt_scan[cl_id][ident]['delta_bsa'] += delta_bsa
                clt_scan[cl_id][ident]['frac_pr'] += 1
    # now average the data
    for cl_id in clt_scan:
        scan_clt_filename = f"scan_clt_{cl_id}.tsv"
        log.info(f"Writing {scan_clt_filename}")
        clt_data = []
        for ident in clt_scan[cl_id]:
            chain = ident.split("-")[0]
            resnum = int(ident.split("-")[1])
            resname = ident.split("-")[2]
            frac_pr = clt_scan[cl_id][ident]["frac_pr"]
            clt_data.append([chain, resnum, resname, ident,
                             clt_scan[cl_id][ident]['delta_score'] / frac_pr,
                             clt_scan[cl_id][ident]['delta_vdw'] / frac_pr,
                             clt_scan[cl_id][ident]['delta_elec'] / frac_pr,
                             clt_scan[cl_id][ident]['delta_desolv'] / frac_pr,
                             clt_scan[cl_id][ident]['delta_bsa'] / frac_pr,
                             clt_scan[cl_id][ident]['frac_pr'] / cl_pops[cl_id]
                             ]
                            )
        df_cols = ['chain', 'resnum', 'resname', 'full_resname', 'delta_score',
                   'delta_vdw', 'delta_elec', 'delta_desolv', 'delta_bsa',
                   'frac_pres']
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
                f.write(f"#######################################################################{os.linesep}")  # noqa E501
                f.write(f"# `alascan` cluster results for cluster {cl_id}{os.linesep}")  # noqa E501
                f.write(f"#{os.linesep}")
                f.write(f"# z_score is calculated with respect to the mean values of all residues{os.linesep}")  # noqa E501
                f.write(f"#######################################################################{os.linesep}")  # noqa E501
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


class ScanJob:
    """A Job dedicated to the parallel alanine scanning of models."""

    def __init__(
            self,
            params,
            scan_obj):

        log.info(f"core {scan_obj.core}, initialising Scan...")
        self.params = params
        self.scan_obj = scan_obj

    def run(self):
        """Run this ScanJob."""
        log.info(f"core {self.scan_obj.core}, running Scan...")
        self.scan_obj.run()
        return


class Scan:
    """Scan class."""

    def __init__(
            self,
            model_list,
            core,
            path,
            **params,
            ):
        """Initialise Scan class."""
        self.model_list = model_list
        self.core = core
        self.path = path
        self.scan_res = params['params']['scan_residue']
        self.int_cutoff = params["params"]["int_cutoff"]
        self.chains = params["params"]["chains"]
        self.output_mutants = params["params"]["output_mutants"]
        # initialising resdic
        if "params" in params.keys():
            self.filter_resdic = {
                key[-1]: value for key, value
                in params["params"].items()
                if key.startswith("resdic")
                }


    def run(self):
        """Run alascan calculations."""
        for native in self.model_list:
            # here we rescore the native model for consistency, as the score
            # attribute could come from any module in principle
            sc_dir = f"haddock3-score-{self.core}"
            n_score, n_vdw, n_elec, n_des, n_bsa = calc_score(native.rel_path,
                                                              run_dir=sc_dir)
            scan_data = []

            # load the coordinates
            atoms = get_atoms(native.rel_path)
            coords, chain_ranges = load_coords(native.rel_path,
                                               atoms,
                                               add_resname=True
                                               )
            
            # check if the user wants to mutate only some residues
            if self.filter_resdic != {'_': []}:
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
                            if res in unique_aas
                            ]
            else:
                interface = CAPRI.identify_interface(
                    native.rel_path,
                    cutoff=self.int_cutoff
                    )
                # in case the user wants to scan only some chains (this is
                # superseded by filter_resdic)
                if self.chains != []:
                    interface = {
                        chain: interface[chain] for chain in self.chains
                        if chain in interface
                        }
            
            resname_dict = {}
            for chain, resid, _atom, resname in coords.keys():
                key = f"{chain}-{resid}"
                if key not in resname_dict:
                    resname_dict[key] = resname
            
            # loop over the interface
            for chain in interface:
                for res in interface[chain]:
                    ori_resname = resname_dict[f"{chain}-{res}"]
                    end_resname = self.scan_res
                    if ori_resname == self.scan_res:
                        # we do not re-score equal residues (e.g. ALA = ALA)
                        c_score = n_score
                        c_vdw = n_vdw
                        c_elec = n_elec
                        c_des = n_des
                        c_bsa = n_bsa
                    else:
                        try:
                            mut_pdb_name = mutate(native.rel_path,
                                                  chain,
                                                  res,
                                                  end_resname)
                        except KeyError:
                            continue
                        # now we score the mutated model
                        c_score, c_vdw, c_elec, c_des, c_bsa = calc_score(
                            mut_pdb_name,
                            run_dir=sc_dir)
                        # now the deltas (wildtype - mutant)
                        delta_score = n_score - c_score
                        delta_vdw = n_vdw - c_vdw
                        delta_elec = n_elec - c_elec
                        delta_desolv = n_des - c_des
                        delta_bsa = n_bsa - c_bsa
 
                        scan_data.append([chain, res, ori_resname, end_resname,
                                          c_score, c_vdw, c_elec, c_des,
                                          c_bsa, delta_score,
                                          delta_vdw, delta_elec, delta_desolv,
                                          delta_bsa])
                        # if self.output_mutants is false, remove the 
                        # mutated pdb file
                        if not self.output_mutants:
                            os.remove(mut_pdb_name)
            # write output
            df_columns = ['chain', 'res', 'ori_resname', 'end_resname',
                          'score', 'vdw', 'elec', 'desolv', 'bsa',
                          'delta_score', 'delta_vdw', 'delta_elec',
                          'delta_desolv', 'delta_bsa']
            self.df_scan = pd.DataFrame(scan_data, columns=df_columns)
            alascan_fname = Path(self.path, f"scan_{native.file_name.removesuffix('.pdb')}.tsv")
            # add zscore
            self.df_scan = add_zscores(self.df_scan, 'delta_score')

            self.df_scan.to_csv(
                alascan_fname,
                index=False,
                float_format='%.2f',
                sep="\t"
                )

            fl_content = open(alascan_fname, 'r').read()
            with open(alascan_fname, 'w') as f:
                f.write(f"##########################################################{os.linesep}")  # noqa E501
                f.write(f"# `alascan` results for {native.file_name}{os.linesep}")  # noqa E501
                f.write(f"#{os.linesep}")
                f.write(f"# native score = {n_score}{os.linesep}")
                f.write(f"#{os.linesep}")
                f.write(f"# z_score is calculated with respect to the other residues")  # noqa E501
                f.write(f"{os.linesep}")
                f.write(f"##########################################################{os.linesep}")  # noqa E501
                f.write(fl_content)
