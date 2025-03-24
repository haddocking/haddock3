"""Module computing contact maps of complexes, alone or grouped by cluster.

Chord diagram functions were adapted from:
https://plotly.com/python/v3/filled-chord-diagram/
"""


import os
import glob
from pathlib import Path

import numpy as np
import plotly.graph_objs as go
from scipy.spatial.distance import pdist, squareform

from haddock import log
from haddock.libs.libontology import PDBFile
from haddock.libs.libpdb import (
    slc_name,
    slc_resname,
    slc_chainid,
    slc_resseq,
    slc_x,
    slc_y,
    slc_z,
    )
from haddock.core.typing import (
    Any,
    NDFloat,
    NDArray,
    Optional,
    Union,
    SupportsRun,
    )
from haddock.libs.libplots import heatmap_plotly, fig_to_html


###############################
# Global variable definitions #
###############################
RESIDUE_POLARITY = {
    "CYS": "polar",
    "HIS": "polar",
    "ASN": "polar",
    "GLN": "polar",
    "SER": "polar",
    "THR": "polar",
    "TYR": "polar",
    "TRP": "polar",
    "ALA": "apolar",
    "PHE": "apolar",
    "GLY": "apolar",
    "ILE": "apolar",
    "VAL": "apolar",
    "MET": "apolar",
    "PRO": "apolar",
    "LEU": "apolar",
    "GLU": "negative",
    "ASP": "negative",
    "LYS": "positive",
    "ARG": "positive",
    }
DNA_RNA_POLARITY = {
    "A": "negative",
    "DA": "negative",
    "T": "negative",
    "DT": "negative",
    "U": "negative",
    "DU": "negative",
    "C": "negative",
    "DC": "negative",
    "G": "negative",
    "DG": "negative",
    }
RESIDUE_POLARITY.update(DNA_RNA_POLARITY)

PI = np.pi

# Define interaction types colors
CONNECT_COLORS = {
    "polar-polar": (153, 255, 153),
    "polar-apolar": (255, 204, 204),
    "polar-negative": (255, 204, 153),
    "polar-positive": (153, 204, 255),
    "apolar-apolar": (255, 255, 0),
    "apolar-negative": (255, 229, 204),
    "apolar-positive": (204, 229, 255),
    "negative-negative": (255, 127, 0),
    "negative-positive": (0, 204, 0),
    "positive-positive": (255, 127, 0),
    }
# Also add reversed keys order
REVERSED_CONNECT_COLORS_KEYS = {
    '-'.join(k.split('-')[::-1]): v
    for k, v in CONNECT_COLORS.items()
    }
CONNECT_COLORS.update(REVERSED_CONNECT_COLORS_KEYS)

# Colors for each amino-acids
RESIDUES_COLORS = {
    "CYS": "rgba(229, 255, 204, 0.80)",
    "MET": "rgba(229, 255, 204, 0.80)",
    "ASN": "rgba(128, 255, 0, 0.80)",
    "GLN": "rgba(128, 255, 0, 0.80)",
    "SER": "rgba(153, 255, 51, 0.80)",
    "THR": "rgba(153, 255, 51, 0.80)",
    "TYR": "rgba(204, 155, 53, 0.80)",
    "TRP": "rgba(204, 155, 53, 0.80)",
    "HIS": "rgba(204, 155, 53, 0.80)",
    "PHE": "rgba(255, 255, 51, 0.80)",
    "ALA": "rgba(255, 255, 0, 0.80)",
    "ILE": "rgba(255, 255, 0, 0.80)",
    "VAL": "rgba(255, 255, 0, 0.80)",
    "PRO": "rgba(255, 255, 0 0.80)",
    "LEU": "rgba(255, 255, 0, 0.80)",
    "GLY": "rgba(255, 255, 255, 0.80)",
    "GLU": "rgba(255, 0, 0, 0.80)",
    "ASP": "rgba(255, 0, 0, 0.80)",
    "LYS": "rgba(0, 0, 255, 0.80)",
    "ARG": "rgba(0, 0, 255, 0.80)",
    }

# Colors for DNA / RNA
# Note: Based on WebLogo color scheme
DNARNA_COLORS = {
    "A": "rgba(51, 204, 51, 0.80)",
    "T": "rgba(204, 0, 0, 0.80)",
    "U": "rgba(204, 0, 0, 0.80)",
    "C": "rgba(51, 102, 255, 0.80)",
    "G": "rgba(255, 163, 26, 0.80)",
    }
FULL_DNARNA_COLORS = {f"D{k}": rgba for k, rgba in DNARNA_COLORS.items()}

# Combine all of them
AA_DNA_RNA_COLORS = {}
AA_DNA_RNA_COLORS.update(RESIDUES_COLORS)
AA_DNA_RNA_COLORS.update(DNARNA_COLORS)
AA_DNA_RNA_COLORS.update(FULL_DNARNA_COLORS)

# Chain colors
CHAIN_COLORS = [
    'rgba(51, 255, 51, 0.85)',
    'rgba(51, 153, 255, 0.85)',
    'rgba(255, 153, 51, 0.85)',
    'rgba(255, 255, 51, 0.85)',
    'rgba(255, 0, 0, 0.85)',
    'rgba(255, 0, 127, 0.85)',
    'rgba(0, 255, 0, 0.85)',
    'rgba(0, 0, 255, 0.85)',
    'rgba(0, 153, 0, 0.85)',
    ]
CHAIN_COLORS = CHAIN_COLORS[::-1]


##################
# Define classes #
##################
class ContactsMap(SupportsRun):
    """ContactMap analysis for single structure."""

    def __init__(
            self,
            model: Path,
            output: Path,
            params: dict,
            ) -> None:
        self.model = model
        self.output = output
        self.params = params
        self.files: dict[str, Union[str, Path]] = {}

    def run(self):
        """Process analysis of contacts of a PDB structure."""
        # Load pdb
        pdb_dt = extract_pdb_dt(self.model)
        # Extract all cordinates
        all_coords, resid_keys, resid_dt = get_ordered_coords(pdb_dt)
        # Compute distance matrix
        full_dist_matrix = compute_distance_matrix(all_coords)

        res_res_contacts = []
        all_heavy_interchain_contacts = []
        # First loop over residues
        for ri, reskey_1 in enumerate(resid_keys):
            # Second loop over residues (half matrix only)
            for _rj, reskey_2 in enumerate(resid_keys[ri + 1:], start=ri + 1):
                # Extract data
                contact_dt = gen_contact_dt(
                    full_dist_matrix,
                    resid_dt,
                    reskey_1,
                    reskey_2,
                    )
                res_res_contacts.append(contact_dt)

                # Extract interchain heavy atoms data
                if reskey_2.split('-')[0] == reskey_1.split('-')[0]:
                    continue
                heavy_atoms_contacts = extract_heavyatom_contacts(
                    full_dist_matrix,
                    resid_dt,
                    reskey_1,
                    reskey_2,
                    contact_distance=self.params['shortest_dist_threshold'],
                    )
                all_heavy_interchain_contacts += heavy_atoms_contacts

        # generate outputs for single models
        if self.params['single_model_analysis']:
            self.generate_output(
                res_res_contacts, all_heavy_interchain_contacts,
                )

        return res_res_contacts, all_heavy_interchain_contacts
    
    def generate_output(
            self,
            res_res_contacts: list[dict],
            all_heavy_interchain_contacts: list[dict],
            ) -> None:
        """Generate several outputs based on contacts.

        Parameters
        ----------
        res_res_contacts : list[dict]
            List of residue-residue contacts
        all_heavy_interchain_contacts : list[dict]
            List of heavy atoms interchain contacts
        """
        # write contacts tsv files
        header = ['res1', 'res2']
        header += [
            v for v in sorted(res_res_contacts[0])
            if v not in header
            ]
        fpath = write_res_contacts(
            res_res_contacts,
            header,
            f'{self.output}_contacts.tsv',
            interchain_data={
                'path': f'{self.output}_interchain_contacts.tsv',
                'data_key': 'ca-ca-dist',
                'contact_threshold': self.params['ca_ca_dist_threshold'],
                }
            )
        log.info(f'Generated contacts file: {fpath}')
        self.files['res-res-contacts'] = fpath

        # Genreate corresponding heatmap
        if self.params['generate_heatmap']:
            heatmap = tsv_to_heatmap(
                fpath,
                data_key='ca-ca-dist',
                contact_threshold=self.params['ca_ca_dist_threshold'],
                colorscale=self.params['color_ramp'],
                output_fname=f'{self.output}_heatmap.html',
                offline=self.params["offline"],
                )
            log.info(f'Generated single model heatmap file: {heatmap}')
            self.files['res-res-contactmap'] = heatmap

        # Generate corresponding chord chart
        if self.params['generate_chordchart']:
            # find theshold type
            if self.params['chordchart_datatype'] == 'ca-ca-dist':
                threshold = self.params['ca_ca_dist_threshold']
            else:
                threshold = self.params['shortest_dist_threshold']
            chordp = tsv_to_chordchart(
                fpath,
                data_key=self.params['chordchart_datatype'],
                contact_threshold=threshold,
                output_fname=f'{self.output}_chordchart.html',
                filter_intermolecular_contacts=True,
                title=Path(self.output).stem.replace('_', ' '),
                offline=self.params["offline"],
                )
            log.info(f'Generated single model chordchart file: {chordp}')
            self.files['res-res-chordchart'] = chordp

        # Write interchain heavy atoms contacts tsv file
        header2 = ['atom1', 'atom2', 'dist']
        fpath2 = write_res_contacts(
            all_heavy_interchain_contacts,
            header2,
            f'{self.output}_heavyatoms_interchain_contacts.tsv',
            )
        log.info(f'Generated contacts file: {fpath2}')
        self.files['atom-atom-interchain-contacts'] = fpath2


class ClusteredContactMap(SupportsRun):
    """ContactMap analysis for set of clustered structures."""

    def __init__(
            self,
            models: list[Path],
            output: Path,
            params: dict,
            ) -> None:
        self.models = models
        self.output = output
        self.params = params
        self.files: dict[str, Union[str, Path]] = {}
        self.terminated = False

    @staticmethod
    def aggregate_contacts(
            contacts_holder: dict,
            contact_keys: list[str],
            contacts: list[dict],
            key1: str,
            key2: str,
            ) -> None:
        """Aggregate single models data belonging to a cluster.

        Parameters
        ----------
        contacts_holder : dict
            Dictionnary holding list of contact data
        contact_keys : list[str]
            Order of the keys to access the dictionnary
        contacts : list[dict]
            Singel model contact data.
        key1 : str
            Name of the key to access first entry in data.
        key2 : str
            Name of the key to access second entry in data.
        """
        # Parse outputs to aggregate contacts in `clusters_contacts`
        for cont in contacts:
            # Check key
            combined_key = f'{cont[key2]}/{cont[key1]}'  # resversed
            if combined_key not in contacts_holder.keys():
                combined_key = f'{cont[key1]}/{cont[key2]}'  # normal
                if combined_key not in contacts_holder.keys():
                    # Add key order
                    contact_keys.append(combined_key)
                    # Initiate key
                    contacts_holder[combined_key] = {
                        k: []
                        for k in cont.keys()
                        if k not in [key1, key2]
                        }
            # Add data
            for dtk in contacts_holder[combined_key].keys():
                contacts_holder[combined_key][dtk].append(cont[dtk])

    def run(self):
        """Process analysis of contacts of a set of PDB structures."""
        # initiate holding variables
        clusters_contacts = {}  # Residue-residue contacts
        resres_keys_list = []  # Ordered residue-residue contacts keys
        clusters_heavyatm_contacts = {}  # Interchain atom-atom contacts
        atat_keys_list = []  # Ordered interchain atom-atom contacts keys

        # loop over models/structures
        for pdb_path in self.models:
            # initiate object
            contact_map_obj = ContactsMap(
                pdb_path,
                f'{self.output}_{pdb_path.stem}',
                self.params,
                )
            # Run it
            pdb_contacts, interchain_heavy_contacts = contact_map_obj.run()

            # Parse outputs to aggregate contacts in `clusters_contacts`
            self.aggregate_contacts(
                clusters_contacts, resres_keys_list,
                pdb_contacts,
                "res1", "res2",
                )

            # Parse outputs for heavy atoms contacts
            self.aggregate_contacts(
                clusters_heavyatm_contacts, atat_keys_list,
                interchain_heavy_contacts,
                "atom1", "atom2",
                )

        # Initiate heavy atoms contact cluster aggrated data
        heavy_atm_clust_list = []
        for atatk in atat_keys_list:
            at1, at2 = atatk.split('/')
            # point corresponding list of distances
            h_dists = clusters_heavyatm_contacts[atatk]['dist']
            # Summerize it
            heavy_atm_clust_list.append({
                "atom1": at1,
                "atom2": at2,
                "nb_dists": len(h_dists),
                "avg_dist": round(np.mean(h_dists), 2),
                "std_dist": round(np.std(h_dists), 2),
                })
        # write contacts
        header = ['atom1', 'atom2']
        header += [
            v for v in sorted(heavy_atm_clust_list[0])
            if v not in header
            ]
        hfpath = write_res_contacts(
            heavy_atm_clust_list,
            header,
            f'{self.output}_heavyatoms_interchain_contacts.tsv',
            )
        log.info(f'Generated heavy atoms interchain contacts file: {hfpath}')
        
        # Initiate cluster aggregated data holder
        combined_clusters_list = []
        # Loop over ordered keys
        for combined_key in resres_keys_list:
            # point data
            dt = clusters_contacts[combined_key]

            # Compute averages Ca-Ca distances
            avg_ca_ca_dist = np.mean(dt['ca-ca-dist'])
            # Compute nb. times cluster members holds a value under threshold
            ca_ca_under_thresh = [
                v for v in dt['ca-ca-dist']
                if v <= self.params['ca_ca_dist_threshold']
                ]
            nb_under = len(ca_ca_under_thresh)
            # Compute probability
            ca_ca_cont_probability = nb_under / len(dt['ca-ca-dist'])

            # Compute averages for shortest distances
            avg_shortest = np.mean(dt['shortest-dist'])
            # Generate list of shortest distances observed between two residues
            short_under_threshold = [
                v for v in dt['shortest-dist']
                if v <= self.params['shortest_dist_threshold']
                ]
            short_nb_und = len(short_under_threshold)
            # Compute nb. time the cluster members holds a value under threshold
            short_cont_proba = short_nb_und / len(dt['shortest-dist'])

            # Find most representative contact type
            cont_ts = list(set(dt['contact-type']))
            # Decreasing sorting of cluster contact types and pick highest one
            cont_t = sorted(
                cont_ts,
                key=lambda k: cont_ts.count(k),
                reverse=True,
                )[0]

            # Split key to recover resiudes names
            res1, res2 = combined_key.split('/')

            # Hold summary data for cluster
            combined_clusters_list.append({
                'res1': res1,
                'res2': res2,
                'ca-ca-dist': round(avg_ca_ca_dist, 1),
                'ca-ca-cont-probability': round(ca_ca_cont_probability, 2),
                'shortest-dist': round(avg_shortest, 1),
                'shortest-cont-probability': round(short_cont_proba, 2),
                'contact-type': cont_t,
                })
        
        # write contacts
        header = ['res1', 'res2']
        header += [
            v for v in sorted(combined_clusters_list[0])
            if v not in header
            ]
        fpath = write_res_contacts(
            combined_clusters_list,
            header,
            f'{self.output}_contacts.tsv',
            interchain_data={
                'path': f'{self.output}_interchain_contacts.tsv',
                'data_key': 'ca-ca-dist',
                'contact_threshold': self.params['ca_ca_dist_threshold'],
                }
            )
        log.info(f'Generated contacts file: {fpath}')
        self.files['res-res-contacts'] = fpath
        self.files['atom-atom-interchain-contacts'] = f'{self.output}_interchain_contacts.tsv'  # noqa : E501

        # Generate corresponding heatmap
        if self.params['generate_heatmap']:
            heatmap_path = tsv_to_heatmap(
                fpath,
                data_key=self.params['cluster_heatmap_datatype'],
                contact_threshold=1,
                colorscale=self.params['color_ramp'],
                output_fname=f'{self.output}_heatmap.html',
                offline=self.params["offline"],
                )
            log.info(f'Generated cluster contacts heatmap: {heatmap_path}')
            self.files['res-res-contactmap'] = heatmap_path

        # Generate corresponding chord chart
        if self.params['generate_chordchart']:
            # find theshold type
            if self.params['chordchart_datatype'] == 'ca-ca-dist':
                threshold = self.params['ca_ca_dist_threshold']
            else:
                threshold = self.params['shortest_dist_threshold']
            chordp = tsv_to_chordchart(
                fpath,
                data_key=self.params['chordchart_datatype'],
                contact_threshold=threshold,
                output_fname=f'{self.output}_chordchart.html',
                filter_intermolecular_contacts=True,
                title=Path(self.output).stem.replace('_', ' '),
                offline=self.params["offline"],
                )
            log.info(f'Generated cluster contacts chordchart file: {chordp}')
            self.files['res-res-chordchart'] = chordp

        self.terminated = True


def make_contactmap_report(
        contactmap_jobs: list[Union[ContactsMap, ClusteredContactMap]],
        outputpath: Union[str, Path],
        ) -> Union[str, Path]:
    """Generate a HTML navigation page holding all generated files.

    Parameters
    ----------
    contact_jobs : list[Union[ClusteredContactMap, ContactsMap]]
        All the terminated jobs
    outputpath : Union[str, Path]
        Output filepath where to write the report.

    Returns
    -------
    outputpath: Union[str, Path]
        Path to the generated report.
    """
    ordered_files = []
    # Loop over terminated jobs
    for job in contactmap_jobs:
        basepath = f"{job.output}_"
        # Gather all files generated by this job
        job_files = glob.glob(f"{basepath}*")

        # Sort them by file extension
        ext_names: dict[str, Union[str, Path]] = {}
        for fpath in job_files:
            _fname, ext = os.path.splitext(fpath)
            if ext not in ext_names.keys():
                ext_names[ext] = []
            ext_names[ext].append(fpath)
        # Sort each keys filepaths
        for ext in ext_names.keys():
            ext_names[ext] = sorted(
                ext_names[ext],
                key=lambda k: k.replace(basepath, ""),
                )
        # Get final list order
        sorted_jobfiles = [
            fpath
            for ext in sorted(ext_names)
            for fpath in ext_names[ext]
            ]

        # Initiate html links holding list
        job_list: list[str] = []
        # Loop over generated files
        for fpath in sorted_jobfiles:
            # Generate html link
            shortname = fpath.replace(basepath, "")
            html_string = f'<a href="{fpath}" target="_blank">{shortname}</a>'
            job_list.append(html_string)
        # Combine all links in one string
        job_list_combined = ', '.join(job_list)
        # Create final string
        job_access = f"<b>{job.output}:</b> {job_list_combined}"
        # Hold that guy
        ordered_files.append(job_access)

    # Combine all jobs outputs as a list
    all_access = '</li>\n            <li>'.join(ordered_files)
    # Generate small html file
    htmldt = f"""
    <div id="contactmap_report">
        <ul>
            <li>
            {all_access}
            </li>
        </ul>
    </div>
"""

    # Write it
    with open(outputpath, 'w') as reportout:
        reportout.write(htmldt)
    log.info(f'Generated report file: {outputpath}')
    # Return generate outputfilepath
    return outputpath


def get_clusters_sets(
        models: list[PDBFile],
        ) -> dict[tuple[Optional[int], Optional[int]], list[PDBFile]]:
    """Split models by clusters ids.

    Parameters
    ----------
    models : list
        List of pdb models/complexes.

    Return
    ------
    clust_sets : dict[tuple[Optional[int], Optional[int]], list[PDBFile]]
        Dictionary of (cluster ids, cluster rank) keys containing their
        respective models as list of PDBFiles.
    """
    clust_sets: dict[tuple[Optional[int], Optional[int]], list[PDBFile]] = {}
    # Loop over models
    for model in models:
        # Set key
        cluster_key = (model.clt_id, model.clt_rank)
        # Create/Gather the corresponding list
        cluster_list = clust_sets.setdefault(cluster_key, [])
        # Add model
        cluster_list.append(model)
    return clust_sets


def topX_models(models: list[PDBFile], topX: int = 10) -> list[Any]:
    """Sort and return subset of top X best models.

    Parameters
    ----------
    models : list
        List of pdb models/complexes.
    topX : int
        Number of models to return after sorting.

    Return
    ------
    subset_bests : list
        List of top `X` best models.
    """
    try:
        sorted_models = sorted(models, key=lambda m: m.score)
    except AttributeError:
        sorted_models = models
    finally:
        subset_bests = sorted_models[:topX]
    return subset_bests


####################
# Define functions #
####################
def extract_pdb_dt(path: Path) -> dict:
    """Read and extract ATOM/HETATM records from a pdb file.

    Parameters
    ----------
    path : Path
        Path to a pdb file.

    Return
    ------
    pdb_chains : dict
        A dictionary of the pdb file accesible using chains as keys.
    """
    pdb_chains: dict = {'chain_order': []}
    # Read file
    with open(path, 'r') as f:
        # Loop over lines
        for _ in f:
            # Skip non ATOM / HETATM lines
            if not any([
               _.startswith('ATOM'),
               _.startswith('HETATM'),
               ]):
                continue

            # Extract residue name
            resname = _[slc_resname]
            # Extract chain id
            chainid = _[slc_chainid]
            # Extract resid
            resid = _[slc_resseq].strip()

            # Check if chain already parsed
            if chainid not in pdb_chains.keys():
                # Add to ordered chains
                pdb_chains['chain_order'].append(chainid)
                # Initiate new chain holder
                pdb_chains[chainid] = {'order': []}

            # Check if new resid id
            if resid not in pdb_chains[chainid].keys():
                # Add to oredered resids
                pdb_chains[chainid]['order'].append(resid)
                # Initiate new residue holder
                pdb_chains[chainid][resid] = {
                    'index': len(pdb_chains[chainid]['order']) - 1,
                    'resname': resname,
                    'chainid': chainid,
                    'resid': resid,
                    'position': len(pdb_chains[chainid]['order']),
                    'atoms_order': [],
                    'atoms': {},
                    }
            # extract atome name
            atname = _[slc_name].strip()
            # check if not an hydrogen
            if atname[0] == 'H':
                continue

            # extact atome coordinates
            coords = extract_pdb_coords(_)
            pdb_chains[chainid][resid]['atoms_order'].append(atname)
            pdb_chains[chainid][resid]['atoms'][atname] = coords

    return pdb_chains


def extract_pdb_coords(line: str) -> list[float]:
    """Extract coordinated from a PDB line.

    Parameters
    ----------
    line : str
        A strandard ATOM/HETATM pdb record.

    Return
    ------
    coords : list[float]
        List of the X, Y and Z coordinate of this atom.
    """
    x = float(line[slc_x].strip())
    y = float(line[slc_y].strip())
    z = float(line[slc_z].strip())
    coords = [x, y, z]
    return coords


def get_ordered_coords(
        pdb_chains: dict,
        ) -> tuple[list[list[float]], list[str], dict]:
    """Generate list of all atom coordinates.

    Parameters
    ----------
    pdb_chains : dict
        A dictionary of the pdb file accesible using chains as keys,
         as provided by the `extract_pdb_dt()` function.

    Return
    ------
    all_coords : list[list[float]]
        All atomic coordinates in a single list.
    resid_keys : list[str]
        Ordered list of residues keys.
    resid_dt : dict
        Dictionary of coordinates indices for each residue.
    """
    # Define holders
    all_coords = []
    resid_keys = []
    resid_dt = {}
    i = 0
    # Loop over chains
    for chainid in pdb_chains['chain_order']:
        # Loop over residues of this chain
        for resid in pdb_chains[chainid]['order']:
            # create a resdiue key
            resname = pdb_chains[chainid][resid]['resname']
            reskey = f'{chainid}-{resid}-{resname}'
            resdt = {
                'atoms_indices': [],
                'resname': resname,
                'atoms_order': pdb_chains[chainid][resid]['atoms_order'],
                }
            # Loop over atoms of this residue
            for atname in pdb_chains[chainid][resid]['atoms_order']:
                # list of internal indices
                resdt['atoms_indices'].append(i)
                # index of a CA
                if atname == 'CA':
                    resdt['CA'] = i
                # Point atome submatrix coordinates index
                all_coords.append(pdb_chains[chainid][resid]['atoms'][atname])
                # increment atom index
                i += 1
            resid_dt[reskey] = resdt
            resid_keys.append(reskey)
    return (all_coords, resid_keys, resid_dt)


def compute_distance_matrix(all_atm_coords: list[list[float]]) -> NDFloat:
    """Compute all vs all distance matrix.

    Paramaters
    ----------
    all_atm_coords : list[list[float]]
        List of atomic coordinates.

    Return
    ------
    dist_matrix : NDFloat
        N*N distance matrix between all coordinates.
    """
    dist_matrix = squareform(pdist(all_atm_coords))
    return dist_matrix


def extract_submatrix(
        matrix: NDFloat,
        indices: list[int],
        indices2: Optional[list[int]] = None,
        ) -> NDFloat:
    """Extract submatrix based on desired indices.

    Paramaters
    ----------
    matrix : NDFloat
        A N*N matrix.
    indices : list[int]
        List of `row` indices to extract from this matrix
    indices2 : list[int]
        List of `columns` indices to extract from this matrix.
         if unspecified, indices2 == indices and symetric matrix
         is extracted.

    Return
    ------
    submat : NDFloat
        The extracted submatrix.
    """
    # Set second set of indices (columns) to first if not defined
    if not indices2:
        indices2 = indices
    # extract submatrix
    submat = matrix[np.ix_(indices, indices2)]
    return submat


def gen_contact_dt(
        matrix: NDFloat,
        resdt: dict,
        res1_key: str,
        res2_key: str,
        ) -> dict:
    """Generate contacts data.

    Parameters
    ----------
    matrix : NDFloat
        The distance matrix.
    resdt : dict
        Residues data with atom indices as returned by `get_ordered_coords()`.
    res1_key : str
        First residue of interest.
    res2_key : str
        Second residue of interest

    Return
    ------
    cont_dt : dict
        Dictionary holding contact data
    """
    # point residues data
    res1_dt = resdt[res1_key]
    res2_dt = resdt[res2_key]
    # point ca-ca dist
    try:
        ca_ca_dist = matrix[res1_dt['CA'], res2_dt['CA']]
    except KeyError:
        ca_ca_dist = 9999
    # obtain submatrix
    res1_res2_atm_submat = extract_submatrix(
        matrix,
        res1_dt['atoms_indices'],
        res2_dt['atoms_indices'],
        )
    # obtain clostest contact
    clostest_contact = min_dist(res1_res2_atm_submat)
    # contact type
    cont_type = get_cont_type(res1_dt['resname'], res2_dt['resname'])
    # set return variable
    cont_dt = {
        'res1': res1_key,
        'res2': res2_key,
        'ca-ca-dist': round(ca_ca_dist, 1),
        'shortest-dist': round(clostest_contact, 1),
        'contact-type': cont_type,
        }
    return cont_dt


def extract_heavyatom_contacts(
        matrix: NDFloat,
        resdt: dict,
        res1_key: str,
        res2_key: str,
        contact_distance: float = 4.5,
        ) -> list[dict[str, Union[float, str]]]:
    """Generate contacts data.

    Parameters
    ----------
    matrix : NDFloat
        The distance matrix.
    resdt : dict
        Residues data with atom indices as returned by `get_ordered_coords()`.
    res1_key : str
        First residue of interest.
    res2_key : str
        Second residue of interest.
    contact_distance : float
        Distance defining a contact.

    Return
    ------
    all_contacts : list[dict[str, Union[float, str]]]
        List holding contact data
    """
    all_contacts: list[dict[str, Union[float, str]]] = []
    # point data for first residue
    res1_indices = resdt[res1_key]['atoms_indices']
    res1_atnames = resdt[res1_key]['atoms_order']
    # point data for second residue
    res2_indices = resdt[res2_key]['atoms_indices']
    res2_atnames = resdt[res2_key]['atoms_order']
    # Loop over res1 atoms / indices
    for r1_atname, r1_atindex in zip(res1_atnames, res1_indices):
        # Loop over res2 atoms / indices
        for r2_atname, r2_atindex in zip(res2_atnames, res2_indices):
            # Point corresponding distance in matrix
            r1_r2_dist = matrix[r1_atindex, r2_atindex]
            # Check if distance <= threshold
            if r1_r2_dist <= contact_distance:
                # Hold data
                contactdt = {
                    'atom1': f'{res1_key}-{r1_atname}',
                    'atom2': f'{res2_key}-{r2_atname}',
                    'dist': r1_r2_dist,
                    }
                all_contacts.append(contactdt)
    return all_contacts


def get_cont_type(resn1: str, resn2: str) -> str:
    """Generate polarity key between two residues.

    Parameters
    ----------
    resn1 : str
       3 letters code of fist residue.
    resn2 : str
       3 letters code of second residue.

    Return
    ------
    pol_key : str
        Combined residues polarities
    """
    pol_keys: list[str] = []
    for resn in [resn1, resn2]:
        if resn.strip() in RESIDUE_POLARITY.keys():
            pol_keys.append(RESIDUE_POLARITY[resn.strip()])
        else:
            pol_keys.append('unknow')
    pol_key = '-'.join(pol_keys)
    return pol_key


def min_dist(matrix: NDFloat) -> float:
    """Find minimum value in a matrix."""
    return np.min(matrix)


def write_res_contacts(
        res_res_contacts: list[dict],
        header: list[str],
        path: Union[Path, str],
        sep: str = '\t',
        interchain_data: Union[bool, dict] = None,
        ) -> Path:
    """Write a tsv file based on residues-residues contacts data.

    Parameters
    ----------
    res_res_contacts : list[dict]
        List of dict holding data for each residue-residue contacts.
    header : list[str]
        Ordered list of keys to access in the dicts.
    path : Path
        Path to the output file to generate.
    sep : str
        Character used to separate data within a line.

    Return
    ------
    path : Path
        Path to the generated file.
    """
    # define README data type content
    dttype_info = {
        'res1': 'Chain-Resname-ResID key identifying first residue',
        'res2': 'Chain-Resname-ResID key identifying second residue',
        'ca-ca-dist': 'Observed distances between the two carbon alpha (Ca) atoms',
        'ca-ca-cont-probability': 'Fraction of times a contact is observed under the ca-ca-dist threshold over all analysed models of the same cluster',  # noqa : E501
        'shortest-dist': 'Observed shortest distance between the two residues',
        'shortest-cont-probability': 'Fraction of times a contact is observed under the shortest-dist threshold over all analysed models of the same cluster',  # noqa : E501
        'contact-type': 'ResidueType - ResidueType contact name',
        'atom1': 'Chain-Resname-ResID-Atome key identifying first atom',
        'atom2': 'Chain-Resname-ResID-Atome key identifying second atom',
        'dist': 'Observed distance between two atoms',
        'nb_dists': 'Total number of observed distances',
        'avg_dist': 'Cluster average distance',
        'std_dist': 'Cluster distance standard deviation',
        }

    # Check for inter chain contacts
    gen_interchain_tsv: bool = False
    if interchain_data and type(interchain_data) == dict:
        expected_keys = ('path', 'contact_threshold', 'data_key', )
        if all([k in interchain_data.keys() for k in expected_keys]):
            gen_interchain_tsv = True
            interchain_tsvdt: list[list[str]] = [header]
        else:
            raise KeyError

    # initiate file content
    tsvdt: list[list[str]] = [header]
    for res_res_cont in res_res_contacts:
        tsvdt.append([str(res_res_cont[h]) for h in header])
        if gen_interchain_tsv:
            chain1 = res_res_cont['res1'].split('-')[0]
            chain2 = res_res_cont['res2'].split('-')[0]
            if chain1 != chain2:
                dist = res_res_cont[interchain_data['data_key']]
                if dist < interchain_data['contact_threshold']:
                    interchain_tsvdt.append(tsvdt[-1])
    tsv_str = '\n'.join([sep.join(_) for _ in tsvdt])

    # generate commented lines to be placed on top of file
    readme = [
        '#' * 80,
        '# This file contains extracted contacts half-matrix information',
        '#' * 80,
        '',
        ]
    for head in header[::-1]:
        readme.insert(2, f'# {head}: {dttype_info[head]}')

    # Write file
    with open(path, 'w') as tsvout:
        tsvout.write('\n'.join(readme))
        tsvout.write(tsv_str)

    # Write inter chain file
    if gen_interchain_tsv:
        # Modify readme
        readme[1] = readme[1].replace(
            'contacts half-matrix',
            'interchain contacts',
            )
        # Write file
        with open(interchain_data['path'], 'w') as f:
            f.write('\n'.join(readme))
            # Write data string
            f.write('\n'.join([sep.join(_) for _ in interchain_tsvdt]))

    return path


def tsv_to_heatmap(
        tsv_path: Path,
        sep: str = '\t',
        data_key: str = 'ca-ca-dist',
        contact_threshold: float = 7.5,
        colorscale: str = 'Greys',
        output_fname: Union[Path, str] = 'contacts.html',
        offline: bool = False,
        ) -> Union[Path, str]:
    """Read a tsv file and generate a heatmap from it.

    Paramters
    ---------
    tsv_path : Path
        Path a the .tsv file containing contact data.
    sep : str
        Separator character used to split data in each line.
    data_key : str
        Data key used to draw the plot.
    contact_threshold : float
        Upper boundary of maximum value to be plotted.
         any value above it will be set to this value.
    output_fname : Path
        Path to the generated graph.

    Return
    ------
    output_filepath : Union[Path, str]
        Path to the generated file.
    """
    half_matrix: list[float] = []
    labels: list[str] = []
    header: Union[bool, list[str]] = None
    with open(tsv_path, 'r') as f:
        for line in f:
            # skip commented lines
            if line.startswith('#'):
                continue
            # split line
            s_ = line.strip().split(sep)
            # gather header
            if not header:
                header = s_
                continue
            # point labels
            label1 = s_[header.index('res1')]
            label2 = s_[header.index('res2')]
            # Add them to set of labels
            if label1 not in labels:
                labels.append(label1)
            if label2 not in labels:
                labels.append(label2)

            # point data
            value = float(s_[header.index(data_key)])
            # bound data to contact_threshold
            bounded_value = min(value, contact_threshold)
            # add it to matrix
            half_matrix.append(bounded_value)

    # Genereate full matrix
    matrix = squareform(half_matrix)

    # set data label
    color_scale = datakey_to_colorscale(data_key, color_scale=colorscale)
    if 'probability' in data_key:
        data_label = 'probability'
        np.fill_diagonal(matrix, 1)
    else:
        data_label = 'distance'
    
    # Compute chains length
    chains_length: dict[str, int] = {}
    ordered_chains: list[str] = []
    for label in labels:
        chainid = label.split('-')[0]
        if chainid not in chains_length.keys():
            chains_length[chainid] = 0
            ordered_chains.append(chainid)
        chains_length[chainid] += 1
    # Compute chains delineations positions
    del_posi = [0]
    for chainid in ordered_chains:
        del_posi.append(del_posi[-1] + chains_length[chainid])
    # Compute chains delineations lines
    chains_limits: list[dict[str, float]] = []
    for delpos in del_posi:
        # Vertical lines
        chains_limits.append({
            "x0": delpos - 0.5,
            "x1": delpos - 0.5,
            "y0": -0.5,
            "y1": len(labels) - 0.5,
            })
        # Horizontal lines
        chains_limits.append({
            "y0": delpos - 0.5,
            "y1": delpos - 0.5,
            "x0": -0.5,
            "x1": len(labels) - 0.5,
            })
    # Generate hover template
    hovertemplate = (
        ' %{y}   &#8621;   %{x} <br>'
        f' Contact {data_label}: %{{z}}'
        '<extra></extra>'
        )

    # Generate heatmap
    output_filepath = heatmap_plotly(
        matrix,
        labels={'color': data_label},
        xlabels=labels,
        ylabels=labels,
        color_scale=color_scale,
        output_fname=output_fname,
        offline=offline,
        delineation_traces=chains_limits,
        hovertemplate=hovertemplate,
        )

    return output_filepath


def datakey_to_colorscale(data_key: str, color_scale: str = 'Greys') -> str:
    """Convert color scale into reverse if data implies to do it.

    data_key : str
        A dictionary key pointing to data type.
    color_scale : str
        Name of a base plotpy color_scale.

    Return
    ------
    color_scale : str
        Possibly the reverse name of the color_scale.
    """
    return f'{color_scale}_r' if 'probability' not in data_key else color_scale


######################################
# Start of the chord chart functions #
######################################
def moduloAB(val: float, lb: float, ub: float) -> float:
    """Map a real number onto the unit circle.

     The unit circle is identified with the interval [lb, ub), ub-lb=2*PI.

    Parameters
    ----------
    val : float
        The value to be mapped into the unit circle.
    lb : float
        The lower boundary.
    ub : float
        The upper boundary

    Return
    ------
    moduloab : float
        The modulo of val between lb and ub
    """
    if lb >= ub:
        raise ValueError('Incorrect interval ends')
    y = (val - lb) % (ub - lb)
    moduloab = y + ub if y < 0 else y + lb
    return moduloab


def within_2PI(val: float) -> bool:
    """Check if float value is within unit circle value range.

    Parameters
    ----------
    val : float
        The value to be tested.
    """
    return 0 <= val < 2 * PI


def check_square_matrix(data_matrix: NDArray) -> int:
    """Check if the matrix is a square one.

    Parameters
    ----------
    data_matrix : NDArray (2DArray)
        The matrix to be checked.

    Return
    ------
    nb_rows : int
        Number of rows in this matrix.
    """
    matrixshape = data_matrix.shape
    nb_rows = matrixshape[0]
    if len(matrixshape) > 2:
        raise ValueError('Data array must have only two dimensions')
    if nb_rows != matrixshape[1]:
        raise ValueError('Data array must have (n,n) shape')
    return nb_rows


def get_ideogram_ends(
        ideogram_len: NDFloat,
        gap: float,
        ) -> list[tuple[float, float]]:
    """Generate ideogram ends.

    Paramaters
    ----------
    ideogram_len : NDArray
        Length of each ideograms.
    gap : float
        Gap to add in between each ideogram.

    Return
    ------
    ideo_ends : list[tuple[float]]
        List of start and end position for each ideograms.
    """
    ideo_ends: list[tuple[float, float]] = []
    start = 0.0
    for k in range(len(ideogram_len)):
        end = float(start + ideogram_len[k])
        ideo_ends.append((start, end))
        # Increment new start by gap for next origin
        start = end + gap
    return ideo_ends


def make_ideogram_arc(
        radius: float,
        _phi: tuple[float, float],
        nb_points: float = 50,
        ) -> NDFloat:
    """Generate ideogran arc.

    Parameters
    ----------
    radius : float
        The circle radius.
    phi : tuple[float, float]
        Tuple of ends angle coordinates of an arc.
    nb_points : float
        Parameter that controls the number of points to be evaluated on an arc

    Return
    ------
    arc_positions : NDArray
        Array of 2D coorinates defining an arc.
    """
    if not within_2PI(_phi[0]) or not within_2PI(_phi[1]):
        phi = [moduloAB(t, 0, 2 * PI) for t in _phi]
    else:
        phi = [t for t in _phi]
    length = (phi[1] - phi[0]) % (2 * PI)
    nr = 5 if length <= (PI / 4) else int((nb_points * length) / PI)
    if phi[0] < phi[1]:
        theta = np.linspace(phi[0], phi[1], nr)
    else:
        theta = np.linspace(
            moduloAB(phi[0], -PI, PI),
            moduloAB(phi[1], -PI, PI),
            nr,
            )
    arc_positions = radius * np.exp(1j * theta)
    return arc_positions


def make_ribbon_ends(
        matrix: NDArray,
        row_sum: list[int],
        ideo_ends: list[tuple[float, float]],
        L: int,
        ) -> list[list[tuple[float, float]]]:
    """Generate all connecting ribbons coordinates.

    Parameters
    ----------
    matrix : NDArray
        The data matrix.
    row_sum : list[int]
        Number of connexions in each row.
    ideo_ends : list[tuple[float, float]]
        List of start and end position for each ideograms.

    Returns
    -------
    ribbon_boundary : list[list[tuple[float, float]]]
        Matrix of per residue ribbons start and end positions.
    """
    ribbon_boundary: list[list[tuple[float, float]]] = []
    for k, ideo_end in enumerate(ideo_ends):
        # Point stating coordinates of this residue ideo
        start = float(ideo_end[0])
        # No ribbon to be formed
        if row_sum[k] == 0:
            # Add empty set of ribbons
            ribbon_boundary.append([(0., 0.) for i in range(len(ideo_ends))])
            continue
        # Initiate row ribbons
        row_ribbon_ends: list[tuple[float, float]] = []
        # Compute increment
        increment = (ideo_end[1] - start) / row_sum[k]
        # Loop over positions
        for j in range(1, L + 1):
            # Skip if no ribbon to add for this k, j pair
            if matrix[k][j - 1] == 0:
                row_ribbon_ends.append((0., 0.))
                continue
            # Define end
            end = float(start + increment)
            # Hold data
            row_ribbon_ends.append((start, end))
            # Set next start to current end
            start = end
        # Add full row
        ribbon_boundary.append(row_ribbon_ends)
    return ribbon_boundary


def control_pts(
        angle: list[float],
        radius: float,
        ) -> list[tuple[float, float]]:
    """Generate control points to draw a SVGpath.

    Parameters
    ----------
    angle : list[float]
        A list containing angular coordinates of the control points b0, b1, b2.
    radius : float
        The distance from b1 to the origin O(0,0)

    Returns
    -------
    control_points : list[tuple[float, float]]
        The set of control points.

    Raises
    ------
    ValueError
        Raised if the number of angular coordinates is not equal to 3.
    """
    # Check number of angular coordinates
    if len(angle) != 3:
        raise ValueError('angle must have len = 3')
    b_cplx = np.array([np.exp(1j * angle[k]) for k in range(3)])
    # Give it its size
    b_cplx[1] = radius * b_cplx[1]
    # Generate control points as a list for two values
    control_points = list(zip(b_cplx.real, b_cplx.imag))
    return control_points


def ctrl_rib_chords(
        side1: tuple[float, float],
        side2: tuple[float, float],
        radius: float,
        ) -> list[list[tuple[float, float]]]:
    """Generate poligons points aiming at drawing ribbons.

    Parameters
    ----------
    side1 : tuple[float, float]
        List of angular variables of the ribbon arc ends defining
         the ribbon starting (ending) arc
    side2 : tuple[float, float]
        List of angular variables of the ribbon arc ends defining
         the ribbon starting (ending) arc
    radius : float, optional
        Circle radius size

    Returns
    -------
    list[list[tuple[float, float]]]
        _description_
    """
    if len(side1) != 2 or len(side2) != 2:
        raise ValueError('the arc ends must be elements in a list of len 2')
    poligons = [
        control_pts(
            [side1[j], (side1[j] + side2[j]) / 2, side2[j]],
            radius,
            )
        for j in range(2)
        ]
    return poligons


def make_q_bezier(control_points: list[tuple[float, float]]) -> str:
    """Define the Plotly SVG path for a quadratic Bezier curve.

        defined by the list of its control points.

    Parameters
    ----------
    control_points : list[tuple[float, float]]
        List of control points

    Return
    ------
    svgpath : str
        An SVG path
    """
    if len(control_points) != 3:
        raise ValueError('control poligon must have 3 points')
    _a, _b, _c = control_points
    svgpath = 'M ' + str(_a[0]) + ',' + str(_a[1]) + ' Q ' +\
        str(_b[0]) + ', ' + str(_b[1]) + ' ' +\
        str(_c[0]) + ', ' + str(_c[1])
    return svgpath


def make_ribbon_arc(theta0: float, theta1: float) -> str:
    """Generate a SVGpath to draw a ribbon arc.

    Parameters
    ----------
    theta0 : float
        Starting angle value
    theta1 : float
        Ending angle value

    Returns
    -------
    string_arc : str
        A string representing the SVGpath of the ribbon arc.

    Raises
    ------
    ValueError
        If provided theta0 and theta1 angles are incorrect for a ribbon.
    ValueError
        If the angle coordinates for an arc side of a ribbon are not
         in the appropriate range [0, 2*pi]
    """
    if within_2PI(theta0) and within_2PI(theta1):
        if theta0 < theta1:
            theta0 = moduloAB(theta0, -PI, PI)
            theta1 = moduloAB(theta1, -PI, PI)
            if theta0 * theta1 > 0:
                raise ValueError('incorrect angle coordinates for ribbon')

        nr = int(40 * (theta0 - theta1) / PI)
        if nr <= 2:
            nr = 3
        theta = np.linspace(theta0, theta1, nr)
        pts = np.exp(1j * theta)  # points on arc in polar complex form

        string_arc: str = ''
        for k in range(len(theta)):
            string_arc += f'L {str(pts.real[k])}, {str(pts.imag[k])} '
        return string_arc
    else:
        raise ValueError('the angle coordinates for an arc side of a ribbon '
                         'must be in [0, 2*pi]')


def make_layout(
        title: str,
        plot_size: float,
        layout_shapes: list[dict],
        ) -> go.Layout:
    """Generate the chart layout.

    Parameters
    ----------
    title : str
        Title to be given to the chart.
    plot_size : float
        Size of the chart.
    layout_shapes : list[dict]
        Shapes to be drawn.

    Returns
    -------
    layout : go.Layout
        The plotly layout.
    """
    # Set axis parameters to hide axis line, grid, ticklabels and title
    axis = {
        'showline': False,
        'zeroline': False,
        'showgrid': False,
        'showticklabels': False,
        'title': '',
        }
    # Getenate the layout
    layout = go.Layout(
        title=title,
        xaxis=axis,
        yaxis=axis,
        showlegend=True,  # Important to show legend
        # legend={'font': {'size': 10}},  # Lower font size
        width=plot_size + 150,  # +150 to accomodate legend / keep circle round
        height=plot_size,
        margin={"t": 25, "b": 25, "l": 25, "r": 25},
        hovermode='closest',
        shapes=layout_shapes,
        )
    return layout


def make_ideo_shape(
        path: str,
        line_color: str,
        fill_color: str,
        ) -> dict:
    """Generate data to draw a ideogram shape.

    Parameters
    ----------
    path : str
        A SVGPath to be drawn.
    line_color : str
        Color of the shape boundary.
    fill_color : str
        Shape filling color fr the ribbon shape.

    Returns
    -------
    dict
        Data enabling to draw a ideogram shape in layout.
    """
    return {
        "line": {"color": line_color, "width": 0.45},
        "path": path,
        "type": 'path',
        "fillcolor": fill_color,
        "layer": 'below',
        }


def make_ribbon(
        side1: tuple[float, float],
        side2: tuple[float, float],
        line_color: str,
        fill_color: str,
        radius: float = 0.2,
        ) -> dict:
    """Generate data to draw a ribbon.

    Parameters
    ----------
    side1 : list[float]
        List of angular variables of first ribbon arc ends defining
         the ribbon starting (ending) arc.
    side2 : list[float]
        List of angular variables of the other ribbon arc ends defining
         the ribbon starting (ending) arc.
    line_color : str
        Color of the shape boundary.
    fill_color : str
        Shape filling color fr the ribbon shape.
    radius : float, optional
        Circle radius size, by default 0.2.

    Returns
    -------
    dict
        Data enabling to draw a ribbon in layout.
    """
    poligon = ctrl_rib_chords(side1, side2, radius)
    _b, _c = poligon
    # Generate the SVGpath
    path = make_q_bezier(_b)
    path += make_ribbon_arc(side2[0], side2[1])
    path += make_q_bezier(_c[::-1])
    path += make_ribbon_arc(side1[1], side1[0])

    return {
        'line': {'color': line_color, 'width': 0.5},
        'path': path,
        'type': 'path',
        'fillcolor': fill_color,
        'layer': 'below',
        }


def invPerm(perm: list[int]) -> list[int]:
    """Generate the inverse of a permutation.

    Parameters
    ----------
    perm : _type_
        A permutation.

    Returns
    -------
    inv : list[int]
        Inverse of a permutation.
    """
    # Fill with zeros
    inv = [0] * len(perm)
    for i, s in enumerate(perm):
        inv[s] = i
    return inv


def get_chains_ideograms_ends(
        chains: dict[str, list[str]],
        gap: float = 2 * PI * 0.005,
        ) -> tuple[list[tuple[float, float]], NDFloat]:
    """Build ideogram ends to represent protein chains.

    Parameters
    ----------
    chains : dict[str, list[str]]
        Dictionary mapping chains with their respective set of residues labels.
    gap : float, optional
        Gap between two ideograms, by default 2*PI*0.005

    Returns
    -------
    chain_ideo_ends : list[tuple[float, float]]
        Ideogram ends to represent protein chains.
    chain_ideogram_length : NDFloat

    """
    chain_row_sum = [len(chains[chain])
                     for chain in sorted(chains, reverse=True)]
    chain_ideogram_length = 2 * PI * np.asarray(chain_row_sum)
    chain_ideogram_length /= sum(chain_row_sum)
    chain_ideogram_length -= gap * np.ones(len(chain_row_sum))
    chain_ideo_ends = get_ideogram_ends(chain_ideogram_length, gap)
    return chain_ideo_ends, chain_ideogram_length


def get_all_ideograms_ends(
        chains: dict,
        gap: float = 2 * PI * 0.005,
        ) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    """Generate both chain and residues ideograms ends.

    Parameters
    ----------
    chains : dict
        Dictionary mapping to list of residues labels.
    gap : float, optional
        Gap distance used to separate two ideograms, by default 2*PI*0.005

    Returns
    -------
    tuple[ideo_ends, chain_ideo_ends]
        A tuple containing residues ideo ends and chains ideo ends.

    ideo_ends : list[tuple[float, float]]
        List of residues ideograms start and ending positions.
    chain_ideo_ends : list[tuple[float, float]]
        List of chain ideograms start and ending positions.
    """
    chain_ideo_ends, chain_ideogram_length = get_chains_ideograms_ends(
        chains,
        gap=gap,
        )

    ideo_ends: list[tuple[float, float]] = []
    left = 0
    for ind, chain in enumerate(sorted(chains, reverse=True)):
        chain_labels = chains[chain]
        for _label in chain_labels:
            right = left + (chain_ideogram_length[ind] / len(chain_labels))
            ideo_ends.append((left, right))
            left = right
        left = right + gap
    return ideo_ends, chain_ideo_ends


def split_labels_by_chains(labels: list[str]) -> dict[str, list[str]]:
    """Map each label to its chain.

    Parameters
    ----------
    labels : list[str]
        List of residues keys. e.g.: A-SER-123 (chain A, serine 123)

    Returns
    -------
    chains : dict[str, list[str]]
        Dictionary mapping chains with their respective set of residues labels.
    """
    chains: dict[str, list] = {}
    for lab in labels:
        chain, resname, resid = lab.split('-')
        if chain not in chains.keys():
            chains[chain] = []
        chains[chain].append(lab)
    return chains


def contacts_to_connect_matrix(
        matrix: NDFloat,
        labels: list[str],
        ) -> list[list[int]]:
    """.

    Parameters
    ----------
    matrix : NDFloat
        A square contact matrix.
    labels : list[str]
        List of labels corresponding row & columns entries.

    Returns
    -------
    connect_matrix : list[list[Union[str, int]]]
        The connectivity matrix without self contacts.
    """
    connect_matrix: list[list[int]] = []
    for ri, label_i in enumerate(labels):
        # Point label chain name
        chain1 = label_i.split('-')[0]
        new_contact_mat_row: list[int] = []
        # Loop over columns
        for ci, label_j in enumerate(labels):
            # Point label chain name
            chain2 = label_j.split('-')[0]
            if chain1 != chain2 and matrix[ri, ci] == 1:
                interaction = 1
            else:
                interaction = 0
            new_contact_mat_row.append(interaction)
        # Hold row
        connect_matrix.append(new_contact_mat_row)
    return connect_matrix


def to_nice_label(label: str) -> str:
    """Convert a label into a user friendly label.

    Parameters
    ----------
    label : str
        Label name as found in csv

    Returns
    -------
    nicelabel : str
        User friendly description of the label.
    """
    slabel = label.split('-')
    nicelabel = f"Chain {slabel[0]}, residue {slabel[2]} {slabel[1]}"
    return nicelabel


def to_color_weight(
        distance: float,
        max_dist: float,
        min_dist: float = 2.,
        min_weight: float = 0.2,
        max_weight: float = 0.90,
        ) -> float:
    """Compute color weight based on distance.

    Parameters
    ----------
    distance : float
        The distance to weight.
    max_dist : float
        The max distance observed in the dataset.
    min_dist : float, optional
        The minumu, distance observed in the dataset, by default 2.
    min_weight : float, optional
        Color wight for the maximum distance, by default 0.2
    max_weight : float, optional
        Color wight for the minimum distance, by default 0.90

    Returns
    -------
    weight : float
        The color weight. in range [min_weight, max_weight]
    """
    # Scale dist into minimum
    dist = max(distance, min_dist)
    # Compute probability
    probability_dist = (dist - min_dist) / (max_dist - min_dist)
    # Obtain corresponding weight
    weight = ((min_weight - max_weight) * probability_dist) + max_weight
    # Return rounded value of weight
    return round(weight, 2)


def to_rgba_color_string(
        connect_color: tuple[int, int, int],
        alpha: float,
        ) -> str:
    """Generate a rgba string from list of colors and alpha.

    Parameters
    ----------
    connect_color : list[int]
        A 3-values list of integers defining the red, green and blue colors.
    alpha : float
        color_weight

    Returns
    -------
    rgba_color : str
        The html like rgba colors. e.g.: 'rgba(123, 123, 123, 0.5)'
    """
    colors_str = ",".join([str(v) for v in connect_color])
    rgba_color = f'rgba({colors_str},{alpha})'
    return rgba_color


def to_full_matrix(
        half_matrix: list[Union[int, float, str]],
        diag_val: Union[int, float, str],
        ) -> NDArray:
    """Generate a full matrix from a half matrix.

    Parameters
    ----------
    half_matrix : list[Any]
        Values of the N*(N-1)/2 half matrix.
    diag_val : Any
        Value to be placed in diagonal of the full matrix.

    Returns
    -------
    matrix : NDArry
        The reconstituted full matrix.
    """
    # Genereate full matrix from N*(N-1)/2 vector
    matrix = squareform(half_matrix)
    # Update diagonal with data
    np.fill_diagonal(matrix, diag_val)
    return matrix


def make_chordchart(
        _contact_matrix: list[list[int]],
        _dist_matrix: list[list[float]],
        _interttype_matrix: list[list[str]],
        _labels: list[str],
        gap: float = 2 * PI * 0.005,
        output_fpath: Union[str, Path] = 'chordchart.html',
        title: str = 'Chord diagram',
        offline: bool = False,
        ) -> Union[str, Path]:
    """Generate a plotly chordchart graph.

    Parameters
    ----------
    _contact_matrix : list[list[int]]
        The contact matrix
    _dist_matrix : list[list[float]]
        The distance matrix
    _interttype_matrix : list[list[str]]
        The interaction type matrix
    _labels : list[str]
        Labels of each matrix rows (and columns as supposed to be symetric)
    gap : float, optional
        Gap between two ideograms, by default 2*PI*0.005
    output_fpath : Union[str, Path], optional
        Path to the output file, by default 'chordchart.html'
    title : str, optional
        Title to give to the diagram, by default 'Chord diagram'

    Returns
    -------
    output_fpath : Union[str, Path]
        Path to the genereated output file.
    """
    # Unpack matrices
    contact_matrix = contacts_to_connect_matrix(_contact_matrix, _labels)

    # Reverse data order so later graph displayed clockwise
    matrix = np.array([ri[::-1] for ri in contact_matrix[::-1]])
    dist_matrix = [ri[::-1] for ri in _dist_matrix[::-1]]
    interttype_matrix = [ri[::-1] for ri in _interttype_matrix[::-1]]
    labels = _labels[::-1]

    # Check matrix shape
    L = check_square_matrix(matrix)

    # Map labels into respective chains
    chains = split_labels_by_chains(labels)
    # Set matrix of indices
    idx_sort = [list(range(L)) for i in range(L)]

    # Compute residues and chain ideograms positions
    ideo_ends, chain_ideo_ends = get_all_ideograms_ends(chains, gap=gap)

    # Compute number of connexion per residues
    row_sum = [np.sum(matrix[k, :]) for k in range(L)]

    # Compute connexion ribbons positions
    ribbon_ends = make_ribbon_ends(matrix, row_sum, ideo_ends, L)

    # Initiate shape holder
    layout_shapes: list[dict] = []
    # Initiate ribbon info holder
    ribbon_info: list[go.Scatter] = []
    # Loop over entries
    for k, label1 in enumerate(labels):
        sigma = idx_sort[k]
        sigma_inv = invPerm(sigma)
        # Half matrix loop to avoid duplicates
        for j in range(k, L):
            # No data to draw
            if matrix[k][j] == 0 and matrix[j][k] == 0:
                continue

            # Obtain ribbon color for this interaction type
            try:
                connect_color = CONNECT_COLORS[interttype_matrix[k][j]]
            except KeyError:
                connect_color = (123, 123, 123)
            color_weight = to_color_weight(dist_matrix[k][j], 9.5)
            rgba_color = to_rgba_color_string(connect_color, color_weight)

            # Point ribbons data
            side1 = ribbon_ends[k][sigma_inv[j]]
            eta = idx_sort[j]
            eta_inv = invPerm(eta)
            side2 = ribbon_ends[j][eta_inv[k]]
            zi = 0.9 * np.exp(1j * (side1[0] + side1[1]) / 2)
            zf = 0.9 * np.exp(1j * (side2[0] + side2[1]) / 2)

            # reverse interaction type for second label
            s_intertype = interttype_matrix[k][j].split('-')
            rev_s_intertype = s_intertype[::-1]
            rev_interttype = '-'.join(rev_s_intertype)
            # Obtain nice labels
            nicelabel1 = to_nice_label(label1)
            nicelabel2 = to_nice_label(labels[j])
            # texti and textf are the strings that will be displayed when
            # hovering the mouse over the two ribbon ends
            texti = f'{nicelabel1} {rev_interttype} with {nicelabel2}'
            textf = f'{nicelabel2} {interttype_matrix[j][k]} with {nicelabel1}'  # noqa : E501
            # Generate interactive labels
            for zv, text in zip([zi, zf], [texti, textf]):
                # Generate ribbon info
                ribbon_info.append(
                    go.Scatter(
                        x=[zv.real],
                        y=[zv.imag],
                        mode='markers',
                        marker={"size": 0.5, "color": rgba_color},
                        text=text,
                        hoverinfo='text',
                        showlegend=False,
                        )
                    )
            # Note: must reverse these arc ends to avoid twisted ribbon
            side2_rev = (side2[1], side2[0])
            # Append the ribbon shape
            layout_shapes.append(
                make_ribbon(
                    side1,
                    side2_rev,
                    'rgba(175,175,175)',
                    rgba_color,
                    )
                )

    ideograms: list[go.Scatter] = []
    # Draw ideograms for residues
    for k, label in enumerate(labels):
        z = make_ideogram_arc(1.1, ideo_ends[k])
        zi = make_ideogram_arc(1.0, ideo_ends[k])

        # Point residue name
        resname = label.split('-')[2]

        # Point corresponding color
        try:
            rescolor = AA_DNA_RNA_COLORS[resname.strip()]
        except KeyError:
            rescolor = 'rgba(123, 123, 123, 0.7)'

        # Build textual info
        text_info = f'{to_nice_label(label)}<br>'
        if row_sum[k] == 0:
            text_info += 'No interaction'
        else:
            text_info += f'Total of {row_sum[k]:d} interaction'
            if row_sum[k] >= 2:
                text_info += 's'
        # Add info
        ideograms.append(
            go.Scatter(
                x=z.real,
                y=z.imag,
                mode='lines',
                line={
                    "color": rescolor,
                    "shape": 'spline',
                    "width": 0.25,
                    },
                text=text_info,
                hoverinfo='text',
                showlegend=False,
                )
            )

        # Build corresponding SVG path
        m = len(z)
        svgpath = 'M '
        for s in range(m):
            svgpath += f'{str(z.real[s])}, {str(z.imag[s])} L '

        Zi = np.array(zi.tolist()[::-1])
        for s in range(m):
            svgpath += f'{str(Zi.real[s])}, {str(Zi.imag[s])} L '
        svgpath += f'{str(z.real[0])}, {str(z.imag[0])}'
        # Hold it
        layout_shapes.append(
            make_ideo_shape(
                svgpath,
                'rgba(150,150,150)',
                rescolor,
                )
            )

    # Draw ideograms for chains
    for k, chainid in enumerate(sorted(chains, reverse=True)):
        z = make_ideogram_arc(1.2, chain_ideo_ends[k])
        zi = make_ideogram_arc(1.11, chain_ideo_ends[k])
        m = len(z)
        ideograms.append(
            go.Scatter(
                x=z.real,
                y=z.imag,
                mode='lines',
                line={
                    "color": CHAIN_COLORS[k],
                    "shape": 'spline',
                    "width": 0.25,
                    },
                text=f'Chain {chainid}',
                hoverinfo='text',
                showlegend=False,
                )
            )

        # Build corresponding SVG path
        svgpath = 'M '
        for s in range(m):
            svgpath += f'{str(z.real[s])}, {str(z.imag[s])} L '
        Zi = np.array(zi.tolist()[::-1])
        for s in range(m):
            svgpath += f'{str(Zi.real[s])}, {str(Zi.imag[s])} L '
        svgpath += f'{str(z.real[0])}, {str(z.imag[0])}'

        layout_shapes.append(
            make_ideo_shape(
                svgpath,
                'rgba(150,150,150)',
                CHAIN_COLORS[k],
                )
            )

    # Compute figure size
    fig_size = 100 * np.log(L * L)
    # Create plotly layout
    layout = make_layout(title, fig_size, layout_shapes)
    # combine all data info
    data = ideograms + ribbon_info
    # Generate the figure
    fig = go.Figure(data=data, layout=layout)
    # Fine tune figure
    fig.update_layout(
        plot_bgcolor='white',
        )
    # Add legend(s)
    add_chordchart_legends(fig)
    # Write it as html file
    fig_to_html(
        fig,
        output_fpath,
        figure_height=fig_size,
        figure_width=fig_size,
        offline=offline,
        )
    return output_fpath


def add_chordchart_legends(fig: go.Figure) -> None:
    """Add custom legend to chordchart.

    Parameters
    ----------
    fig : go.Figure
        A plotly figure.
    """
    # Add connection types legends
    for key_key, color in REVERSED_CONNECT_COLORS_KEYS.items():
        # Create dummy traces
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                legendgroup="connect_color",
                legendgrouptitle_text="Interaction types",
                showlegend=True,
                name=key_key.replace('-', '&#8621;'),
                mode="lines",
                marker={
                    "color": to_rgba_color_string(color, 0.75),
                    "size": 10,
                    "symbol": "line-ew-open",
                    },
                )
            )
    # Add unknown interaction type
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            legendgroup="connect_color",
            legendgrouptitle_text="Interaction types",
            showlegend=True,
            name="Unknown",
            mode="lines",
            marker={
                "color": to_rgba_color_string((111, 111, 111), 0.75),
                "size": 10,
                "symbol": "line-ew-open",
                },
            )
        )

    # Add aa types legend
    for aa, rgba_color in RESIDUES_COLORS.items():
        # Create dummy traces
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                legendgroup="aa_color",
                legendgrouptitle_text="Residues/Bases",
                showlegend=True,
                name=aa,
                mode="lines",
                marker={
                    "color": rgba_color,
                    "size": 10,
                    "symbol": "line-ew-open",
                    },
                )
            )
    # Add nucleobases legend
    for na, rgba_color in DNARNA_COLORS.items():
        # Create dummy traces
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                legendgroup="aa_color",
                legendgrouptitle_text="Residues/Bases",
                showlegend=True,
                name=na,
                mode="lines",
                marker={
                    "color": rgba_color,
                    "size": 10,
                    "symbol": "line-ew-open",
                    },
                )
            )
    # Add unknown type
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            legendgroup="aa_color",
            legendgrouptitle_text="Residues/Bases",
            showlegend=True,
            name="Unknown",
            mode="lines",
            marker={
                "color": to_rgba_color_string((111, 111, 111), 0.75),
                "size": 10,
                "symbol": "line-ew-open",
                },
            )
        )


def tsv_to_chordchart(
        tsv_path: Path,
        sep: str = '\t',
        data_key: str = 'ca-ca-dist',
        contact_threshold: float = 7.5,
        filter_intermolecular_contacts: bool = True,
        output_fname: Union[Path, str] = 'contacts_chordchart.html',
        title: str = 'Chord diagram',
        offline: bool = False,
        ) -> Union[Path, str]:
    """Read a tsv file and generate a chord diagram from it.

    Paramters
    ---------
    tsv_path : Path
        Path a the .tsv file containing contact data.
    sep : str
        Separator character used to split data in each line.
    data_key : str
        Data key used to draw the plot.
    contact_threshold : float
        Upper boundary of maximum value to be plotted.
         any value above it will be set to this value.
    output_fname : Union[Path, str]
        Path where to generate the graph.
    title : str
        Title to give to the Chord diagram

    Return
    ------
    chord_chart_fpath : Union[Path, str]
        Path to the generated graph
    """
    # Initiate holders
    half_contact_matrix: list[int] = []
    half_value_matrix: list[float] = []
    half_intertype_matrix: list[str] = []
    labels: list[str] = []
    header: Union[bool, list[str]] = None
    # Read tsv file
    with open(tsv_path, 'r') as f:
        for line in f:
            # skip commented lines
            if line.startswith('#'):
                continue
            # split line
            s_ = line.strip().split(sep)
            # gather header
            if not header:
                header = s_
                continue
            # point labels
            label1 = s_[header.index('res1')]
            label2 = s_[header.index('res2')]
            # Add them to set of labels
            if label1 not in labels:
                labels.append(label1)
            if label2 not in labels:
                labels.append(label2)

            # point data
            value = float(s_[header.index(data_key)])
            # check if in contact
            contact = 1 if value <= contact_threshold else 0
            # Point interaction type
            inter_type = s_[header.index('contact-type')]
            # add it to matrix
            half_contact_matrix.append(contact)
            half_value_matrix.append(value)
            half_intertype_matrix.append(inter_type)

    # Genereate full matrices
    contact_matrix = to_full_matrix(half_contact_matrix, 1)
    dist_matrix = to_full_matrix(half_value_matrix, 0.)
    intertype_matrix = to_full_matrix(half_intertype_matrix, 'self-self')

    # Check if must get only the intermolecular contacts submatrix
    if filter_intermolecular_contacts:
        # Filter positions where: intermolecular contacts + dist <= threshold
        intmol_cont_labels: list[str] = []
        for ri, label1 in enumerate(labels):
            label1_chain = label1.split('-')[0]
            for ci, label2 in enumerate(labels):
                label2_chain = label2.split('-')[0]
                # Skip same chains
                if label1_chain == label2_chain:
                    continue
                # Check contacts
                if contact_matrix[ri, ci] == 1:
                    intmol_cont_labels += [label1, label2]

        # Obtain sorted subset
        nodoubles_intmol_cont_labels = list(set(intmol_cont_labels))
        sorted_intmol_cont_labels = sorted(
            nodoubles_intmol_cont_labels,
            key=lambda k: labels.index(k),
            )

        # Get indices to be extracted
        sorted_indices = [labels.index(k) for k in sorted_intmol_cont_labels]

        # Get submatrix
        contact_submatrix = extract_submatrix(contact_matrix, sorted_indices)
        dist_submatrix = extract_submatrix(dist_matrix, sorted_indices)
        intert_submatrix = extract_submatrix(intertype_matrix, sorted_indices)
        sublabels = sorted_intmol_cont_labels

    else:
        contact_submatrix = contact_matrix
        dist_submatrix = dist_matrix
        intert_submatrix = intertype_matrix
        sublabels = labels

    # Generate chord chart
    chord_chart_fpath = make_chordchart(
        contact_submatrix,
        dist_submatrix,
        intert_submatrix,
        sublabels,
        output_fpath=output_fname,
        title=title,
        offline=offline,
        )
    
    return chord_chart_fpath
