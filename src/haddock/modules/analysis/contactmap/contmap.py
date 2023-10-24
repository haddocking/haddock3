"""Module computing contact maps of complexes, alone or grouped by cluster."""

from pathlib import Path

import numpy as np
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
from haddock.core.typing import Any, NDFloat, Optional, Union
from haddock.libs.libplots import heatmap_plotly


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
    "GLU": "neg_charged",
    "ASP": "neg_charged",
    "LYS": "pos_charged",
    "ARG": "pos_charged",
    }


##################
# Define classes #
##################
class ContactsMap():
    """ContactMap analysis for single structure."""

    def __init__(self, model: Path, output: Path, params: dict) -> None:
        self.model = model
        self.output = output
        self.params = params

    def run(self):
        """Process analysis of contacts of a PDB structure."""
        # Load pdb
        pdb_dt = extract_pdb_dt(self.model)
        # Extract all cordinates
        all_coords, resid_keys, resid_dt = get_ordered_coords(pdb_dt)
        # Compute distance matrix
        full_dist_matrix = compute_distance_matrix(all_coords)

        res_res_contacts = []
        # First loop over residues
        for ri, reskey_1 in enumerate(resid_keys):
            # Second loop over residues
            for _rj, reskey_2 in enumerate(resid_keys[ri + 1:], start=ri + 1):
                contact_dt = gen_contact_dt(
                    full_dist_matrix,
                    resid_dt,
                    reskey_1,
                    reskey_2,
                    )
                res_res_contacts.append(contact_dt)

        # generate outputs for single models
        if self.params['single_model_analysis']:

            # write contacts
            header = ['res1', 'res2']
            header += [
                v for v in sorted(res_res_contacts[0])
                if v not in header
                ]
            fpath = write_res_contacts(
                res_res_contacts,
                header,
                f'{self.output}_contacts.tsv',
                )
            log.info(f'Generated contacts file: {fpath}')

            # Genreate corresponding heatmap
            if self.params['generate_heatmap']:
                heatmap = tsv_to_heatmap(
                    fpath,
                    data_key='ca-ca-dist',
                    contact_threshold=self.params['ca_ca_dist_threshold'],
                    output_fname=f'{self.output}_heatmap.html',
                    )
                log.info(f'Generated single model heatmap file: {heatmap}')

        return res_res_contacts


class ClusteredContactMap():
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
        self.terminated = False

    def run(self):
        """Process analysis of contacts of a set of PDB structures."""
        # initiate holding variables
        clusters_contacts = {}
        keys_list = []
        # loop over models/structures
        for pdb_path in self.models:
            # initiate object
            contact_map_obj = ContactsMap(
                pdb_path,
                f'{self.output}_{pdb_path.stem}',
                self.params,
                )
            # Run it
            pdb_contacts = contact_map_obj.run()
            # Parse outputs
            for cont in pdb_contacts:
                # Check key
                combined_key = f'{cont["res2"]}/{cont["res1"]}'
                if combined_key not in clusters_contacts.keys():
                    combined_key = f'{cont["res1"]}/{cont["res2"]}'
                    if combined_key not in clusters_contacts.keys():
                        keys_list.append(combined_key)
                        clusters_contacts[combined_key] = {
                            k: []
                            for k in cont.keys()
                            if k not in ["res1", "res2"]
                            }
                # Add data
                for dtk in clusters_contacts[combined_key].keys():
                    clusters_contacts[combined_key][dtk].append(cont[dtk])

        combined_clusters_list = []
        # Loop over keys
        for combined_key in keys_list:
            # point data
            res1, res2 = combined_key.split('/')
            dt = clusters_contacts[combined_key]

            # Compute averages and ratio
            avg_ca_ca_dist = np.mean(dt['ca-ca-dist'])
            ca_ca_above_thresh = [
                v for v in dt['ca-ca-dist']
                if v <= self.params['ca_ca_dist_threshold']
                ]
            ca_ca_cont_ratio = len(ca_ca_above_thresh) / len(dt['ca-ca-dist'])

            avg_shortest = np.mean(dt['shortest-dist'])
            short_ab_t = [
                v for v in dt['shortest-dist']
                if v <= self.params['shortest_dist_threshold']
                ]
            shortest_cont_ratio = len(short_ab_t) / len(dt['shortest-dist'])

            # most representative contact type
            cont_ts = list(set(dt['contact-type']))
            cont_t = sorted(
                cont_ts,
                key=lambda k: cont_ts.count(k),
                reverse=True,
                )[0]
            cont_t_ratio = cont_ts.count(cont_t) / len(dt['contact-type'])

            # Hold summery data for cluster
            combined_clusters_list.append({
                'res1': res1,
                'res2': res2,
                'ca-ca-dist': round(avg_ca_ca_dist, 1),
                'ca-ca-cont-ratio': round(ca_ca_cont_ratio, 2),
                'shortest-dist': round(avg_shortest, 1),
                'shortest-cont-ratio': round(shortest_cont_ratio, 2),
                'contact-type': cont_t,
                'contact-type-ratio': cont_t_ratio,
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
            )
        log.info(f'Generated contacts file: {fpath}')

        # Generate corresponding heatmap
        if self.params['generate_heatmap']:
            heatmap_path = tsv_to_heatmap(
                fpath,
                data_key=self.params['cluster-heatmap-datatype'],
                contact_threshold=1,
                output_fname=f'{self.output}_heatmap.html',
                )
            log.info(f'Generated cluster contacts heatmap: {heatmap_path}')

        self.terminated = True


def get_clusters_sets(models: list[PDBFile]) -> dict:
    """Split models by clusters ids.

    Parameters
    ----------
    models : list
        List of pdb models/complexes.

    Return
    ------
    clusters_sets : dict
        Dictionary of models acccessible by their cluster ids as keys.
    """
    clusters_sets: dict = {}
    for model in models:
        if model.clt_id not in clusters_sets.keys():
            clusters_sets[model.clt_id] = []
        clusters_sets[model.clt_id].append(model)
    return clusters_sets


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
                }
            # Loop over atoms of this residue
            for atname in pdb_chains[chainid][resid]['atoms_order']:
                # list of internal indices
                resdt['atoms_indices'].append(i)
                # index of a CA
                if atname == 'CA':
                    resdt['CA'] = i
                # Point atome cooct_submatrixrdinate
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
        # add other types of data ??
        }
    return cont_dt


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
        if resn in RESIDUE_POLARITY.keys():
            pol_keys.append(RESIDUE_POLARITY[resn])
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
        path: Path,
        sep: str = '\t',
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
    # initiate header
    header = ['res1', 'res2']
    header += [v for v in sorted(res_res_contacts[0]) if v not in header]

    # initiate file content
    tsvdt = [header]
    for res_res_cont in res_res_contacts:
        tsvdt.append([str(res_res_cont[h]) for h in header])
    tsv_str = '\n'.join([sep.join(_) for _ in tsvdt])

    # Write file
    with open(path, 'w') as tsvout:
        tsvout.write(tsv_str)

    return path


def tsv_to_heatmap(
        tsv_path: Path,
        sep: str = '\t',
        data_key: str = 'ca-ca-dist',
        contact_threshold: float = 7.5,
        output_fname: Union[Path, str] = 'contacts.html',
        ) -> None:
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
    """
    half_matrix = []
    labels = []
    with open(tsv_path, 'r') as f:
        header = f.readline().strip().split(sep)
        for line in f:
            s_ = line.strip().split(sep)
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
    color_scale = datakey_to_colorscale(data_key, color_scale='Greys')
    if 'ratio' in data_key:
        data_label = 'ratio'
        np.fill_diagonal(matrix, 1)
    else:
        data_label = 'distance'

    # Generate heatmap
    heatmap_plotly(
        matrix,
        labels={'color': data_label},
        xlabels=labels,
        ylabels=labels,
        color_scale=color_scale,
        output_fname=output_fname,
        )


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
    return f'{color_scale}_r' if 'ratio' not in data_key else color_scale
