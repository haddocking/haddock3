"""Module computing contact maps of complexes, alone or grouped by cluster.

Chord diagram functions were adapted from:
https://plotly.com/python/v3/filled-chord-diagram/
"""

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
from haddock.core.typing import Any, NDFloat, NDArray, Optional, Union
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
    "GLU": "negative",
    "ASP": "negative",
    "LYS": "positive",
    "ARG": "positive",
    }

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
reversed_keys_dict = {'-'.join(k.split('-')[::-1]): v
                      for k, v in CONNECT_COLORS.items()}
CONNECT_COLORS.update(reversed_keys_dict)

# Colors for each residue
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

# Chain colors ( in +- pymol order )
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
                    colorscale=self.params['color_ramp'],
                    output_fname=f'{self.output}_heatmap.html',
                    )
                log.info(f'Generated single model heatmap file: {heatmap}')

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
                    title="",
                    )
                log.info(f'Generated single model chordchart file: {chordp}')

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
        # Loop over ordered keys
        for combined_key in keys_list:
            # point data
            res1, res2 = combined_key.split('/')
            dt = clusters_contacts[combined_key]

            # Compute averages Ca-Ca distances
            avg_ca_ca_dist = np.mean(dt['ca-ca-dist'])
            # Compute number of time the cluster members holds a value under threshold
            ca_ca_above_thresh = [
                v for v in dt['ca-ca-dist']
                if v <= self.params['ca_ca_dist_threshold']
                ]
            ca_ca_cont_ratio = len(ca_ca_above_thresh) / len(dt['ca-ca-dist'])

            # Compute averages for shortest distances
            avg_shortest = np.mean(dt['shortest-dist'])
            # Generate list of shortest distances observed between two residues
            short_ab_t = [
                v for v in dt['shortest-dist']
                if v <= self.params['shortest_dist_threshold']
                ]
            # Compute number of time the cluster members holds a value under threshold
            shortest_cont_ratio = len(short_ab_t) / len(dt['shortest-dist'])

            # Find most representative contact type
            cont_ts = list(set(dt['contact-type']))
            # Decreasing sorting of cluster contact types and pick highest one
            cont_t = sorted(
                cont_ts,
                key=lambda k: cont_ts.count(k),
                reverse=True,
                )[0]
            # Compute ratio for this contact type to be found
            cont_t_ratio = cont_ts.count(cont_t) / len(dt['contact-type'])

            # Hold summary data for cluster
            combined_clusters_list.append({
                'res1': res1,
                'res2': res2,
                'ca-ca-dist': round(avg_ca_ca_dist, 1),
                'ca-ca-cont-ratio': round(ca_ca_cont_ratio, 2),
                'shortest-dist': round(avg_shortest, 1),
                'shortest-cont-ratio': round(shortest_cont_ratio, 2),
                'contact-type': cont_t,
                'contact-type-ratio': round(cont_t_ratio, 2),
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
                data_key=self.params['cluster_heatmap_datatype'],
                contact_threshold=1,
                colorscale=self.params['color_ramp'],
                output_fname=f'{self.output}_heatmap.html',
                )
            log.info(f'Generated cluster contacts heatmap: {heatmap_path}')

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
                title="",
                )
            log.info(f'Generated cluster contacts chordchart file: {chordp}')

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
    # define readme data type content
    dttype_info = {
        'res1': 'Chain-Resname-ResID key identifying first residue',
        'res2': 'Chain-Resname-ResID key identifying second residue',
        'ca-ca-dist': 'observed distances between the two Ca',
        'ca-ca-cont-ratio': 'ratio of times the ca-ca-dist was observed under threshold',  # noqa : E501
        'shortest-dist': 'observed shortest distance between the two residues',
        'shortest-cont-ratio': 'ratio of times the shortest distance was observed under threshold',  # noqa : E501
        'contact-type': 'type of contacts between the two residues',
        'contact-type-ratio': 'ratio of times the type of contacts between the two residues is observed',  # noqa : E501
        }

    # initiate file content
    tsvdt = [header]
    for res_res_cont in res_res_contacts:
        tsvdt.append([str(res_res_cont[h]) for h in header])
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

    return path


def tsv_to_heatmap(
        tsv_path: Path,
        sep: str = '\t',
        data_key: str = 'ca-ca-dist',
        contact_threshold: float = 7.5,
        colorscale: str = 'Greys',
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
    if 'ratio' in data_key:
        data_label = 'ratio'
        np.fill_diagonal(matrix, 1)
    else:
        data_label = 'distance'

    # Generate heatmap
    output_filepath = heatmap_plotly(
        matrix,
        labels={'color': data_label},
        xlabels=labels,
        ylabels=labels,
        color_scale=color_scale,
        output_fname=output_fname,
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
    return f'{color_scale}_r' if 'ratio' not in data_key else color_scale


#############
# Start of the chord chart functions
#############
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


def get_ideogram_ends(ideogram_len: NDFloat, gap: float) -> list[tuple[float, float]]:
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
        showlegend=False,
        width=plot_size,
        height=plot_size,
        margin=dict(t=25, b=25, l=25, r=25),
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
    # Compute ratio
    ratio_dist = (dist - min_dist) / (max_dist - min_dist)
    # Obtain corresponding weight
    weight = ((min_weight - max_weight) * ratio_dist) + max_weight
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
    str
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
            color_weith = to_color_weight(dist_matrix[k][j], 9.5)
            rgba_color = to_rgba_color_string(connect_color, color_weith)

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
            rescolor = RESIDUES_COLORS[resname]
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
    # Write it as html file
    fig.write_html(output_fpath)
    return output_fpath


def tsv_to_chordchart(
        tsv_path: Path,
        sep: str = '\t',
        data_key: str = 'ca-ca-dist',
        contact_threshold: float = 7.5,
        filter_intermolecular_contacts: bool = True,
        output_fname: Union[Path, str] = 'contacts_chordchart.html',
        title: str = 'Chord diagram',
        ) -> None:
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
        )
    
    return chord_chart_fpath
