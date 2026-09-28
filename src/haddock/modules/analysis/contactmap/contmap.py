"""Module computing contact maps of complexes, alone or grouped by cluster.

Conventions used throughout this module:

* Residues are identified by keys of the form ``chain-resid-resname``
  (e.g. ``A-52A-ALA``), where ``resid`` includes the insertion code, if any.
  Always use :func:`parse_reskey` to split them, as ``resid`` can be negative.
* The ``ca-ca-dist`` is computed between one reference atom per residue:
  ``CA`` for amino acids, ``C4'`` for nucleotides and ``C1`` for
  carbohydrates. Residues without a reference atom (ions, ligands, ...)
  get a ``nan`` distance.
* The ``contact-type`` combines the *classes* of the two residues
  (apolar, polar, positive, negative, nucleotide, carbohydrate or unknown).
  It describes which kinds of residues are in contact, it is **not** a
  detected interaction such as a hydrogen bond or a salt bridge. Note the
  following choices: CYS and TRP are considered polar, HIS is polar unless
  named with a protonated variant (HIP, HSP), GLY is apolar, aromatic
  residues are not grouped in a class of their own, and charged sugars
  (e.g. sialic acids) are in the carbohydrate class.

Chord diagram functions were adapted from:
https://plotly.com/python/v3/filled-chord-diagram/
"""

import html
import os
from collections import Counter
from pathlib import Path

import numpy as np
import plotly.graph_objs as go
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist, squareform

from haddock import log
from haddock.core.typing import (
    Any,
    NDArray,
    NDFloat,
    Optional,
    SupportsRun,
    Union,
)
from haddock.libs.libontology import PDBFile
from haddock.libs.libpdb import (
    slc_chainid,
    slc_element,
    slc_icode,
    slc_name,
    slc_resname,
    slc_resseq,
    slc_x,
    slc_y,
    slc_z,
)
from haddock.libs.libplots import fig_to_html, heatmap_plotly

###############################
# Global variable definitions #
###############################
AMINO_ACID_CLASSES = {
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
# Protonation variants and common modified amino acids
MODIFIED_AMINO_ACID_CLASSES = {
    "HID": "polar",
    "HIE": "polar",
    "HSD": "polar",
    "HSE": "polar",
    "HIP": "positive",
    "HSP": "positive",
    "CYX": "polar",
    "CYM": "polar",
    "CYC": "polar",
    "CYF": "polar",
    "CFE": "polar",
    "SEC": "polar",
    "ASH": "polar",
    "GLH": "polar",
    "LYN": "polar",
    "ALY": "polar",
    "PCA": "polar",
    "CIR": "polar",
    "MSE": "apolar",
    "HYP": "apolar",
    "HY3": "apolar",
    "SEP": "negative",
    "TPO": "negative",
    "TOP": "negative",
    "PTR": "negative",
    "TYP": "negative",
    "TYS": "negative",
    "NEP": "negative",
    "CSP": "negative",
    "MLY": "positive",
    "MLZ": "positive",
    "M3L": "positive",
    "HLY": "positive",
}
# Standard and common modified nucleotides
NUCLEOTIDES = (
    "A", "C", "G", "T", "U", "I",
    "DA", "DC", "DG", "DT", "DU", "DI", "DJ",
    "PSU", "5MU", "5MC", "1MA", "2MG", "M2G", "7MG",
    "OMC", "OMG", "OMU", "H2U", "4SU",
)  # fmt: skip
# Carbohydrates supported by HADDOCK (see cns/toppar/carbohydrate.top)
CARBOHYDRATES = (
    "GLC", "BGC", "GLA", "GAL", "MAN", "BMA", "NGA", "NDG", "NAM", "NAA",
    "NAG", "GCS", "MAG", "A2G", "NGM", "FUC", "FUL", "FCA", "FCB", "SIA",
    "SIB", "XYP", "RAM", "GXL", "BDP", "MMA", "XYS", "ABE",
)  # fmt: skip
CARBOHYDRATE_CLASS = "carbohydrate"
RESIDUE_CLASSES = {
    **AMINO_ACID_CLASSES,
    **MODIFIED_AMINO_ACID_CLASSES,
    **{nuc: "nucleotide" for nuc in NUCLEOTIDES},
    **{sugar: CARBOHYDRATE_CLASS for sugar in CARBOHYDRATES},
}
UNKNOWN_CLASS = "unknown"

# Reference atoms used to compute the `ca-ca-dist`, by order of preference
NUCLEOTIDE_REFERENCE_ATOMS = ("C4'", "C4*")
CARBOHYDRATE_REFERENCE_ATOM = "C1"

PI = np.pi

# Colors of residue classes (Okabe-Ito colorblind safe palette)
RESIDUE_CLASS_COLORS = {
    "apolar": (230, 159, 0),
    "polar": (0, 158, 115),
    "positive": (0, 114, 178),
    "negative": (213, 94, 0),
    "nucleotide": (204, 121, 167),
    CARBOHYDRATE_CLASS: (86, 180, 233),
    UNKNOWN_CLASS: (153, 153, 153),
}
# Colors of the most informative residue class pairs: same class pairs use
# the class color, and opposite charges have their own color.
# Keys are sorted class names, see `get_pair_color()`.
CONNECT_COLORS = {
    **{
        f"{resclass}-{resclass}": color
        for resclass, color in RESIDUE_CLASS_COLORS.items()
        if resclass != UNKNOWN_CLASS
    },
    "negative-positive": (106, 61, 154),
}
# Color of all other residue class pairs
OTHER_PAIR_COLOR = (200, 200, 200)

# Chain colors, kept neutral to not be confused with residue classes
CHAIN_COLORS = [
    "rgba(64, 64, 64, 0.85)",
    "rgba(160, 160, 160, 0.85)",
    "rgba(96, 72, 48, 0.85)",
    "rgba(52, 82, 110, 0.85)",
    "rgba(110, 110, 60, 0.85)",
    "rgba(100, 64, 100, 0.85)",
    "rgba(40, 96, 88, 0.85)",
    "rgba(170, 130, 100, 0.85)",
    "rgba(80, 80, 130, 0.85)",
]

# Chord chart sizes (in pixels)
MIN_CHORDCHART_SIZE = 500
CHORDCHART_LEGEND_WIDTH = 150

# Maximum number of atom-atom distances held in memory at once
MAX_DISTANCES_BLOCK = 5_000_000

# Output files headers
RES_CONTACTS_HEADER = ["res1", "res2", "ca-ca-dist", "contact-type", "shortest-dist"]
CLUSTER_RES_CONTACTS_HEADER = [
    "res1",
    "res2",
    "ca-ca-cont-probability",
    "ca-ca-dist",
    "contact-type",
    "shortest-cont-probability",
    "shortest-dist",
]
HEAVY_CONTACTS_HEADER = ["atom1", "atom2", "dist"]
CLUSTER_HEAVY_CONTACTS_HEADER = ["atom1", "atom2", "avg_dist", "nb_dists", "std_dist"]

# Suffixes of the files generated by a job, in the order of the report
CONTACTMAP_OUTPUT_SUFFIXES = (
    "heatmap.html",
    "chordchart.html",
    "contacts.tsv",
    "interchain_contacts.tsv",
    "heavyatoms_interchain_contacts.tsv",
)


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

    def run(self) -> tuple[list[dict], list[dict]]:
        """Process analysis of contacts of a PDB structure."""
        # Load pdb
        pdb_dt = extract_pdb_dt(self.model)
        # Extract all coordinates
        all_coords, resid_keys, resid_dt = get_ordered_coords(pdb_dt)
        # Compute residue-residue distances
        ref_dists, shortest_dists = compute_residue_distances(
            all_coords,
            resid_keys,
            resid_dt,
        )
        res_res_contacts = gen_contacts_dt(
            ref_dists,
            shortest_dists,
            resid_keys,
            resid_dt,
        )
        all_heavy_interchain_contacts = extract_heavyatom_contacts(
            all_coords,
            resid_keys,
            resid_dt,
            contact_distance=self.params["shortest_dist_threshold"],
        )

        # generate outputs for single models
        if self.params["single_model_analysis"]:
            self.generate_output(
                res_res_contacts,
                all_heavy_interchain_contacts,
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
        fpath = write_res_contacts(
            res_res_contacts,
            RES_CONTACTS_HEADER,
            f"{self.output}_contacts.tsv",
            interchain_data=interchain_tsv_data(self.output, self.params),
        )
        log.info(f"Generated contacts file: {fpath}")
        self.files["res-res-contacts"] = fpath
        self.files["res-res-interchain-contacts"] = (
            f"{self.output}_interchain_contacts.tsv"
        )

        # Generate corresponding heatmap and chord chart
        generate_figures(
            fpath,
            self.output,
            self.params,
            heatmap_datatype="ca-ca-dist",
            files=self.files,
        )

        # Write interchain heavy atoms contacts tsv file
        fpath2 = write_res_contacts(
            all_heavy_interchain_contacts,
            HEAVY_CONTACTS_HEADER,
            f"{self.output}_heavyatoms_interchain_contacts.tsv",
        )
        log.info(f"Generated contacts file: {fpath2}")
        self.files["atom-atom-interchain-contacts"] = fpath2


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
            Dictionary holding list of contact data
        contact_keys : list[str]
            Order of the keys to access the dictionary
        contacts : list[dict]
            Single model contact data.
        key1 : str
            Name of the key to access first entry in data.
        key2 : str
            Name of the key to access second entry in data.
        """
        for cont in contacts:
            # Check key
            combined_key = f"{cont[key2]}/{cont[key1]}"  # reversed
            if combined_key not in contacts_holder:
                combined_key = f"{cont[key1]}/{cont[key2]}"  # normal
                if combined_key not in contacts_holder:
                    # Add key order
                    contact_keys.append(combined_key)
                    # Initiate key
                    contacts_holder[combined_key] = {
                        k: [] for k in cont if k not in (key1, key2)
                    }
            # Add data
            for dtk, values in contacts_holder[combined_key].items():
                values.append(cont[dtk])

    def run(self) -> None:
        """Process analysis of contacts of a set of PDB structures."""
        # initiate holding variables
        clusters_contacts: dict = {}  # Residue-residue contacts
        resres_keys_list: list[str] = []  # Ordered residue-residue keys
        clusters_heavyatm_contacts: dict = {}  # Interchain atom-atom contacts
        atat_keys_list: list[str] = []  # Ordered interchain atom-atom keys

        # loop over models/structures
        for pdb_path in self.models:
            contact_map_obj = ContactsMap(
                pdb_path,
                f"{self.output}_{pdb_path.stem}",
                self.params,
            )
            pdb_contacts, interchain_heavy_contacts = contact_map_obj.run()

            # Aggregate residue-residue contacts
            self.aggregate_contacts(
                clusters_contacts,
                resres_keys_list,
                pdb_contacts,
                "res1",
                "res2",
            )
            # Aggregate heavy atoms contacts
            self.aggregate_contacts(
                clusters_heavyatm_contacts,
                atat_keys_list,
                interchain_heavy_contacts,
                "atom1",
                "atom2",
            )

        # Summarize heavy atoms contacts
        heavy_atm_clust_list = []
        for atatk in atat_keys_list:
            at1, at2 = atatk.split("/")
            h_dists = clusters_heavyatm_contacts[atatk]["dist"]
            heavy_atm_clust_list.append(
                {
                    "atom1": at1,
                    "atom2": at2,
                    "nb_dists": len(h_dists),
                    "avg_dist": round(float(np.mean(h_dists)), 2),
                    "std_dist": round(float(np.std(h_dists)), 2),
                }
            )
        if not heavy_atm_clust_list:
            log.info(f"No interchain heavy atoms contacts found for {self.output}")
        hfpath = write_res_contacts(
            heavy_atm_clust_list,
            CLUSTER_HEAVY_CONTACTS_HEADER,
            f"{self.output}_heavyatoms_interchain_contacts.tsv",
        )
        log.info(f"Generated heavy atoms interchain contacts file: {hfpath}")
        self.files["atom-atom-interchain-contacts"] = hfpath

        # Summarize residue-residue contacts
        ca_ca_threshold = self.params["ca_ca_dist_threshold"]
        shortest_threshold = self.params["shortest_dist_threshold"]
        combined_clusters_list = []
        for combined_key in resres_keys_list:
            dt = clusters_contacts[combined_key]
            ca_ca_dists = np.asarray(dt["ca-ca-dist"], dtype=float)
            shortest_dists = np.asarray(dt["shortest-dist"], dtype=float)
            # Most represented contact type
            cont_t = Counter(dt["contact-type"]).most_common(1)[0][0]
            res1, res2 = combined_key.split("/")
            combined_clusters_list.append(
                {
                    "res1": res1,
                    "res2": res2,
                    "ca-ca-dist": round(nanmean(ca_ca_dists), 1),
                    "ca-ca-cont-probability": contact_probability(
                        ca_ca_dists, ca_ca_threshold
                    ),
                    "shortest-dist": round(nanmean(shortest_dists), 1),
                    "shortest-cont-probability": contact_probability(
                        shortest_dists, shortest_threshold
                    ),
                    "contact-type": cont_t,
                }
            )

        # write contacts
        fpath = write_res_contacts(
            combined_clusters_list,
            CLUSTER_RES_CONTACTS_HEADER,
            f"{self.output}_contacts.tsv",
            interchain_data=interchain_tsv_data(self.output, self.params),
        )
        log.info(f"Generated contacts file: {fpath}")
        self.files["res-res-contacts"] = fpath
        self.files["res-res-interchain-contacts"] = (
            f"{self.output}_interchain_contacts.tsv"
        )

        # Generate corresponding heatmap and chord chart
        generate_figures(
            fpath,
            self.output,
            self.params,
            heatmap_datatype=self.params["cluster_heatmap_datatype"],
            files=self.files,
        )

        self.terminated = True


def get_data_threshold(data_key: str, params: dict) -> float:
    """Return the threshold defining a contact for a given data type.

    Parameters
    ----------
    data_key : str
        Name of the data column (e.g. `shortest-dist`).
    params : dict
        Module parameters.

    Return
    ------
    threshold : float
        1 for probabilities, otherwise the corresponding distance threshold.
    """
    if "probability" in data_key:
        return 1.0
    if data_key.startswith("ca-ca"):
        return params["ca_ca_dist_threshold"]
    return params["shortest_dist_threshold"]


def interchain_tsv_data(output: Union[str, Path], params: dict) -> dict:
    """Build the parameters used to write the interchain contacts file.

    The interchain file is filtered with the same data as the chord chart.
    """
    data_key = params["chordchart_datatype"]
    return {
        "path": f"{output}_interchain_contacts.tsv",
        "data_key": data_key,
        "contact_threshold": get_data_threshold(data_key, params),
    }


def generate_figures(
    fpath: Union[str, Path],
    output: Union[str, Path],
    params: dict,
    heatmap_datatype: str,
    files: dict[str, Union[str, Path]],
) -> None:
    """Generate the heatmap and chord chart of a contacts file.

    Parameters
    ----------
    fpath : Union[str, Path]
        Path to the residue-residue contacts tsv file.
    output : Union[str, Path]
        Output basename.
    params : dict
        Module parameters.
    heatmap_datatype : str
        Data key used to draw the heatmap.
    files : dict[str, Union[str, Path]]
        Holder of generated files, updated in place.
    """
    if params["generate_heatmap"]:
        heatmap = tsv_to_heatmap(
            fpath,
            data_key=heatmap_datatype,
            contact_threshold=get_data_threshold(heatmap_datatype, params),
            colorscale=params["color_ramp"],
            output_fname=f"{output}_heatmap.html",
            offline=params["offline"],
        )
        if heatmap:
            log.info(f"Generated contacts heatmap file: {heatmap}")
            files["res-res-contactmap"] = heatmap

    if params["generate_chordchart"]:
        data_key = params["chordchart_datatype"]
        chordp = tsv_to_chordchart(
            fpath,
            data_key=data_key,
            contact_threshold=get_data_threshold(data_key, params),
            output_fname=f"{output}_chordchart.html",
            filter_intermolecular_contacts=True,
            title=Path(output).stem.replace("_", " "),
            offline=params["offline"],
        )
        if chordp:
            log.info(f"Generated contacts chordchart file: {chordp}")
            files["res-res-chordchart"] = chordp


def nanmean(values: NDFloat) -> float:
    """Compute the mean of non-nan values, nan if there are none."""
    finite = values[~np.isnan(values)]
    return float(np.mean(finite)) if finite.size else float("nan")


def contact_probability(values: NDFloat, threshold: float) -> float:
    """Compute the fraction of values under threshold (nan are not contacts)."""
    return round(float(np.sum(values <= threshold)) / len(values), 2)


def list_job_outputs(output: Union[str, Path]) -> list[str]:
    """List existing files generated for a given output basename."""
    candidates = [f"{output}_{suffix}" for suffix in CONTACTMAP_OUTPUT_SUFFIXES]
    return [fpath for fpath in candidates if os.path.exists(fpath)]


def make_contactmap_report(
    contactmap_jobs: list[Union[ContactsMap, ClusteredContactMap]],
    outputpath: Union[str, Path],
) -> Union[str, Path]:
    """Generate a HTML navigation page holding all generated files.

    Parameters
    ----------
    contactmap_jobs : list[Union[ClusteredContactMap, ContactsMap]]
        All the terminated jobs
    outputpath : Union[str, Path]
        Output filepath where to write the report.

    Returns
    -------
    outputpath: Union[str, Path]
        Path to the generated report.
    """
    # List (basename, files) to be reported
    ordered_outputs: list[tuple[str, list[str]]] = []
    for job in contactmap_jobs:
        ordered_outputs.append((str(job.output), list_job_outputs(job.output)))
        # Single model outputs of clustered analyses
        if isinstance(job, ClusteredContactMap) and job.params.get(
            "single_model_analysis"
        ):
            for model in job.models:
                model_output = f"{job.output}_{Path(model).stem}"
                ordered_outputs.append((model_output, list_job_outputs(model_output)))

    ordered_files: list[str] = []
    for output, job_files in ordered_outputs:
        basepath = f"{output}_"
        job_list = [
            (
                f'<a href="{html.escape(fpath)}" target="_blank">'
                f"{html.escape(fpath.replace(basepath, '', 1))}</a>"
            )
            for fpath in job_files
        ]
        ordered_files.append(f"<b>{html.escape(output)}:</b> {', '.join(job_list)}")

    # Combine all jobs outputs as a list
    all_access = "</li>\n            <li>".join(ordered_files)
    htmldt = f"""
    <div id="contactmap_report">
        <ul>
            <li>
            {all_access}
            </li>
        </ul>
    </div>
"""

    with open(outputpath, "w") as reportout:
        reportout.write(htmldt)
    log.info(f"Generated report file: {outputpath}")
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
    for model in models:
        cluster_key = (model.clt_id, model.clt_rank)
        clust_sets.setdefault(cluster_key, []).append(model)
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
        List of top `X` best models. If models cannot be sorted by score
        (no score attribute or undefined scores), the input order is kept.
    """
    try:
        sorted_models = sorted(models, key=lambda m: m.score)
    except (AttributeError, TypeError):
        sorted_models = list(models)
    return sorted_models[:topX]


####################
# Define functions #
####################
def parse_reskey(key: str) -> tuple[str, str, str]:
    """Split a residue key into its chain, residue id and residue name.

    Parameters
    ----------
    key : str
        Residue key of the form `chain-resid-resname`, e.g. `A--5-ALA`.

    Return
    ------
    chain, resid, resname : tuple[str, str, str]
    """
    chain, rest = key.split("-", 1)
    resid, resname = rest.rsplit("-", 1)
    return chain, resid, resname


def get_residue_class(resname: str) -> str:
    """Return the class of a residue (apolar, polar, ..., unknown)."""
    return RESIDUE_CLASSES.get(resname.strip(), UNKNOWN_CLASS)


def is_hydrogen(line: str) -> bool:
    """Check if an ATOM/HETATM record is a hydrogen (or deuterium).

    The element column is used when present. Otherwise the atom name is used,
    ignoring leading digits (e.g. `1HB`), and single atom residues named after
    their atom (e.g. mercury `HG`) are not considered hydrogens.
    """
    element = line[slc_element].strip()
    if element:
        return element.upper() in ("H", "D")
    atname = line[slc_name].strip()
    if atname == line[slc_resname].strip():
        return False
    return atname.lstrip("0123456789")[:1] == "H"


def extract_pdb_dt(path: Path) -> dict:
    """Read and extract ATOM/HETATM records from a pdb file.

    Only the first model of multi-model files is read, hydrogens are skipped
    and only the first alternate location of each atom is kept.

    Parameters
    ----------
    path : Path
        Path to a pdb file.

    Return
    ------
    pdb_chains : dict
        A dictionary of the pdb file accessible using chains as keys.
    """
    pdb_chains: dict = {"chain_order": []}
    with open(path, "r") as f:
        for line in f:
            if line.startswith("ENDMDL"):
                break
            if not line.startswith(("ATOM", "HETATM")):
                continue

            resname = line[slc_resname].strip()
            chainid = line[slc_chainid]
            # Include insertion code in residue identifier
            resid = line[slc_resseq].strip() + line[slc_icode].strip()

            # Check if chain already parsed
            if chainid not in pdb_chains:
                pdb_chains["chain_order"].append(chainid)
                pdb_chains[chainid] = {"order": []}
            chain_dt = pdb_chains[chainid]

            # Check if new resid id
            if resid not in chain_dt:
                chain_dt["order"].append(resid)
                chain_dt[resid] = {
                    "index": len(chain_dt["order"]) - 1,
                    "resname": resname,
                    "chainid": chainid,
                    "resid": resid,
                    "position": len(chain_dt["order"]),
                    "atoms_order": [],
                    "atoms": {},
                }

            if is_hydrogen(line):
                continue
            atname = line[slc_name].strip()
            residue = chain_dt[resid]
            # Keep only the first alternate location of an atom
            if atname in residue["atoms"]:
                continue
            residue["atoms_order"].append(atname)
            residue["atoms"][atname] = extract_pdb_coords(line)

    return pdb_chains


def extract_pdb_coords(line: str) -> list[float]:
    """Extract coordinates from a PDB line.

    Parameters
    ----------
    line : str
        A standard ATOM/HETATM pdb record.

    Return
    ------
    coords : list[float]
        List of the X, Y and Z coordinate of this atom.
    """
    x = float(line[slc_x].strip())
    y = float(line[slc_y].strip())
    z = float(line[slc_z].strip())
    return [x, y, z]


def get_reference_atom(atoms_order: list[str], resname: str = "") -> Optional[str]:
    """Find the atom used to compute the `ca-ca-dist` of a residue.

    `CA` is only considered an alpha carbon if a backbone `N` or `C` atom is
    also present, so that calcium ions are not mistaken for amino acids.
    `C1` is only used for known carbohydrates, as many ligands have one.

    Parameters
    ----------
    atoms_order : list[str]
        Names of the heavy atoms of the residue.
    resname : str
        Name of the residue.

    Return
    ------
    ref_atom : Optional[str]
        Name of the reference atom, None if the residue has none.
    """
    if "CA" in atoms_order and ("N" in atoms_order or "C" in atoms_order):
        return "CA"
    if get_residue_class(resname) == CARBOHYDRATE_CLASS:
        if CARBOHYDRATE_REFERENCE_ATOM in atoms_order:
            return CARBOHYDRATE_REFERENCE_ATOM
        return None
    for atname in NUCLEOTIDE_REFERENCE_ATOMS:
        if atname in atoms_order:
            return atname
    return None


def get_ordered_coords(
    pdb_chains: dict,
) -> tuple[NDFloat, list[str], dict]:
    """Generate array of all atom coordinates.

    Residues without any heavy atom are ignored.

    Parameters
    ----------
    pdb_chains : dict
        A dictionary of the pdb file accessible using chains as keys,
         as provided by the `extract_pdb_dt()` function.

    Return
    ------
    all_coords : NDFloat
        (N, 3) array of all atomic coordinates, contiguous by residue.
    resid_keys : list[str]
        Ordered list of residues keys.
    resid_dt : dict
        Dictionary of coordinates indices for each residue.
    """
    all_coords: list[list[float]] = []
    resid_keys: list[str] = []
    resid_dt: dict = {}
    for chainid in pdb_chains["chain_order"]:
        for resid in pdb_chains[chainid]["order"]:
            residue = pdb_chains[chainid][resid]
            atoms_order = residue["atoms_order"]
            reskey = f"{chainid}-{resid}-{residue['resname']}"
            if not atoms_order:
                log.debug(f"Residue {reskey} has no heavy atom and is ignored")
                continue
            start = len(all_coords)
            all_coords += [residue["atoms"][atname] for atname in atoms_order]
            ref_atom = get_reference_atom(atoms_order, residue["resname"])
            resid_dt[reskey] = {
                "atoms_indices": list(range(start, len(all_coords))),
                "resname": residue["resname"],
                "chainid": chainid,
                "atoms_order": atoms_order,
                "ref": None
                if ref_atom is None
                else start + atoms_order.index(ref_atom),
            }
            resid_keys.append(reskey)
    return np.asarray(all_coords, dtype=float).reshape(-1, 3), resid_keys, resid_dt


def compute_residue_distances(
    all_coords: NDFloat,
    resid_keys: list[str],
    resid_dt: dict,
    max_block: int = MAX_DISTANCES_BLOCK,
) -> tuple[NDFloat, NDFloat]:
    """Compute residue-residue reference atom and shortest distances.

    Atom-atom distances are computed by blocks of residues, so that the full
    atom-atom distance matrix is never held in memory.

    Parameters
    ----------
    all_coords : NDFloat
        (N, 3) array of atomic coordinates, contiguous by residue.
    resid_keys : list[str]
        Ordered list of residues keys.
    resid_dt : dict
        Residues data as returned by `get_ordered_coords()`.
    max_block : int
        Maximum number of atom-atom distances computed at once.

    Return
    ------
    ref_dists : NDFloat
        (R, R) distances between residues reference atoms, nan if missing.
    shortest_dists : NDFloat
        (R, R) shortest heavy atom distances between residues.
    """
    nb_res = len(resid_keys)
    ref_dists = np.full((nb_res, nb_res), np.nan)
    shortest_dists = np.zeros((nb_res, nb_res))
    if nb_res == 0:
        return ref_dists, shortest_dists

    # Reference atoms distances
    ref_indices = np.array(
        [-1 if resid_dt[k]["ref"] is None else resid_dt[k]["ref"] for k in resid_keys]
    )
    has_ref = ref_indices >= 0
    ref_coords = all_coords[ref_indices[has_ref]]
    ref_dists[np.ix_(has_ref, has_ref)] = cdist(ref_coords, ref_coords)

    # Shortest distances, computed by blocks of residues
    bounds = [
        (resid_dt[k]["atoms_indices"][0], resid_dt[k]["atoms_indices"][-1] + 1)
        for k in resid_keys
    ]
    starts = np.array([start for start, _end in bounds])
    block_atoms = max(1, max_block // len(all_coords))
    ri = 0
    while ri < nb_res:
        rj = ri + 1
        while rj < nb_res and bounds[rj][1] - bounds[ri][0] <= block_atoms:
            rj += 1
        first_atom, last_atom = bounds[ri][0], bounds[rj - 1][1]
        block = cdist(all_coords[first_atom:last_atom], all_coords)
        # Minimum over the atoms of each residue (columns then rows)
        per_residue = np.minimum.reduceat(block, starts, axis=1)
        shortest_dists[ri:rj] = np.minimum.reduceat(
            per_residue,
            starts[ri:rj] - first_atom,
            axis=0,
        )
        ri = rj
    return ref_dists, shortest_dists


def extract_submatrix(
    matrix: NDArray,
    indices: list[int],
    indices2: Optional[list[int]] = None,
) -> NDArray:
    """Extract submatrix based on desired indices.

    Parameters
    ----------
    matrix : NDArray
        A N*N matrix.
    indices : list[int]
        List of `row` indices to extract from this matrix
    indices2 : list[int]
        List of `columns` indices to extract from this matrix.
         if unspecified, indices2 == indices and symmetric matrix
         is extracted.

    Return
    ------
    submat : NDArray
        The extracted submatrix.
    """
    if indices2 is None:
        indices2 = indices
    return matrix[np.ix_(indices, indices2)]


def gen_contacts_dt(
    ref_dists: NDFloat,
    shortest_dists: NDFloat,
    resid_keys: list[str],
    resid_dt: dict,
) -> list[dict]:
    """Generate residue-residue contacts data (half matrix).

    Parameters
    ----------
    ref_dists : NDFloat
        (R, R) distances between residues reference atoms.
    shortest_dists : NDFloat
        (R, R) shortest distances between residues.
    resid_keys : list[str]
        Ordered list of residues keys.
    resid_dt : dict
        Residues data as returned by `get_ordered_coords()`.

    Return
    ------
    contacts : list[dict]
        One dictionary per residue pair, in row major half matrix order.
    """
    classes = [get_residue_class(resid_dt[k]["resname"]) for k in resid_keys]
    rows, cols = np.triu_indices(len(resid_keys), k=1)
    ref_values = np.round(ref_dists[rows, cols], 1).tolist()
    shortest_values = np.round(shortest_dists[rows, cols], 1).tolist()
    return [
        {
            "res1": resid_keys[i],
            "res2": resid_keys[j],
            "ca-ca-dist": ref_dist,
            "shortest-dist": shortest_dist,
            "contact-type": f"{classes[i]}-{classes[j]}",
        }
        for i, j, ref_dist, shortest_dist in zip(
            rows.tolist(), cols.tolist(), ref_values, shortest_values
        )
    ]


def extract_heavyatom_contacts(
    all_coords: NDFloat,
    resid_keys: list[str],
    resid_dt: dict,
    contact_distance: float = 4.5,
) -> list[dict[str, Union[float, str]]]:
    """Find interchain heavy atom pairs closer than a distance.

    Parameters
    ----------
    all_coords : NDFloat
        (N, 3) array of atomic coordinates, contiguous by residue.
    resid_keys : list[str]
        Ordered list of residues keys.
    resid_dt : dict
        Residues data as returned by `get_ordered_coords()`.
    contact_distance : float
        Distance defining a contact (inclusive).

    Return
    ------
    all_contacts : list[dict[str, Union[float, str]]]
        Interchain contacts, sorted by residue pair and then atoms.
    """
    if len(all_coords) < 2:
        return []
    atom_labels: list[str] = []
    atom_chains: list[str] = []
    atom_res = np.empty(len(all_coords), dtype=int)
    for ri, reskey in enumerate(resid_keys):
        resdt = resid_dt[reskey]
        atom_res[resdt["atoms_indices"]] = ri
        atom_labels += [f"{reskey}-{atname}" for atname in resdt["atoms_order"]]
        atom_chains += [resdt["chainid"]] * len(resdt["atoms_order"])
    chains = np.array(atom_chains)

    # Pairs (i, j) with i < j and distance <= contact_distance
    pairs = cKDTree(all_coords).query_pairs(r=contact_distance, output_type="ndarray")
    if len(pairs) == 0:
        return []
    ai, aj = pairs[:, 0], pairs[:, 1]
    interchain = chains[ai] != chains[aj]
    ai, aj = ai[interchain], aj[interchain]
    order = np.lexsort((aj, ai, atom_res[aj], atom_res[ai]))
    ai, aj = ai[order], aj[order]
    dists = np.round(np.linalg.norm(all_coords[ai] - all_coords[aj], axis=1), 2)
    return [
        {"atom1": atom_labels[i], "atom2": atom_labels[j], "dist": dist}
        for i, j, dist in zip(ai.tolist(), aj.tolist(), dists.tolist())
    ]


def get_cont_type(resn1: str, resn2: str) -> str:
    """Generate the residue class pair of two residues.

    Parameters
    ----------
    resn1 : str
       3 letters code of first residue.
    resn2 : str
       3 letters code of second residue.

    Return
    ------
    pol_key : str
        Combined residues classes, e.g. `polar-negative`.
    """
    return f"{get_residue_class(resn1)}-{get_residue_class(resn2)}"


def write_res_contacts(
    res_res_contacts: list[dict],
    header: list[str],
    path: Union[Path, str],
    sep: str = "\t",
    interchain_data: Optional[dict] = None,
) -> Union[Path, str]:
    """Write a tsv file based on residues-residues contacts data.

    Parameters
    ----------
    res_res_contacts : list[dict]
        List of dict holding data for each residue-residue contacts.
    header : list[str]
        Ordered list of keys to access in the dicts.
    path : Union[Path, str]
        Path to the output file to generate.
    sep : str
        Character used to separate data within a line.
    interchain_data : Optional[dict]
        If provided, also write the interchain contacts for which the
        `data_key` value is <= `contact_threshold` in the file `path`.

    Return
    ------
    path : Union[Path, str]
        Path to the generated file.
    """
    dttype_info = {
        "res1": "Chain-ResID-Resname key identifying first residue (ResID includes the insertion code, if any)",
        "res2": "Chain-ResID-Resname key identifying second residue (ResID includes the insertion code, if any)",
        "ca-ca-dist": "Distance between the reference atoms of the two residues (CA for amino acids, C4' for nucleotides, C1 for carbohydrates), nan if one of them has no reference atom",
        "ca-ca-cont-probability": "Fraction of times a contact is observed under the ca-ca-dist threshold over all analysed models of the same cluster",
        "shortest-dist": "Observed shortest distance between the heavy atoms of the two residues",
        "shortest-cont-probability": "Fraction of times a contact is observed under the shortest-dist threshold over all analysed models of the same cluster",
        "contact-type": "Classes of the two residues (apolar, polar, positive, negative, nucleotide, carbohydrate or unknown), not a detected interaction type",
        "atom1": "Chain-ResID-Resname-AtomName key identifying first atom",
        "atom2": "Chain-ResID-Resname-AtomName key identifying second atom",
        "dist": "Observed distance between two atoms",
        "nb_dists": "Total number of observed distances",
        "avg_dist": "Cluster average distance",
        "std_dist": "Cluster distance standard deviation",
    }

    # Check for inter chain contacts
    gen_interchain_tsv = False
    if isinstance(interchain_data, dict):
        expected_keys = ("path", "contact_threshold", "data_key")
        missing = [k for k in expected_keys if k not in interchain_data]
        if missing:
            raise KeyError(f"Missing keys in interchain_data: {missing}")
        gen_interchain_tsv = True
        interchain_tsvdt: list[list[str]] = [header]

    # initiate file content
    tsvdt: list[list[str]] = [header]
    for res_res_cont in res_res_contacts:
        tsvdt.append([str(res_res_cont[h]) for h in header])
        if gen_interchain_tsv:
            chain1 = parse_reskey(res_res_cont["res1"])[0]
            chain2 = parse_reskey(res_res_cont["res2"])[0]
            if chain1 != chain2:
                value = res_res_cont[interchain_data["data_key"]]
                if value <= interchain_data["contact_threshold"]:
                    interchain_tsvdt.append(tsvdt[-1])

    # generate commented lines to be placed on top of file
    readme = [
        "#" * 80,
        "# This file contains extracted contacts half-matrix information",
        "#" * 80,
        "",
    ]
    for head in header[::-1]:
        readme.insert(2, f"# {head}: {dttype_info.get(head, 'No description')}")

    with open(path, "w") as tsvout:
        tsvout.write("\n".join(readme))
        tsvout.write("\n".join([sep.join(_) for _ in tsvdt]) + "\n")

    # Write inter chain file
    if gen_interchain_tsv:
        readme[1] = readme[1].replace("contacts half-matrix", "interchain contacts")
        with open(interchain_data["path"], "w") as f:
            f.write("\n".join(readme))
            f.write("\n".join([sep.join(_) for _ in interchain_tsvdt]) + "\n")

    return path


def read_contacts_tsv(
    tsv_path: Union[Path, str],
    sep: str = "\t",
) -> tuple[list[str], list[list[str]]]:
    """Read the header and data rows of a contacts tsv file."""
    header: list[str] = []
    rows: list[list[str]] = []
    with open(tsv_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            s_ = line.strip("\n").split(sep)
            if not header:
                header = s_
            else:
                rows.append(s_)
    return header, rows


def tsv_to_heatmap(
    tsv_path: Union[Path, str],
    sep: str = "\t",
    data_key: str = "ca-ca-dist",
    contact_threshold: float = 7.5,
    colorscale: str = "Greys",
    output_fname: Union[Path, str] = "contacts.html",
    offline: bool = False,
) -> Optional[Union[Path, str]]:
    """Read a tsv file and generate a heatmap from it.

    Parameters
    ----------
    tsv_path : Union[Path, str]
        Path a the .tsv file containing contact data.
    sep : str
        Separator character used to split data in each line.
    data_key : str
        Data key used to draw the plot.
    contact_threshold : float
        Upper boundary of maximum value to be plotted.
         any value above it will be set to this value.
    output_fname : Union[Path, str]
        Path to the generated graph.

    Return
    ------
    output_filepath : Optional[Union[Path, str]]
        Path to the generated file, None if there is not enough data.
    """
    header, rows = read_contacts_tsv(tsv_path, sep=sep)
    half_matrix: list[float] = []
    labels: list[str] = []
    seen: set[str] = set()
    idx1, idx2, idxv = (header.index(k) for k in ("res1", "res2", data_key))
    for s_ in rows:
        for label in (s_[idx1], s_[idx2]):
            if label not in seen:
                seen.add(label)
                labels.append(label)
        # bound data to contact_threshold (nan values are kept)
        half_matrix.append(min(float(s_[idxv]), contact_threshold))

    if len(labels) < 2:
        log.warning(f"Not enough residues in {tsv_path} to generate a heatmap")
        return None

    matrix = squareform(half_matrix)

    # set data label
    color_scale = datakey_to_colorscale(data_key, color_scale=colorscale)
    if "probability" in data_key:
        data_label = "probability"
        np.fill_diagonal(matrix, 1)
    else:
        data_label = "distance"

    # Compute chains length
    chains_length: dict[str, int] = {}
    for label in labels:
        chainid = parse_reskey(label)[0]
        chains_length[chainid] = chains_length.get(chainid, 0) + 1
    # Compute chains delineations positions
    del_posi = [0]
    for length in chains_length.values():
        del_posi.append(del_posi[-1] + length)
    # Compute chains delineations lines
    chains_limits: list[dict[str, float]] = []
    for delpos in del_posi:
        # Vertical lines
        chains_limits.append(
            {
                "x0": delpos - 0.5,
                "x1": delpos - 0.5,
                "y0": -0.5,
                "y1": len(labels) - 0.5,
            }
        )
        # Horizontal lines
        chains_limits.append(
            {
                "y0": delpos - 0.5,
                "y1": delpos - 0.5,
                "x0": -0.5,
                "x1": len(labels) - 0.5,
            }
        )
    hovertemplate = (
        f" %{{y}}   &#8621;   %{{x}} <br> Contact {data_label}: %{{z}}<extra></extra>"
    )

    output_filepath = heatmap_plotly(
        matrix,
        labels={"color": data_label},
        xlabels=labels,
        ylabels=labels,
        color_scale=color_scale,
        output_fname=output_fname,
        offline=offline,
        delineation_traces=chains_limits,
        hovertemplate=hovertemplate,
    )

    return output_filepath


def datakey_to_colorscale(data_key: str, color_scale: str = "Greys") -> str:
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
    return f"{color_scale}_r" if "probability" not in data_key else color_scale


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
        raise ValueError("Incorrect interval ends")
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
        raise ValueError("Data array must have only two dimensions")
    if nb_rows != matrixshape[1]:
        raise ValueError("Data array must have (n,n) shape")
    return nb_rows


def get_ideogram_ends(
    ideogram_len: NDFloat,
    gap: float,
) -> list[tuple[float, float]]:
    """Generate ideogram ends.

    Parameters
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
    """Generate ideogram arc.

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
        Array of 2D coordinates defining an arc.
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
        # Point starting coordinates of this residue ideo
        start = float(ideo_end[0])
        # No ribbon to be formed
        if row_sum[k] == 0:
            ribbon_boundary.append([(0.0, 0.0) for i in range(len(ideo_ends))])
            continue
        row_ribbon_ends: list[tuple[float, float]] = []
        increment = (ideo_end[1] - start) / row_sum[k]
        for j in range(1, L + 1):
            # Skip if no ribbon to add for this k, j pair
            if matrix[k][j - 1] == 0:
                row_ribbon_ends.append((0.0, 0.0))
                continue
            end = float(start + increment)
            row_ribbon_ends.append((start, end))
            start = end
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
    if len(angle) != 3:
        raise ValueError("angle must have len = 3")
    b_cplx = np.array([np.exp(1j * angle[k]) for k in range(3)])
    b_cplx[1] = radius * b_cplx[1]
    control_points = list(zip(b_cplx.real, b_cplx.imag))
    return control_points


def ctrl_rib_chords(
    side1: tuple[float, float],
    side2: tuple[float, float],
    radius: float,
) -> list[list[tuple[float, float]]]:
    """Generate polygons points aiming at drawing ribbons.

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
    polygons : list[list[tuple[float, float]]]
        Control points of the two ribbon sides.
    """
    if len(side1) != 2 or len(side2) != 2:
        raise ValueError("the arc ends must be elements in a list of len 2")
    polygons = [
        control_pts(
            [side1[j], (side1[j] + side2[j]) / 2, side2[j]],
            radius,
        )
        for j in range(2)
    ]
    return polygons


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
        raise ValueError("control polygon must have 3 points")
    _a, _b, _c = control_points
    return f"M {_a[0]},{_a[1]} Q {_b[0]}, {_b[1]} {_c[0]}, {_c[1]}"


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
                raise ValueError("incorrect angle coordinates for ribbon")

        nr = int(40 * (theta0 - theta1) / PI)
        if nr <= 2:
            nr = 3
        theta = np.linspace(theta0, theta1, nr)
        pts = np.exp(1j * theta)  # points on arc in polar complex form

        string_arc: str = ""
        for k in range(len(theta)):
            string_arc += f"L {pts.real[k]!s}, {pts.imag[k]!s} "
        return string_arc
    else:
        raise ValueError(
            "the angle coordinates for an arc side of a ribbon must be in [0, 2*pi]"
        )


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
        "showline": False,
        "zeroline": False,
        "showgrid": False,
        "showticklabels": False,
        "title": "",
    }
    layout = go.Layout(
        title=title,
        xaxis=axis,
        yaxis=axis,
        showlegend=True,
        # Extra width accommodates the legend and keeps the circle round
        width=plot_size + CHORDCHART_LEGEND_WIDTH,
        height=plot_size,
        margin={"t": 25, "b": 25, "l": 25, "r": 25},
        hovermode="closest",
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
        "type": "path",
        "fillcolor": fill_color,
        "layer": "below",
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
    polygon = ctrl_rib_chords(side1, side2, radius)
    _b, _c = polygon
    path = make_q_bezier(_b)
    path += make_ribbon_arc(side2[0], side2[1])
    path += make_q_bezier(_c[::-1])
    path += make_ribbon_arc(side1[1], side1[0])

    return {
        "line": {"color": line_color, "width": 0.5},
        "path": path,
        "type": "path",
        "fillcolor": fill_color,
        "layer": "below",
    }


def get_chains_ideograms_ends(
    chains: dict[str, list[str]],
    gap: float = 2 * PI * 0.005,
) -> tuple[list[tuple[float, float]], NDFloat]:
    """Build ideogram ends to represent protein chains.

    Chains are drawn in the insertion order of the `chains` dictionary.

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
        Angular length of each chain ideogram.
    """
    chain_row_sum = [len(labels) for labels in chains.values()]
    chain_ideogram_length = 2 * PI * np.asarray(chain_row_sum)
    chain_ideogram_length /= sum(chain_row_sum)
    chain_ideogram_length -= gap * np.ones(len(chain_row_sum))
    chain_ideo_ends = get_ideogram_ends(chain_ideogram_length, gap)
    return chain_ideo_ends, chain_ideogram_length


def get_all_ideograms_ends(
    chains: dict[str, list[str]],
    gap: float = 2 * PI * 0.005,
) -> tuple[list[tuple[float, float]], list[tuple[float, float]]]:
    """Generate both chain and residues ideograms ends.

    Residues of each chain must be contiguous and in the same order as in
    the labels used to build the matrices.

    Parameters
    ----------
    chains : dict[str, list[str]]
        Dictionary mapping chains to list of residues labels.
    gap : float, optional
        Gap distance used to separate two ideograms, by default 2*PI*0.005

    Returns
    -------
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
    left = 0.0
    for ind, chain_labels in enumerate(chains.values()):
        right = left
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
        List of residues keys. e.g.: A-123-SER (chain A, serine 123)

    Returns
    -------
    chains : dict[str, list[str]]
        Dictionary mapping chains, in order of first appearance, with their
        respective set of residues labels.
    """
    chains: dict[str, list[str]] = {}
    for lab in labels:
        chains.setdefault(parse_reskey(lab)[0], []).append(lab)
    return chains


def group_indices_by_chain(labels: list[str]) -> list[int]:
    """Order label indices so that residues of each chain are contiguous.

    Chains are kept in order of first appearance, as well as residues within
    a chain.
    """
    label_chains = [parse_reskey(lab)[0] for lab in labels]
    chain_rank: dict[str, int] = {}
    for chain in label_chains:
        chain_rank.setdefault(chain, len(chain_rank))
    return sorted(range(len(labels)), key=lambda i: chain_rank[label_chains[i]])


def contacts_to_connect_matrix(
    matrix: NDArray,
    labels: list[str],
) -> NDArray:
    """Keep only interchain contacts of a contact matrix.

    Parameters
    ----------
    matrix : NDArray
        A square contact matrix (1 for contacts).
    labels : list[str]
        List of labels corresponding row & columns entries.

    Returns
    -------
    connect_matrix : NDArray
        The interchain connectivity matrix, without self contacts.
    """
    chains = np.array([parse_reskey(lab)[0] for lab in labels])
    interchain = chains[:, None] != chains[None, :]
    return (interchain & (np.asarray(matrix) == 1)).astype(int)


def to_nice_label(label: str) -> str:
    """Convert a label into a user friendly label.

    Parameters
    ----------
    label : str
        Label name as found in tsv

    Returns
    -------
    nicelabel : str
        User friendly description of the label.
    """
    chain, resid, resname = parse_reskey(label)
    return f"Chain {chain}, residue {resname} {resid}"


def to_color_weight(
    distance: float,
    max_dist: float,
    min_dist: float = 2.0,
    min_weight: float = 0.35,
    max_weight: float = 0.90,
) -> float:
    """Compute color weight based on distance.

    Parameters
    ----------
    distance : float
        The distance to weight.
    max_dist : float
        The maximum distance, usually the contact threshold.
    min_dist : float, optional
        The minimum distance, by default 2.
    min_weight : float, optional
        Color weight for the maximum distance, by default 0.35
    max_weight : float, optional
        Color weight for the minimum distance, by default 0.90

    Returns
    -------
    weight : float
        The color weight, in range [min_weight, max_weight]
    """
    if max_dist <= min_dist:
        return max_weight
    # Relative position of the distance in [min_dist, max_dist]
    relative_dist = (distance - min_dist) / (max_dist - min_dist)
    relative_dist = min(max(relative_dist, 0.0), 1.0)
    weight = ((min_weight - max_weight) * relative_dist) + max_weight
    return round(weight, 2)


def to_rgba_color_string(
    connect_color: tuple[int, int, int],
    alpha: float,
) -> str:
    """Generate a rgba string from list of colors and alpha.

    Parameters
    ----------
    connect_color : tuple[int, int, int]
        A 3-values tuple of integers defining the red, green and blue colors.
    alpha : float
        color_weight

    Returns
    -------
    rgba_color : str
        The html like rgba colors. e.g.: 'rgba(123,123,123,0.5)'
    """
    colors_str = ",".join([str(v) for v in connect_color])
    return f"rgba({colors_str},{alpha})"


def get_pair_color(cont_type: str) -> tuple[int, int, int]:
    """Return the ribbon color of a residue class pair, in any order."""
    pair_key = "-".join(sorted(str(cont_type).split("-")))
    return CONNECT_COLORS.get(pair_key, OTHER_PAIR_COLOR)


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
    matrix : NDArray
        The reconstituted full matrix.
    """
    matrix = squareform(half_matrix)
    np.fill_diagonal(matrix, diag_val)
    return matrix


def make_arc_svgpath(outer: NDArray, inner: NDArray) -> str:
    """Build the closed SVG path between an outer and an inner arc."""
    svgpath = "M "
    for point in outer:
        svgpath += f"{point.real}, {point.imag} L "
    for point in inner[::-1]:
        svgpath += f"{point.real}, {point.imag} L "
    svgpath += f"{outer[0].real}, {outer[0].imag}"
    return svgpath


def make_chordchart(
    _contact_matrix: NDArray,
    _dist_matrix: NDArray,
    _interttype_matrix: NDArray,
    _labels: list[str],
    gap: float = 2 * PI * 0.005,
    output_fpath: Union[str, Path] = "chordchart.html",
    title: str = "Chord diagram",
    offline: bool = False,
    contact_threshold: float = 9.5,
) -> Optional[Union[str, Path]]:
    """Generate a plotly chordchart graph.

    Parameters
    ----------
    _contact_matrix : NDArray
        The contact matrix
    _dist_matrix : NDArray
        The distance matrix
    _interttype_matrix : NDArray
        The residue class pair matrix
    _labels : list[str]
        Labels of each matrix rows (and columns as supposed to be symmetric)
    gap : float, optional
        Gap between two ideograms, by default 2*PI*0.005
    output_fpath : Union[str, Path], optional
        Path to the output file, by default 'chordchart.html'
    title : str, optional
        Title to give to the diagram, by default 'Chord diagram'
    contact_threshold : float, optional
        Distance threshold used to define contacts, sets the range of
        the ribbons transparency.

    Returns
    -------
    output_fpath : Optional[Union[str, Path]]
        Path to the generated output file, None if there are no labels.
    """
    L = check_square_matrix(np.asarray(_contact_matrix))
    if L == 0:
        log.warning("No residue to draw, chord chart not generated")
        return None

    # Group residues by chain, then reverse order so graph displays clockwise
    order = group_indices_by_chain(_labels)[::-1]
    grid = np.ix_(order, order)
    matrix = contacts_to_connect_matrix(_contact_matrix, _labels)[grid]
    dist_matrix = np.asarray(_dist_matrix, dtype=float)[grid]
    interttype_matrix = np.asarray(_interttype_matrix)[grid]
    labels = [_labels[i] for i in order]

    # Map labels into respective chains
    chains = split_labels_by_chains(labels)

    # Compute residues and chain ideograms positions
    ideo_ends, chain_ideo_ends = get_all_ideograms_ends(chains, gap=gap)

    # Compute number of connexion per residues
    row_sum = matrix.sum(axis=1).tolist()

    # Compute connexion ribbons positions
    ribbon_ends = make_ribbon_ends(matrix, row_sum, ideo_ends, L)

    classes = [get_residue_class(parse_reskey(lab)[2]) for lab in labels]
    nicelabels = [
        f"{to_nice_label(lab)} ({resclass})" for lab, resclass in zip(labels, classes)
    ]

    layout_shapes: list[dict] = []
    ribbon_info: list[go.Scatter] = []
    for k in range(L):
        # Half matrix loop to avoid duplicates
        for j in range(k + 1, L):
            if matrix[k, j] == 0 and matrix[j, k] == 0:
                continue

            connect_color = get_pair_color(interttype_matrix[k, j])
            color_weight = to_color_weight(dist_matrix[k, j], contact_threshold)
            rgba_color = to_rgba_color_string(connect_color, color_weight)

            side1 = ribbon_ends[k][j]
            side2 = ribbon_ends[j][k]
            zi = 0.9 * np.exp(1j * (side1[0] + side1[1]) / 2)
            zf = 0.9 * np.exp(1j * (side2[0] + side2[1]) / 2)

            # Strings displayed when hovering the two ribbon ends
            dist_text = f"Distance: {dist_matrix[k, j]:.1f} &#8491;"
            texti = f"{nicelabels[k]} &#8621; {nicelabels[j]}<br>{dist_text}"
            textf = f"{nicelabels[j]} &#8621; {nicelabels[k]}<br>{dist_text}"
            for zv, text in zip([zi, zf], [texti, textf]):
                ribbon_info.append(
                    go.Scatter(
                        x=[zv.real],
                        y=[zv.imag],
                        mode="markers",
                        marker={"size": 0.5, "color": rgba_color},
                        text=text,
                        hoverinfo="text",
                        showlegend=False,
                    )
                )
            # Note: must reverse these arc ends to avoid twisted ribbon
            side2_rev = (side2[1], side2[0])
            layout_shapes.append(
                make_ribbon(
                    side1,
                    side2_rev,
                    "rgba(175,175,175)",
                    rgba_color,
                )
            )

    ideograms: list[go.Scatter] = []
    # Draw ideograms for residues, colored by residue class
    for k in range(L):
        z = make_ideogram_arc(1.1, ideo_ends[k])
        zi = make_ideogram_arc(1.0, ideo_ends[k])
        rescolor = to_rgba_color_string(RESIDUE_CLASS_COLORS[classes[k]], 0.8)

        text_info = f"{nicelabels[k]}<br>"
        if row_sum[k] == 0:
            text_info += "No contact"
        else:
            text_info += f"Total of {row_sum[k]:d} contact"
            if row_sum[k] >= 2:
                text_info += "s"
        ideograms.append(
            go.Scatter(
                x=z.real,
                y=z.imag,
                mode="lines",
                line={
                    "color": rescolor,
                    "shape": "spline",
                    "width": 0.25,
                },
                text=text_info,
                hoverinfo="text",
                showlegend=False,
            )
        )
        layout_shapes.append(
            make_ideo_shape(
                make_arc_svgpath(z, zi),
                "rgba(150,150,150)",
                rescolor,
            )
        )

    # Draw ideograms for chains
    for k, chainid in enumerate(chains):
        chain_color = CHAIN_COLORS[k % len(CHAIN_COLORS)]
        z = make_ideogram_arc(1.2, chain_ideo_ends[k])
        zi = make_ideogram_arc(1.11, chain_ideo_ends[k])
        ideograms.append(
            go.Scatter(
                x=z.real,
                y=z.imag,
                mode="lines",
                line={
                    "color": chain_color,
                    "shape": "spline",
                    "width": 0.25,
                },
                text=f"Chain {chainid}",
                hoverinfo="text",
                showlegend=False,
            )
        )
        layout_shapes.append(
            make_ideo_shape(
                make_arc_svgpath(z, zi),
                "rgba(150,150,150)",
                chain_color,
            )
        )

    fig_size = max(MIN_CHORDCHART_SIZE, 100 * np.log(L * L))
    layout = make_layout(title, fig_size, layout_shapes)
    fig = go.Figure(data=ideograms + ribbon_info, layout=layout)
    fig.update_layout(plot_bgcolor="white")
    add_chordchart_legends(fig)
    fig_to_html(
        fig,
        output_fpath,
        figure_height=fig_size,
        figure_width=fig_size + CHORDCHART_LEGEND_WIDTH,
        offline=offline,
    )
    return output_fpath


def add_legend_entry(
    fig: go.Figure,
    name: str,
    color: tuple[int, int, int],
    group: str,
    group_title: str,
) -> None:
    """Add a dummy trace to the figure to display a legend entry."""
    fig.add_trace(
        go.Scatter(
            x=[None],
            y=[None],
            legendgroup=group,
            legendgrouptitle_text=group_title,
            showlegend=True,
            name=name,
            mode="lines",
            line={"color": to_rgba_color_string(color, 0.9), "width": 6},
        )
    )


def add_chordchart_legends(fig: go.Figure) -> None:
    """Add custom legend to chordchart.

    Parameters
    ----------
    fig : go.Figure
        A plotly figure.
    """
    pair_title = "Residue class pair"
    for pair_key, color in CONNECT_COLORS.items():
        add_legend_entry(
            fig,
            pair_key.replace("-", " &#8621; "),
            color,
            "connect_color",
            pair_title,
        )
    add_legend_entry(fig, "other pairs", OTHER_PAIR_COLOR, "connect_color", pair_title)

    for resclass, color in RESIDUE_CLASS_COLORS.items():
        add_legend_entry(fig, resclass, color, "class_color", "Residue class")


def tsv_to_chordchart(
    tsv_path: Union[Path, str],
    sep: str = "\t",
    data_key: str = "ca-ca-dist",
    contact_threshold: float = 7.5,
    filter_intermolecular_contacts: bool = True,
    output_fname: Union[Path, str] = "contacts_chordchart.html",
    title: str = "Chord diagram",
    offline: bool = False,
) -> Optional[Union[Path, str]]:
    """Read a tsv file and generate a chord diagram from it.

    Parameters
    ----------
    tsv_path : Union[Path, str]
        Path a the .tsv file containing contact data.
    sep : str
        Separator character used to split data in each line.
    data_key : str
        Data key used to draw the plot.
    contact_threshold : float
        Values <= to this threshold are considered as contacts.
    filter_intermolecular_contacts : bool
        Only draw residues involved in interchain contacts.
    output_fname : Union[Path, str]
        Path where to generate the graph.
    title : str
        Title to give to the Chord diagram

    Return
    ------
    chord_chart_fpath : Optional[Union[Path, str]]
        Path to the generated graph, None if there is no contact to draw.
    """
    header, rows = read_contacts_tsv(tsv_path, sep=sep)
    half_contact_matrix: list[int] = []
    half_value_matrix: list[float] = []
    half_intertype_matrix: list[str] = []
    labels: list[str] = []
    seen: set[str] = set()
    idx1, idx2, idxv, idxt = (
        header.index(k) for k in ("res1", "res2", data_key, "contact-type")
    )
    for s_ in rows:
        for label in (s_[idx1], s_[idx2]):
            if label not in seen:
                seen.add(label)
                labels.append(label)
        value = float(s_[idxv])
        half_contact_matrix.append(1 if value <= contact_threshold else 0)
        half_value_matrix.append(value)
        half_intertype_matrix.append(s_[idxt])

    if len(labels) < 2:
        log.warning(f"Not enough residues in {tsv_path} to generate a chord chart")
        return None

    contact_matrix = to_full_matrix(half_contact_matrix, 1)
    dist_matrix = to_full_matrix(half_value_matrix, 0.0)
    intertype_matrix = to_full_matrix(half_intertype_matrix, "self-self")

    if filter_intermolecular_contacts:
        # Keep residues involved in at least one interchain contact
        connect_matrix = contacts_to_connect_matrix(contact_matrix, labels)
        sorted_indices = np.flatnonzero(connect_matrix.any(axis=1)).tolist()
        if not sorted_indices:
            log.info(
                f"No interchain contact under threshold in {tsv_path}, "
                "chord chart not generated"
            )
            return None
        contact_matrix = extract_submatrix(contact_matrix, sorted_indices)
        dist_matrix = extract_submatrix(dist_matrix, sorted_indices)
        intertype_matrix = extract_submatrix(intertype_matrix, sorted_indices)
        labels = [labels[i] for i in sorted_indices]

    return make_chordchart(
        contact_matrix,
        dist_matrix,
        intertype_matrix,
        labels,
        output_fpath=output_fname,
        title=title,
        offline=offline,
        contact_threshold=contact_threshold,
    )
