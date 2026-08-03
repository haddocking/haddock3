"""dnascan module.

Core logic for the DNA base-pair scan.

Unlike the RNA scan, where each interface nucleotide is mutated independently,
DNA is double stranded and its bases are engaged in Watson-Crick base pairs
(A:T and G:C). Mutating a single base in isolation would break the base pair
and is not physically meaningful. Therefore every mutation performed by this
module is a *double* mutation: the selected interface nucleotide is mutated to a
target base and its base-pairing partner on the complementary strand is
mutated to the complementary base so that a valid Watson-Crick pair is
preserved (e.g. A:T -> G:C).

Same-ring-type base-pair mutations are applied simultaneously and scored in a
single CNS call. Cross-ring-type base-pair mutations (where the two nucleotides
swap ring type, one purine -> pyrimidine and its partner pyrimidine -> purine)
are instead performed in two sequential, CNS-regularised steps: first the
purine -> pyrimidine mutation is applied and energy-minimised by CNS, then the
pyrimidine -> purine mutation is applied to that minimised intermediate and
scored by CNS. The score and energies of the second (final) CNS call are the
ones reported and plotted.

To keep the comparison consistent, each mutant is compared against a wild-type
baseline that went through the same number of CNS minimisation passes: a
one-pass baseline for same-ring-type mutants and a two-pass baseline (the wild
type minimised and then re-scored) for cross-ring-type mutants.
"""

import os
import shutil
from pathlib import Path
from typing import List, Tuple, Dict

import numpy as np


from haddock import log
from haddock.core.exceptions import ConfigurationError
from haddock.core.typing import Any, Optional, Union
from haddock.libs.libalign import get_atoms, load_coords
from haddock.libs.libontology import PDBFile
from haddock.libs.libscan import (
    MutationResult,
    calc_score,
    AddDeltaBFactor as _AddDeltaBFactor,
    ClusterOutputer as _ClusterOutputer,
    write_scan_out as _write_scan_out,
    group_scan_by_cluster as _group_scan_by_cluster,
)
from haddock.libs.libcapri import CAPRI

# Heavy atoms of the deoxyribose-phosphate backbone that are common to all DNA
# nucleotides. These atoms are always preserved when mutating a base. Unlike
# RNA, DNA has a 2'-deoxyribose, so there is no O2' backbone atom. Both the
# modern (OP1/OP2) and legacy (O1P/O2P) phosphate oxygen namings are accepted
# so that models from either convention are handled correctly.
BACKBONE_ATOMS = [
    "P",
    "OP1",
    "OP2",
    "O1P",
    "O2P",
    "O5'",
    "C5'",
    "C4'",
    "O4'",
    "C1'",
    "C2'",
    "C3'",
    "O3'",
]

# Base ring atoms shared by both purines (DA and DG). Preserving them when
# mutating one purine into the other keeps the base plane and glycosidic
# orientation, so CNS only has to rebuild the differing substituents.
PURINE_BASE_ATOMS = ["N9", "C8", "N7", "C5", "C4", "N3", "C2", "N1", "C6"]

# Base ring atoms shared by both pyrimidines (DC and DT). Same rationale as for
# the purines above.
PYRIMIDINE_BASE_ATOMS = ["N1", "C2", "O2", "N3", "C4", "C5", "C6"]

# DNA residues that can be scanned/mutated by this module.
DNA_RESIDUES = ("DA", "DC", "DG", "DT")

# Ring-type classification of the DNA bases.
PURINES = ("DA", "DG")
PYRIMIDINES = ("DC", "DT")

# Glycosidic-region anchor atoms kept (and renamed) on cross-type mutations
# (purine <-> pyrimidine). Because the purine and pyrimidine ring systems do
# not share atom names, these three atoms are preserved and renamed to their
# counterpart in the target ring so that CNS rebuilds the new base with the
# correct glycosidic orientation (and syn/anti conformation). The
# correspondence is:
#   pyrimidine N1 <-> purine N9   (glycosidic nitrogen)
#   pyrimidine C2 <-> purine C4
#   pyrimidine C6 <-> purine C8
PYRIMIDINE_TO_PURINE_ANCHORS = {"N1": "N9", "C2": "C4", "C6": "C8"}
PURINE_TO_PYRIMIDINE_ANCHORS = {"N9": "N1", "C4": "C2", "C8": "C6"}

# Watson-Crick base-pair complementarity. Mutating a base to a target
# implies mutating its partner to the complementary base to keep a valid pair.
COMPLEMENT = {"DA": "DT", "DT": "DA", "DG": "DC", "DC": "DG"}

# Mapping of the two-letter DNA residue name to the short one-letter code used
# when building mutation identifiers and output file names.
RES_CODES = dict(
    [
        # DNA (two-letter name -> short one-letter code used in identifiers)
        ("DA", "A"),
        ("DC", "C"),
        ("DG", "G"),
        ("DT", "T"),
    ]
)

# Distance cutoff (Å) used when detecting Watson-Crick base pairs from the
# distance between the hydrogen-bonding ring nitrogens (N1 of the purine and
# N3 of the pyrimidine). A canonical WC pair has an N1...N3 distance of ~2.8 Å;
# a small margin is allowed to tolerate energy-minimised / slightly distorted
# geometries. Overridden at runtime by the ``bp_cutoff`` value from the module's
# defaults.yaml.
BP_CUTOFF = 3.5

# Default distance cutoff (Å) used to define interface contacts between two
# interacting molecules. Overridden at runtime by the ``int_cutoff`` value from
# the module's defaults.yaml.
INT_CUTOFF = 5.0


def wc_atom(resname: str) -> str:
    """Return the Watson-Crick hydrogen-bonding ring nitrogen for a base.

    Purines (DA, DG) pair through their N1, pyrimidines (DC, DT) through their
    N3. The distance between these two atoms is used to identify base pairs.

    Parameters
    ----------
    resname : str
        DNA residue name.

    Returns
    -------
    str
        Atom name of the Watson-Crick ring nitrogen (``N1`` or ``N3``).
    """
    return "N1" if resname in PURINES else "N3"


def get_atoms_to_keep(ori_resname: str, target_resname: str) -> Dict[str, str]:
    """Return the atoms to preserve when mutating one base into another.

    The result maps each atom name to keep (as found in the original residue)
    to the atom name it must be written with in the mutated residue. Backbone
    atoms and same-ring-type base atoms keep their name (mapped to themselves).

    The deoxyribose-phosphate backbone is always kept. In addition, when the
    original and target bases are of the same ring type (both purines or both
    pyrimidines), the base ring atoms common to that type are kept as well, so
    that the base orientation is preserved and CNS only rebuilds the differing
    substituents. For cross-type mutations (purine <-> pyrimidine) the three
    glycosidic-region anchor atoms are kept and renamed to their counterpart in
    the target ring system (pyrimidine N1/C2/C6 <-> purine N9/C4/C8), so that the
    base orientation is preserved and CNS rebuilds the rest of the new base.

    Parameters
    ----------
    ori_resname : str
        Original (wild-type) residue name.
    target_resname : str
        Target base residue name.

    Returns
    -------
    dict of str -> str
        Mapping of original atom name -> atom name to write for the mutated
        nucleotide.
    """
    # Backbone is always kept, with unchanged atom names.
    atoms_to_keep: Dict[str, str] = {atom: atom for atom in BACKBONE_ATOMS}
    if ori_resname in PURINES and target_resname in PURINES:
        atoms_to_keep.update({atom: atom for atom in PURINE_BASE_ATOMS})
    elif ori_resname in PYRIMIDINES and target_resname in PYRIMIDINES:
        atoms_to_keep.update({atom: atom for atom in PYRIMIDINE_BASE_ATOMS})
    elif ori_resname in PYRIMIDINES and target_resname in PURINES:
        atoms_to_keep.update(PYRIMIDINE_TO_PURINE_ANCHORS)
    elif ori_resname in PURINES and target_resname in PYRIMIDINES:
        atoms_to_keep.update(PURINE_TO_PYRIMIDINE_ANCHORS)
    return atoms_to_keep


# Default set of target bases tested for each selected interface nucleotide.
DEFAULT_SCAN_BASES = ["DA", "DC", "DG", "DT"]


def validate_scan_bases(scan_bases: List[str]) -> List[str]:
    """Validate and normalise the list of target DNA bases.

    Only the canonical two-letter DNA residue names (``DA``, ``DC``, ``DG``,
    ``DT``) are accepted, in a case-insensitive manner. The one-letter names
    (``A``, ``C``, ``G``, ``U``) denote RNA bases in HADDOCK and are therefore
    rejected by this DNA-specific module.

    Parameters
    ----------
    scan_bases : list of str
        Target bases requested by the user.

    Returns
    -------
    list of str
        Normalised, de-duplicated list of canonical DNA residue names, in the
        order first seen.

    Raises
    ------
    ConfigurationError
        If ``scan_bases`` is empty or contains an unrecognised base.
    """
    if not scan_bases:
        raise ConfigurationError(
            "'scan_bases' must contain at least one DNA base "
            f"among {', '.join(DNA_RESIDUES)}."
        )
    normalised: List[str] = []
    for base in scan_bases:
        canonical = str(base).strip().upper()
        if canonical not in DNA_RESIDUES:
            raise ConfigurationError(
                f"Invalid 'scan_bases' entry {base!r}. Allowed values are the "
                f"two-letter DNA residue names {', '.join(DNA_RESIDUES)} "
                "(one-letter names such as A, C, G, U denote RNA bases and are "
                "not supported)."
            )
        if canonical not in normalised:
            normalised.append(canonical)
    return normalised


def _norm_atom_name(atom_name: str) -> str:
    """Normalise atom names so that primes are always written as `'`."""
    return atom_name.replace("*", "'")


def is_cross_type(ori_resname: str, target_resname: str) -> bool:
    """Return True if a mutation changes the base ring type.

    A cross-type mutation converts a purine into a pyrimidine or vice versa.

    Parameters
    ----------
    ori_resname : str
        Original (wild-type) residue name.
    target_resname : str
        Target base residue name.
    """
    return (ori_resname in PURINES) != (target_resname in PURINES)


def _mutate_residues(pdb_f, mutations: Dict[Tuple[str, int], str]) -> Path:
    """Apply one or more nucleotide mutations to a PDB file.

    For each targeted nucleotide the deoxyribose-phosphate backbone is always
    kept; when the original and target bases share a ring type the common base
    ring atoms are kept as well to preserve the base orientation (see
    ``get_atoms_to_keep``). For cross-type mutations (purine <-> pyrimidine) the
    three glycosidic-region anchor atoms are kept and renamed to their
    counterpart in the target ring (pyrimidine N1/C2/C6 <-> purine N9/C4/C8). The
    remaining base atoms are dropped and rebuilt by CNS during scoring.

    Parameters
    ----------
    pdb_f : Path
        Path to the pdb file.
    mutations : dict
        Mapping of ``(chain, resid)`` -> target residue name. The output file
        name and mutation identifier follow the insertion order of this dict.

    Returns
    -------
    mut_pdb_fname : Path
        Path to the mutated pdb file.
    """
    # Map (chain, resid) -> (target base name, discovered original name)
    targets = {key: {"mut": mut, "ori": ""} for key, mut in mutations.items()}
    mut_pdb_l = []
    with open(pdb_f, "r") as fh:
        for line in fh.readlines():
            if line.startswith("ATOM"):
                chain = line[21]
                resid = int(line[22:26])
                atom_name = _norm_atom_name(line[12:16].strip())
                key = (chain, resid)
                if key in targets:
                    entry = targets[key]
                    if not entry["ori"]:
                        entry["ori"] = line[17:20].strip()
                    atoms_to_keep = get_atoms_to_keep(entry["ori"], entry["mut"])
                    if atom_name in atoms_to_keep:
                        # DNA residue names are shorter than 3 characters, so
                        # they must be right-justified to keep the PDB columns
                        # (18-20) aligned.
                        resname_field = entry["mut"].rjust(3)
                        line = line[:17] + resname_field + line[20:]
                        # rename the atom if it maps to a different name in the
                        # target ring system (cross-type anchor atoms)
                        new_atom_name = atoms_to_keep[atom_name]
                        if new_atom_name != atom_name:
                            new_field = line[12:16].replace(atom_name, new_atom_name, 1)
                            line = line[:12] + new_field + line[16:]
                        mut_pdb_l.append(line)
                else:
                    mut_pdb_l.append(line)

    # Build the mutation identifier, one "<ori><resid><target>" token per
    # mutated nucleotide, in the order the mutations were provided.
    id_tokens = []
    for (chain, resid), entry in targets.items():
        ori = entry["ori"]
        try:
            id_tokens.append(f"{RES_CODES[ori]}{resid}{RES_CODES[entry['mut']]}")
        except KeyError:
            raise KeyError(f"Could not mutate {ori} into {entry['mut']}.")
    mut_id = "-".join(id_tokens)
    first_chain = next(iter(targets))[0]
    mut_pdb_fname = Path(pdb_f.name.replace(".pdb", f"-{first_chain}_{mut_id}.pdb"))
    with open(mut_pdb_fname, "w") as fh:
        fh.write("".join(mut_pdb_l))
    return mut_pdb_fname


def mutate(
    pdb_f,
    target_chain,
    target_resid,
    mut_resname,
    partner_chain,
    partner_resid,
    partner_mut_resname,
):
    """
    Perform a Watson-Crick double mutation of a base pair in a PDB file.

    Both the selected nucleotide and its base-pairing partner are mutated at
    once so that a valid Watson-Crick pair is preserved. See
    ``_mutate_residues`` for the details of how atoms are kept and renamed.

    Parameters
    ----------
    pdb_f : str
        Path to the pdb file.
    target_chain : str
        Chain of the primary nucleotide to be mutated.
    target_resid : int
        Residue number of the primary nucleotide to be mutated.
    mut_resname : str
        Residue name of the target base for the primary nucleotide
        (e.g. ``DA``, ``DC``, ``DG``, ``DT``).
    partner_chain : str
        Chain of the base-pairing partner nucleotide.
    partner_resid : int
        Residue number of the base-pairing partner nucleotide.
    partner_mut_resname : str
        Residue name of the target base for the partner nucleotide (the
        Watson-Crick complement of ``mut_resname``).

    Returns
    -------
    mut_pdb_fname : Path
        Path to the mutated pdb file.
    """
    return _mutate_residues(
        pdb_f,
        {
            (target_chain, target_resid): mut_resname,
            (partner_chain, partner_resid): partner_mut_resname,
        },
    )


def find_base_pairs(
    coords: Dict[Tuple[str, int, str, str], Any],
    cutoff: float = BP_CUTOFF,
) -> Dict[Tuple[str, int], Tuple[str, int, str]]:
    """Detect Watson-Crick base pairs from atomic coordinates.

    For every DNA nucleotide, the Watson-Crick hydrogen-bonding ring nitrogen
    (N1 for purines, N3 for pyrimidines) is located and the closest nucleotide
    of the complementary ring type (purine <-> pyrimidine) whose corresponding
    nitrogen lies within ``cutoff`` is taken as its base-pairing partner. Base
    pairing is not restricted to a single chain: partners on the same strand
    (e.g. hairpin duplexes) or on a different chain are both detected.

    Parameters
    ----------
    coords : dict
        Coordinate dictionary as returned by
        :func:`haddock.libs.libalign.load_coords` with ``add_resname=True``;
        keys are ``(chain, resid, atom_name, resname)`` tuples.
    cutoff : float, optional
        Maximum N1...N3 distance (Å) for two bases to be considered paired.

    Returns
    -------
    dict
        Mapping ``(chain, resid) -> (partner_chain, partner_resid,
        partner_resname)`` for every nucleotide for which a partner was found.
    """
    # Gather the Watson-Crick nitrogen coordinate of every DNA nucleotide.
    nucleotides: Dict[Tuple[str, int], Tuple[str, Any]] = {}
    for (chain, resid, atom, resname), xyz in coords.items():
        if resname not in DNA_RESIDUES:
            continue
        if atom == wc_atom(resname):
            nucleotides[(chain, resid)] = (resname, np.asarray(xyz))

    pairs: Dict[Tuple[str, int], Tuple[str, int, str]] = {}
    keys = list(nucleotides.keys())
    for key_i in keys:
        resname_i, xyz_i = nucleotides[key_i]
        best_key = None
        best_dist = cutoff
        for key_j in keys:
            if key_j == key_i:
                continue
            resname_j, xyz_j = nucleotides[key_j]
            # Partners must be of complementary ring type (purine <-> pyrimidine)
            if (resname_i in PURINES) == (resname_j in PURINES):
                continue
            dist = float(np.linalg.norm(xyz_i - xyz_j))
            if dist < best_dist:
                best_dist = dist
                best_key = key_j
        if best_key is not None:
            pairs[key_i] = (best_key[0], best_key[1], nucleotides[best_key][0])
    return pairs


class ClusterOutputer(_ClusterOutputer):
    """Manage the generation of dnascan outputs for cluster-based analysis."""

    module_name = "dnascan"
    default_scan_residue = "DNA base pair"
    sort_columns = ["chain", "resid", "target_resname"]
    zscore_reference = "mutations"

    def _identity_columns(self):
        return [
            "chain",
            "resid",
            "resname",
            "target_resname",
            "partner_chain",
            "partner_resid",
            "partner_resname",
            "partner_target_resname",
            "full_resname",
        ]

    def _extra_header_lines(self):
        return ["each row is a Watson-Crick double mutation (base pair)"]

    def _identity_row(self, ident, clt_res_dt):
        # Build a compact human-readable label for the base-pair mutation
        full_resname = (
            f"{clt_res_dt['chain']}{clt_res_dt['resid']}"
            f"{clt_res_dt['ori_resname']}>{clt_res_dt['target_resname']}"
            f"/{clt_res_dt['partner_chain']}{clt_res_dt['partner_resid']}"
            f"{clt_res_dt['partner_ori_resname']}>"
            f"{clt_res_dt['partner_target_resname']}"
        )
        return [
            clt_res_dt["chain"],
            clt_res_dt["resid"],
            clt_res_dt["ori_resname"],
            clt_res_dt["target_resname"],
            clt_res_dt["partner_chain"],
            clt_res_dt["partner_resid"],
            clt_res_dt["partner_ori_resname"],
            clt_res_dt["partner_target_resname"],
            full_resname,
        ]


class AddDeltaBFactor(_AddDeltaBFactor):
    """Add dnascan delta score in the b-factor column of a PDB.

    The delta score of a base-pair mutation is attributed to both nucleotides
    of the pair so that both light up when colouring by b-factor.
    """

    module_name = "dnascan"

    def _residue_keys(self, mut_result):
        return [
            (mut_result.chain, mut_result.resid),
            (mut_result.partner_chain, mut_result.partner_resid),
        ]


_DNA_SCAN_COLUMNS = [
    "chain",
    "res",
    "ori_resname",
    "end_resname",
    "partner_chain",
    "partner_res",
    "partner_ori_resname",
    "partner_end_resname",
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


def _dna_scan_row(result):
    """Build a per-model TSV row for a base-pair (double) mutation result."""
    m_score, m_vdw, m_elec, m_des, m_bsa = result.mutant_scores
    d_score, d_vdw, d_elec, d_des, d_bsa = result.delta_scores
    return [
        result.chain,
        result.resid,
        result.ori_resname,
        result.target_resname,
        result.partner_chain,
        result.partner_resid,
        result.partner_ori_resname,
        result.partner_target_resname,
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


def _dna_cluster_metadata(result):
    """Static per-mutation metadata stored for each base-pair identifier."""
    return {
        "chain": result.chain,
        "resid": result.resid,
        "ori_resname": result.ori_resname,
        "target_resname": result.target_resname,
        "partner_chain": result.partner_chain,
        "partner_resid": result.partner_resid,
        "partner_ori_resname": result.partner_ori_resname,
        "partner_target_resname": result.partner_target_resname,
    }


def group_scan_by_cluster(models, results_by_model):
    """Group dnascan data per cluster, keyed by base-pair mutation."""
    return _group_scan_by_cluster(
        models,
        results_by_model,
        ident_builder=lambda r: (
            f"{r.chain}-{r.resid}-{r.ori_resname}-{r.target_resname}"
        ),
        metadata_builder=_dna_cluster_metadata,
    )


def write_scan_out(results, model_id):
    """Save dnascan base-pair mutation results for one model to a tsv file."""
    _write_scan_out(
        results,
        model_id,
        module_name="dnascan",
        sort_columns=["chain", "res", "end_resname"],
        row_builder=_dna_scan_row,
        columns=_DNA_SCAN_COLUMNS,
        zscore_reference="mutations",
        extra_header_lines=["each row is a Watson-Crick double mutation (base pair)"],
    )


class InterfaceScanner:
    """Scan interface of a model to get target base pairs and create
    corresponding double-mutation jobs.
    """

    def __init__(
        self,
        model: Union[str, Path, Any],
        scan_bases: Optional[List[str]] = None,
        params: Optional[Dict[str, Any]] = None,
    ) -> None:
        """
        Initialize InterfaceScanner for a single model.

        Parameters
        ----------
        model : str, Path, or model object
            HADDOCK ``PDBFile`` model object or a path to a PDB file.
        scan_bases : list of str, optional
            Target bases tested for each interface nucleotide
            (default: DA, DC, DG, DT)
        params : dict, optional
            Additional parameters for interface detection
            (list on top of dnascan/__init__.py)
        """
        self.model = model
        self.scan_bases = list(scan_bases) if scan_bases else list(DEFAULT_SCAN_BASES)
        self.params = params or {}
        self.point_mutations_jobs = []
        self.ligand_param_fname = self.params.get("ligand_param_fname", "")
        self.ligand_top_fname = self.params.get("ligand_top_fname", "")
        self.filter_resdic = {
            key[-1]: value
            for key, value in self.params.items()
            if key.startswith("resdic")
        }
        if isinstance(model, PDBFile):
            self.model_path = model.rel_path
            self.model_id = model.file_name.removesuffix(".pdb")
        else:
            # model given as a plain path
            self.model_path = Path(model)
            self.model_id = self.model_path.stem

    def _compute_native_baselines(
        self,
    ) -> Tuple[
        Tuple[float, float, float, float, float],
        Tuple[float, float, float, float, float],
    ]:
        """Return the one-pass and two-pass wild-type score baselines.

        The one-pass baseline is the score of the wild-type model (matching the
        single CNS call used for same-ring-type mutants). The two-pass baseline
        scores the CNS energy-minimised wild type a second time, so that it
        matches the two sequential CNS steps used to score cross-ring-type
        base-pair mutations. Both are computed once per model.

        Returns
        -------
        tuple
            ``(native_scores, native_scores_2step)``.
        """
        sc_dir_1 = f"haddock3-score-{self.model_id}-{os.getpid()}-nat1"
        sc_dir_2 = f"haddock3-score-{self.model_id}-{os.getpid()}-nat2"
        min_wt = Path(f"{Path(self.model_path).stem}_hs.pdb")
        try:
            # one CNS pass: score the WT and write its energy-minimised structure
            native_scores = calc_score(
                self.model_path,
                run_dir=sc_dir_1,
                outputpdb=True,
                ligand_param_fname=self.ligand_param_fname,
                ligand_top_fname=self.ligand_top_fname,
            )
            # two CNS passes: score the energy-minimised WT a second time
            if min_wt.exists():
                native_scores_2step = calc_score(
                    min_wt,
                    run_dir=sc_dir_2,
                    outputpdb=False,
                    ligand_param_fname=self.ligand_param_fname,
                    ligand_top_fname=self.ligand_top_fname,
                )
            else:
                # the minimised structure was not produced (e.g. calc_score is
                # mocked in tests); fall back to the one-pass baseline
                native_scores_2step = native_scores
        finally:
            if min_wt.exists():
                os.remove(min_wt)
            for d in (sc_dir_1, sc_dir_2):
                if os.path.exists(d):
                    shutil.rmtree(d)
        return native_scores, native_scores_2step

    @staticmethod
    def _iter_mutations(interface, resname_dict, base_pairs, scan_bases):
        """Yield each valid Watson-Crick base-pair (double) mutation to perform.

        Flattens the interface (chain -> residues) and the requested target
        bases into a single stream of fully-resolved mutation tuples. Skips
        non-DNA residues, nucleotides without a Watson-Crick partner, no-op
        mutations (target base equal to the original) and deduplicates each base
        pair so it is scanned once (not once per nucleotide).

        Parameters
        ----------
        interface : dict
            Mapping of chain id to the list of interface residue numbers.
        resname_dict : dict
            Mapping of ``"{chain}-{resid}"`` to the residue name.
        base_pairs : dict
            Mapping of ``(chain, resid)`` to its Watson-Crick partner tuple
            ``(partner_chain, partner_resid, partner_ori_resname)``.
        scan_bases : list of str
            Target bases to scan each nucleotide into.

        Yields
        ------
        tuple
            ``(chain, resid, ori_resname, target_resname, partner_chain,
            partner_resid, partner_ori_resname, partner_target_resname)`` for
            each base-pair mutation.
        """
        # Keep track of base pairs already scheduled so that a pair is not
        # scanned twice (once from each of its two nucleotides).
        scheduled_pairs = set()
        for chain, residues in interface.items():
            for res in residues:
                ori_resname = resname_dict[f"{chain}-{res}"]
                # Only scan DNA nucleotides, skip protein/other residues
                if ori_resname not in DNA_RESIDUES:
                    continue
                # DNA mutations must be done as base pairs: locate partner
                partner = base_pairs.get((chain, res))
                if partner is None:
                    log.warning(
                        f"No Watson-Crick partner found for {chain}:{res} "
                        f"{ori_resname}; skipping it (dnascan only mutates "
                        "base pairs)."
                    )
                    continue
                partner_chain, partner_resid, partner_ori = partner
                # Deduplicate on the unordered base pair
                pair_key = frozenset({(chain, res), (partner_chain, partner_resid)})
                if pair_key in scheduled_pairs:
                    continue
                scheduled_pairs.add(pair_key)
                for end_resname in scan_bases:
                    # Skip no-op mutation (e.g. DA -> DA)
                    if ori_resname == end_resname:
                        continue
                    # Partner is mutated to the Watson-Crick complement to keep
                    # a valid base pair
                    partner_target = COMPLEMENT[end_resname]
                    yield (
                        chain,
                        res,
                        ori_resname,
                        end_resname,
                        partner_chain,
                        partner_resid,
                        partner_ori,
                        partner_target,
                    )

    def run(self):
        """
        Get interface base pairs and create the double-mutation jobs.

        The jobs are returned (not executed): the caller hands them to a haddock
        Engine so that all mutations are scheduled together.

        Returns
        -------
        List[ModelBasePairMutation]
            The base-pair mutation jobs to perform for this model.
        """
        try:
            # Calculate the wild-type score baselines. Two baselines are used
            # so that each mutant is compared against a wild type that has gone
            # through the same number of CNS minimisation passes:
            #   * native_scores       : one CNS pass  (same-ring-type mutants)
            #   * native_scores_2step : two CNS passes (cross-ring-type mutants,
            #     which are themselves scored in two sequential CNS steps)
            native_scores, native_scores_2step = self._compute_native_baselines()

            # Load coordinates
            atoms = get_atoms(self.model_path)
            coords, _chain_ranges = load_coords(
                self.model_path,
                atoms,
                add_resname=True,
            )

            # Detect Watson-Crick base pairs across the whole model
            bp_cutoff = self.params.get("bp_cutoff", BP_CUTOFF)
            base_pairs = find_base_pairs(coords, cutoff=bp_cutoff)

            # Determine target nucleotides: get interface, then apply user filters
            # Get all interface residues
            cutoff = self.params.get("int_cutoff", INT_CUTOFF)
            interface = CAPRI.identify_interface(self.model_path, cutoff=cutoff)

            # get user_chains for the check down the line
            user_chains = self.params.get("chains", [])

            # if user defined target residues, check they are in the interface
            if self.filter_resdic != {"_": []}:
                filtered_interface = {}
                for chain in self.filter_resdic:
                    if chain in interface:
                        # Search for the intersection of user queried residues and interface residues
                        user_res_valid = list(
                            set(self.filter_resdic[chain]).intersection(
                                set(interface[chain])
                            )
                        )
                        # If at least one residue must be analyzed, add it to residues to be scanned
                        if user_res_valid:
                            filtered_interface[chain] = user_res_valid
                interface = filtered_interface

            # if (user defined target chains) & (no user target residues) - do use user chains
            elif user_chains:
                interface = {
                    chain: res
                    for chain, res in interface.items()
                    if chain in user_chains
                }

            # get all atoms of the model to verify residue type down the line
            resname_dict = {}
            for chain, resid, _atom, resname in coords.keys():
                key = f"{chain}-{resid}"
                if key not in resname_dict:
                    resname_dict[key] = resname

            # Create mutations
            output_mutants = self.params.get("output_mutants", False)
            for mut in self._iter_mutations(
                interface, resname_dict, base_pairs, self.scan_bases
            ):
                (
                    chain,
                    res,
                    ori_resname,
                    end_resname,
                    partner_chain,
                    partner_resid,
                    partner_ori,
                    partner_target,
                ) = mut
                job = ModelBasePairMutation(
                    model_path=self.model_path,
                    model_id=self.model_id,
                    chain=chain,
                    resid=res,
                    ori_resname=ori_resname,
                    target_resname=end_resname,
                    partner_chain=partner_chain,
                    partner_resid=partner_resid,
                    partner_ori_resname=partner_ori,
                    partner_target_resname=partner_target,
                    native_scores=native_scores,
                    native_scores_2step=native_scores_2step,
                    output_mutants=output_mutants,
                    ligand_param_fname=self.ligand_param_fname,
                    ligand_top_fname=self.ligand_top_fname,
                )
                self.point_mutations_jobs.append(job)

            return self.point_mutations_jobs

        except Exception as e:
            log.error(f"Failed to scan model {self.model_id}: {e}")
            raise


class ModelBasePairMutation:
    """Executes a single Watson-Crick base-pair (double) mutation."""

    def __init__(
        self,
        model_path: Path,
        model_id: str,
        chain: str,
        resid: int,
        ori_resname: str,
        target_resname: str,
        partner_chain: str,
        partner_resid: int,
        partner_ori_resname: str,
        partner_target_resname: str,
        native_scores: Tuple[float, float, float, float, float],
        native_scores_2step: Optional[Tuple[float, float, float, float, float]] = None,
        output_mutants: bool = False,
        ligand_param_fname: Union[Path, str] = "",
        ligand_top_fname: Union[Path, str] = "",
    ) -> None:
        """
        Initialize a single base-pair mutation job.

        Parameters
        ----------
        model_path : Path
            Path to the PDB file
        model_id : str
            Identifier for the model
        chain : str
            Chain identifier of the primary nucleotide
        resid : int
            Residue number of the primary nucleotide
        ori_resname : str
            Original residue name of the primary nucleotide
        target_resname : str
            Target base name for the primary nucleotide
        partner_chain : str
            Chain identifier of the base-pairing partner
        partner_resid : int
            Residue number of the base-pairing partner
        partner_ori_resname : str
            Original residue name of the base-pairing partner
        partner_target_resname : str
            Target base name for the base-pairing partner (WC complement)
        native_scores : tuple
            Wild-type model scores from a single CNS pass (score, vdw, elec,
            desolv, bsa); used as the baseline for same-ring-type mutants.
        native_scores_2step : tuple, optional
            Wild-type model scores from two CNS passes; used as the baseline for
            cross-ring-type mutants (which are scored in two sequential CNS
            steps). Defaults to ``native_scores`` when not provided.
        output_mutants : bool
            Whether to keep mutant PDB files
        ligand_param_fname : Union[Path, str]
            Path to additional parameter file used by CNS
        ligand_top_fname : Union[Path, str]
            Path to additional topology file used by CNS
        """
        self.model_path = Path(model_path)
        self.model_id = model_id
        self.chain = chain
        self.resid = resid
        self.ori_resname = ori_resname
        self.target_resname = target_resname
        self.partner_chain = partner_chain
        self.partner_resid = partner_resid
        self.partner_ori_resname = partner_ori_resname
        self.partner_target_resname = partner_target_resname
        self.native_scores = native_scores
        self.native_scores_2step = (
            native_scores_2step if native_scores_2step is not None else native_scores
        )
        self.output_mutants = output_mutants
        self.ligand_param_fname = ligand_param_fname
        self.ligand_top_fname = ligand_top_fname

    def _canonical_mutant_name(self) -> Path:
        """Return the canonical output name for the fully mutated base pair."""
        mut_id = (
            f"{RES_CODES[self.ori_resname]}{self.resid}"
            f"{RES_CODES[self.target_resname]}"
            f"-{RES_CODES[self.partner_ori_resname]}{self.partner_resid}"
            f"{RES_CODES[self.partner_target_resname]}"
        )
        return Path(self.model_path.name.replace(".pdb", f"-{self.chain}_{mut_id}.pdb"))

    def _score_same_type(self, mutation_id):
        """Score a same-ring-type base-pair mutation with a single CNS call.

        Both nucleotides keep their ring type (purine <-> purine and
        pyrimidine <-> pyrimidine), so the two mutations can be applied
        simultaneously and rebuilt by CNS in one pass.
        """
        sc_dir = f"haddock3-score-{mutation_id}"
        os.makedirs(sc_dir, exist_ok=True)

        # Perform the double (base-pair) mutation on the pdb file
        mut_pdb = mutate(
            self.model_path,
            self.chain,
            self.resid,
            self.target_resname,
            self.partner_chain,
            self.partner_resid,
            self.partner_target_resname,
        )

        mutant_scores = calc_score(
            mut_pdb,
            run_dir=sc_dir,
            outputpdb=self.output_mutants,
            ligand_param_fname=self.ligand_param_fname,
            ligand_top_fname=self.ligand_top_fname,
        )

        # Handle output files
        em_mut_pdb = Path(f"{mut_pdb.stem}_hs.pdb")
        if not self.output_mutants:
            for f in (mut_pdb, em_mut_pdb):
                if os.path.exists(f):
                    os.remove(f)
        elif os.path.exists(em_mut_pdb):
            # keep the energy-minimized mutant under the canonical name
            shutil.move(em_mut_pdb, mut_pdb)
        if os.path.exists(sc_dir):
            shutil.rmtree(sc_dir)

        return mutant_scores

    def _score_cross_type(self, mutation_id):
        """Score a cross-ring-type base-pair mutation with two CNS calls.

        A cross-type base-pair mutation swaps the ring types of the two
        nucleotides (one purine -> pyrimidine and its partner
        pyrimidine -> purine). Building a whole purine ring from its three
        anchor atoms in the same pass as the partner mutation is unreliable, so
        the mutation is performed in two sequential, CNS-regularised steps:

        1. the purine -> pyrimidine mutation is applied and energy-minimised by
           CNS, and
        2. the pyrimidine -> purine mutation is then applied to that minimised
           intermediate and scored by CNS.

        The scores and energies of the *second* (final) CNS call are returned
        and used for the plots.
        """
        # Identify the purine->pyrimidine nucleotide (step 1) and the
        # pyrimidine->purine one (step 2).
        if self.ori_resname in PURINES:
            step1_key, step1_target = (self.chain, self.resid), self.target_resname
            step2_key, step2_target = (
                (self.partner_chain, self.partner_resid),
                self.partner_target_resname,
            )
        else:
            step1_key, step1_target = (
                (self.partner_chain, self.partner_resid),
                self.partner_target_resname,
            )
            step2_key, step2_target = (self.chain, self.resid), self.target_resname

        sc_dir_1 = f"haddock3-score-{mutation_id}-1"
        sc_dir_2 = f"haddock3-score-{mutation_id}-2"
        os.makedirs(sc_dir_1, exist_ok=True)
        os.makedirs(sc_dir_2, exist_ok=True)

        intermediates = []
        try:
            # Step 1: purine -> pyrimidine, energy-minimised by CNS
            step1_pdb = _mutate_residues(self.model_path, {step1_key: step1_target})
            intermediates.append(step1_pdb)
            calc_score(
                step1_pdb,
                run_dir=sc_dir_1,
                outputpdb=True,
                ligand_param_fname=self.ligand_param_fname,
                ligand_top_fname=self.ligand_top_fname,
            )
            step1_min = Path(f"{step1_pdb.stem}_hs.pdb")
            intermediates.append(step1_min)
            if not step1_min.exists():
                raise FileNotFoundError(
                    f"CNS did not produce the energy-minimized intermediate "
                    f"{step1_min} for the first mutation step."
                )

            # Step 2: pyrimidine -> purine on the minimised intermediate
            step2_pdb = _mutate_residues(step1_min, {step2_key: step2_target})
            intermediates.append(step2_pdb)
            mutant_scores = calc_score(
                step2_pdb,
                run_dir=sc_dir_2,
                outputpdb=self.output_mutants,
                ligand_param_fname=self.ligand_param_fname,
                ligand_top_fname=self.ligand_top_fname,
            )
            step2_min = Path(f"{step2_pdb.stem}_hs.pdb")

            # Handle output files
            if self.output_mutants and step2_min.exists():
                # keep the final energy-minimized mutant under the canonical name
                shutil.move(step2_min, self._canonical_mutant_name())
            elif step2_min.exists():
                intermediates.append(step2_min)
        finally:
            for f in intermediates:
                if os.path.exists(f):
                    os.remove(f)
            for d in (sc_dir_1, sc_dir_2):
                if os.path.exists(d):
                    shutil.rmtree(d)

        return mutant_scores

    def run(self):
        """Execute the base-pair (double) mutation."""
        mutation_id = (
            f"{self.model_id}_{self.chain}{self.resid}{self.target_resname}"
            f"_{self.partner_chain}{self.partner_resid}{self.partner_target_resname}"
        )

        try:
            # Cross-type base-pair mutations (a purine and a pyrimidine swapping
            # ring types) are performed in two sequential CNS-regularised steps;
            # same-type mutations can be done in a single CNS call. Each mutant
            # is compared against the wild-type baseline that went through the
            # same number of CNS minimisation passes.
            if is_cross_type(self.ori_resname, self.target_resname):
                mutant_scores = self._score_cross_type(mutation_id)
                baseline_scores = self.native_scores_2step
            else:
                mutant_scores = self._score_same_type(mutation_id)
                baseline_scores = self.native_scores

            # Calculate deltas (native - mutant)
            n_score, n_vdw, n_elec, n_des, n_bsa = baseline_scores
            m_score, m_vdw, m_elec, m_des, m_bsa = mutant_scores
            delta_scores = (
                n_score - m_score,
                n_vdw - m_vdw,
                n_elec - m_elec,
                n_des - m_des,
                n_bsa - m_bsa,
            )

            return MutationResult(
                model_id=self.model_id,
                chain=self.chain,
                resid=self.resid,
                ori_resname=self.ori_resname,
                target_resname=self.target_resname,
                partner_chain=self.partner_chain,
                partner_resid=self.partner_resid,
                partner_ori_resname=self.partner_ori_resname,
                partner_target_resname=self.partner_target_resname,
                mutant_scores=mutant_scores,
                delta_scores=delta_scores,
                success=True,
            )

        except Exception as e:
            return MutationResult(
                model_id=self.model_id,
                chain=self.chain,
                resid=self.resid,
                ori_resname=self.ori_resname,
                target_resname=self.target_resname,
                partner_chain=self.partner_chain,
                partner_resid=self.partner_resid,
                partner_ori_resname=self.partner_ori_resname,
                partner_target_resname=self.partner_target_resname,
                mutant_scores=(0, 0, 0, 0, 0),
                delta_scores=(0, 0, 0, 0, 0),
                success=False,
                error_msg=str(e),
            )
