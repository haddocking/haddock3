"""Ligand topology utilities."""

import math
import re
import subprocess
import tempfile
from pathlib import Path
from string import ascii_uppercase, digits

from haddock.core.supported_molecules import supported_residues
from haddock.core.defaults import prodrg_exec, prodrg_param
from haddock.libs.libpdb import (
    format_atom_name,
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
from haddock import log
import shutil

from haddock.core.typing import FilePath, Optional


def identify_unknown_hetatms(pdb_file: FilePath) -> list[str]:
    """Return residue names in a PDB that are not in the supported residues.

    Parameters
    ----------
    pdb_file : FilePath
        Path to the PDB file to inspect.

    Returns
    -------
    list[str]
        Unique residue names not found in the supported residues set,
        in order of first appearance.
    """
    seen: list[str] = []
    with open(pdb_file) as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[17:20].strip()
            if resname and resname not in supported_residues and resname not in seen:
                seen.append(resname)
    return seen


# Two-letter elements prodrg is able to parametrise (halides). Every other
# supported element (H, C, N, O, S, P, F, I) is a single letter.
_PRODRG_TWO_LETTER_ELEMENTS = ("CL", "BR")


def _unsupported_metal_symbols() -> set[str]:
    """Two-letter metal/ion element symbols prodrg cannot handle.

    CNS can mislabel a ligand atom's element as a two-letter metal: a
    beta-phosphorus named ``PB`` ends up with element ``PB`` (lead), which
    prodrg rejects. The set of such symbols is read from ``ion.top`` (its
    residue names are the ion element symbols, optionally suffixed with a
    charge). Isolated ions are never passed to prodrg, so any of these symbols
    appearing on a ligand atom is a mislabelling. The two-letter halides prodrg
    *does* support (``Cl``, ``Br``) are excluded so genuine halogen atoms are
    left untouched.

    Returns
    -------
    set[str]
        Upper-cased two-letter element symbols to treat as mislabellings.
    """
    from haddock import toppar_path

    ion_top = Path(toppar_path, "ion.top")
    metals: set[str] = set()
    for line in ion_top.read_text().splitlines():
        match = re.match(r"\s*RESI\w*\s+([A-Za-z]+)", line)
        if match:
            symbol = match.group(1).upper()
            if len(symbol) == 2 and symbol not in _PRODRG_TWO_LETTER_ELEMENTS:
                metals.add(symbol)
    return metals


def _demetalise_atom(line: str, metals: set[str]) -> str:
    """Fix a ligand atom line that CNS mislabelled as a two-letter metal.

    When the element column holds an unsupported two-letter metal symbol (e.g.
    ``PB`` for a phosphate phosphorus), the element is collapsed to its first
    character and the atom name is re-justified to the single-letter-element
    convention. Lines that are not mislabelled are returned unchanged.

    Parameters
    ----------
    line : str
        A PDB ``ATOM``/``HETATM`` line.
    metals : set[str]
        Two-letter symbols to treat as mislabellings (see
        :func:`_unsupported_metal_symbols`).

    Returns
    -------
    str
        The (possibly corrected) PDB line, newline-terminated.
    """
    element = line[slc_element].strip().upper()
    if element not in metals:
        return line
    stripped = line.rstrip("\n").ljust(78)
    atom_name = stripped[slc_name].strip()
    # organic ligand atoms only carry single-letter elements here, so the real
    # element is the first character of the mislabelled two-letter symbol.
    true_element = element[0]
    name_field = format_atom_name(atom_name, true_element)
    fixed = (
        stripped[:12]
        + name_field
        + stripped[16:76]
        + true_element.rjust(2)
        + stripped[78:]
    )
    return fixed + "\n"


SUPPORTED_DOC_URL = (
    "https://www.bonvinlab.org/haddock3-user-manual/supported_particles.html"
)
ION_DOC_URL = f"{SUPPORTED_DOC_URL}#supported-ions"


def _charged_ions() -> dict[str, list[tuple[str, str]]]:
    """Map each ion element to its charged residues in ``ion.top``.

    Only single-atom residues whose name carries a charge (e.g. ``ZN2``) are
    kept, which leaves out the neutral atoms (e.g. ``ZN``) and multi-atom
    ions (e.g. ``SO4``).

    Returns
    -------
    dict[str, list[tuple[str, str]]]
        Element symbol mapped to the (residue name, atom name) pairs, e.g.
        ``{"CL": [("CL1", "CL-1")], "FE": [("FE2", "FE+2"), ("FE3", "FE+3")]}``.
    """
    from haddock import toppar_path

    ions: dict[str, list[tuple[str, str]]] = {}
    resname = None
    atoms: list[str] = []
    for line in Path(toppar_path, "ion.top").read_text().splitlines():
        words = line.split()
        if not words:
            continue
        keyword = words[0].upper()
        if keyword.startswith("RESI") and len(words) > 1:
            resname, atoms = words[1].upper(), []
        elif keyword == "ATOM" and resname and len(words) > 1:
            atoms.append(words[1].upper())
        elif keyword.startswith("END") and resname:
            match = re.fullmatch(r"([A-Z]+)(\d)", resname)
            if len(atoms) == 1 and match:
                ions.setdefault(match.group(1), []).append((resname, atoms[0]))
            resname = None
    return ions


def ion_naming_hint(pdb_file: FilePath, resnames: list[str]) -> str:
    """Explain the ion naming for unknown residues that look like ions.

    A residue looks like an ion when it consists of a single atom, or when its
    name is an ion element followed by a charge that is not defined (e.g.
    ``CU3``, or ``LI`` without a charge).

    Parameters
    ----------
    pdb_file : FilePath
        PDB file containing the residues.
    resnames : list[str]
        Names of the unknown residues.

    Returns
    -------
    str
        The explanation, or an empty string if no residue looks like an ion.
    """
    charged = _charged_ions()
    # atom elements of the first copy of each residue
    elements: dict[str, list[str]] = {}
    first_copy: dict[str, str] = {}
    with open(pdb_file) as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[slc_resname].strip()
            if resname in resnames and (
                first_copy.setdefault(resname, line[21:27]) == line[21:27]
            ):
                element = line[slc_element].strip().upper()
                elements.setdefault(resname, []).append(element)

    hints = []
    for resname in resnames:
        atoms = elements.get(resname, [])
        candidates = [resname.rstrip("0123456789+-").upper()]
        if len(atoms) == 1:
            candidates.append(atoms[0])
        known = next((charged[c] for c in candidates if c in charged), None)
        if known:
            names = " or ".join(
                f"residue {res} with atom {atom}" for res, atom in known
            )
            hints.append(f"- {resname}: if this is an ion, use {names}")
        elif len(atoms) == 1:
            hints.append(f"- {resname}: single atom, not a supported ion")
    if not hints:
        return ""
    return (
        "Ions need the absolute charge in the residue name and the exact charge, "
        "positive or negative, in the atom name (e.g. residue ZN2 with atom ZN+2, "
        "residue CL1 with atom CL-1):\n" + "\n".join(hints) + f"\nSee {ION_DOC_URL}"
    )


# Atoms bonding consecutive residues of a polymer chain: (previous, next)
_CHAIN_LINKS = (("C", "N"), ("O3'", "P"))


def find_chain_residues(
    pdb_file: FilePath,
    resnames: list[str],
    cutoff: float = 2.0,
) -> list[str]:
    """Find the residues covalently linked to their neighbours in a chain.

    Such residues (e.g. a modified amino acid within a protein) cannot be
    parameterised by prodrg, which treats every residue as a free molecule
    and does not create the bonds to the rest of the chain.

    A residue is linked when its N (or P) is within ``cutoff`` of the C (or
    O3') of the previous residue in the same chain, or its C (or O3') within
    ``cutoff`` of the N (or P) of the next one.

    Parameters
    ----------
    pdb_file : FilePath
        PDB file containing the residues.
    resnames : list[str]
        Names of the residues to check.
    cutoff : float
        Maximum bond length, in Angstrom.

    Returns
    -------
    list[str]
        The linked residues, as ``"<resname> (chain <chain>, residue <resid>)"``.
    """
    residues = _read_residues(pdb_file)

    def bonded(previous: dict, following: dict) -> bool:
        return any(
            atom1 in previous
            and atom2 in following
            and math.dist(previous[atom1], following[atom2]) < cutoff
            for atom1, atom2 in _CHAIN_LINKS
        )

    def neighbour(index: int, chain: str) -> dict:
        """Atoms of the residue at ``index`` if it is in ``chain``."""
        if 0 <= index < len(keys) and keys[index][0] == chain:
            return residues[keys[index]][1]
        return {}

    keys = list(residues)
    linked = []
    for i, key in enumerate(keys):
        resname, atoms = residues[key]
        if resname not in resnames:
            continue
        if bonded(neighbour(i - 1, key[0]), atoms) or bonded(
            atoms, neighbour(i + 1, key[0])
        ):
            linked.append(f"{resname} ({_residue_label(key)})")
    return linked


# Ring atoms of the sugars and their anomeric carbon
_SUGAR_RINGS = (
    ({"C1", "C2", "C3", "C4", "C5", "O5"}, "C1"),  # aldopyranoses
    ({"C1", "C2", "C3", "C4", "O4"}, "C1"),  # aldofuranoses
    ({"C2", "C3", "C4", "C5", "C6", "O6"}, "C2"),  # ketopyranoses, e.g. SIA
    ({"C2", "C3", "C4", "C5", "O5"}, "C2"),  # ketofuranoses
)


def find_linked_sugars(
    pdb_file: FilePath,
    resnames: list[str],
    cutoff: float = 2.0,
) -> list[str]:
    """Find the sugars linked to another residue by a glycosidic bond.

    Like residues within a chain, such sugars (e.g. an unsupported glycan
    unit bound to a protein) cannot be parameterised by prodrg, which does
    not create the glycosidic bonds.

    A residue is a sugar when it has the ring atoms of a pyranose or a
    furanose (see ``_SUGAR_RINGS``). It is linked when its anomeric carbon
    (C1, or C2 for ketoses such as sialic acids) is within ``cutoff`` of an O,
    N or S atom of another residue (e.g. ASN ND2, or O4 of another sugar), or
    when one of its O atoms is within ``cutoff`` of the anomeric carbon of
    another sugar.

    Parameters
    ----------
    pdb_file : FilePath
        PDB file containing the residues.
    resnames : list[str]
        Names of the residues to check.
    cutoff : float
        Maximum bond length, in Angstrom.

    Returns
    -------
    list[str]
        The linked sugars and their partners, as
        ``"<resname> (chain <chain>, residue <resid>) linked to <partner>"``.
    """
    residues = _read_residues(pdb_file)
    # possible anomeric carbons of each sugar; a wrong candidate (e.g. C2 of
    # an aldose) is harmless, as only bonds to other residues are considered
    anomeric = {
        key: list(dict.fromkeys(c for ring, c in _SUGAR_RINGS if ring <= set(atoms)))
        for key, (_, atoms) in residues.items()
    }
    anomeric = {key: carbons for key, carbons in anomeric.items() if carbons}
    linked = []
    for key, carbons in anomeric.items():
        resname, atoms = residues[key]
        if resname not in resnames:
            continue
        partners = []
        for other, (other_resname, other_atoms) in residues.items():
            if other == key:
                continue
            label = f"{other_resname} {other[1].strip()}"
            # glycosidic bond from the anomeric carbon of this sugar...
            partners.extend(
                f"{label} {name}"
                for carbon in carbons
                for name, coords in other_atoms.items()
                if name[0] in "ONS" and math.dist(atoms[carbon], coords) < cutoff
            )
            # ...or to the anomeric carbon of another sugar
            partners.extend(
                f"{label} {carbon}"
                for carbon in anomeric.get(other, [])
                for name, coords in atoms.items()
                if name[0] == "O" and math.dist(other_atoms[carbon], coords) < cutoff
            )
        if partners:
            linked.append(
                f"{resname} ({_residue_label(key)}) linked to {', '.join(partners)}"
            )
    return sorted(linked)


def _read_residues(
    pdb_file: FilePath,
) -> dict[tuple[str, str], tuple[str, dict[str, tuple[float, float, float]]]]:
    """Read the atom coordinates of each residue, in the order of the file.

    Returns
    -------
    dict
        ``(chain, resid)`` mapped to ``(resname, {atom name: coordinates})``.
    """
    residues: dict = {}
    with open(pdb_file) as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            key = (line[slc_chainid], line[slc_resseq] + line[slc_icode])
            atoms = residues.setdefault(key, (line[slc_resname].strip(), {}))[1]
            atoms[line[slc_name].strip()] = (
                float(line[slc_x]),
                float(line[slc_y]),
                float(line[slc_z]),
            )
    return residues


def _residue_label(key: tuple[str, str]) -> str:
    """Format a ``(chain, resid)`` key as ``"chain A, residue 12"``."""
    return f"chain {key[0].strip() or '-'}, residue {key[1].strip()}"


def extract_ligand(
    pdb_file: FilePath,
    resnames: list[str],
    dest: FilePath,
) -> Path:
    """Write a PDB with a single copy of each of the given residue names.

    prodrg can only handle a single small molecule. When the ligand is part of
    a larger system (e.g. a protein/ligand complex), the surrounding
    macromolecule must be stripped out before prodrg is run, otherwise it fails
    to generate topology and parameter files.

    When a ligand is present in several copies (same residue name, different
    residue numbers/chains), all copies share the same topology, so only the
    first copy of each residue name is kept to avoid feeding prodrg redundant
    (and potentially too many) atoms.

    Atoms that CNS mislabelled as a two-letter metal element (e.g. a phosphate
    ``PB`` written with element ``Pb``) are corrected back to their real
    single-letter element so prodrg can parametrise them.

    Parameters
    ----------
    pdb_file : FilePath
        Path to the input PDB file (may contain a full system).
    resnames : list[str]
        Residue names to keep (typically the unknown ligand residues).
    dest : FilePath
        Path to the PDB file to write with only the selected residues.

    Returns
    -------
    Path
        Path to the written ligand-only PDB file (``dest``).
    """
    keep = set(resnames)
    metals = _unsupported_metal_symbols()
    # residue identity (chain + residue number + insertion code) of the first
    # copy seen for each residue name; only atoms of that copy are written out.
    first_copy: dict[str, str] = {}
    dest = Path(dest)
    with open(pdb_file) as fh, open(dest, "w") as out:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            resname = line[slc_resname].strip()
            if resname not in keep:
                continue
            # chainID (col 22) + resSeq (cols 23-26) + iCode (col 27)
            res_id = line[21:27]
            if first_copy.setdefault(resname, res_id) == res_id:
                out.write(_demetalise_atom(line, metals))
        out.write("END" + "\n")
    return dest


def _run_prodrg_single(
    pdb_file: Path,
    cnssep: Optional[str] = None,
) -> tuple[str, str]:
    """Run prodrg on a single small-molecule PDB and return its CNS content.

    prodrg writes its output to fixed filenames in the current working
    directory, so the process is executed inside a managed temporary directory
    that is cleaned up automatically.

    Parameters
    ----------
    pdb_file : Path
        Path to a PDB file containing a single small molecule.
    cnssep : Optional[str]
        Single character passed to prodrg as ``CNSSEP=<x>``. prodrg embeds it
        in every generated atom type name (e.g. ``CA82``/``HA82`` for ``A``),
        so giving each ligand a unique letter keeps their atom types from
        clashing when several topologies are concatenated.

    Returns
    -------
    tuple[str, str]
        Sanitized topology and parameter file contents.

    Raises
    ------
    RuntimeError
        If prodrg exits with a non-zero return code or the expected output
        files are not created.
    """
    natoms = sum(
        line.startswith(("ATOM  ", "HETATM"))
        for line in pdb_file.read_text().splitlines()
    )
    if natoms < 3:
        # prodrg rejects them, e.g. ions
        raise RuntimeError(
            f"prodrg cannot handle molecules of less than three atoms (found {natoms})"
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        # NOTE: PRODRG has a max of 79 chars as input, make sure to pass relative
        #  paths so we need to copy everything to this temporary dir
        shutil.copy(prodrg_param, Path(tmpdir, prodrg_param.name))
        shutil.copy(pdb_file, Path(tmpdir, pdb_file.name))

        # NOTE: We need this `PDBELEM` flag here
        prodrg_args = [
            str(prodrg_exec),
            str(pdb_file.name),
            str(prodrg_param.name),
            "PDBELEM",
        ]
        if cnssep:
            prodrg_args.append(f"CNSSEP={cnssep}")

        try:
            result = subprocess.run(
                prodrg_args,
                cwd=tmpdir,
                capture_output=True,
                text=True,
            )
        except OSError as err:
            raise RuntimeError(f"could not execute {prodrg_exec}: {err}") from err

        messages = _prodrg_messages(f"{result.stdout}\n{result.stderr}")
        if result.returncode != 0:
            raise RuntimeError(
                f"prodrg failed with return code {result.returncode}:\n{messages}"
            )

        tmp_top = Path(tmpdir) / "DRGCNS.TOP"
        tmp_par = Path(tmpdir) / "DRGCNS.PAR"

        # TODO: Check if the atom names have been changed!
        # tmp_pdb = Path(tmpdir) / "<???>.pdb"

        # prodrg exits with 0 even when it fails, without writing its output
        if not tmp_top.exists() or not tmp_par.exists():
            log.debug(f"prodrg output:\n{result.stdout}\n{result.stderr}")
            raise RuntimeError(f"prodrg did not generate any output:\n{messages}")

        top_content = _sanitize_atom_names(tmp_top.read_text())
        par_content = _sanitize_atom_names(_remove_nbonds(tmp_par.read_text()))

    return top_content, par_content


def _prodrg_messages(output: str, nlines: int = 5) -> str:
    """Return the part of the prodrg output that explains a failure.

    prodrg prints its errors on stdout, prefixed with ``ERRDRG>``, after a
    long banner. The line just before the first error, which usually names
    the offending atom, is kept as well. Without any error line, the last
    ``nlines`` lines are returned.

    Parameters
    ----------
    output : str
        Standard output (and error) of prodrg.
    nlines : int
        Number of lines returned when no ``ERRDRG>`` line is found.

    Returns
    -------
    str
        The relevant prodrg output lines.
    """
    lines = [
        line.strip()
        for line in output.splitlines()
        if line.strip() not in ("", "PRODRG>")
    ]
    first_error = next(
        (i for i, line in enumerate(lines) if line.startswith("ERRDRG>")),
        None,
    )
    if first_error is None:
        return "\n".join(lines[-nlines:])
    return "\n".join(lines[max(first_error - 1, 0) :])


def _used_cofactor_separators() -> set[str]:
    """Return the CNSSEP separator characters already used by the cofactors.

    prodrg encodes its ``CNSSEP`` character as the second character of every
    atom type name it generates (e.g. ``CA82`` for ``A``). The built-in
    ``cofactors.top`` topology was generated the same way, so its atom type
    names occupy some of these separators. Reusing one of them for an
    auto-generated ligand could produce atom types that clash with the
    cofactor ones, hence they must be avoided.

    Returns
    -------
    set[str]
        Upper-cased second characters of every ``MASS`` atom type name defined
        in ``cofactors.top``.
    """
    from haddock import toppar_path

    cofactors_top = Path(toppar_path, "cofactors.top")
    separators: set[str] = set()
    for line in cofactors_top.read_text().splitlines():
        match = re.match(r"\s*MASS\s+(\S+)", line)
        if match and len(match.group(1)) >= 2:
            separators.add(match.group(1)[1].upper())
    return separators


def run_prodrg(
    pdb_file: FilePath,
    output_dir: FilePath,
    ligand_resnames: Optional[list[str]] = None,
) -> tuple[Path, Path]:
    """Run prodrg on a ligand PDB and write CNS topology and parameter files.

    The resulting files are written to ``output_dir`` named after the input
    PDB stem.

    prodrg can only process a single small molecule at a time. When
    ``ligand_resnames`` lists several distinct ligands, prodrg is run once per
    ligand (on a single extracted copy of each) and the resulting topology and
    parameter files are concatenated into one ``.top`` and one ``.param``. Each
    ligand is given a unique ``CNSSEP`` character so that prodrg generates
    non-overlapping atom type names across the concatenated files; the
    characters used are also chosen to avoid the separators already present in
    the built-in ``cofactors.top`` topology, so a single auto-generated ligand
    likewise cannot clash with the cofactor atom types.

    Parameters
    ----------
    pdb_file : FilePath
        Path to the ligand PDB file. It may also be a full system (e.g. a
        protein/ligand complex), in which case ``ligand_resnames`` must be
        provided so that only the ligand residues are passed to prodrg.
    output_dir : FilePath
        Directory where the named ``.top`` and ``.param`` files are written.
    ligand_resnames : Optional[list[str]]
        Residue names of the ligand(s) to extract before running prodrg. When
        provided, each distinct ligand is extracted (a single copy) and passed
        to prodrg separately, stripping out the rest of the system (e.g. the
        protein). When ``None``, the whole ``pdb_file`` is passed to prodrg
        unchanged.

    Returns
    -------
    tuple[Path, Path]
        Paths to the written ``<stem>_prodrg.top`` and ``<stem>_prodrg.param`` files.

    Raises
    ------
     RuntimeError
        If prodrg exits with a non-zero return code or the expected output
        files are not created.
    """
    if prodrg_exec is None or prodrg_param is None:
        raise RuntimeError(
            "prodrg is not available on this platform. "
            "Provide the binary path via the PRODRG_EXEC environment variable."
        )

    pdb_file = Path(pdb_file).resolve()
    output_dir = Path(output_dir).resolve()

    top_path = output_dir / f"{pdb_file.stem}_prodrg.top"
    par_path = output_dir / f"{pdb_file.stem}_prodrg.param"

    if ligand_resnames:
        top_parts: list[str] = []
        par_parts: list[str] = []
        # PRODRG only handles one molecule at a time, so run it once per distinct
        # ligand (on a single extracted copy) and concatenate the results.
        distinct_resnames = list(dict.fromkeys(ligand_resnames))
        # Give each ligand a unique CNSSEP character so their prodrg atom type
        # names do not overlap once concatenated, avoiding the separators
        # already used by the built-in cofactors topology.
        excluded = _used_cofactor_separators()
        separators = [c for c in ascii_uppercase + digits if c not in excluded]
        if len(distinct_resnames) > len(separators):
            raise RuntimeError(
                "Cannot auto-generate topologies for more than "
                f"{len(separators)} distinct ligands."
            )
        for idx, resname in enumerate(distinct_resnames):
            cnssep = separators[idx]
            with tempfile.TemporaryDirectory() as tmpdir:
                ligand_pdb = extract_ligand(
                    pdb_file, [resname], Path(tmpdir, f"{resname}.pdb")
                )
                try:
                    top_content, par_content = _run_prodrg_single(
                        ligand_pdb, cnssep=cnssep
                    )
                except RuntimeError as err:
                    raise RuntimeError(f"residue {resname}: {err}") from err
            top_parts.append(top_content)
            par_parts.append(par_content)
        top_path.write_text("\n".join(top_parts))
        par_path.write_text("\n".join(par_parts))
    else:
        top_content, par_content = _run_prodrg_single(pdb_file)
        top_path.write_text(top_content)
        par_path.write_text(par_content)

    return top_path, par_path


def _sanitize_atom_names(content: str) -> str:
    """Remove colons from atom type names in prodrg CNS output.

    prodrg may generate atom type names containing colons (e.g. ``HT:A``)
    and this is not compatible with CNS so they must be removed.

    Parameters
    ----------
    content : str
        Contents of a prodrg-generated CNS file.

    Returns
    -------
    str
        Content with colons stripped from non-comment lines.
    """
    lines = []
    for line in content.splitlines(keepends=True):
        if line.lstrip().startswith("!"):
            lines.append(line)
        else:
            lines.append(line.replace(":", ""))
    return "".join(lines)


def _remove_nbonds(par_content: str) -> str:
    """Remove the NBONds...END block from a prodrg CNS parameter string.

    PRODRG definition of NBONds might interfere with HADDOCK's internal parameters,
    so we must remove the NBONds lines from the param generated by HADDOCK

    Parameters
    ----------
    par_content : str
        Contents of the ``DRGCNS.PAR`` file.

    Returns
    -------
    str
        Parameter content with the NBONds block removed.
    """
    return re.sub(r"(?s)NBONds.*?END", "", par_content)
