"""Ligand topology utilities."""

import re
import subprocess
import tempfile
from pathlib import Path
from string import ascii_uppercase, digits

from haddock.core.supported_molecules import supported_residues
from haddock.core.defaults import prodrg_exec, prodrg_param
from haddock.libs.libpdb import format_atom_name
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
    element = line[76:78].strip().upper()
    if element not in metals:
        return line
    stripped = line.rstrip("\n").ljust(78)
    atom_name = stripped[12:16].strip()
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
            resname = line[17:20].strip()
            if resname not in keep:
                continue
            # chainID (col 22) + resSeq (cols 23-26) + iCode (col 27)
            res_id = line[21:27]
            if resname not in first_copy:
                first_copy[resname] = res_id
            if first_copy[resname] == res_id:
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

        result = subprocess.run(
            prodrg_args,
            cwd=tmpdir,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"prodrg failed with return code {result.returncode}:\n{result.stderr}"
            )

        tmp_top = Path(tmpdir) / "DRGCNS.TOP"
        tmp_par = Path(tmpdir) / "DRGCNS.PAR"

        # TODO: Check if the atom names have been changed!
        # tmp_pdb = Path(tmpdir) / "<???>.pdb"

        if not tmp_top.exists() or not tmp_par.exists():
            ls = list(Path(tmpdir).iterdir())
            prodrg_err = Path(tmpdir) / "DRGDRG.ERR"
            prodrg_log = Path(tmpdir) / "DRGDRG.LOG"
            log.error(f"DRGDRG.log: {prodrg_log.read_text()}")
            log.error(f"DRGDRG.err: {prodrg_err.read_text()}")
            log.debug(f"ls: {ls}")
            raise RuntimeError(
                f"prodrg finished but expected output files are missing in {tmpdir} "
            )

        top_content = _sanitize_atom_names(tmp_top.read_text())
        par_content = _sanitize_atom_names(_remove_nbonds(tmp_par.read_text()))

    return top_content, par_content


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
                top_content, par_content = _run_prodrg_single(ligand_pdb, cnssep=cnssep)
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
