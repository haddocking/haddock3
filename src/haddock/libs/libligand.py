"""Ligand topology utilities."""

import re
import subprocess
import tempfile
from pathlib import Path

from haddock.core.supported_molecules import supported_residues
from haddock.core.defaults import prodrg_exec, prodrg_param
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
                out.write(line)
        out.write("END" + "\n")
    return dest


def _run_prodrg_single(pdb_file: Path) -> tuple[str, str]:
    """Run prodrg on a single small-molecule PDB and return its CNS content.

    prodrg writes its output to fixed filenames in the current working
    directory, so the process is executed inside a managed temporary directory
    that is cleaned up automatically.

    Parameters
    ----------
    pdb_file : Path
        Path to a PDB file containing a single small molecule.

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

        result = subprocess.run(
            # NOTE: We need this `PDBELEM` flag here
            [str(prodrg_exec), str(pdb_file.name), str(prodrg_param.name), "PDBELEM"],
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
    parameter files are concatenated into one ``.top`` and one ``.param``.

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
        for resname in dict.fromkeys(ligand_resnames):
            with tempfile.TemporaryDirectory() as tmpdir:
                ligand_pdb = extract_ligand(
                    pdb_file, [resname], Path(tmpdir, f"{resname}.pdb")
                )
                top_content, par_content = _run_prodrg_single(ligand_pdb)
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
