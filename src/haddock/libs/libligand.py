"""Ligand topology utilities."""

import re
import subprocess
import tempfile
from pathlib import Path

from haddock.core.supported_molecules import supported_residues
from haddock.core.defaults import prodrg_exec, prodrg_param
from haddock import log
import shutil

from haddock.core.typing import FilePath


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


def run_prodrg(
    pdb_file: FilePath,
    output_dir: FilePath,
) -> tuple[Path, Path]:
    """Run prodrg on a ligand PDB and write CNS topology and parameter files.

    prodrg writes its output to fixed filenames in the current working
    directory, so the process is executed inside a managed temporary directory
    that is cleaned up automatically. The resulting files are written to
    ``output_dir`` named after the input PDB stem.

    Parameters
    ----------
    pdb_file : FilePath
        Path to the ligand PDB file.
    output_dir : FilePath
        Directory where the named ``.top`` and ``.param`` files are written.

    Returns
    -------
    tuple[Path, Path]
        Paths to the written ``<stem>_prodrg.top`` and ``<stem>_prodg.param`` files.

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

    # PRODRG will write its output in its `pwd`, so run it in a temporary one
    with tempfile.TemporaryDirectory() as tmpdir:
        # NOTE: PRODRG has a max of 79 chars as input, make sure to pass relative paths
        #  so we need to copy everything to this temporary dir
        dst = Path(tmpdir, prodrg_param.name)
        shutil.copy(prodrg_param, dst)

        dst = Path(tmpdir, pdb_file.name)
        shutil.copy(pdb_file, dst)

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

        top_path = output_dir / f"{pdb_file.stem}_prodrg.top"
        par_path = output_dir / f"{pdb_file.stem}_prodrg.param"

        top_content = _sanitize_atom_names(tmp_top.read_text())
        par_content = _sanitize_atom_names(_remove_nbonds(tmp_par.read_text()))

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
