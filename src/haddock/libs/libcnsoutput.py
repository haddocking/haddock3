"""Helpers for normalizing CNS-generated output artifacts."""

import gzip
import os
import uuid
from pathlib import Path

from haddock.core.typing import FilePath


CNS_PDB_VOLATILE_PREFIXES = (
    b"REMARK FILENAME=",
    b"REMARK initial structure ",
    b"REMARK DATE:",
    b"REMARK HADDOCK stats for ",
)
CNS_PSF_STABLE_TITLE = b"; HADDOCK3 normalized topology"


def _is_volatile_psf_line(line: bytes) -> bool:
    """Whether a PSF line is the wall-clock stamp CNS writes into its title."""
    stripped = line.strip()
    return stripped.startswith(b"DATE:") and b"created by user:" in stripped


def _normalize_psf_line(line: bytes) -> bytes:
    """Return stable PSF title content."""
    if line.strip().startswith(b"; FILENAME="):
        return CNS_PSF_STABLE_TITLE + _line_ending(line)
    if _is_volatile_psf_line(line):
        return b""
    return line


def _rewrite_atomically(path: Path, stable: bytes) -> None:
    """Replace a file's bytes without writing into it in place."""
    temporary = path.with_name(f".{path.name}.normalize-{uuid.uuid4().hex}")
    try:
        temporary.write_bytes(stable)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def normalize_cns_pdb(path: FilePath) -> bool:
    """Remove run-volatile CNS header lines from a PDB file.

    Returns ``True`` when the file bytes were changed.
    """
    pdb_path = Path(path)
    if not pdb_path.exists():
        return False

    original = pdb_path.read_bytes()
    logical = gzip.decompress(original) if pdb_path.name.endswith(".gz") else original
    stable = normalize_cns_pdb_bytes(logical)

    if stable == logical:
        return False

    _rewrite_atomically(
        pdb_path,
        gzip.compress(stable, mtime=0) if pdb_path.name.endswith(".gz") else stable,
    )
    return True


def normalize_cns_pdb_bytes(pdb_bytes: bytes) -> bytes:
    """Return PDB bytes without CNS run-volatile header lines."""
    stable_lines = [
        line
        for line in _split_lf_records(pdb_bytes)
        if not line.startswith(CNS_PDB_VOLATILE_PREFIXES)
    ]
    return b"".join(stable_lines)


def is_normalized_cns_pdb(path: FilePath) -> bool:
    """Return whether a PDB, compressed or not, has stable CNS headers."""
    return _is_normalized(path, normalize_cns_pdb_bytes)


def normalize_cns_psf(path: FilePath) -> bool:
    """Remove the run-volatile CNS date stamp from a PSF file.

    Returns ``True`` when the file bytes were changed.
    """
    psf_path = Path(path)
    if not psf_path.exists():
        return False

    original = psf_path.read_bytes()
    logical = gzip.decompress(original) if psf_path.name.endswith(".gz") else original
    stable = normalize_cns_psf_bytes(logical)

    if stable == logical:
        return False

    _rewrite_atomically(
        psf_path,
        gzip.compress(stable, mtime=0) if psf_path.name.endswith(".gz") else stable,
    )
    return True


def normalize_cns_psf_bytes(psf_bytes: bytes) -> bytes:
    """Return PSF bytes without the CNS date stamp."""
    stable_lines = [_normalize_psf_line(line) for line in _split_lf_records(psf_bytes)]
    return b"".join(stable_lines)


def is_normalized_cns_psf(path: FilePath) -> bool:
    """Return whether a PSF, compressed or not, has a stable CNS title."""
    return _is_normalized(path, normalize_cns_psf_bytes)


def _is_normalized(path: FilePath, normalize) -> bool:
    """Whether a file, compressed or not, is unchanged by ``normalize``."""
    artifact_bytes = Path(path).read_bytes()
    if Path(path).name.endswith(".gz"):
        artifact_bytes = gzip.decompress(artifact_bytes)
    return normalize(artifact_bytes) == artifact_bytes


def _split_lf_records(data: bytes) -> list[bytes]:
    """Split bytes into LF-terminated records without interpreting text."""
    records: list[bytes] = []
    start = 0
    for idx, byte in enumerate(data):
        if byte == 0x0A:
            records.append(data[start : idx + 1])
            start = idx + 1
    if start < len(data):
        records.append(data[start:])
    return records


def _line_ending(line: bytes) -> bytes:
    """Return the existing line ending for an LF-split record."""
    if line.endswith(b"\r\n"):
        return b"\r\n"
    if line.endswith(b"\n"):
        return b"\n"
    return b""
