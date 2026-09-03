"""Pure, standalone normalization of CNS-generated output artifacts.

This module deliberately imports only the Python standard library. Stage 4
materializes this exact source file as a transformation input and executes it
directly, while HADDOCK imports the same functions for production CNS jobs.
"""

from __future__ import annotations

import gzip
import os
import sys
import uuid
from pathlib import Path
from typing import Callable, Sequence


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


def normalize_cns_pdb(path: str | Path) -> bool:
    """Remove run-volatile CNS header lines from a PDB file."""
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
    return b"".join(
        line
        for line in _split_lf_records(pdb_bytes)
        if not line.startswith(CNS_PDB_VOLATILE_PREFIXES)
    )


def is_normalized_cns_pdb(path: str | Path) -> bool:
    """Return whether a PDB, compressed or not, has stable CNS headers."""
    return _is_normalized(path, normalize_cns_pdb_bytes)


def normalize_cns_psf(path: str | Path) -> bool:
    """Remove the run-volatile CNS date stamp from a PSF file."""
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
    return b"".join(_normalize_psf_line(line) for line in _split_lf_records(psf_bytes))


def is_normalized_cns_psf(path: str | Path) -> bool:
    """Return whether a PSF, compressed or not, has a stable CNS title."""
    return _is_normalized(path, normalize_cns_psf_bytes)


def _is_normalized(
    path: str | Path,
    normalize: Callable[[bytes], bytes],
) -> bool:
    """Whether a file, compressed or not, is unchanged by ``normalize``."""
    artifact_path = Path(path)
    artifact_bytes = artifact_path.read_bytes()
    if artifact_path.name.endswith(".gz"):
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


def is_normalized_cns_artifact(
    path: str | Path,
    logical_name: str | None = None,
) -> bool:
    """Return whether a CNS output artifact is free of run-volatile content."""
    name = logical_name if logical_name is not None else Path(path).name
    if name.endswith(".gz"):
        name = name[:-3]
    suffix = Path(name).suffix.lower()
    if suffix == ".pdb":
        return is_normalized_cns_pdb(path)
    if suffix == ".psf":
        return is_normalized_cns_psf(path)
    return True


def main(argv: Sequence[str] | None = None) -> int:
    """Normalize each named PDB or PSF, for a staged CNS transformation."""
    paths = list(sys.argv[1:] if argv is None else argv)
    if not paths:
        print("usage: normalize-cns-output.py FILE.pdb [FILE.psf ...]", file=sys.stderr)
        return 2
    for raw_path in paths:
        path = Path(raw_path)
        if not path.is_file():
            print(f"CNS did not produce required output: {path}", file=sys.stderr)
            return 1
        suffix = Path(path.name.removesuffix(".gz")).suffix.lower()
        if suffix == ".pdb":
            normalize_cns_pdb(path)
        elif suffix == ".psf":
            normalize_cns_psf(path)
        else:
            print(f"unsupported CNS output type: {path}", file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
