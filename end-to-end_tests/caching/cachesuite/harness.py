"""Run ``haddock3`` and observe the result directory.

Everything here is black-box.  A test invokes ``haddock3 file.cfg --cache ...``
as a subprocess and then looks at inodes, bytes and the clock.

Gate 1 -- inode identity -- is the primary, per-job, discriminating observable:
a cached output is a hit iff ``(st_dev, st_ino)`` equals that of a specific
file in a specific source run.  Gate 2 -- wall clock -- is aggregate and
coarse, and catches the one thing Gate 1 cannot: a "hit" that ran CNS anyway.
Gate 3 -- content equality -- is weak and used only to confirm that a copy
delivered the right bytes.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import signal
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path

from .config import CACHEABLE_MODULES, TOPOLOGY_MODULES
from .runio import step_outputs

#: Compression suffixes a cache may hardlink an artifact under.
COMPRESSED_SUFFIXES = (".gz", ".zst")


def haddock3_executable() -> str:
    """The ``haddock3`` under test."""
    explicit = os.environ.get("HADDOCK3_CACHE_EXECUTABLE")
    if explicit:
        return explicit
    found = shutil.which("haddock3")
    if found is None:
        raise RuntimeError(
            "haddock3 is not on PATH; set HADDOCK3_CACHE_EXECUTABLE to point at it"
        )
    return found


@dataclass
class RunResult:
    """What one ``haddock3`` invocation produced."""

    run_dir: Path
    returncode: int | None
    duration: float
    stdout: str
    #: True only when the process was alive at the deadline and had to be
    #: killed.  A self-exit before the deadline is never a timeout.
    killed: bool = False
    command: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return self.returncode == 0

    def tail(self, lines: int = 25) -> str:
        return "\n".join(self.stdout.splitlines()[-lines:])


def run_haddock3(
    config: Path,
    *,
    cache: tuple[Path, ...] | list[Path] = (),
    env: dict[str, str] | None = None,
    cwd: Path | None = None,
    timeout: float | None = None,
    install: Path | None = None,
    executable: str | None = None,
) -> RunResult:
    """Invoke ``haddock3`` as a subprocess in its own process group.

    ``timeout`` kills the whole process group, not just the ``haddock3``
    parent.  One orphaned CNS process per killed run would accumulate fast in a
    suite that deliberately kills hundreds of runs.
    """
    config = Path(config)
    cwd = Path(cwd) if cwd is not None else config.parent
    command = [executable or haddock3_executable(), str(config)]
    for source in cache:
        command += ["--cache", str(source)]

    process_env = dict(os.environ)
    process_env.pop("HADDOCK_CACHE_HARDLINK", None)
    if install is not None:
        existing = process_env.get("PYTHONPATH", "")
        process_env["PYTHONPATH"] = (
            f"{install}{os.pathsep}{existing}" if existing else str(install)
        )
    if env:
        for key, value in env.items():
            if value is None:
                process_env.pop(key, None)
            else:
                process_env[key] = str(value)

    started = time.monotonic()
    process = subprocess.Popen(
        command,
        cwd=str(cwd),
        env=process_env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
        start_new_session=True,
    )
    killed = False
    try:
        stdout, _ = process.communicate(timeout=timeout)
        returncode = process.returncode
    except subprocess.TimeoutExpired:
        killed = True
        _kill_process_group(process)
        try:
            stdout, _ = process.communicate(timeout=60)
        except subprocess.TimeoutExpired:  # pragma: no cover - defensive
            process.kill()
            stdout, _ = process.communicate()
        returncode = None
    duration = time.monotonic() - started

    run_dir = _resolve_run_dir(config, cwd)
    return RunResult(
        run_dir=run_dir,
        returncode=returncode,
        duration=duration,
        stdout=stdout or "",
        killed=killed,
        command=command,
    )


def _kill_process_group(process: subprocess.Popen) -> None:
    try:
        group = os.getpgid(process.pid)
    except ProcessLookupError:  # pragma: no cover - race with a clean exit
        return
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(group, sig)
        except ProcessLookupError:
            return
        deadline = time.monotonic() + 10.0
        while time.monotonic() < deadline:
            if process.poll() is not None:
                return
            time.sleep(0.05)


def _resolve_run_dir(config: Path, cwd: Path) -> Path:
    for line in config.read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if stripped.startswith("run_dir"):
            _, _, value = stripped.partition("=")
            name = value.strip().strip('"').strip("'")
            candidate = Path(name)
            return candidate if candidate.is_absolute() else (cwd / candidate)
    raise ValueError(f"{config} declares no run_dir")


# --------------------------------------------------------------------------
# Gate 1 -- inode identity
# --------------------------------------------------------------------------


def inode(path: Path) -> tuple[int, int]:
    """``(st_dev, st_ino)``.  Both fields; ``st_ino`` alone is not identity."""
    stat = path.stat()
    return (stat.st_dev, stat.st_ino)


def artifact_variants(path: Path) -> list[Path]:
    """``path`` plus the compressed spellings a cache may deliver it under."""
    return [path] + [Path(f"{path}{suffix}") for suffix in COMPRESSED_SUFFIXES]


def link_source(path: Path, sources: list[Path]) -> Path | None:
    """Return the file in ``sources`` sharing an inode with ``path``.

    ``st_nlink > 1`` proves *a* link exists but not *to what*, so the whole
    source corpus is indexed by ``(st_dev, st_ino)`` and matched exactly.  That
    is what makes a key collision -- a hit from the *wrong* entry, the
    catastrophic failure -- visible rather than merely plausible.
    """
    index = source_inode_index(sources)
    for variant in artifact_variants(path):
        if not variant.exists():
            continue
        found = index.get(inode(variant))
        if found is not None:
            return found
    return None


_INODE_INDEX_CACHE: dict[tuple[str, ...], dict[tuple[int, int], Path]] = {}


def source_inode_index(sources: list[Path]) -> dict[tuple[int, int], Path]:
    """Map ``(st_dev, st_ino)`` to a path, over every file in every source.

    Sources are read-only by contract during the suite, so this is cached per
    source set.  Inode numbers are never stored in a manifest -- they are
    collected at assertion time, because they go stale the moment a fixture is
    regenerated.
    """
    key = tuple(sorted(str(Path(source).resolve()) for source in sources))
    cached = _INODE_INDEX_CACHE.get(key)
    if cached is not None:
        return cached
    index: dict[tuple[int, int], Path] = {}
    for source in sources:
        root = Path(source).resolve()
        for path in root.rglob("*"):
            if not path.is_file() or path.is_symlink():
                continue
            index.setdefault(inode(path), path)
    _INODE_INDEX_CACHE[key] = index
    return index


def forget_inode_index() -> None:
    """Drop the cached index; call after any fixture is regenerated."""
    _INODE_INDEX_CACHE.clear()


# --------------------------------------------------------------------------
# Declared outputs of a run
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class Artifact:
    """One cacheable output file, addressed run-relative."""

    relative: str
    module: str
    #: Index of this step among the run's steps.
    step_index: int
    #: Index of this step among the run's steps *of the same module*.
    occurrence: int
    kind: str  # "pdb" or "psf"

    def path(self, run_dir: Path) -> Path:
        return run_dir / self.relative


def step_folders(run_dir: Path) -> list[tuple[int, str, Path]]:
    """``(index, module, path)`` for every numbered step folder, in order."""
    folders = []
    for entry in sorted(Path(run_dir).iterdir()):
        if not entry.is_dir():
            continue
        head, _, module = entry.name.partition("_")
        if not head.isdigit() or not module:
            continue
        folders.append((int(head), module, entry))
    folders.sort(key=lambda item: item[0])
    return folders


def cacheable_artifacts(run_dir: Path) -> list[Artifact]:
    """Every ``.pdb`` (and topology ``.psf``) a CNS module of this run declared.

    Enumerated from each step's public ``io.json`` rather than by globbing:
    ``topoaa`` writes its split ensemble *inputs* into its own folder next to
    its outputs, and globbing would declare those as cacheable results.
    """
    artifacts: list[Artifact] = []
    seen_modules: dict[str, int] = {}
    for index, module, folder in step_folders(run_dir):
        occurrence = seen_modules.get(module, 0)
        seen_modules[module] = occurrence + 1
        if module not in CACHEABLE_MODULES:
            continue
        for output in step_outputs(folder):
            relative = _relative(run_dir, folder, output.file_name)
            artifacts.append(Artifact(relative, module, index, occurrence, "pdb"))
            if module in TOPOLOGY_MODULES and output.psf_name:
                artifacts.append(
                    Artifact(
                        _relative(run_dir, folder, output.psf_name),
                        module,
                        index,
                        occurrence,
                        "psf",
                    )
                )
    return artifacts


def _relative(run_dir: Path, folder: Path, name: str) -> str:
    return (folder / name).relative_to(Path(run_dir)).as_posix()


# --------------------------------------------------------------------------
# Gate 3 -- content equality
# --------------------------------------------------------------------------


def content_checksum(path: Path) -> str:
    """SHA-256 of a file, transparently through ``.gz``/``.zst`` storage."""
    for variant in artifact_variants(path):
        if not variant.exists():
            continue
        data = variant.read_bytes()
        if variant.suffix == ".gz":
            import gzip

            data = gzip.decompress(data)
        elif variant.suffix == ".zst":
            import zstandard

            data = zstandard.ZstdDecompressor().decompress(data, max_output_size=1 << 30)
        return hashlib.sha256(data).hexdigest()
    raise FileNotFoundError(path)


def exists_any(path: Path) -> bool:
    """Whether the artifact is present under any storage spelling."""
    return any(variant.exists() for variant in artifact_variants(path))
