"""Virtual canonical representations for CNS job identity.

The representation is checksum-side only: ordinary HADDOCK runs execute their
generated inputs in the normal step layout. Cache-key construction will become
the production consumer in the caching stage; executable canonical-workspace
validation is deliberately deferred to the later audit stage.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path

from haddock.core.typing import Optional, Sequence, Union


_ASSIGNMENT_PATTERNS = (
    re.compile(
        r"eval(?:uate)?\s*\(\s*\$(?P<name>[A-Za-z0-9_]+)\s*=\s*"
        r"(?P<value>.*)\s*\)"
    ),
    re.compile(r"\{===>\}\s*(?P<name>[A-Za-z0-9_]+)\s*=\s*(?P<value>.*?)\s*;"),
)

_REFERENCE_PATTERN = re.compile(
    r"@@?(?P<target>"
    r"MODULE:[^\s;]+|TOPPAR:[^\s;]+|"
    r"MODULE/[^\s;]+|TOPPAR/[^\s;]+|"
    r"[$&][A-Za-z0-9_]*_[$&][A-Za-z0-9_]+|"
    r"[$&][A-Za-z0-9_]+|[^\s;]+)"
)
_INDEXED_VARIABLE_REFERENCE_PATTERN = re.compile(
    r"(?P<prefix>[$&][A-Za-z0-9_]*_)(?P<index>[$&][A-Za-z0-9_]+)$"
)
_DYNAMIC_TOPPAR_PREFIX_PATTERN = re.compile(
    r'"?(?P<prefix>TOPPAR(?::|/)[^"\s]+?)"?\s*\+\s*encode\('
)
_COUNT_SUFFIX_REFERENCE_PATTERN = re.compile(
    r"(?P<base>[$&][A-Za-z0-9_]+|[^+\s]+)"
    r'\s*\+\s*"(?P<separator>[^"]*)"\s*\+\s*'
    r"encode\(\s*(?P<count>[$&][A-Za-z0-9_]+|\d+)\s*\)"
)
_LOCATOR_COUNT_ASSIGNMENT_PATTERN = re.compile(
    r"(?im)(?P<prefix>"
    r"(?:eval(?:uate)?\s*\(\s*\$|\{===>\}\s*)"
    r"count\s*=\s*"
    r")(?P<value>[^);\r\n]*)(?P<suffix>[);])"
)


#: What the ``canonical-cns`` pin is bound to.
#:
#: The CNS executable is not an input to the computation; it is the machine
#: that evaluates it.  The computation a job declares is its script together
#: with the data that script reads, and the binary is the interpreter of that
#: declaration.  So the pin exists -- the executable is a named position in
#: the mapping, not something the representation forgets about -- but it is
#: bound to a constant rather than to the executable's own bytes.
#:
#: Binding it to the bytes would make an identity that no two installations
#: can ever share, because no two of them compile or download the same binary.
#: That is not a corner case: sharing identities across installations is the
#: principal use of ever computing one -- a lab-wide store, a store published
#: alongside a paper, a workstation result reused on a cluster.
#:
#: Nor would the bytes buy safety.  They cannot detect a CNS build that
#: computes *different* results; they only prevent two builds that agree from
#: being recognised as agreeing, which is the overwhelmingly common case.  A
#: build that genuinely disagrees is a reproducibility problem this
#: representation cannot fix and should not pretend to: pointing HADDOCK3 at a
#: different engine and expecting different numbers is a user error that no
#: identity scheme catches for free.  The honest answer is to record which
#: executable produced a result as *provenance*, where a mixture of builds is
#: visible and auditable without being part of identity.
_CNS_EXECUTABLE_POLICY_PIN = hashlib.sha256(
    b"haddock3-canonical-cns-executable-is-not-an-input"
).hexdigest()


@dataclass(frozen=True)
class CanonicalDependency:
    """One immutable CNS input in its location-independent role."""

    original_path: Path
    canonical_name: str
    checksum: str


@dataclass(frozen=True)
class CanonicalMapping:
    """The complete, immutable identity mapping of a CNS job."""

    canonical_script: str
    dependencies: tuple[CanonicalDependency, ...]
    cns_exec: Path
    output_paths: tuple[Path, ...]
    canonical_output_names: tuple[str, ...]
    output_shape: str
    invariant_dependencies: tuple[str, ...]
    work_dir: Path
    module_dir: Path
    toppar_dir: Path

    @property
    def checksums(self) -> dict[str, str]:
        """Return the canonical checksum tree for this CNS job."""
        return {
            "canonical-cns": _CNS_EXECUTABLE_POLICY_PIN,
            "canonical.inp": _checksum_bytes(self.canonical_script.encode()),
            **{dep.canonical_name: dep.checksum for dep in self.dependencies},
        }

    @property
    def dependency_paths(self) -> dict[Path, str]:
        """Map resolved original paths to their canonical names."""
        return {
            dependency.original_path: dependency.canonical_name
            for dependency in self.dependencies
        }


@dataclass
class CNSDependencyScan:
    """Resolved CNS read dependencies for one job."""

    read_files: list[Path]
    unresolved_reads: list[str]


class _IgnoredReference:
    """Sentinel for optional references that do not resolve to a file."""


IGNORED_REFERENCE = _IgnoredReference()


def build_canonical_mapping(
    input_file: Path | str,
    *,
    envvars: dict[str, str],
    cns_exec: Path,
    output_files: Sequence[Path] = (),
    output_pdb_files: Sequence[Path] = (),
    work_dir: Path | None = None,
) -> CanonicalMapping:
    """Resolve every CNS read and rewrite the job into canonical names."""
    work_dir = (work_dir or Path.cwd()).resolve()
    module_dir = _resolve_env_path(envvars["MODULE"], work_dir)
    toppar_dir = _resolve_env_path(envvars["TOPPAR"], work_dir)
    script = _cns_job_script(input_file, work_dir)
    scan = _scan_cns_job_dependencies(input_file, envvars, work_dir)
    if scan.unresolved_reads:
        unresolved = ", ".join(scan.unresolved_reads)
        raise ValueError(f"Canonical CNS input has unresolved reads: {unresolved}")

    paths = list(dict.fromkeys(path.resolve() for path in scan.read_files))
    variables = _script_path_variables(script, work_dir, module_dir, toppar_dir)
    canonical_names = _canonical_dependency_names(
        paths,
        module_dir,
        toppar_dir,
        variables,
    )
    dependencies = tuple(
        CanonicalDependency(
            path,
            canonical_names[path],
            compression_transparent_checksum(path),
        )
        for path in paths
    )

    normalized_outputs = tuple(_absolute_path(path, work_dir) for path in output_files)
    pdb_outputs = {_absolute_path(path, work_dir) for path in output_pdb_files}
    for output in normalized_outputs:
        if output.suffix.lower() == ".pdb":
            pdb_outputs.add(output)
    if len(pdb_outputs) != 1 or len(normalized_outputs) not in (1, 2):
        raise ValueError(
            "A canonical CNS job must declare exactly one PDB and at most one "
            "PSF output."
        )
    psf_outputs = [
        output for output in normalized_outputs if output.suffix.lower() == ".psf"
    ]
    if len(psf_outputs) > 1:
        raise ValueError("A canonical CNS job may declare at most one PSF output.")

    pdb_output = next(iter(pdb_outputs))
    canonical_output_names = ("canonical-output.pdb",) + (
        ("canonical-output.psf",) if psf_outputs else ()
    )
    output_name_map = {pdb_output: "canonical-output.pdb"}
    if psf_outputs:
        output_name_map[psf_outputs[0]] = "canonical-output.psf"

    canonical_script = _rewrite_canonical_script(
        script,
        work_dir,
        module_dir,
        toppar_dir,
        {
            dependency.original_path: dependency.canonical_name
            for dependency in dependencies
        },
        output_name_map,
        text_names=_script_dependency_aliases(
            script,
            work_dir,
            module_dir,
            toppar_dir,
            canonical_names,
        ),
    )
    _assert_canonical_script(
        canonical_script,
        work_dir,
        module_dir,
        toppar_dir,
        [
            *[
                dependency.original_path
                for dependency in dependencies
                if not dependency.canonical_name.startswith(("module/", "toppar/"))
            ],
            *normalized_outputs,
        ],
        [
            *[dependency.canonical_name for dependency in dependencies],
            *output_name_map.values(),
        ],
        canonical_output_names,
    )
    cns_exec = _absolute_path(cns_exec, work_dir)
    return CanonicalMapping(
        canonical_script=canonical_script,
        dependencies=dependencies,
        cns_exec=cns_exec,
        output_paths=(pdb_output, *psf_outputs),
        canonical_output_names=canonical_output_names,
        output_shape="pdb+psf" if psf_outputs else "pdb",
        invariant_dependencies=tuple(
            sorted(
                ["canonical-cns"]
                + [
                    dependency.canonical_name
                    for dependency in dependencies
                    if dependency.canonical_name.startswith(("module/", "toppar/"))
                ]
            )
        ),
        work_dir=work_dir,
        module_dir=module_dir,
        toppar_dir=toppar_dir,
    )


def canonical_mapping_for_job(job) -> CanonicalMapping:
    """Build the canonical mapping for a CNSJob without importing its class."""
    return build_canonical_mapping(
        job.input_file,
        envvars=job.envvars,
        cns_exec=Path(job.cns_exec),
        output_files=job.output_files,
        output_pdb_files=job.output_pdb_files,
        work_dir=job.work_dir,
    )


def scan_cns_dependencies(
    input_file: Path,
    envvars: dict[str, str],
) -> CNSDependencyScan:
    """Resolve explicit CNS read dependencies."""
    input_file = input_file.resolve()
    workdir = input_file.parent
    module_dir = _resolve_env_path(envvars["MODULE"], workdir)
    toppar_dir = _resolve_env_path(envvars["TOPPAR"], workdir)

    read_files: list[Path] = []
    seen_read_files: set[Path] = set()
    unresolved_reads: list[str] = []
    visited: set[Path] = set()

    def add_read_file(path: Path) -> None:
        path = path.resolve()
        if path in seen_read_files:
            return
        seen_read_files.add(path)
        read_files.append(path)

    def scan_file(path: Path, variables: dict[str, str]) -> None:
        path = path.resolve()
        if path in visited:
            return
        visited.add(path)

        if not path.exists():
            unresolved_reads.append(str(path))
            return

        add_read_file(path)
        local_vars = dict(variables)
        guard_stack: list[tuple[str, bool]] = []

        for raw_line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
            line = raw_line.split("!", 1)[0].strip()
            if not line:
                continue

            lowered = line.lower()
            if lowered.startswith("if ("):
                guard_var = _extract_nonempty_guard_variable(line)
                if guard_var is not None:
                    guard_stack.append((guard_var, True))
            elif lowered.startswith("end if") and guard_stack:
                guard_stack.pop()

            for pattern in _ASSIGNMENT_PATTERNS:
                match = pattern.search(line)
                if match:
                    value = _normalize_assignment_value(match.group("value"))
                    local_vars[match.group("name")] = value
                    dynamic_prefix = _extract_dynamic_toppar_prefix(value)
                    if dynamic_prefix is not None:
                        local_vars[f"__dynamic_prefix__{match.group('name')}"] = (
                            dynamic_prefix
                        )
                    break

            for match in _REFERENCE_PATTERN.finditer(line):
                token = match.group("target")
                resolved = _resolve_reference(
                    token=token,
                    current_file=path,
                    workdir=workdir,
                    module_dir=module_dir,
                    toppar_dir=toppar_dir,
                    variables=local_vars,
                )
                if resolved is IGNORED_REFERENCE:
                    continue
                if resolved is None:
                    if _is_guarded_optional_reference(token, guard_stack):
                        continue
                    unresolved_reads.append(token)
                    continue
                if isinstance(resolved, list):
                    for resolved_path in resolved:
                        if not resolved_path.exists():
                            unresolved_reads.append(str(resolved_path))
                            continue
                        add_read_file(resolved_path)
                        if resolved_path.suffix.lower() == ".cns":
                            scan_file(resolved_path, local_vars)
                    continue
                if not isinstance(resolved, Path):
                    continue
                if not resolved.exists():
                    unresolved_reads.append(str(resolved))
                    continue

                add_read_file(resolved)
                if resolved.suffix.lower() == ".cns":
                    scan_file(resolved, local_vars)

    scan_file(input_file, {})
    return CNSDependencyScan(
        read_files=read_files,
        unresolved_reads=sorted(set(unresolved_reads)),
    )


def compression_transparent_checksum(path: Path) -> str:
    """Checksum a file by its logical bytes, transparently to gzip storage."""
    path = path.resolve()
    data = path.read_bytes()
    if path.name.endswith(".gz"):
        data = gzip.decompress(data)
    return _checksum_bytes(data)


def _checksum_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _cns_job_script(input_file: Path | str, work_dir: Path) -> str:
    """Read a materialized CNS input or return its in-memory contents."""
    if isinstance(input_file, Path):
        return _absolute_path(input_file, work_dir).read_text(encoding="utf-8")
    return input_file


def _scan_cns_job_dependencies(
    input_file: Path | str,
    envvars: dict[str, str],
    work_dir: Path,
) -> CNSDependencyScan:
    """Resolve the files read by a materialized or in-memory CNS input."""
    if isinstance(input_file, Path):
        script_path = _absolute_path(input_file, work_dir)
        scan = scan_cns_dependencies(script_path, envvars)
        return CNSDependencyScan(
            [path for path in scan.read_files if path != script_path],
            scan.unresolved_reads,
        )

    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=work_dir,
        suffix=".inp",
        delete=False,
    ) as handle:
        handle.write(input_file)
        temporary_input = Path(handle.name)
    try:
        scan = scan_cns_dependencies(temporary_input, envvars)
        return CNSDependencyScan(
            [path for path in scan.read_files if path != temporary_input.resolve()],
            scan.unresolved_reads,
        )
    finally:
        temporary_input.unlink(missing_ok=True)


def _script_path_variables(
    script: str,
    work_dir: Path,
    module_dir: Path,
    toppar_dir: Path,
) -> dict[Path, str]:
    """Associate assignment variables with resolved paths for named roles."""
    result: dict[Path, str] = {}
    variables: dict[str, str] = {}
    for line in script.splitlines():
        for pattern in _ASSIGNMENT_PATTERNS:
            match = pattern.search(line)
            if match is None:
                continue
            variable = match.group("name")
            value = _normalize_assignment_value(match.group("value"))
            variables[variable] = value
            resolved = _resolve_count_suffix_reference(
                token=value,
                current_file=work_dir / "canonical.inp",
                workdir=work_dir,
                module_dir=module_dir,
                toppar_dir=toppar_dir,
                variables=variables,
                seen_variables=set(),
            )
            if resolved is None:
                resolved = _resolve_reference(
                    token=value,
                    current_file=work_dir / "canonical.inp",
                    workdir=work_dir,
                    module_dir=module_dir,
                    toppar_dir=toppar_dir,
                    variables=variables,
                )
            if isinstance(resolved, Path) and resolved.exists():
                result[resolved.resolve()] = _canonical_role_variable(
                    variable,
                    value,
                )
            break
    return result


def _script_dependency_aliases(
    script: str,
    work_dir: Path,
    module_dir: Path,
    toppar_dir: Path,
    canonical_names: dict[Path, str],
) -> dict[str, str]:
    """Map textual dependency aliases to the resolved canonical role name."""
    aliases: dict[str, str] = {}
    variables: dict[str, str] = {}
    for line in script.splitlines():
        for pattern in _ASSIGNMENT_PATTERNS:
            match = pattern.search(line)
            if match is None:
                continue
            variable = match.group("name")
            value = _normalize_assignment_value(match.group("value"))
            variables[variable] = value
            dynamic_prefix = _extract_dynamic_toppar_prefix(value)
            if dynamic_prefix is not None:
                variables[f"__dynamic_prefix__{variable}"] = dynamic_prefix

            resolved = _resolve_count_suffix_reference(
                token=value,
                current_file=work_dir / "canonical.inp",
                workdir=work_dir,
                module_dir=module_dir,
                toppar_dir=toppar_dir,
                variables=variables,
                seen_variables=set(),
            )
            if not isinstance(resolved, Path):
                resolved = _resolve_reference(
                    token=value,
                    current_file=work_dir / "canonical.inp",
                    workdir=work_dir,
                    module_dir=module_dir,
                    toppar_dir=toppar_dir,
                    variables=variables,
                )
            if not isinstance(resolved, Path):
                break

            canonical_name = canonical_names.get(resolved.resolve())
            if canonical_name is None:
                break
            if canonical_name.startswith(("module/", "toppar/")):
                # install-tree files keep their MODULE:/TOPPAR: spelling
                break

            if not value.startswith(("$", "&")):
                # a variable *reference* is already location-independent; rewriting it
                # would replace the CNS symbol name itself rather than a filename
                aliases[value] = canonical_name
            count_match = _COUNT_SUFFIX_REFERENCE_PATTERN.fullmatch(value.strip())
            if count_match is not None:
                base_value = _resolve_variable_text(
                    count_match.group("base"),
                    variables,
                )
                if base_value:
                    aliases[base_value] = canonical_name
            break
    return aliases


def _resolve_variable_text(token: str, variables: dict[str, str]) -> Optional[str]:
    """Resolve a textual CNS variable reference without resolving filesystem paths."""
    seen: set[str] = set()
    while token.startswith(("$", "&")):
        variable_name = token[1:]
        if variable_name in seen or variable_name not in variables:
            return None
        seen.add(variable_name)
        token = variables[variable_name]
    return token


def _canonical_role_variable(variable: str, value: str) -> str:
    """Return the variable name that carries a dependency's stable role."""
    match = _COUNT_SUFFIX_REFERENCE_PATTERN.fullmatch(value.strip())
    if match is not None and match.group("base").startswith(("$", "&")):
        return match.group("base")[1:].lower()
    return variable.lower()


def _canonical_dependency_names(
    paths: Sequence[Path],
    module_dir: Path,
    toppar_dir: Path,
    variables: dict[Path, str],
) -> dict[Path, str]:
    """Assign stable roles in dependency first-reference order."""
    names: dict[Path, str] = {}
    counters = {"pdb": 0, "psf": 0, "generic": 0}
    named_roles = (
        ("unambig", "canonical-unambig.tbl"),
        ("hbond", "canonical-hbond.tbl"),
        ("ambig", "canonical-ambig.tbl"),
        ("dihe", "canonical-dihe.tbl"),
        ("sym", "canonical-symmetry.tbl"),
        ("tensor", "canonical-tensor.tbl"),
        ("ligand_top", "canonical-ligand.top"),
        ("ligand_param", "canonical-ligand.param"),
    )
    for path in paths:
        if path.is_relative_to(module_dir):
            names[path] = f"module/{path.relative_to(module_dir).as_posix()}"
            continue
        if path.is_relative_to(toppar_dir):
            names[path] = f"toppar/{path.relative_to(toppar_dir).as_posix()}"
            continue
        variable = variables.get(path, "")
        named = next(
            (
                name
                for marker, name in named_roles
                if re.search(rf"(?:^|_){re.escape(marker)}(?:_|$)", variable)
            ),
            None,
        )
        if named is not None and named not in names.values():
            names[path] = named
            continue
        logical_name = _strip_known_compression_suffix(path.name)
        suffix = Path(logical_name).suffix.lower()
        if suffix in (".pdb", ".psf"):
            counters[suffix[1:]] += 1
            names[path] = f"canonical-input-{counters[suffix[1:]]}{suffix}"
        elif suffix == ".tbl" and variable.startswith(("input_aa_", "input_cgtbl_")):
            counters["generic"] += 1
            names[path] = f"canonical-cg-input-{counters['generic']}.tbl"
        else:
            counters["generic"] += 1
            names[path] = f"canonical-input-{counters['generic']}{suffix}"
    return names


def _canonicalize_install_paths(
    script: str,
    module_dir: Path,
    toppar_dir: Path,
) -> str:
    """Give absolute install-tree paths their environment-relative CNS spelling.

    ``MODULE:`` and ``TOPPAR:`` are resolved by CNS through environment variables,
    so a reference spelled that way is already independent of where HADDOCK3 is
    installed. Paths written out in full are not, and reach the canonical form even
    when the recipe never reads them -- which is how the install path used to enter
    the identity of every topology job.
    """

    def replace(match: re.Match[str]) -> str:
        candidate = Path(match.group("path")).resolve()
        for root, prefix in ((module_dir, "MODULE:"), (toppar_dir, "TOPPAR:")):
            if candidate.is_relative_to(root):
                return f"{prefix}{candidate.relative_to(root).as_posix()}"
        return match.group(0)

    return re.sub(r'(?<![\w.-])(?P<path>/[^\s"\';]+)', replace, script)


def _rewrite_canonical_script(
    script: str,
    work_dir: Path,
    module_dir: Path,
    toppar_dir: Path,
    dependency_names: dict[Path, str],
    output_names: dict[Path, str],
    text_names: Optional[dict[str, str]] = None,
) -> str:
    """Replace only resolved job-specific path spellings in CNS text.

    Two rules keep the result both canonical and executable:

    - install-tree references keep their ``MODULE:``/``TOPPAR:`` spelling. The
      ``module/`` and ``toppar/`` names are pin names for the checksum tree; writing
      them into the script instead would leave the top-level input and the module
      scripts it includes disagreeing about what ``MODULE`` points at.
    - only *path* spellings are rewritten, never a bare basename, so that an input
      and an output sharing a basename cannot consume each other's binding.
    """
    result = _canonicalize_install_paths(script, module_dir, toppar_dir)
    job_paths = {
        path: name
        for path, name in dependency_names.items()
        if not name.startswith(("module/", "toppar/"))
    }
    for path, name in {**job_paths, **output_names}.items():
        candidates = {str(path), path.as_posix()}
        candidates.add(os.path.relpath(path, start=work_dir))
        for match in re.finditer(
            r'(?<![\w.-])(?P<path>(?:\.\.?/|/)[^\s"\';]+)',
            result,
        ):
            candidate = match.group("path")
            if _absolute_path(Path(candidate), work_dir) == path:
                candidates.add(candidate)
        for candidate in sorted(candidates, key=len, reverse=True):
            if candidate:
                result = _replace_script_token(result, candidate, name)
    for candidate, name in sorted(
        (text_names or {}).items(),
        key=lambda item: len(item[0]),
        reverse=True,
    ):
        if candidate:
            result = _replace_script_token(result, candidate, name)
    return _canonicalize_logging_level(_canonicalize_locator_count_assignment(result))


def _replace_script_token(script: str, token: str, replacement: str) -> str:
    """Replace a path-like token without matching inside another path token."""
    pattern = re.compile(rf"(?<![A-Za-z0-9_./-]){re.escape(token)}(?![A-Za-z0-9_./-])")
    return pattern.sub(replacement, script)


def _canonicalize_locator_count_assignment(script: str) -> str:
    """Erase generated structure index values from canonical CNS input."""
    return _LOCATOR_COUNT_ASSIGNMENT_PATTERN.sub(
        r"\g<prefix>canonical-count\g<suffix>",
        script,
    )


def _canonicalize_logging_level(script: str) -> str:
    """Normalize CNS interpreter verbosity, which cannot affect job artifacts."""
    return re.sub(
        r'(?P<prefix>\$log_level\s*=\s*)"[^"]*"',
        r'\g<prefix>"canonical-log-level"',
        script,
    )


def _assert_output_bindings(
    script: str,
    canonical_output_names: Sequence[str],
) -> None:
    """The script must write exactly the outputs the mapping declares.

    Output shape is part of the transformation, so a script whose
    ``$output_pdb_filename`` says something other than the declared canonical name
    is describing a different computation from the one being checksummed.
    """
    for variable, canonical in (
        ("output_pdb_filename", "canonical-output.pdb"),
        ("output_psf_filename", "canonical-output.psf"),
    ):
        match = re.search(rf'\${variable}\s*=\s*"(?P<name>[^"]*)"', script)
        if match is None:
            continue
        if canonical not in canonical_output_names:
            raise ValueError(
                f"Canonical CNS script writes ${variable} but no {canonical!r} "
                "output was declared"
            )
        if match.group("name") != canonical:
            raise ValueError(
                f"Canonical CNS script binds ${variable} to "
                f"{match.group('name')!r} instead of the declared {canonical!r}"
            )


def _assert_canonical_script(
    script: str,
    work_dir: Path,
    module_dir: Path,
    toppar_dir: Path,
    paths: Sequence[Path],
    canonical_names: Sequence[str] = (),
    canonical_output_names: Sequence[str] = (),
) -> None:
    """Reject location-dependent canonical scripts before they become identity."""
    _assert_output_bindings(script, canonical_output_names)
    scanned = script
    for name in sorted(canonical_names, key=len, reverse=True):
        if name:
            scanned = scanned.replace(name, "\x00")
    leaked = [str(work_dir), str(module_dir), str(toppar_dir)]
    leaked.extend(path.name for path in paths if path.name)
    for token in leaked:
        if token and token in scanned:
            raise ValueError(
                f"Canonical CNS script leaked {token!r} from job in {work_dir}"
            )
    match = re.search(r"(?:^|/)\d+_[A-Za-z][A-Za-z0-9_]*", scanned)
    if match:
        raise ValueError(
            f"Canonical CNS script leaked step-folder token {match.group(0)!r}"
        )


def _absolute_path(path: Path, work_dir: Path) -> Path:
    return path.resolve() if path.is_absolute() else (work_dir / path).resolve()


def _resolve_env_path(raw_path: str, workdir: Path) -> Path:
    path = Path(raw_path)
    if path.is_absolute():
        return path.resolve()
    return (workdir / path).resolve()


def _resolve_reference(
    *,
    token: str,
    current_file: Path,
    workdir: Path,
    module_dir: Path,
    toppar_dir: Path,
    variables: dict[str, str],
    seen_variables: Optional[set[str]] = None,
) -> Optional[Union[Path, list[Path], _IgnoredReference]]:
    seen_variables = seen_variables or set()
    indexed_reference = _INDEXED_VARIABLE_REFERENCE_PATTERN.fullmatch(token)
    if indexed_reference is not None:
        prefix = indexed_reference.group("prefix")[1:]
        indexed_variables = sorted(
            variable
            for variable in variables
            if re.fullmatch(rf"{re.escape(prefix)}\d+", variable)
        )
        if not indexed_variables:
            return None
        resolved_paths: list[Path] = []
        for variable in indexed_variables:
            resolved = _resolve_reference(
                token=f"${variable}",
                current_file=current_file,
                workdir=workdir,
                module_dir=module_dir,
                toppar_dir=toppar_dir,
                variables=variables,
                seen_variables=set(seen_variables),
            )
            if not isinstance(resolved, Path):
                return None
            resolved_paths.append(resolved)
        return resolved_paths
    resolved_from_variable = False
    if token.startswith(("$", "&")):
        variable_name = token[1:]
        if variable_name in seen_variables:
            return None
        seen_variables.add(variable_name)
        dynamic_prefix = variables.get(f"__dynamic_prefix__{variable_name}")
        if variable_name not in variables:
            if dynamic_prefix is not None:
                prefix_path = _resolve_reference(
                    token=dynamic_prefix,
                    current_file=current_file,
                    workdir=workdir,
                    module_dir=module_dir,
                    toppar_dir=toppar_dir,
                    variables=variables,
                    seen_variables=seen_variables,
                )
                if isinstance(prefix_path, Path):
                    return sorted(prefix_path.parent.glob(f"{prefix_path.name}*"))
            return None
        token = variables[variable_name]
        if token == "":
            return IGNORED_REFERENCE
        resolved_from_variable = True
        dynamic_count_path = _resolve_count_suffix_reference(
            token=token,
            current_file=current_file,
            workdir=workdir,
            module_dir=module_dir,
            toppar_dir=toppar_dir,
            variables=variables,
            seen_variables=seen_variables,
        )
        if dynamic_count_path is not None:
            return dynamic_count_path
        if token.startswith(("$", "&")):
            return _resolve_reference(
                token=token,
                current_file=current_file,
                workdir=workdir,
                module_dir=module_dir,
                toppar_dir=toppar_dir,
                variables=variables,
                seen_variables=seen_variables,
            )
        if dynamic_prefix is not None and any(
            fragment in token for fragment in ("+", "encode(")
        ):
            prefix_path = _resolve_reference(
                token=dynamic_prefix,
                current_file=current_file,
                workdir=workdir,
                module_dir=module_dir,
                toppar_dir=toppar_dir,
                variables=variables,
                seen_variables=seen_variables,
            )
            if isinstance(prefix_path, Path):
                return sorted(prefix_path.parent.glob(f"{prefix_path.name}*"))

    if not token or any(fragment in token for fragment in ("+", "encode(", "$", "&")):
        return None

    if token.startswith("MODULE:"):
        return (module_dir / token.split(":", 1)[1].lstrip("/")).resolve()

    if token.startswith("TOPPAR:"):
        return (toppar_dir / token.split(":", 1)[1].lstrip("/")).resolve()

    if token.startswith("MODULE/"):
        return (module_dir / token.removeprefix("MODULE/").lstrip("/")).resolve()

    if token.startswith("TOPPAR/"):
        return (toppar_dir / token.removeprefix("TOPPAR/").lstrip("/")).resolve()

    candidate = Path(token)
    if candidate.is_absolute():
        return candidate.resolve()
    relative_base = workdir if resolved_from_variable else current_file.parent
    return (relative_base / candidate).resolve()


def _extract_nonempty_guard_variable(line: str) -> Optional[str]:
    match = re.search(r'\$([A-Za-z0-9_]+)\s*#\s*""', line)
    if match is not None:
        return match.group(1)

    match = re.search(r"&BLANK%([A-Za-z0-9_]+)\s*=\s*false", line, re.IGNORECASE)
    if match is not None:
        return match.group(1)

    return None


def _resolve_count_suffix_reference(
    *,
    token: str,
    current_file: Path,
    workdir: Path,
    module_dir: Path,
    toppar_dir: Path,
    variables: dict[str, str],
    seen_variables: set[str],
) -> Optional[Path]:
    match = _COUNT_SUFFIX_REFERENCE_PATTERN.fullmatch(token.strip())
    if match is None:
        return None

    base = _resolve_reference(
        token=match.group("base"),
        current_file=current_file,
        workdir=workdir,
        module_dir=module_dir,
        toppar_dir=toppar_dir,
        variables=variables,
        seen_variables=set(seen_variables),
    )
    count = _resolve_literal_or_variable(
        match.group("count"),
        variables=variables,
    )
    if not isinstance(base, Path) or count is None:
        return None

    counted = Path(f"{base}{match.group('separator')}{count}")
    if counted.exists():
        return counted.resolve()
    if base.exists():
        return base.resolve()
    return None


def _resolve_literal_or_variable(
    token: str,
    *,
    variables: dict[str, str],
) -> Optional[str]:
    if token.isdigit():
        return token
    if not token.startswith(("$", "&")):
        return None
    variable = variables.get(token[1:])
    return variable if variable and variable.isdigit() else None


def _is_guarded_optional_reference(
    token: str,
    guard_stack: list[tuple[str, bool]],
) -> bool:
    if not token.startswith(("$", "&")):
        return False
    variable_name = token[1:]
    return any(
        guard_var == variable_name and is_optional
        for guard_var, is_optional in guard_stack
    )


def _extract_dynamic_toppar_prefix(value: str) -> Optional[str]:
    match = _DYNAMIC_TOPPAR_PREFIX_PATTERN.search(value)
    if match is None:
        return None
    return match.group("prefix")


def _normalize_assignment_value(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        return value[1:-1]
    return value


def _strip_known_compression_suffix(name: str) -> str:
    if name.endswith(".gz"):
        return name[:-3]
    return name
