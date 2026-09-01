"""The perturbation vocabulary.

A test case is a **pair of runs**, A then B, where B differs from A by exactly
one controlled perturbation.  A is a corpus base run; B is built here, by
applying one operation to A's configuration, inputs, or environment.

Operations are structural rather than textual: ``module`` changes a parameter
in a named block, ``insert`` adds a step at a named position.  A regular
expression over config text can silently match nothing and produce a case that
passes because it perturbed *nothing*; these cannot.

The three questions that classify any perturbation, in order:

1. does it change the content of anything in the declared read-set?
   -> MUST-MISS, localised to the jobs that read the changed thing
2. does it change a *binding* -- which pin an input occupies, or the output
   shape? -> MUST-MISS
3. does it change only a *locator* -- a path, a filename, a step ordinal, a
   rank, a run directory, a version label? -> MUST-HIT
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

from .config import Config, Step


@dataclass
class Context:
    """Everything a perturbation may act on."""

    config: Config
    #: A private copy of the system's input files, safe to edit.
    data: dict[str, Path]
    data_dir: Path
    tmp: Path
    env: dict[str, str] = field(default_factory=dict)
    #: ``PYTHONPATH`` entry for a shadowed HADDOCK3 install, if one was made.
    install: Path | None = None
    cwd: Path | None = None
    notes: list[str] = field(default_factory=list)

    def expand(self, value):
        if isinstance(value, str):
            return value.format(tmp=self.tmp, data=self.data_dir)
        if isinstance(value, list):
            return [self.expand(item) for item in value]
        if isinstance(value, dict):
            return {key: self.expand(item) for key, item in value.items()}
        return value


class PerturbationError(RuntimeError):
    """A case asks for something the vocabulary cannot express."""


def apply_all(operations: list[dict], context: Context) -> None:
    for operation in operations or []:
        apply_one(dict(operation), context)


def apply_one(operation: dict, context: Context) -> None:
    name = operation.pop("op")
    handler = _OPS.get(name)
    if handler is None:
        raise PerturbationError(f"unknown perturbation op {name!r}")
    handler(operation, context)


# -- configuration --------------------------------------------------------


def _op_top(spec: dict, context: Context) -> None:
    """Set or remove top-level configuration parameters."""
    for key, value in (spec.get("set") or {}).items():
        context.config.top[key] = context.expand(value)
    for key in spec.get("unset") or []:
        context.config.top.pop(key, None)


def _op_module(spec: dict, context: Context) -> None:
    """Set or remove parameters of one ``[module]`` block."""
    step = context.config.find(spec["module"], spec.get("occurrence", 0))
    for key, value in (spec.get("set") or {}).items():
        step.params[key] = context.expand(value)
    for key in spec.get("unset") or []:
        step.params.pop(key, None)
    for name, params in (spec.get("subsection") or {}).items():
        step.subsections.setdefault(name, {}).update(context.expand(params))


def _op_insert(spec: dict, context: Context) -> None:
    """Insert a step, by position or relative to an existing module."""
    step = Step(spec["module"], context.expand(spec.get("params") or {}))
    context.config.steps.insert(_position(spec, context), step)


def _op_remove(spec: dict, context: Context) -> None:
    index = context.config.index_of(spec["module"], spec.get("occurrence", 0))
    context.config.steps.pop(index)


def _op_duplicate(spec: dict, context: Context) -> None:
    """Repeat a module, to check a job's identity does not absorb its
    occurrence index (Axis 2.6)."""
    index = context.config.index_of(spec["module"], spec.get("occurrence", 0))
    step = context.config.steps[index]
    context.config.steps.insert(index + 1, Step(step.module, dict(step.params)))


def _op_swap_steps(spec: dict, context: Context) -> None:
    first = context.config.index_of(spec["a"], spec.get("a_occurrence", 0))
    second = context.config.index_of(spec["b"], spec.get("b_occurrence", 0))
    steps = context.config.steps
    steps[first], steps[second] = steps[second], steps[first]


def _position(spec: dict, context: Context) -> int:
    if "at" in spec:
        return int(spec["at"])
    if "before" in spec:
        return context.config.index_of(spec["before"], spec.get("occurrence", 0))
    if "after" in spec:
        return context.config.index_of(spec["after"], spec.get("occurrence", 0)) + 1
    raise PerturbationError("insert needs one of: at, before, after")


def _op_molecules(spec: dict, context: Context) -> None:
    """Reorder, replace or swap the input molecules.

    Swapping two molecules rebinds them to different pins.  That is genuine
    science, not renaming, and it must miss -- it is the probe that keeps
    locator-erasure from over-shooting into "names never matter".
    """
    molecules = list(context.config.top["molecules"])
    if "swap" in spec:
        first, second = spec["swap"]
        molecules[first], molecules[second] = molecules[second], molecules[first]
    elif spec.get("reverse"):
        molecules.reverse()
    elif "set" in spec:
        molecules = [str(context.expand(item)) for item in spec["set"]]
    else:
        raise PerturbationError("molecules needs one of: swap, reverse, set")
    context.config.top["molecules"] = molecules


# -- input files ----------------------------------------------------------


def _target(spec: dict, context: Context) -> Path:
    name = spec["file"]
    path = context.data.get(name)
    if path is None:
        raise PerturbationError(f"no input file {name!r} in this system")
    return path


def _op_edit_input(spec: dict, context: Context) -> None:
    """Change an input file's *content*.  Always MUST-MISS for its readers."""
    path = _target(spec, context)
    how = spec["how"]
    editor = _EDITORS.get(how)
    if editor is None:
        raise PerturbationError(f"unknown input edit {how!r}")
    editor(path, spec, context)


def _edit_coordinate(path: Path, spec: dict, context: Context) -> None:
    """Move one atom.  The canonical MUST-MISS probe for Axis 6.1."""
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    for index, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")) and len(line) > 54:
            x = float(line[30:38])
            lines[index] = f"{line[:30]}{x + 1.0:8.3f}{line[38:]}"
            break
    else:
        raise PerturbationError(f"{path} has no ATOM record to move")
    path.write_text("".join(lines), encoding="utf-8")


def _edit_remark(path: Path, spec: dict, context: Context) -> None:
    """Add a REMARK line: a semantically-null edit that still changes bytes."""
    text = path.read_text(encoding="utf-8")
    path.write_text("REMARK   a comment the science does not read\n" + text, "utf-8")


def _edit_whitespace(path: Path, spec: dict, context: Context) -> None:
    """Trailing whitespace and CRLF line endings; bytes change, meaning does not."""
    text = path.read_text(encoding="utf-8")
    path.write_text(text.replace("\n", " \r\n"), encoding="utf-8")


def _edit_segid(path: Path, spec: dict, context: Context) -> None:
    """Rename a chain/segid consistently throughout."""
    old, new = spec.get("from", "B"), spec.get("to", "Z")
    lines = []
    for line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith(("ATOM", "HETATM")) and len(line) > 22 and line[21] == old:
            line = f"{line[:21]}{new}{line[22:]}"
        lines.append(line)
    path.write_text("".join(lines), encoding="utf-8")


def _edit_restraint(path: Path, spec: dict, context: Context) -> None:
    """Change one distance in a restraint table."""
    text = path.read_text(encoding="utf-8")
    match = re.search(r"(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)", text)
    if match is None:
        raise PerturbationError(f"{path} has no distance triple to change")
    replaced = f"{float(match.group(1)) + 0.5:.1f} {match.group(2)} {match.group(3)}"
    path.write_text(text[: match.start()] + replaced + text[match.end() :], "utf-8")


def _edit_gzip(path: Path, spec: dict, context: Context) -> None:
    """Recompress: different bytes, same content (Axis 6.15)."""
    data = path.read_bytes()
    archive = Path(f"{path}.gz")
    with gzip.open(archive, "wb", mtime=0) as handle:
        handle.write(data)
    path.unlink()
    context.data[path.name] = archive
    _retarget(context, str(path), str(archive))


def _edit_models(path: Path, spec: dict, context: Context) -> None:
    """Add, remove or reorder members of a multi-model ensemble PDB.

    ``remove`` drops a member from the *middle*. Dropping the last one would
    leave every remaining member's name aligned with its content, so the case
    would pass without ever exercising the renumbering it claims to test.

    ``add`` and ``add_distinct`` are deliberately different perturbations and
    have opposite verdicts. ``add`` appends a byte-identical copy of an
    existing member, which is not new work at all -- the copy's topology job
    reads the same bytes as the original's and must be served from it.
    ``add_distinct`` appends a member that is genuinely a different conformer,
    which is new work and must be computed.
    """
    text = path.read_text(encoding="utf-8")
    blocks = re.findall(r"MODEL.*?ENDMDL\s*\n", text, flags=re.S)
    if len(blocks) < 3:
        raise PerturbationError(f"{path} has too few models to edit in the middle")
    how = spec.get("models", "reorder")
    if how == "reorder":
        blocks = [blocks[1], blocks[0]] + blocks[2:]
    elif how == "remove":
        middle = len(blocks) // 2
        context.notes.append(f"removed ensemble member {middle + 1} of {len(blocks)}")
        blocks = blocks[:middle] + blocks[middle + 1 :]
    elif how == "add":
        blocks = blocks + [blocks[0]]
        context.notes.append("appended a byte-identical copy of member 1")
    elif how == "add_distinct":
        blocks = blocks + [_moved_atom(blocks[0])]
        context.notes.append("appended a new conformer, one atom moved")
    else:
        raise PerturbationError(f"unknown ensemble edit {how!r}")
    head = text[: text.index(blocks_start(text))]
    path.write_text(head + "".join(blocks), encoding="utf-8")


def _moved_atom(block: str) -> str:
    """One ensemble member with a single atom displaced by 1 angstrom."""
    lines = block.splitlines(keepends=True)
    for index, line in enumerate(lines):
        if line.startswith(("ATOM", "HETATM")) and len(line) > 54:
            x = float(line[30:38])
            lines[index] = f"{line[:30]}{x + 1.0:8.3f}{line[38:]}"
            return "".join(lines)
    raise PerturbationError("ensemble member has no ATOM record to move")


def blocks_start(text: str) -> str:
    match = re.search(r"MODEL", text)
    if match is None:
        raise PerturbationError("no MODEL record")
    return text[match.start() : match.start() + 5]


def _op_copy_input(spec: dict, context: Context) -> None:
    """Copy an input to a new name, byte-identical, and use the copy.

    A path is a locator.  The same content under another name is the same
    input, so this is MUST-HIT -- and it is the probe that must sit next to
    every content edit, because a key that ignores content passes those.
    """
    path = _target(spec, context)
    destination = context.data_dir / spec["as"]
    shutil.copy2(path, destination)
    context.data[spec["as"]] = destination
    _retarget(context, str(path), str(destination))


def _op_move_input(spec: dict, context: Context) -> None:
    """Move an input elsewhere on disk, content unchanged (Axis 1.6)."""
    path = _target(spec, context)
    destination = Path(context.expand(spec["to"]))
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(path, destination)
    # Register under both spellings: later operations in the same case refer
    # to the file by whichever name they know it as, and a case that moves a
    # file and then edits it -- the adjacent MUST-MISS probe for Axis 1.6 --
    # needs the new one.
    context.data[spec["file"]] = destination
    context.data[destination.name] = destination
    _retarget(context, str(path), str(destination))


def _retarget(context: Context, old: str, new: str) -> None:
    """Point every configuration reference at a moved or renamed input."""
    molecules = context.config.top.get("molecules")
    if isinstance(molecules, list):
        context.config.top["molecules"] = [
            new if str(item) == old else item for item in molecules
        ]
    for step in context.config.steps:
        for key, value in list(step.params.items()):
            if isinstance(value, str) and value == old:
                step.params[key] = new


# -- environment and installation ----------------------------------------


def _op_env(spec: dict, context: Context) -> None:
    """Perturb the process environment (Axis 8's user-level half)."""
    for key, value in (spec.get("set") or {}).items():
        context.env[key] = str(context.expand(value))
    for key in spec.get("unset") or []:
        context.env[key] = None  # type: ignore[assignment]


def shadow_install(tmp: Path, repo_src: Path) -> Path:
    """Copy the HADDOCK3 package tree so B runs from another install path.

    Used both for Axis 1.3 (a different install path must not change any key)
    and for Axis 6.8/6.9 (a toppar file or CNS template genuinely changed).
    The copy is put on ``PYTHONPATH``, which precedes the editable install's
    ``.pth`` entry, so the subprocess imports this tree instead.
    """
    target = tmp / "install" / "src"
    target.mkdir(parents=True, exist_ok=True)
    package = target / "haddock"
    if not package.exists():
        shutil.copytree(repo_src / "haddock", package, symlinks=True)
    return target


def _op_install(spec: dict, context: Context) -> None:
    """Shadow the installation, optionally editing a file inside it.

    An edit is either the name of one of the editors below, or a
    ``{replace: [old, new]}`` mapping that names the substitution in the case
    itself. The second form exists so that a case which claims to change what
    a recipe *computes* can be read, and disbelieved, without opening this
    file.
    """
    from .systems import REPO_ROOT

    context.install = shadow_install(context.tmp, REPO_ROOT / "src")
    for relative, how in (spec.get("edit") or {}).items():
        path = context.install / "haddock" / relative
        if not path.is_file():
            raise PerturbationError(f"no installed file at {relative}")
        if isinstance(how, dict):
            wanted, replacement = how["replace"]
            _install_replace(path, wanted, replacement)
            context.notes.append(f"edited installed {relative}: {wanted!r} -> {replacement!r}")
            continue
        _INSTALL_EDITORS[how](path)
        context.notes.append(f"edited installed {relative} ({how})")


def _install_replace(path: Path, wanted: str, replacement: str) -> None:
    """Substitute one exact string, and fail if it is not there exactly once.

    A substitution that silently matches nothing would turn a case claiming a
    substantive change into a second copy of the case claiming an inert one --
    asserting the opposite verdict, and passing for the wrong reason.
    """
    text = path.read_text(encoding="utf-8")
    occurrences = text.count(wanted)
    if occurrences != 1:
        raise PerturbationError(
            f"{path.name} contains {occurrences} occurrences of {wanted!r}, "
            "so this edit is not the one the case describes"
        )
    path.write_text(text.replace(wanted, replacement, 1), encoding="utf-8")


def _install_append_comment(path: Path) -> None:
    """Append a CNS comment.  Changes the file's content, so every job that
    reads it must miss -- and every job that does not must still hit."""
    with open(path, "a", encoding="utf-8") as handle:
        handle.write("\n{ perturbed by the caching contract suite }\n")


def _install_change_parameter(path: Path) -> None:
    """Change a force-field number, which genuinely changes the science."""
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"(\d+\.\d{3})", text)
    if match is None:
        _install_append_comment(path)
        return
    value = f"{float(match.group(1)) + 0.001:.3f}"
    path.write_text(text[: match.start()] + value + text[match.end() :], "utf-8")


def _op_cns_wrapper(spec: dict, context: Context) -> None:
    """Point ``cns_exec`` at a shell wrapper around the real binary.

    Axis 6.10 needs a *different executable*, and only one CNS binary runs on
    this platform.  A wrapper is a different executable by every test the
    filesystem can apply -- different path, different bytes, different inode --
    and it must cost nothing, because the executable is the machine that
    evaluates a computation and not part of the computation itself.  It stands
    in for the case that cannot be built here: a second, independently
    installed HADDOCK3, whose binary shares no bytes with this one.
    """
    real_cns = _resolve_cns_exec(context)
    wrapper = context.tmp / "cns-wrapper.sh"
    wrapper.write_text(
        f'#!/bin/sh\nexec "{real_cns}" "$@"\n',
        encoding="utf-8",
    )
    wrapper.chmod(0o755)
    context.config.top["cns_exec"] = str(wrapper)


def _resolve_cns_exec(context: Context) -> Path:
    """Ask the installation under test where its CNS binary is.

    A subprocess rather than an import, for two reasons.  The suite never
    loads the package it is testing into its own process; and when a case has
    shadowed the installation, the binary that matters is *that* copy's, which
    an import here would not see.

    ``cns_exec`` is a documented configuration parameter, but its default
    resolves at install time and ``haddock3-cfg`` reports it as empty, so
    there is no public command that prints the path.
    """
    environment = dict(os.environ)
    if context.install is not None:
        existing = environment.get("PYTHONPATH", "")
        environment["PYTHONPATH"] = (
            f"{context.install}{os.pathsep}{existing}" if existing else str(context.install)
        )
    finished = subprocess.run(
        [sys.executable, "-c", "from haddock.core.defaults import cns_exec; print(cns_exec)"],
        capture_output=True,
        text=True,
        env=environment,
        check=False,
    )
    path = Path(finished.stdout.strip())
    if finished.returncode != 0 or not path.is_file():
        raise PerturbationError(
            "could not locate the CNS executable of the installation under "
            f"test: {finished.stderr.strip() or finished.stdout.strip()!r}"
        )
    return path


def _op_cwd(spec: dict, context: Context) -> None:
    """Invoke haddock3 from a different working directory (Axis 1.4, 8.7)."""
    target = Path(context.expand(spec["path"]))
    target.mkdir(parents=True, exist_ok=True)
    context.cwd = target


def _op_note(spec: dict, context: Context) -> None:
    context.notes.append(str(spec.get("text", "")))


def _op_cosmetic(spec: dict, context: Context) -> None:
    """Reorder top-level parameters, which changes only the config's spelling.

    Comments and blank lines are added at write time by :func:`cosmetic_text`;
    parameter order is changed here, because the renderer emits ``top`` in
    insertion order.
    """
    context.config.top = dict(reversed(list(context.config.top.items())))
    context.notes.append("top-level parameters reordered; comments added")


def cosmetic_text(text: str) -> str:
    """Add comments and blank lines to a rendered config."""
    lines = ["# a comment, which the science does not read", ""]
    for line in text.splitlines():
        lines.append(line)
        if line.startswith("["):
            lines.append("    # another comment")
    return "\n".join(lines) + "\n\n\n"


def _op_explicit_default(spec: dict, context: Context) -> None:
    """Write a parameter out explicitly at its own documented default value.

    The config text changes; the computation does not.  A key built by copying
    the module parameters and deleting known-irrelevant ones would notice;
    a key built by explicit inclusion cannot.
    """
    import yaml  # noqa: PLC0415

    from .systems import REPO_ROOT  # noqa: PLC0415

    module, parameter = spec["module"], spec["param"]
    matches = list((REPO_ROOT / "src" / "haddock" / "modules").glob(f"*/{module}/defaults.yaml"))
    if not matches:
        raise PerturbationError(f"no defaults.yaml for module {module!r}")
    defaults = yaml.safe_load(matches[0].read_text(encoding="utf-8"))
    if parameter not in defaults:
        raise PerturbationError(f"{module} has no parameter {parameter!r}")
    value = defaults[parameter]["default"]
    context.config.find(module, spec.get("occurrence", 0)).params[parameter] = value
    context.notes.append(f"{module}.{parameter} written explicitly as {value!r}")


def _op_restraints_archive(spec: dict, context: Context) -> None:
    """Give a refinement module one restraint table *per model*.

    HADDOCK3 assigns the members of a ``.tgz`` to models **by position**.  That
    makes this the sharpest probe on Axis 5:

    * ``identical`` members -- every model reads the same bytes, so a rank
      shift changes no read-set and the job MUST-HIT.  A miss here is a
      locator (the index) leaking into the key.
    * ``distinct`` members -- a rank shift changes which table a model is
      refined against, so the read-set genuinely differs and the job
      MUST-MISS.

    Nothing has to detect intent: the restraint assignment is not metadata
    *about* the job, it is one of the job's inputs, and it enters the key as
    content.
    """
    import tarfile  # noqa: PLC0415

    source = _target(spec, context)
    count = int(spec.get("count", 8))
    flavour = spec.get("members", "identical")
    folder = context.data_dir / f"per-model-{flavour}"
    folder.mkdir(parents=True, exist_ok=True)
    text = source.read_text(encoding="utf-8")
    for index in range(count):
        member = folder / f"restraints_{index + 1}.tbl"
        if flavour == "identical":
            member.write_text(text, encoding="utf-8")
        elif flavour == "distinct":
            member.write_text(_shift_first_distance(text, index), encoding="utf-8")
        else:
            raise PerturbationError(f"unknown archive flavour {flavour!r}")
    archive = context.data_dir / f"per-model-{flavour}.tgz"
    with tarfile.open(archive, "w:gz") as handle:
        for member in sorted(folder.iterdir()):
            handle.add(member, arcname=member.name)
    context.data[archive.name] = archive
    for module in spec.get("modules", [spec.get("module", "flexref")]):
        context.config.find(module).params["ambig_fname"] = str(archive)
    context.notes.append(f"{count} {flavour} per-model restraint tables")


def _shift_first_distance(text: str, index: int) -> str:
    match = re.search(r"(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)", text)
    if match is None:
        raise PerturbationError("restraint table has no distance triple")
    replaced = (
        f"{float(match.group(1)) + 0.01 * index:.3f} "
        f"{match.group(2)} {match.group(3)}"
    )
    return text[: match.start()] + replaced + text[match.end() :]


_EDITORS = {
    "coordinate": _edit_coordinate,
    "remark": _edit_remark,
    "whitespace": _edit_whitespace,
    "segid": _edit_segid,
    "restraint": _edit_restraint,
    "gzip": _edit_gzip,
    "models": _edit_models,
}

_INSTALL_EDITORS = {
    "append-comment": _install_append_comment,
    "change-parameter": _install_change_parameter,
}

_OPS = {
    "top": _op_top,
    "module": _op_module,
    "insert": _op_insert,
    "remove": _op_remove,
    "duplicate": _op_duplicate,
    "swap_steps": _op_swap_steps,
    "molecules": _op_molecules,
    "edit_input": _op_edit_input,
    "copy_input": _op_copy_input,
    "move_input": _op_move_input,
    "env": _op_env,
    "install": _op_install,
    "cns_wrapper": _op_cns_wrapper,
    "note": _op_note,
    "cwd": _op_cwd,
    "cosmetic": _op_cosmetic,
    "explicit_default": _op_explicit_default,
    "restraints_archive": _op_restraints_archive,
}
