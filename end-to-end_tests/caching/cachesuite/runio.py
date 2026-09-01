"""Read HADDOCK3's public per-step ``io.json``.

This is the only HADDOCK3-internal format the suite reads, and it is
deliberately not part of the caching implementation: ``io.json`` is the
documented data-flow record between modules, and it is what
``haddock3-traceback`` consumes.  Whatever bookkeeping the cache keeps for
itself is never read, because it is the thing under test.

It is read for one purpose: to learn what each job read, so that the expected
source of its output can be identified by the *content* of those inputs rather
than by the output's rank or filename.  Three fields carry that:

* ``ori_name`` -- the input a per-model job consumed, and the split ensemble
  member a topology job was built from;
* ``topology`` -- which topologies a sampling job docked, and therefore which
  models the combination was made of;
* ``seed`` -- the job's random seed, which separates two repeats of one job
  from each other.

There is no way to derive any of this from the outside.

``io.json`` is jsonpickle output.  It is walked as plain JSON rather than
deserialized, so the suite does not depend on the ``haddock`` package it is
testing being importable at the same version.  One jsonpickle detail has to be
handled rather than ignored -- repeated objects are written as back-references
-- and :func:`_shared_topologies` explains why, and checks itself.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator

_PDBFILE = "haddock.libs.libontology.PDBFile"


@dataclass(frozen=True)
class ModelRef:
    """One ``PDBFile`` entry as recorded in a step's ``io.json``."""

    file_name: str
    directory: str
    psf_name: str | None
    psf_directory: str | None
    ori_name: str | None
    score: float | None
    seed: int | None
    #: Every topology this model carries, in the order recorded.  A model
    #: built from one molecule carries one; a docked complex carries one per
    #: component, and that tuple is what says which combination produced it.
    topology_names: tuple[str, ...] = ()

    @property
    def path(self) -> Path:
        return Path(self.directory) / self.file_name

    @property
    def stem(self) -> str:
        return Path(self.file_name).stem


def _walk(node: Any) -> Iterator[dict]:
    """Yield every ``PDBFile`` dict, in document order."""
    if isinstance(node, dict):
        if node.get("py/object") == _PDBFILE:
            yield node
            return
        for value in node.values():
            yield from _walk(value)
    elif isinstance(node, list):
        for item in node:
            yield from _walk(item)


def _as_ref(node: dict, shared: dict[int, dict] | None = None) -> ModelRef:
    topology = node.get("topology")
    psf_name = psf_directory = None
    topologies = _topology_entries(topology, shared or {})
    if len(topologies) == 1:
        psf_name = topologies[0].get("file_name")
        psf_directory = topologies[0].get("path")
    score = node.get("score")
    if isinstance(score, str) or score != score:  # NaN or a jsonpickle marker
        score = None
    return ModelRef(
        file_name=node.get("file_name", ""),
        directory=node.get("path", ""),
        psf_name=psf_name,
        psf_directory=psf_directory,
        ori_name=node.get("ori_name"),
        score=score,
        seed=node.get("seed"),
        topology_names=tuple(
            str(entry.get("file_name", "")) for entry in topologies
        ),
    )


def _topology_entries(topology: Any, shared: dict[int, dict] | None = None) -> list[dict]:
    """The topology records of one model, in the order recorded.

    ``topology`` is a single object for a model built from one molecule and a
    list for a docked complex, and jsonpickle wraps a list as
    ``{"py/tuple": [...]}`` or similar.  Both spellings are walked the same
    way, because what matters is which topologies the model carries and in
    what order.
    """
    shared = shared or {}
    found: list[dict] = []
    stack = [topology]
    while stack:
        node = stack.pop(0)
        if isinstance(node, dict):
            if "file_name" in node:
                found.append(node)
            elif "py/id" in node:
                resolved = shared.get(node["py/id"])
                if resolved is not None:
                    found.append(resolved)
            else:
                stack = list(node.values()) + stack
        elif isinstance(node, list):
            stack = list(node) + stack
    return found


def _read(step_dir: Path, key: str) -> list[ModelRef]:
    path = Path(step_dir) / "io.json"
    if not path.is_file():
        return []
    # ``io.json`` may contain bare NaN, which strict JSON forbids but the
    # standard library accepts.
    document = json.loads(path.read_text(encoding="utf-8"))
    shared = _shared_topologies(document)
    return [_as_ref(node, shared) for node in _walk(document.get(key, []))]


def _shared_topologies(document: Any) -> dict[int, dict]:
    """Resolve the back-references jsonpickle writes for repeated objects.

    A docked complex carries the very same topology objects its inputs carry,
    and jsonpickle writes the second and later appearances as ``{"py/id": n}``
    rather than repeating them.  So the one field that says which models a
    docking job combined is written as a set of pointers, and reading it means
    resolving them.

    The pointers are resolved without reproducing jsonpickle's numbering,
    which the file does not preserve -- the document is written with sorted
    keys, and the numbers were handed out in attribute order.  What does
    survive is that they are handed out in *encoding* order, so a topology
    appearing earlier in the document has a lower number than one appearing
    later.  Sorting the pointers and pairing them with the topologies written
    out in full, in document order, therefore recovers the mapping.

    It also checks itself: the two must be the same length.  A document where
    they are not is one this reasoning does not hold for, and the caller is
    told it cannot say what a job combined rather than being handed a guess.
    """
    inlined = [node for node in _walk_topologies(document) if "file_name" in node]
    referenced = sorted(
        {
            node["py/id"]
            for node in _walk_topologies(document)
            if "py/id" in node and "file_name" not in node
        }
    )
    if not referenced:
        return {}
    if len(referenced) != len(inlined):
        return {}
    return dict(zip(referenced, inlined))


def _walk_topologies(node: Any, inside: bool = False) -> Iterator[dict]:
    """Every value written under a ``topology`` key, in document order."""
    if isinstance(node, dict):
        if inside and ("file_name" in node or "py/id" in node):
            yield node
            return
        for key, value in node.items():
            yield from _walk_topologies(value, inside or key == "topology")
    elif isinstance(node, list):
        for item in node:
            yield from _walk_topologies(item, inside)


def step_outputs(step_dir: Path) -> list[ModelRef]:
    """Models this step produced, in the order it recorded them."""
    return _read(step_dir, "output")


def step_inputs(step_dir: Path) -> list[ModelRef]:
    """Models this step consumed, in the order it recorded them."""
    return _read(step_dir, "input")
