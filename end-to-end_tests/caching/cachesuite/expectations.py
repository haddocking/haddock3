"""Turn a case's declared verdicts into a per-artifact source mapping.

Every case specifies a **mapping**, not a verdict word::

    B output path  ->  expected source path,  or  None (must miss)

That is strictly more informative than an executed-vs-reused signal, because
it identifies *which* source entry was used.  A name-keyed or position-keyed
cache fails it loudly, and a key collision -- a hit from the *wrong* entry --
is caught by the same check.

**A source entry is identified by content, never by name.**  That is the whole
point of the axis this suite is built around, so the oracle may not take the
shortcut the feature is forbidden to take.  Asking "which file in the source
run has this name, in a step with this module and occurrence?" would identify
a job by its name; the feature identifies a job by what it reads, and a suite
that asked the first question would pass a cache that answered it -- and fail
a correct one.

So the expected source of an output is found by computing, for the job that
produced it, an identity made of:

* the **content** of everything the job started from -- the split ensemble
  member behind a topology job, the model behind a refinement or scoring job,
  the models of the combination behind a sampling job, in order; and
* the job's recorded **seed**, which is what separates two repeats of one job
  from each other and is therefore part of what the job *is*.

and then looking for a job in each source run with the same identity.  Three
consequences are worth stating, because each is a defect this oracle used to
have:

* a job whose inputs changed has **no** source entry, so a content edit comes
  out as a miss without the case having to declare one artifact at a time;
* two jobs that are byte-identical resolve to the same source entry, so a
  duplicate is a hit -- an added ensemble member that copies an existing one
  is not new work, and declaring it new work would be a test only an
  incorrect implementation could pass;
* a source may serve several targets, since after a duplicate is introduced
  two outputs correctly hardlink to one source file.

Nothing here reads the cache's own records, and nothing here re-derives what
CNS reads.  It reads ``io.json``, which is HADDOCK3's public data-flow record.
"""

from __future__ import annotations

import fnmatch
from dataclasses import dataclass
from pathlib import Path

from .config import CACHEABLE_MODULES, TOPOLOGY_MODULES
from .harness import (
    Artifact,
    cacheable_artifacts,
    content_checksum,
    exists_any,
    step_folders,
)
from .runio import ModelRef, step_inputs, step_outputs

#: Modules whose jobs consume exactly one model, named by the output's
#: ``ori_name``.
PER_MODEL_MODULES = frozenset(
    {"flexref", "emref", "mdref", "emscoring", "mdscoring", "cgtoaa"}
)

#: Modules whose jobs consume a *combination* of models -- one per component
#: of the complex -- so the job's identity is the combination, in order.
COMBINATION_MODULES = frozenset({"rigidbody"})


@dataclass(frozen=True)
class Expectation:
    """What one output of B must be."""

    relative: str
    module: str
    #: The source files this output may be a hardlink to.  Empty means the
    #: output must be recomputed.
    #:
    #: Usually there is exactly one, and then the assertion is as strong as it
    #: can be: *this* output came from *that* entry.  It holds more than one
    #: only when several sources genuinely hold the same job -- two caches
    #: that overlap and agree -- where reuse from either is correct and
    #: insisting on one would be testing an ordering the design does not
    #: promise.
    sources: tuple[Path, ...]
    #: How the expectation was derived, for the failure message.
    origin: str

    @property
    def must_miss(self) -> bool:
        return not self.sources

    @property
    def source(self) -> Path | None:
        """The single expected source, when there is only one."""
        return self.sources[0] if len(self.sources) == 1 else None

    def describe(self) -> str:
        if not self.sources:
            return "nothing (must be recomputed)"
        if len(self.sources) == 1:
            return str(self.sources[0])
        return " or ".join(str(source) for source in self.sources)


class ExpectationError(RuntimeError):
    """A case declares an expectation that cannot be resolved."""


def resolve(
    run_dir: Path,
    sources: dict[str, Path],
    spec: dict,
) -> list[Expectation]:
    """Build the expected mapping for every cacheable artifact in ``run_dir``."""
    default = spec.get("default", "hit")
    exact: dict[str, str | None] = spec.get("paths", {}) or {}
    patterns: list[dict] = spec.get("patterns", []) or []
    module_expect: dict[str, str] = spec.get("modules", {}) or {}

    expectations: list[Expectation] = []
    for artifact in cacheable_artifacts(run_dir):
        if artifact.relative in exact:
            declared = exact[artifact.relative]
            if declared is None:
                expectations.append(
                    Expectation(artifact.relative, artifact.module, (), "declared miss")
                )
            else:
                expectations.append(
                    Expectation(
                        artifact.relative,
                        artifact.module,
                        (_named_source(declared, sources),),
                        f"declared hit from {declared}",
                    )
                )
            continue

        verdict = _verdict_for(artifact, default, patterns, module_expect)
        if verdict == "ignore":
            # Deliberately unasserted. Used where a case's boundary is real
            # but this suite cannot know on which side a particular job falls
            # without re-deriving what CNS reads -- which would make the test
            # a restatement of the implementation. Recorded, not hidden.
            continue
        if verdict == "miss":
            expectations.append(
                Expectation(artifact.relative, artifact.module, (), "declared miss")
            )
            continue

        found = _find_sources(artifact, run_dir, sources)
        if verdict == "auto":
            # The source either holds a usable entry for this job or it does
            # not, and which it is follows from content, not from the case
            # author's bookkeeping.  This is what makes Axis 5 and Axis 9
            # expressible: the expected mapping for a reordered selection, or
            # for an interrupted fixture, must be *derived from what the
            # source actually contains* rather than declared a priori.
            expectations.append(
                Expectation(
                    artifact.relative,
                    artifact.module,
                    found,
                    "auto: a source holds this job" if found else "auto: none does",
                )
            )
            continue
        if not found:
            raise ExpectationError(
                f"case declares {artifact.relative} a hit, but no source run "
                "holds a job with the same inputs and seed; the case is "
                "under-specified, not the implementation at fault"
            )
        expectations.append(
            Expectation(artifact.relative, artifact.module, found, "same job by content")
        )
    return expectations


#: Provenance a finished HADDOCK3 run strips from every artifact it publishes.
#:
#: Their presence means the opposite of what it looks like.  An artifact
#: carrying them was written by CNS and never published: the run died between
#: CNS closing the file and HADDOCK3 normalizing it.  It is raw output, not a
#: result -- it still names its own filename, the step folders it read from,
#: and the wall-clock minute it was written -- so serving it would inject a
#: timestamp and a pair of step paths into the run that reused it.
_UNPUBLISHED_MARKERS = (b"\nREMARK FILENAME=", b"\nREMARK DATE:")


def plausible(path: Path) -> bool:
    """Whether a source artifact looks complete enough to be reusable.

    Used only by the ``auto`` verdict, and only to keep an interrupted or
    damaged fixture from being *expected* to serve a job it visibly cannot.
    A torn PDB is present on disk but has no terminating record; expecting a
    hit from it would turn correct MUST-DEGRADE behaviour into a test failure.

    Completeness is not enough on its own.  An interrupted run can leave a
    *whole* artifact that was never published, and that one is the more
    interesting case precisely because nothing about its length or its final
    record gives it away -- it is detected by the provenance a published
    artifact no longer has.
    """
    for candidate in (path, Path(f"{path}.gz"), Path(f"{path}.zst")):
        if not candidate.exists():
            continue
        try:
            if candidate.stat().st_size == 0:
                return False
            if candidate.suffix == ".pdb":
                data = candidate.read_bytes()
                if any(marker in data for marker in _UNPUBLISHED_MARKERS):
                    return False
                return b"END" in data[-4096:]
            return True
        except OSError:
            return False
    return False


def _verdict_for(
    artifact: Artifact,
    default: str,
    patterns: list[dict],
    module_expect: dict[str, str],
) -> str:
    for rule in patterns:
        if fnmatch.fnmatch(artifact.relative, rule["match"]):
            return rule["expect"]
    if artifact.module in module_expect:
        return module_expect[artifact.module]
    return default


def _named_source(declared: str, sources: dict[str, Path]) -> Path:
    name, _, relative = declared.partition(":")
    if not relative:
        if len(sources) != 1:
            raise ExpectationError(
                f"{declared!r} does not name a source, and this case has "
                f"{len(sources)} of them"
            )
        (root,) = sources.values()
        return root / name
    if name not in sources:
        raise ExpectationError(f"{declared!r} names an unknown source {name!r}")
    return sources[name] / relative


def _find_sources(
    artifact: Artifact,
    run_dir: Path,
    sources: dict[str, Path],
) -> tuple[Path, ...]:
    """The source entries that hold the job which produced ``artifact``.

    Empty when none of them does -- which is the answer whenever the job's
    inputs or seed changed, and is what makes a content edit come out as a
    miss without the case enumerating artifacts.
    """
    folder = run_dir / artifact.relative.rsplit("/", 1)[0]
    basename = artifact.relative.rsplit("/", 1)[-1]
    outputs = step_outputs(folder)
    position = _output_position(outputs, basename)
    if position is None:
        raise ExpectationError(
            f"{artifact.relative} is not declared in its step's io.json, so "
            "the oracle cannot say what job produced it"
        )
    wanted = _job_identity(run_dir, folder, artifact.module, outputs, position)

    found = []
    for root in sources.values():
        source_folder = _step_folder(root, artifact.module, artifact.occurrence)
        if source_folder is None:
            continue
        source_outputs = step_outputs(source_folder)
        for index in range(len(source_outputs)):
            try:
                identity = _job_identity(
                    root, source_folder, artifact.module, source_outputs, index
                )
            except (ExpectationError, OSError):
                # A damaged or truncated source cannot answer for this job.
                # That is a miss, not an oracle failure: Axis 9 and Axis 12
                # build exactly such sources on purpose.
                continue
            if identity != wanted:
                continue
            name = (
                source_outputs[index].file_name
                if artifact.kind == "pdb"
                else source_outputs[index].psf_name
            )
            if not name:
                continue
            candidate = source_folder / name
            if not plausible(candidate):
                # This source holds the right job but not a usable result --
                # an interrupted or damaged fixture. Keep looking, in case it
                # holds the same job twice, and otherwise leave this source
                # out: it cannot serve what it does not have.
                continue
            found.append(candidate)
            break
    return tuple(found)


def _step_folder(run_dir: Path, module: str, occurrence: int) -> Path | None:
    """The step of ``run_dir`` running ``module`` for the ``occurrence``-th time.

    The step *ordinal* is deliberately not part of the match: Axis 2 moves it,
    and a step is identified by which module it runs and which of that
    module's occurrences it is, never by where it sits in the workflow.
    """
    seen = 0
    for _index, name, folder in step_folders(run_dir):
        if name != module:
            continue
        if seen == occurrence:
            return folder
        seen += 1
    return None


def _output_position(outputs: list[ModelRef], basename: str) -> int | None:
    """Which declared output of the step this file is, coordinates or topology."""
    for index, model in enumerate(outputs):
        if basename in (model.file_name, model.psf_name):
            return index
    return None


def _job_identity(
    run_dir: Path,
    folder: Path,
    module: str,
    outputs: list[ModelRef],
    position: int,
) -> tuple[tuple[str, ...], int | None]:
    """What the job that produced ``outputs[position]`` was.

    The content of everything it started from, in order, together with its
    seed.  Two jobs with the same identity are the same computation, whatever
    they are called and wherever they sit.

    A job whose inputs cannot be located raises, rather than silently
    identifying itself by less than it reads: an oracle that quietly compared
    fewer inputs would report hits it cannot justify.
    """
    model = outputs[position]
    if module in TOPOLOGY_MODULES:
        reads = (_checksum(_topology_input(run_dir, folder, model)),)
    elif module in COMBINATION_MODULES:
        reads = tuple(
            _checksum(path) for path in _combination_inputs(run_dir, folder, model)
        )
    elif module in PER_MODEL_MODULES:
        reads = (_checksum(_per_model_input(run_dir, folder, model)),)
    else:
        raise ExpectationError(
            f"{module} is a cacheable module whose job shape this oracle does "
            "not know; classify it before asserting anything about it"
        )
    return (reads, model.seed)


def _checksum(path: Path) -> str:
    try:
        return content_checksum(path)
    except FileNotFoundError as error:
        raise ExpectationError(f"cannot read job input {path}") from error


def _topology_input(run_dir: Path, folder: Path, model: ModelRef) -> Path:
    """The single structure a topology job was built from.

    ``ori_name`` records it. Where it sits depends on how the job was reached:
    a topology step that is not the run's first one consumes models declared
    in its own ``io.json``; the run's first one consumes a molecule, which
    HADDOCK3 either splits into the step folder -- one file per ensemble
    member, which is exactly what makes a per-member verdict possible -- or,
    for a single-model molecule, leaves in the run's ``data`` folder.
    """
    if not model.ori_name:
        raise ExpectationError(
            f"{model.file_name} does not record the structure it was built "
            "from, so the oracle cannot identify its job"
        )
    for declared in step_inputs(folder):
        if declared.stem == model.ori_name:
            return _rebase(declared.path, run_dir)
    candidate = folder / f"{model.ori_name}.pdb"
    if exists_any(candidate):
        return candidate
    for found in sorted(Path(run_dir, "data").rglob(f"{model.ori_name}.pdb")):
        return found
    raise ExpectationError(
        f"cannot find the structure {model.ori_name!r} that {model.file_name} "
        "was built from"
    )


def _per_model_input(run_dir: Path, folder: Path, model: ModelRef) -> Path:
    """The one model a refinement, scoring or conversion job consumed.

    Named by ``ori_name`` and resolved to a file through the step's own input
    list, so the name is only ever a join key inside one step; what enters the
    identity is the file's content.
    """
    declared = {entry.file_name: entry for entry in step_inputs(folder)}
    entry = declared.get(model.ori_name or "")
    if entry is None:
        raise ExpectationError(
            f"{model.file_name} records no input model this step declares "
            f"(ori_name={model.ori_name!r})"
        )
    return _rebase(entry.path, run_dir)


def _combination_inputs(run_dir: Path, folder: Path, model: ModelRef) -> list[Path]:
    """The models a sampling job docked, in the order they were bound.

    A docked complex carries one topology per component, and the step's input
    list pairs each topology with the structure it belongs to.  Order matters:
    the same two molecules bound the other way round is a different docking
    job, not a renaming of this one.
    """
    by_topology = {
        entry.psf_name: entry for entry in step_inputs(folder) if entry.psf_name
    }
    paths = []
    for psf_name in model.topology_names:
        entry = by_topology.get(psf_name)
        if entry is None:
            raise ExpectationError(
                f"{model.file_name} was docked from a topology "
                f"{psf_name!r} that this step does not declare as an input"
            )
        paths.append(_rebase(entry.path, run_dir))
    if not paths:
        raise ExpectationError(
            f"{model.file_name} records no topology, so the oracle cannot say "
            "which models were docked to produce it"
        )
    return paths


def _rebase(recorded: Path, run_dir: Path) -> Path:
    """Re-root a path ``io.json`` recorded, onto the run it was found in.

    ``io.json`` stores absolute paths.  A fixture that was *copied* -- and the
    corpus copies several -- therefore records paths pointing back at the run
    it was copied from, which may since have been damaged or removed.  The
    step folder and file name are what identify the file; the prefix is a
    locator, exactly as everywhere else in this suite.
    """
    if recorded.is_absolute() and not recorded.exists():
        rebased = Path(run_dir) / recorded.parent.name / recorded.name
        if rebased.exists():
            return rebased
    return recorded


def summarize(expectations: list[Expectation]) -> tuple[dict[str, int], int]:
    """``(misses per module, number of hits)`` -- the Gate 2 bound's inputs."""
    misses: dict[str, int] = {}
    hits = 0
    for expectation in expectations:
        if expectation.must_miss:
            misses[expectation.module] = misses.get(expectation.module, 0) + 1
        else:
            hits += 1
    return misses, hits


def assert_cacheable_modules_present(run_dir: Path) -> None:
    """Fail loudly if a run produced no cacheable job at all."""
    if not any(
        module in CACHEABLE_MODULES for _index, module, _folder in step_folders(run_dir)
    ):
        raise ExpectationError(f"{run_dir} contains no cacheable CNS step")
