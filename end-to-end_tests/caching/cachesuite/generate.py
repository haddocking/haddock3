"""Phase 1 -- build the ``OLD-RUN-DIR`` corpus.

**Scripts, not tests.**  Expensive (real CNS), run once, reused by the whole
suite.  What they do is deliberately ordinary: invoke ``haddock3`` on a config,
then act on the resulting directory from the outside.  There is no
cache-population mode to write, because every ordinary run is already a usable
cache source -- which is what makes a mock-free suite possible at all.
"""

from __future__ import annotations

import os
import shutil
import signal
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from . import damage, record_format
from .corpus import CorpusNotBuilt, Fixture, Manifest, corpus_root, freeze, thaw
from .harness import cacheable_artifacts, haddock3_executable, run_haddock3
from .systems import BASE_RUNS, SYSTEMS, build_config


_STARTED = time.monotonic()


#: Fixtures the rest of the corpus is derived from. If one of these fails to
#: build there is nothing to continue with.
REQUIRED_BASE_RUNS = frozenset({"tiny"})


def log(message: str) -> None:
    elapsed = time.monotonic() - _STARTED
    minutes, seconds = divmod(int(elapsed), 60)
    hours, minutes = divmod(minutes, 60)
    print(f"[corpus {hours:d}:{minutes:02d}:{seconds:02d}] {message}", flush=True)


def prepare_systems(root: Path) -> dict[str, dict[str, Path]]:
    """Copy every system's input files into the corpus.

    Perturbation cases edit inputs (Axis 6), and they must never edit the
    repository's ``examples/`` tree.
    """
    data: dict[str, dict[str, Path]] = {}
    for system in SYSTEMS.values():
        target = root / "systems" / system.name
        target.mkdir(parents=True, exist_ok=True)
        files = {}
        for name in system.files:
            source = system.source / name
            if not source.is_file():
                raise FileNotFoundError(f"{system.name}: missing input {source}")
            destination = target / name
            if not destination.exists():
                shutil.copy2(source, destination)
            files[name] = destination.resolve()
        data[system.name] = files
    return data


def build_base_run(
    root: Path,
    name: str,
    data: dict[str, dict[str, Path]],
) -> Fixture:
    """Run one base workflow and calibrate it against itself."""
    spec = next(entry for entry in BASE_RUNS if entry.name == name)
    run_dir = root / "base" / name
    if run_dir.exists():
        thaw(run_dir)
        shutil.rmtree(run_dir)
    config_path = root / "configs" / f"{name}.cfg"
    config = build_config(name, run_dir, data[spec.system])
    config.write(config_path)

    log(f"base run {name} ({spec.cost}) ...")
    result = run_haddock3(config_path, cwd=root)
    if not result.ok:
        raise RuntimeError(
            f"base run {name} failed with exit code {result.returncode}\n"
            f"{result.tail(40)}"
        )
    outputs = [artifact.relative for artifact in cacheable_artifacts(run_dir)]
    log(f"base run {name}: {result.duration:.1f} s, {len(outputs)} cacheable outputs")

    allhit = _calibrate(root, name, config, config_path, run_dir)
    freeze(run_dir)
    return Fixture(
        name=name,
        kind="base",
        path=str(run_dir.relative_to(root)),
        system=spec.system,
        purpose=spec.purpose,
        config=str(config_path.relative_to(root)),
        duration=result.duration,
        allhit=allhit,
        outputs=outputs,
    )


def _calibrate(
    root: Path,
    name: str,
    config,
    config_path: Path,
    run_dir: Path,
) -> float:
    """Measure ``overhead(config)`` by rerunning the config against itself.

    Every job hits, so what is left is exactly the uncacheable remainder:
    startup, analysis and clustering modules.  This is the only measured
    quantity in the timing model -- ``t_miss`` and ``t_hit`` are declared.

    It doubles as Phase 1's self-check: a base run that cannot serve itself is
    not a usable fixture, and saying so here is much clearer than a hundred
    Phase 2 failures later.
    """
    calibrate_dir = root / "calibrate" / name
    if calibrate_dir.exists():
        thaw(calibrate_dir)
        shutil.rmtree(calibrate_dir)
    calibrate_config = config.copy()
    calibrate_config.top["run_dir"] = str(calibrate_dir)
    calibrate_path = root / "configs" / f"{name}--allhit.cfg"
    calibrate_config.write(calibrate_path)

    result = run_haddock3(calibrate_path, cache=[run_dir], cwd=root)
    if not result.ok:
        raise RuntimeError(
            f"all-hit calibration of {name} failed: exit {result.returncode}\n"
            f"{result.tail(40)}"
        )
    log(f"  calibration (all-hit): {result.duration:.1f} s")
    freeze(calibrate_dir)
    return result.duration


# --------------------------------------------------------------------------
# Derived-by-damage fixtures
# --------------------------------------------------------------------------

#: ``(fixture name, base run, what it is for, mutation)``.
DAMAGE_RECIPES: tuple[tuple[str, str, str, Callable[[Path], None]], ...] = (
    (
        "deleted",
        "tiny",
        "12.1 -- one rigidbody artifact deleted; MUST-DEGRADE for that job only",
        lambda run: damage.delete(damage.pick(run, "rigidbody", 0)),
    ),
    (
        "modified",
        "tiny",
        "12.2 -- one artifact modified in place, same length",
        lambda run: damage.modify_in_place(damage.pick(run, "rigidbody", 0)),
    ),
    (
        "truncated",
        "tiny",
        "12.3 / 9.3 -- one artifact truncated, as a torn write would leave it",
        lambda run: damage.truncate(damage.pick(run, "rigidbody", 0)),
    ),
    (
        "same-size",
        "tiny",
        "12.4 -- one artifact replaced by a same-size file of other bytes",
        lambda run: damage.replace_same_size(damage.pick(run, "rigidbody", 0)),
    ),
    (
        "replaced-dir",
        "tiny",
        "12.5 -- one artifact replaced by a directory",
        lambda run: damage.replace_with_directory(damage.pick(run, "rigidbody", 0)),
    ),
    (
        "replaced-symlink",
        "tiny",
        "12.5 -- one artifact replaced by a symlink to a sibling artifact",
        lambda run: damage.replace_with_symlink(
            damage.pick(run, "rigidbody", 0), damage.pick(run, "rigidbody", 1)
        ),
    ),
    (
        "dangling-symlink",
        "tiny",
        "12.5 -- one artifact replaced by a dangling symlink",
        lambda run: damage.replace_with_dangling_symlink(
            damage.pick(run, "rigidbody", 0)
        ),
    ),
    (
        "unreadable",
        "tiny",
        "12.6 -- one artifact chmod 000",
        lambda run: damage.make_unreadable(damage.pick(run, "rigidbody", 0)),
    ),
    (
        "compressed",
        "tiny",
        "1.7 / 12.7 -- every artifact gzipped, as clean/postprocess would leave them",
        damage.compress_run,
    ),
    (
        "no-store",
        "tiny",
        "11.9 / 12.11 -- records complete, artifact store empty",
        damage.strip_artifacts,
    ),
    (
        "poisoned",
        "tiny",
        "12.10 -- every record present, every artifact holding the wrong bytes",
        damage.poison_store,
    ),
    (
        "rigidbody-only",
        "tiny",
        "11.2 -- half of a disjoint pair: the topology outputs are gone, so "
        "this run can serve the sampling jobs and nothing else",
        lambda run: damage.strip_module_artifacts(run, "topoaa"),
    ),
    (
        "records-truncated",
        "tiny",
        "11.11 -- the record store ends mid-record",
        lambda run: damage.corrupt_records(run, "truncated"),
    ),
    (
        "records-blank",
        "tiny",
        "11.11 -- the record store contains a blank line",
        lambda run: damage.corrupt_records(run, "blank-line"),
    ),
    (
        "records-arity",
        "tiny",
        "11.11 -- a record has the wrong number of fields",
        lambda run: damage.corrupt_records(run, "wrong-arity"),
    ),
    (
        "records-nonchecksum",
        "tiny",
        "11.12 -- a record key is not checksum-shaped",
        lambda run: damage.corrupt_records(run, "non-checksum-key"),
    ),
    (
        "records-duplicate",
        "tiny",
        "10.5 / 11.13 -- a duplicated record that agrees; must be accepted silently",
        lambda run: damage.corrupt_records(run, "duplicate-agreeing"),
    ),
    (
        "records-conflict",
        "tiny",
        "10.x -- two records claiming one key with different results",
        lambda run: damage.corrupt_records(run, "duplicate-conflicting"),
    ),
)


def build_damaged(root: Path, manifest: Manifest, skip=None) -> list[Fixture]:
    """Copy base runs and damage them.  Cheap: no CNS."""
    fixtures = []
    for name, base, purpose, mutate in DAMAGE_RECIPES:
        if skip is not None and skip(name):
            log(f"damaged fixture {name}: already built, skipping")
            continue
        base_fixture = manifest.require(base)
        target = root / "damaged" / name
        if target.exists():
            thaw(target)
        log(f"damaged fixture {name}")
        damage.copy_run(base_fixture.run_dir(root), target)
        mutate(target)
        outputs = [artifact.relative for artifact in cacheable_artifacts(target)]
        freeze(target)
        fixtures.append(
            Fixture(
                name=name,
                kind="damaged",
                path=str(target.relative_to(root)),
                system=base_fixture.system,
                purpose=purpose,
                outputs=outputs,
            )
        )
    return fixtures


def build_relocated(root: Path, manifest: Manifest, other_fs: Path | None) -> list[Fixture]:
    """Fixtures that move a run rather than damage it (Axes 1.2, 1.8, 12.12)."""
    fixtures = []
    base_fixture = manifest.require("tiny")

    moved = root / "moved" / "renamed-after-the-fact" / "deeply" / "nested"
    if moved.exists():
        thaw(moved)
    log("relocated fixture moved")
    damage.copy_run(base_fixture.run_dir(root), moved)
    freeze(moved)
    fixtures.append(
        Fixture(
            name="moved",
            kind="damaged",
            path=str(moved.relative_to(root)),
            system=base_fixture.system,
            purpose="1.2 / 1.5 -- source relocated and renamed after A finished",
            outputs=list(base_fixture.outputs),
        )
    )

    if other_fs is not None:
        elsewhere = other_fs / "haddock3-cache-suite" / "crossfs"
        if elsewhere.exists():
            thaw(elsewhere)
        log(f"relocated fixture crossfs on {other_fs}")
        elsewhere.parent.mkdir(parents=True, exist_ok=True)
        damage.copy_run(base_fixture.run_dir(root), elsewhere)
        freeze(elsewhere)
        fixtures.append(
            Fixture(
                name="crossfs",
                kind="damaged",
                path=str(elsewhere),
                system=base_fixture.system,
                purpose="1.8 / 12.12 -- source on another filesystem; forces copy",
                outputs=list(base_fixture.outputs),
                notes=["absolute path: outside the corpus root"],
            )
        )
    return fixtures


# --------------------------------------------------------------------------
# Interrupted runs
# --------------------------------------------------------------------------


def run_until_and_kill(
    config_path: Path,
    cwd: Path,
    run_dir: Path,
    ready: Callable[[Path], bool],
    sig: int,
    limit: float = 900.0,
) -> str:
    """Start ``haddock3``, wait for ``ready``, then signal the process group.

    Genuinely requires running and killing a real process, which makes this the
    least deterministic generator in the corpus and the likeliest source of
    flakes.  The kill point is reported so the fixture's actual content can be
    inspected afterwards -- the expected mapping for Axis 9 must be derived
    from what the interrupted fixture *contains*, not declared a priori.
    """
    command = [haddock3_executable(), str(config_path)]
    process = subprocess.Popen(
        command,
        cwd=str(cwd),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )
    deadline = time.monotonic() + limit
    reached = False
    try:
        while time.monotonic() < deadline:
            if process.poll() is not None:
                return "run finished before the kill point was reached"
            if ready(run_dir):
                reached = True
                break
            time.sleep(0.05)
        if not reached:
            return "kill point was never reached"
        os.killpg(os.getpgid(process.pid), sig)
        try:
            process.wait(timeout=120)
        except subprocess.TimeoutExpired:
            os.killpg(os.getpgid(process.pid), signal.SIGKILL)
            process.wait(timeout=60)
        return "killed at the requested point"
    finally:
        if process.poll() is None:  # pragma: no cover - defensive
            try:
                os.killpg(os.getpgid(process.pid), signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait(timeout=60)


def _after_topoaa(run_dir: Path) -> bool:
    return any(
        entry.name.partition("_")[2] == "rigidbody" for entry in _dirs(run_dir)
    )


def _some_rigidbody_models(run_dir: Path) -> bool:
    for entry in _dirs(run_dir):
        if entry.name.partition("_")[2] != "rigidbody":
            continue
        if len(list(entry.glob("rigidbody_*.pdb"))) >= 2:
            return True
    return False


def _dirs(run_dir: Path):
    if not run_dir.is_dir():
        return []
    return [entry for entry in run_dir.iterdir() if entry.is_dir()]


INTERRUPT_RECIPES = (
    (
        "interrupted-between-steps-int",
        "tiny",
        "9.1 / 9.7 -- SIGINT between steps, cleanly",
        _after_topoaa,
        signal.SIGINT,
    ),
    (
        "interrupted-between-steps-kill",
        "tiny",
        "9.1 / 9.7 -- SIGKILL between steps; less flushing than SIGINT",
        _after_topoaa,
        signal.SIGKILL,
    ),
    (
        "interrupted-mid-step-int",
        "tiny",
        "9.2 -- SIGINT mid-step, some jobs complete; recovery must be per-job",
        _some_rigidbody_models,
        signal.SIGINT,
    ),
    (
        "interrupted-mid-step-kill",
        "tiny",
        "9.2 / 9.6 -- SIGKILL mid-step; the record write may itself be torn",
        _some_rigidbody_models,
        signal.SIGKILL,
    ),
)


def build_interrupted(
    root: Path,
    data: dict[str, dict[str, Path]],
    skip=None,
) -> list[Fixture]:
    """Start real runs and kill them at chosen points."""
    fixtures = []
    for name, base, purpose, ready, sig in INTERRUPT_RECIPES:
        if skip is not None and skip(name):
            log(f"interrupted fixture {name}: already built, skipping")
            continue
        spec = next(entry for entry in BASE_RUNS if entry.name == base)
        run_dir = root / "interrupted" / name
        if run_dir.exists():
            thaw(run_dir)
            shutil.rmtree(run_dir)
        config_path = root / "configs" / f"{name}.cfg"
        config = build_config(base, run_dir, data[spec.system])
        # A bigger schedule makes the mid-step kill point reachable.
        if "mid-step" in name:
            config.find("rigidbody").params["sampling"] = 12
        config.write(config_path)
        log(f"interrupted fixture {name}")
        note = run_until_and_kill(config_path, root, run_dir, ready, sig)
        log(f"  {note}")
        if not run_dir.exists():
            fixtures.append(
                Fixture(
                    name=name,
                    kind="interrupted",
                    path="",
                    purpose=purpose,
                    notes=[note, "fixture unusable: no run directory was created"],
                )
            )
            continue
        outputs = [artifact.relative for artifact in cacheable_artifacts(run_dir)]
        records_present = record_format.available(run_dir)
        freeze(run_dir)
        fixtures.append(
            Fixture(
                name=name,
                kind="interrupted",
                path=str(run_dir.relative_to(root)),
                system=spec.system,
                purpose=purpose,
                config=str(config_path.relative_to(root)),
                outputs=outputs,
                notes=[note, f"records present: {records_present}"],
            )
        )
    return fixtures


# --------------------------------------------------------------------------
# Entry point
# --------------------------------------------------------------------------


def find_other_filesystem(root: Path) -> Path | None:
    """A writable directory on a different filesystem from ``root``.

    Without root privileges ``/dev/shm`` is usually a separate tmpfs.  Failing
    to find one is a **coverage hole**, not a neutral skip: it leaves the
    cross-filesystem cases unvalidated.
    """
    device = root.stat().st_dev
    for candidate in (Path("/dev/shm"), Path("/run/user") / str(os.getuid()), Path("/tmp")):
        try:
            if candidate.is_dir() and os.access(candidate, os.W_OK):
                if candidate.stat().st_dev != device:
                    return candidate
        except OSError:
            continue
    return None


def build(
    root: Path | None = None,
    *,
    only: tuple[str, ...] | None = None,
    skip_interrupted: bool = False,
    resume: bool = False,
) -> Manifest:
    """Build the whole corpus.

    The manifest is written after **every** fixture, not once at the end, so
    an interrupted build leaves a usable record of what it already has.  With
    ``resume``, a fixture already recorded and still present on disk is
    skipped -- which is what makes an eight-hour build restartable rather than
    all-or-nothing.
    """
    root = (root or corpus_root()).resolve()
    root.mkdir(parents=True, exist_ok=True)
    (root / ".gitignore").write_text("*\n", encoding="utf-8")

    data = prepare_systems(root)
    manifest = Manifest(
        generated_at=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        haddock3=haddock3_executable(),
    )
    # Always merge whatever is already recorded, so a partial build never
    # silently discards fixtures a previous invocation paid CNS time for.
    try:
        existing = Manifest.read(root)
    except CorpusNotBuilt:
        existing = None
    if existing is not None:
        manifest.fixtures.update(existing.fixtures)
        manifest.notes = list(existing.notes)
        log(f"{len(existing.fixtures)} fixtures already recorded")

    def done(name: str) -> bool:
        if not resume:
            return False
        fixture = manifest.fixtures.get(name)
        if fixture is None:
            return False
        if not fixture.path:
            return True  # recorded as unusable; do not retry automatically
        return fixture.run_dir(root).is_dir()

    total = len(BASE_RUNS)
    for index, spec in enumerate(BASE_RUNS, start=1):
        if only and spec.name not in only:
            continue
        if done(spec.name):
            log(f"[{index}/{total}] base run {spec.name}: already built, skipping")
            continue
        log(f"[{index}/{total}] base run {spec.name}")
        try:
            manifest.fixtures[spec.name] = build_base_run(root, spec.name, data)
        except Exception as error:  # noqa: BLE001 - the build must not stop here
            if spec.name in REQUIRED_BASE_RUNS:
                raise
            # A base run that will not build is a real finding, but it is a
            # finding about one job shape. Recording it and carrying on is
            # what keeps a long unattended build from losing everything after
            # it, and the cases that need it fail loudly with this reason
            # rather than mysteriously.
            log(f"  FAILED: {error}")
            manifest.notes.append(
                f"COVERAGE HOLE: base run {spec.name!r} could not be built, so "
                f"the job shapes it carries ({', '.join(spec.tags)}) are "
                "untested. First line of the error: "
                + str(error).splitlines()[0]
            )
            manifest.fixtures[spec.name] = Fixture(
                name=spec.name,
                kind="base",
                path="",
                system=spec.system,
                purpose=spec.purpose,
                notes=[f"build failed: {error}"],
            )
        manifest.write(root)

    if not only:
        for fixture in build_damaged(root, manifest, skip=done):
            manifest.fixtures[fixture.name] = fixture
            manifest.write(root)
        other_fs = find_other_filesystem(root)
        if other_fs is None:
            note = (
                "COVERAGE HOLE: no second filesystem found; Axes 1.8 and 12.12 "
                "are unvalidated."
            )
            if note not in manifest.notes:
                manifest.notes.append(note)
        for fixture in build_relocated(root, manifest, other_fs):
            manifest.fixtures[fixture.name] = fixture
            manifest.write(root)
        if not skip_interrupted:
            for fixture in build_interrupted(root, data, skip=done):
                manifest.fixtures[fixture.name] = fixture
                manifest.write(root)

    # Corpus hygiene, enforced once at the end as well as per fixture: every
    # source must be read-only for the whole of Phase 2. A test that writes
    # into a fixture invalidates every inode assertion that follows it.
    for fixture in manifest.fixtures.values():
        if not fixture.path:
            continue
        directory = fixture.run_dir(root)
        if directory.is_dir():
            freeze(directory)

    manifest.write(root)
    log(f"done: {len(manifest.fixtures)} fixtures")
    return manifest
