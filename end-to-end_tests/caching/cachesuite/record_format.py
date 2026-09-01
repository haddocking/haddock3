"""The only file in this suite that knows how a run records its results.

**The contract this suite tests is two things and nothing else:**

* ``haddock3 B.cfg --cache OLD-RUN-DIR``, repeatable, and
* the ``HADDOCK_CACHE_HARDLINK`` environment variable.

Everything else is observed from the outside: what files a run directory
contains, which of them share an inode with a file in a source run, what bytes
they hold, and how long the run took.  That is deliberate.  The implementation
is expected to be rewritten, and a suite phrased in terms of a record format,
a file layout or an internal API dies with it.

Four items in the taxonomy cannot be expressed that way, because they are
*about* the record format rather than about caching behaviour:

* 11.11 -- a truncated, blank or wrong-arity record
* 11.12 -- a record key that is not checksum-shaped
* 11.13 / 10.5 -- a duplicated record that agrees
* 10.6 -- two records claiming one key with different results

Each is a parsing question with nothing to do with artifacts, and none can be
constructed without writing a malformed record.  So the coupling is confined
here, and it is **opt-in**: if this adapter no longer matches the
implementation, :func:`available` returns ``False``, those fixtures are
recorded as unavailable with a pointer to this file, and every other case in
the suite is unaffected.

**To port this suite to a different implementation, update this file. Nothing
else in the suite needs to change.**

The format described below is the one the current implementation writes: a
``CACHE`` file at the root of a run directory, one tab-separated record per
line, four fields -- job checksum, result checksum, PDB path, PSF path -- with
the two paths relative to the run directory.
"""

from __future__ import annotations

import re
from pathlib import Path

#: Where a run records what it computed, relative to the run directory.
STORE_NAME = "CACHE"

#: How one record is spelled.
SEPARATOR = "\t"
FIELD_COUNT = 4
KEY_FIELD = 0
RESULT_FIELD = 1
#: Fields holding a run-relative artifact path, in output order.
PATH_FIELDS = (2, 3)

_CHECKSUM = re.compile(r"^[0-9a-f]{64}$")


def store_path(run_dir: Path) -> Path:
    """The record store of a run."""
    return Path(run_dir) / STORE_NAME


def available(run_dir: Path) -> bool:
    """Whether this adapter still recognises how ``run_dir`` records results.

    Checked rather than assumed, so that a rewritten implementation produces a
    clear "the record-format adapter is stale" message instead of a fixture
    that is silently wrong.
    """
    path = store_path(run_dir)
    if not path.is_file():
        return False
    try:
        lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
    except OSError:
        return False
    if not lines:
        return False
    return all(
        len(line.split(SEPARATOR)) == FIELD_COUNT
        and _CHECKSUM.fullmatch(line.split(SEPARATOR)[KEY_FIELD])
        for line in lines
    )


def why_unavailable(run_dir: Path) -> str:
    """An actionable message for a fixture this adapter can no longer build."""
    return (
        f"the record-format adapter does not recognise {store_path(run_dir)}. "
        "This is a suite-maintenance issue, not a finding: update "
        "cachesuite/record_format.py to describe how this implementation "
        "records results. Only the malformed-record fixtures depend on it."
    )


def read(run_dir: Path) -> list[list[str]]:
    """Every record, as a list of fields."""
    text = store_path(run_dir).read_text(encoding="utf-8")
    return [line.split(SEPARATOR) for line in text.splitlines() if line]


def write(run_dir: Path, records: list[list[str]]) -> None:
    """Replace the record store with ``records``."""
    store_path(run_dir).write_text(
        "".join(SEPARATOR.join(fields) + "\n" for fields in records), encoding="utf-8"
    )


def artifact_step(record: list[str]) -> str:
    """The step folder a record's artifacts live in, or ``""``."""
    if len(record) <= PATH_FIELDS[0]:
        return ""
    return record[PATH_FIELDS[0]].split("/")[0]


#: The malformed states Axis 11.11-11.13 and 10.6 require.
CORRUPTIONS = (
    "truncated",
    "blank-line",
    "wrong-arity",
    "non-checksum-key",
    "duplicate-agreeing",
    "duplicate-conflicting",
)


def corrupt(run_dir: Path, how: str) -> None:
    """Make the record store malformed in a specific, named way."""
    path = store_path(run_dir)
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()
    if how == "truncated":
        path.write_text(
            "\n".join(lines[:-1]) + "\n" + lines[-1][:20], encoding="utf-8"
        )
    elif how == "blank-line":
        path.write_text("\n".join(lines[:1] + [""] + lines[1:]) + "\n", encoding="utf-8")
    elif how == "wrong-arity":
        broken = SEPARATOR.join(lines[0].split(SEPARATOR)[: FIELD_COUNT - 2])
        path.write_text("\n".join([broken] + lines[1:]) + "\n", encoding="utf-8")
    elif how == "non-checksum-key":
        fields = lines[0].split(SEPARATOR)
        fields[KEY_FIELD] = "not-a-checksum"
        path.write_text(
            "\n".join([SEPARATOR.join(fields)] + lines[1:]) + "\n", encoding="utf-8"
        )
    elif how == "duplicate-agreeing":
        path.write_text(text + lines[0] + "\n", encoding="utf-8")
    elif how == "duplicate-conflicting":
        fields = lines[0].split(SEPARATOR)
        other = "0" * 63 + "1"
        fields[RESULT_FIELD] = other if fields[RESULT_FIELD] != other else "0" * 64
        path.write_text(text + SEPARATOR.join(fields) + "\n", encoding="utf-8")
    else:  # pragma: no cover - programming error
        raise ValueError(f"unknown corruption {how!r}")
