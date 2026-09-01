"""Shared library for the CNS caching contract-compliance suite.

**The contract under test is exactly two things:**

* ``haddock3 B.cfg --cache OLD-RUN-DIR``, repeatable, and
* the ``HADDOCK_CACHE_HARDLINK`` environment variable.

The suite observes the cache only through those, and through what an ordinary
run leaves on disk: which files exist, which of them share an inode with a
file in a source run, what bytes they hold, and how long the run took.  It
does not import the caching implementation, monkeypatch it, call into it, or
read the records it keeps.  The implementation is expected to be rewritten,
and a suite phrased in terms of a record format or an internal API dies with
it.

Two dependencies sit outside that line and are declared rather than hidden:

* **HADDOCK3 itself** -- module names, step-folder naming, the installed tree
  layout that Axis 1.3 and Axis 6.8/6.9 perturb, and the public per-step
  ``io.json`` that ``runio.py`` reads to learn which model a job consumed.
  ``io.json`` is HADDOCK3's documented data-flow record between modules and is
  what ``haddock3-traceback`` consumes; it is not part of the cache.
* **``record_format.py``** -- the one file that knows how results are
  recorded, needed only to *construct* the malformed-record fixtures that Axis
  11.11-11.13 and 10.6 are about.  Nothing reads it to decide a verdict, and
  if it goes stale those fixtures report themselves unavailable while the rest
  of the suite carries on.
"""
