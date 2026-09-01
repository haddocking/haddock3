#!/usr/bin/env python3
"""Print which results in a run were reused from an older run.

    python check_reuse.py OLD-RUN NEW-RUN

No test framework, no dependencies beyond Python itself.  Run it after any

    haddock3 B.cfg --cache OLD-RUN

and it will tell you, file by file, whether the result was taken from the old
run or computed again.

How it knows
------------
When HADDOCK3 reuses a result it does not copy the file -- it makes a
**hardlink**, which is a second name for bytes that are already on disk.  Two
names for the same bytes have the same "inode number", and different bytes
always have different inode numbers.  So:

    reused      the file in the new run has the same inode number as a file
                in the old run
    recomputed  it does not

You can see exactly the same thing with ordinary shell tools:

    ls -i  old-run/1_rigidbody/rigidbody_1.pdb  new-run/1_rigidbody/rigidbody_1.pdb
    find old-run -samefile new-run/1_rigidbody/rigidbody_1.pdb

Why not just compare the file contents?  Because identical contents prove
nothing.  HADDOCK3 is meant to be reproducible, so a *recomputed* result has
the same contents as the original.  Only the inode number tells you whether
the work was actually skipped.

If you set HADDOCK_CACHE_HARDLINK=0, HADDOCK3 copies instead of linking, and
this script will report everything as recomputed -- correctly, because with
copies there is no way to tell from the files alone.
"""

import sys
from pathlib import Path

# HADDOCK3 modules that run CNS.  Only these produce cacheable results;
# everything else (caprieval, clustfcc, seletop, ...) is cheap analysis and is
# always recomputed, by design.
CNS_MODULES = {
    "topoaa", "topocg", "rigidbody", "flexref", "emref", "mdref",
    "emscoring", "mdscoring", "cgtoaa",
}


def results_of(run_dir):
    """Yield the result files a run produced, as (module, path)."""
    for folder in sorted(run_dir.iterdir()):
        if not folder.is_dir() or "_" not in folder.name:
            continue
        number, _, module = folder.name.partition("_")
        if not number.isdigit() or module not in CNS_MODULES:
            continue
        for path in sorted(folder.iterdir()):
            # Strip any .gz the run may have been compressed with.
            name = path.name[:-3] if path.name.endswith(".gz") else path.name
            if not name.endswith((".pdb", ".psf")):
                continue
            # A topology step writes the *input* ensemble members into its own
            # folder alongside its outputs. Only the "_haddock" files are
            # results.
            if module in ("topoaa", "topocg") and "_haddock" not in name:
                continue
            yield module, path


def inode(path):
    """The identity of the bytes: device plus inode number, never one alone."""
    info = path.stat()
    return (info.st_dev, info.st_ino)


def main(argv):
    if len(argv) != 3:
        print(__doc__)
        return 2
    old_run, new_run = Path(argv[1]).resolve(), Path(argv[2]).resolve()
    for run in (old_run, new_run):
        if not run.is_dir():
            print(f"not a run directory: {run}")
            return 2

    # Index every file in the old run by its inode, so a reused file can be
    # traced back to *which* old result it came from -- not merely that it
    # came from somewhere. That matters: a cache serving the right file is
    # correct, and a cache serving the wrong file is the worst thing it can do.
    old_by_inode = {}
    for path in old_run.rglob("*"):
        if path.is_file() and not path.is_symlink():
            old_by_inode.setdefault(inode(path), path)

    reused = recomputed = 0
    print(f"old run: {old_run}")
    print(f"new run: {new_run}\n")
    for module, path in results_of(new_run):
        source = old_by_inode.get(inode(path))
        shown = path.relative_to(new_run)
        if source is None:
            recomputed += 1
            print(f"  RECOMPUTED  {shown}")
        else:
            reused += 1
            origin = source.relative_to(old_run)
            note = "" if origin == shown else f"   <- {origin}"
            print(f"  reused      {shown}{note}")

    total = reused + recomputed
    print(f"\n{reused} of {total} results reused, {recomputed} recomputed.")
    if reused and recomputed:
        print("A mixture is normal: the cache is meant to reuse exactly the")
        print("jobs whose inputs did not change, and no others.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
