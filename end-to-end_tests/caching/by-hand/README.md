# Try the cache by hand

Six things a HADDOCK3 user actually does, each as two ordinary runs you drive
yourself. No test framework, no fixtures, no harness — just `haddock3`, and a
small script that prints which results were reused.

Every one of these is also an automated case in [`../cases/`](../cases/). If
you want to check a claim the automated suite makes, this is where to check it
by hand.

Run everything from **this directory**. The configs reference the standard
protein–protein example data by relative path, so nothing needs copying.

```bash
cd end-to-end_tests/caching/by-hand
```

## The one thing to understand first

When HADDOCK3 reuses a result it does **not** copy the file. It makes a
*hardlink*: a second name for bytes that are already on disk. That is how you
can tell reuse from repetition:

```bash
ls -i run-a/1_rigidbody/rigidbody_1.pdb run-b/1_rigidbody/rigidbody_1.pdb
```

Same number on both lines → the same bytes → the result was reused.
Different numbers → it was computed again.

Comparing the *contents* would tell you nothing, because HADDOCK3 is
reproducible: a recomputed model has the same contents as the original. Only
the inode number distinguishes "skipped the work" from "did the work again".

`check_reuse.py` does exactly this for every result file at once:

```bash
python check_reuse.py OLD-RUN NEW-RUN
```

It is about eighty lines, all of them commented. Read it — it is the whole
measurement.

---

## 1. I renamed my run directory

*Nothing about the science changed. Everything should be reused.*

```bash
haddock3 01-renamed-run/a.cfg
haddock3 01-renamed-run/b.cfg --cache run-a
python check_reuse.py run-a a-completely-different-name
```

**Expect:** every result reused, and the second run finishing in seconds
rather than the minute the first one took.

A directory name is a label, not part of the calculation. If anything is
recomputed here, HADDOCK3 has put a path inside the thing it uses to
recognise work it has already done — and then *no* other reuse in this
document is possible either.

---

## 2. I want more models

*The most common reason to rerun anything.*

```bash
haddock3 02-more-models/a.cfg            # sampling = 10
haddock3 02-more-models/b.cfg --cache run-a   # sampling = 15
python check_reuse.py run-a run-b
```

**Expect:** `rigidbody_1` … `rigidbody_10` reused, `rigidbody_11` …
`rigidbody_15` recomputed. The topologies reused as well.

This is the property worth having: asking for more models must not throw away
the ones you already have. If all fifteen are recomputed, raising `sampling`
costs full price every time.

---

## 3. I edited my restraints

*Something real changed. Exactly what reads it should be recomputed.*

First make an edited copy of the restraints:

```bash
cp ../../../examples/docking-protein-protein/data/e2a-hpr_air.tbl edited_air.tbl
# now change one distance in edited_air.tbl -- any number on any "assign" line
haddock3 03-changed-restraints/a.cfg
haddock3 03-changed-restraints/b.cfg --cache run-a
python check_reuse.py run-a run-b
```

**Expect:** the topologies (`0_topoaa/…_haddock.pdb` and `.psf`) reused,
because building a topology does not read the restraints; every
`1_rigidbody/rigidbody_*.pdb` recomputed, because docking does.

This is the boundary that matters in both directions. If the topologies are
recomputed, the cache is being needlessly cautious. If the docking is
**reused**, it is serving you models that were computed against your old
restraints while telling you they came from your new ones — and nothing in the
output would say so.

That is the failure this whole test suite exists to prevent, and it is the one
worth checking yourself.

---

## 4. I changed the clustering

*The case where caching is worth the most, and the hardest one to get right.*

```bash
haddock3 04-different-clustering/a.cfg          # min_population = 2
haddock3 04-different-clustering/b.cfg --cache run-a   # min_population = 1
python check_reuse.py run-a run-b
```

**Expect:** most of the `flexref` results reused — and look closely at the
`<-` arrows in the output. A model refined as `flexref_3.pdb` in the new run
may have been `flexref_1.pdb` in the old one.

That is the point. Clustering *chooses* structures; it does not change them. A
model that survives into refinement is the same model whatever cluster it
lands in, whatever rank it gets, and whatever the file ends up being called.
Refinement is the expensive step, so reusing it across a clustering change is
where most of the saving is.

If everything is recomputed here, the cache is recognising work by file name
or by position rather than by what the file actually contains — which means
touching a clustering parameter costs you the entire refinement stage.

---

## 5. My run died overnight

```bash
haddock3 05-crashed-run/a.cfg
#   ... wait until it is partway through the rigidbody step,
#   then press Ctrl-C
haddock3 05-crashed-run/b.cfg --cache run-a
python check_reuse.py run-a run-b
```

**Expect:** the topologies and however many `rigidbody_*.pdb` files the first
run had finished are reused; the rest are computed.

Two things are being claimed at once. Recovery is per **job**, not per step:
a step that was half finished still contains completed, correct models, and
they are reusable. And the job that was *in flight* when you interrupted it —
whose file may exist but be half written — must **not** be reused. Look for it
in the `RECOMPUTED` list.

---

## 6. Something happened to my old run

*What if the cache is damaged?*

```bash
haddock3 06-damaged-cache/a.cfg
cp -r run-a run-a-damaged
# break one result, in whichever way you like:
truncate -s 500 run-a-damaged/1_rigidbody/rigidbody_1.pdb
haddock3 06-damaged-cache/b.cfg --cache run-a-damaged
python check_reuse.py run-a-damaged run-b
```

**Expect:** `rigidbody_1.pdb` recomputed, everything else reused. The run
should succeed normally — a damaged cache entry is a reason to redo one job,
not a reason to fail.

Try the other ways of breaking it too; they should all behave the same way:

```bash
rm run-a-damaged/1_rigidbody/rigidbody_1.pdb                  # deleted
chmod 000 run-a-damaged/1_rigidbody/rigidbody_1.pdb           # unreadable
gzip run-a-damaged/1_rigidbody/rigidbody_1.pdb                # compressed
```

The last one is different from the others: compressing a file changes its
bytes but not its contents, so it should still be **reused**. Compression is
a storage decision, not a scientific one.

---

## Changing these examples

They are ordinary config files. Edit `sampling`, swap in your own molecules,
add a `[caprieval]`, point them at a different system — then rerun the pair
and read `check_reuse.py`'s output. The question it answers is always the
same:

> given what I changed, is exactly the right work being redone?

Two answers are wrong, and only one of them is obvious. Too much recomputed
wastes your time and you will notice. Too much reused gives you results
computed from inputs you no longer have, and **nothing in the output will tell
you** — which is why the automated suite spends most of its effort on that
direction.

## If nothing is ever reused

Check, in this order:

1. Was the old run produced by a HADDOCK3 that records its results? Every
   completed run does, with no extra flags — but a run made before the feature
   existed, or one that died before finishing its first job, has nothing to
   offer.
2. Is `mode = "local"` in both configs? `--cache` is refused in batch modes
   rather than silently ignored.
3. Are both runs on the **same filesystem**? A hardlink cannot cross one, so
   results are copied instead, and copies are indistinguishable from
   recomputation to `check_reuse.py`. Put both run directories in the same
   place and try again.
4. Set `HADDOCK_CACHE_HARDLINK=1` before running. Then any failure to link is
   reported as an error instead of quietly falling back to copying.

   *Note:* `HADDOCK_CACHE_HARDLINK` is part of what the test suite requires,
   and at the time of writing HADDOCK3 does not act on it yet — setting it
   changes nothing. `end-to-end_tests/caching/test_phase0_hardlink.py` is the
   check for that, and it currently fails. It matters for the automated
   suite, which needs linking to be guaranteed rather than best-effort; for
   working through this document by hand, keeping both runs on the same
   filesystem is enough.
