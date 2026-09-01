# The corpus: the "previous runs" the tests reuse from

Every test in this set is a *pair* of runs: an old one, and a new one that
tries to reuse the old one's results. The old runs are the corpus. They are
built once, by

```bash
python end-to-end_tests/caching/build_corpus.py
```

and then reused by all ~120 live cases. That is why the test set is expensive
once and cheap thereafter.

The corpus is a **build artifact**: regenerable, not committed to git, and
required before the tests will run. A fresh checkout cannot run Phase 2 until
this script has been run. It lives in `end-to-end_tests/caching/corpus/`, or
wherever `HADDOCK3_CACHE_CORPUS` points.

```bash
python build_corpus.py --list              # describe it without building
python build_corpus.py --resume            # continue an interrupted build
python build_corpus.py --only tiny         # just one base run
python build_corpus.py --skip-interrupted  # skip the least reliable fixtures
```

`--resume` matters: the manifest is written after *every* fixture, so a build
that is interrupted — or one you extend later by adding a base run — continues
rather than starting over.

---

## There is nothing clever in here

The generator does exactly two things:

1. runs `haddock3` on a config, and
2. edits the resulting directory afterwards.

That is the whole mechanism. There is no cache-population mode, because **every
ordinary HADDOCK3 run is already a usable cache source**, whether or not anyone
ever passes `--cache`. And the damaged and interrupted fixtures are made by
acting on real run directories from outside, which is what keeps the whole test
set free of mocks and stubs.

The six malformed-record fixtures are the one exception: they have to write a
record the implementation will reject, so they have to know how records are
spelled. Everything they know is in
[`cachesuite/record_format.py`](cachesuite/record_format.py), the only
implementation-coupled file in the suite. If it goes stale, those six report
themselves unavailable and the other twenty-six carry on.

Once built, every fixture is made **read-only**. A test that wrote into one
would invalidate every measurement after it.

---

## The base runs

Eight ordinary HADDOCK3 runs. Their configs are written out to
`corpus/configs/*.cfg` and are perfectly readable; the workflows themselves are
defined in [`cachesuite/systems.py`](cachesuite/systems.py).

| Name | System | Workflow | Why it exists |
|---|---|---|---|
| `tiny` | protein–glycan | `topoaa`, `rigidbody` (4 models) | The workhorse. Small enough that a change invalidating *everything* costs eight recomputations rather than hundreds — which is what makes the global cases (force field changed, seed changed) affordable at all. |
| `tiny-wide` | protein–glycan | same, 40 models | Identical to `tiny` except for the number of jobs. Subtracting the two all-hit rerun times cancels startup and leaves the per-job cost of a cache hit, which is how the performance requirement gets measured instead of guessed. |
| `tiny-topoaa` | protein–glycan | `topoaa` only | Half of a pair of caches with no overlap; also a workflow that is a strict prefix of another. |
| `refine` | protein–glycan | `topoaa`, `rigidbody`, `caprieval`, `rmsdmatrix`, `flexref`, `emref`, `mdref` | The refinement job shapes, on a small system, so they do not have to be paid for at protein–protein scale. Its two analysis modules are what make "swap two independent modules" expressible. |
| `pp` | protein–protein (e2A–HPr) | the standard docking workflow via `seletop`, 40 models | Ten steps — deliberately one below the point where HADDOCK3 renames `9_x` to `09_x`, so a single added module crosses that boundary. |
| `pp-cluster` | protein–protein (e2A–HPr) | the same 40 models, selected via `clustfcc` + `seletopclusts` | The workhorse for selection and ranking. Shares its sampling stage with `pp` exactly, so a structure selected by both routes is the same refinement job reached two different ways. |
| `scoring` | two complexes | `topoaa`, `emscoring`, `mdscoring` | Scoring rather than sampling: no seed, no restraints. A shape with fewer moving parts is where a shortcut is most likely to hide. |
| `cg` | protein–protein + shape | `topoaa`, `topocg`, `rigidbody`, `seletop`, `flexref`, `cgtoaa`, `emref` | The coarse-grained shapes, reachable no other way. |

Each base run is also **rerun immediately against itself**, with every job
hitting. That rerun does double duty: it is the measurement of how long a run
takes when nothing is computed (which every timing check is built on), and it
is the build's own sanity check — a base run that cannot serve itself is not a
usable fixture, and saying so here is far clearer than a hundred confusing
failures later.

---

## The damaged fixtures

Copies of `tiny` — content copies, so they share no files with the original —
with one thing broken. No CNS runs, so these are the cheapest fixtures in the
corpus and cover most of what can go wrong with an old run.

Each lesion targets a *single* file and leaves the rest intact, so each case
asserts a boundary: one result must be recomputed, the other seven must still
be reused. "It refused to use anything" would be safe and useless.

| Fixture | What was done to it | What the new run must do |
|---|---|---|
| `deleted` | one output file removed | recompute that one |
| `modified` | one output changed in place, same length | recompute that one |
| `truncated` | one output cut in half | recompute that one |
| `same-size` | one output overwritten with junk of the same size | recompute that one |
| `replaced-dir` | one output replaced by a directory | recompute that one |
| `replaced-symlink` | one output replaced by a link to a *different* result | recompute that one — this is the nastiest: the link resolves, the bytes are a valid PDB, and they are the wrong model |
| `dangling-symlink` | one output replaced by a broken link | recompute that one |
| `unreadable` | one output `chmod 000` | recompute that one |
| `compressed` | every output gzipped, as `clean = true` would leave them | reuse everything — compression is packaging, not identity |
| `no-store` | records kept, every output deleted | recompute everything, without erroring |
| `poisoned` | records kept, every output's bytes replaced | recompute everything, and not be poisoned by it |
| `rigidbody-only` | the topology half removed entirely | supply the other half from a second cache |
| `records-*` | the record file truncated, blanked, given a short line, given a bad key, duplicated, made self-contradictory | refuse clearly, rather than guess |
| `moved` | the whole run copied to a deeper path under a new name | reuse everything |
| `crossfs` | the whole run copied to another filesystem | reuse everything, by copying rather than linking |

---

## The interrupted fixtures

Real runs, started and then killed at a chosen moment: after the topology step
finished, or partway through the sampling step, with `SIGINT` (the polite
signal, which lets things flush) and with `SIGKILL` (which does not).

These are the least reproducible things in the corpus, because the kill lands
where it lands. So the tests that use them do not declare in advance what
should be reused — they read what the interrupted run *actually* contains and
require exactly that much to be reusable, and no more. A half-written file must
never be served.

This is the everyday case behind the whole axis: a run died overnight and you
restart it in the morning. Resume and reuse are the same operation.

---

## The manifest

`corpus/CORPUS.json` records, for every fixture: where it is, which workflow
made it, what it is for, every cacheable output it holds, how long it took to
build, and how long an all-hit rerun took. Nothing in it is an inode number —
those go stale the moment a fixture is rebuilt, so the tests always ask the
filesystem fresh.

If something could not be built, it is recorded as unusable with the reason,
and the cases that need it fail with that reason rather than with a missing
path.

---

## Cost

On an ordinary workstation the whole corpus is roughly twenty minutes to half
an hour of real CNS time; on a slow or busy machine, longer. The `--only` and
`--resume` flags exist so you never have to pay it twice.
