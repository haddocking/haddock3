# The CNS caching test set

HADDOCK3 can reuse CNS results from a previous run:

```bash
haddock3 my-workflow.cfg --cache some-previous-run
```

Anything the previous run already computed, and that your new workflow needs
unchanged, is taken from it instead of being computed again. Everything else
runs normally. This directory is the test set for that feature.

**If you are reviewing this and have half an hour, start in
[`by-hand/`](by-hand/).** It is six things a HADDOCK3 user actually does, each
as two ordinary runs you drive yourself, with no test framework anywhere in
sight. Everything else in this directory is the automated version of the same
claims.

---

## The one idea the whole test set is built on

A cached result is delivered as a **hardlink**: the new run directory and the
old one end up pointing at *the same file on disk*. So the question "was this
reused or recomputed?" is answered by asking the filesystem whether two names
are the same file. That is all the measurement is.

```bash
python by-hand/check_reuse.py run1 run2
```

```
  reused      0_topoaa/e2aP_1F3G_haddock.pdb
  reused      0_topoaa/e2aP_1F3G_haddock.psf
  reused      1_rigidbody/rigidbody_1.pdb   <- 1_rigidbody/rigidbody_4.pdb
  RECOMPUTED  1_rigidbody/rigidbody_5.pdb
```

`by-hand/check_reuse.py` imports nothing but the Python standard library and
works on any two HADDOCK3 run directories, whether or not they came from this
test set. It is about eighty commented lines, and it is worth reading: it is
the entire measurement.

If you would rather not take a script's word for it, `ls -i` prints inode
numbers, and two names with the same inode number are the same file:

```bash
ls -i run1/*/rigidbody_1.pdb run2/*/rigidbody_1.pdb
```

Same number, reused. Different numbers, recomputed.

Notice that the report says **which** old file each result came from, not just
that it came from somewhere. That matters more than it looks: the dangerous
failure is not "it recomputed something it could have reused" (that only wastes
time) but "it reused *the wrong model*" (that silently gives you wrong
science). Naming the source is how the test set catches the second kind.

---

## What must be reused, and what must not

Three questions settle any case, in this order.

1. **Did the content of anything CNS reads change?** Then it must be
   recomputed — the coordinates, the restraints, the force field, a module's
   CNS script.
2. **Did a *role* change** — the same file now used as a different kind of
   input, or two molecules swapped so each takes the other's place? Then it
   must be recomputed. Same bytes in a different role is a different
   calculation.
3. **Did only a *name* change** — a directory, a filename, a step number, a
   model's rank, where HADDOCK3 is installed? Then it must be reused. None of
   those are part of what a result *is*.

There is one deliberate exception, and it is a policy rather than a
consequence: **the CNS executable itself is not part of what a result is.** It
is the machine that evaluates the calculation, not an input to it. Keying on
it would mean no two HADDOCK3 installations could ever share a result — nobody
compiles or downloads the same binary twice — and sharing between
installations is the principal use: a lab-wide store, results published with a
paper, a workstation result reused on a cluster. It would not buy safety
either: it cannot detect a CNS build that computes *different* answers, only
stop two builds that agree from being recognised as agreeing.

The cost is real and worth knowing: if you point `cns_exec` at a genuinely
different engine and expect different numbers, this feature will not notice.
Case `axis6.10-cns-executable-changed` states that policy and asserts it.
Whether it should be possible to opt into a stricter rule is a fair question
for you to raise.

The third category is where the value of the feature lives, and the second is
what stops the third from going too far. "Names never matter" is nearly right
and dangerously wrong: a model's rank is just a name, but which molecule is
molecule 1 is not.

The two failure modes are not symmetric, and the test set is weighted
accordingly:

| Failure | Costs you |
|---|---|
| something was recomputed that could have been reused | time |
| something was reused that should have been recomputed | **wrong results** |

---

## Reading a test case

Every case is a small YAML entry in [`cases/`](cases/), written to be read.
Here is one, complete:

```yaml
- case: axis4.1-sampling-increased
  taxonomy: "4.1"
  title: sampling raised from 4 to 6
  why: >
    The first four jobs are the same four jobs. Only the two new ones are new
    work. If all six miss, the schedule is not prefix-stable and raising
    sampling throws away every model already computed -- the single most
    expensive way for this feature to fail.
  base: tiny
  sources: [tiny]
  perturb:
    - {op: module, module: rigidbody, set: {sampling: 6}}
  expect:
    default: hit
    paths:
      "1_rigidbody/rigidbody_5.pdb": null
      "1_rigidbody/rigidbody_6.pdb": null
```

* `base` — which workflow the "before" run used.
* `sources` — which previous run(s) are offered to `--cache`.
* `perturb` — what is different about the "after" run. Usually one line.
* `expect` — `hit` means "must be reused", `miss` (written as `null` for a
  single file) means "must be recomputed", `auto` means "reused wherever the
  old run genuinely holds that job".

To see any case in plain English:

```bash
python try_case.py axis4.1-sampling-increased
```

To read the whole matrix at once, without opening a single YAML file:

```bash
python list_cases.py               # everything
python list_cases.py --axis 5      # one axis
python list_cases.py --markdown    # as a document
```

---

## Running any case, by hand, without any of this

The six scenarios in [`by-hand/`](by-hand/) are already written out this way.
For any *other* case, `try_case.py` will write one out for you:

```bash
python try_case.py axis4.1-sampling-increased --export ~/cache-demo
cd ~/cache-demo && ./run.sh
```

The exported folder contains `A.cfg`, `B.cfg`, a `data/` folder and a
`run.sh` — nothing else, and nothing that imports this test suite. `run.sh`
runs A, runs B with `--cache runA`, and prints the reuse table.

Then edit `B.cfg` and run it again. That is the intended way to use this: change
a parameter and watch which files stop being reused. If you can predict it, the
feature is behaving; if you cannot, either the feature or the documentation
needs work, and both are worth reporting.

Start with the everyday cases:

```bash
python try_case.py --list --common
```

| Case | What a user would be doing |
|---|---|
| `axis1.1-different-run-dir-name` | rerunning the same workflow in a new directory |
| `axis2.1-insert-non-cns-module-upstream` | adding a `[caprieval]` step to an existing workflow |
| `axis3.1-ncores` | running on a machine with a different number of cores |
| `axis4.1-sampling-increased` | "give me more models" |
| `axis4.2-sampling-decreased` | "that was too many models" |
| `axis5.2-clustering-cutoff-moved` | retuning the clustering and re-refining |
| `axis6.1-one-molecule-changed` | fixing something in an input structure |
| `axis6.6-restraint-content-changed` | changing the restraints |
| `axis7.1-parameter-read-by-one-module` | changing a scoring weight |
| `axis9.2a-sigint-mid-step` | the run died overnight; restarting it |
| `axis11.2-two-caches-disjoint-coverage` | drawing on two earlier runs at once |
| `axis12.1-artifact-deleted` | a file went missing from the old run |

---

## Running the whole set

```bash
# once, and it takes a while: build the previous runs the tests reuse from
python build_corpus.py

# then, and this part is fast, because almost everything is reused
pytest end-to-end_tests/caching/
```

Three groups run, in this order:

1. **Phase 0** (`test_phase0_*.py`) checks the *instrument*, not the feature:
   that a reused file really is detectable as one, that each perturbation
   really makes the change its case describes, and that the search for "which
   old result should this one have come from?" goes by content rather than by
   filename. If it fails, none of the rest means anything, so Phase 2 refuses
   to run. This part is written for developers and is not worth your time as a
   reviewer.
2. **Phase 1** (`build_corpus.py`) produces the previous runs everything else
   reuses from. It is a script, not a test. See [corpus.md](corpus.md).
3. **Phase 2** (`test_phase2_cases.py`) is the test set proper — the
   `cases/*.yaml` files. This is the part worth reviewing.

Useful flags:

```bash
pytest end-to-end_tests/caching/ -k axis5        # one axis
pytest end-to-end_tests/caching/ --no-timing     # skip the stopwatch checks
pytest end-to-end_tests/caching/ --ignore-phase0 # run anyway, unvalidated
```

`--no-timing` is worth knowing about. Most of what the suite checks is pure
filesystem identity and is true on any machine at any speed. A few checks also
time the run, and those assume a reasonably quick, unloaded machine with an
SSD. **A timing failure on a busy or slow machine may mean nothing. A reuse
failure means something wherever it appears.**

---

## What is deliberately not tested here

Twenty-four cases are recorded with a `skip:` and a reason instead of being
run. They are not oversights and they are not passes — they are the boundary of
what you can find out by running `haddock3` and looking at the results.

Most of them need to see the *key*: the internal description of a job that the
cache compares. Two jobs that must differ but produce the same key is the worst
possible bug, and from outside you can only catch it when it happens to change
an answer. A companion test set that inspects keys directly is the intended
place for those; each skipped case says so and says why.

```bash
grep -A3 'skip:' cases/*.yaml     # read them all
```

---

## What this suite is coupled to

The implementation is expected to be rewritten. A suite phrased in terms of a
record format or an internal API would die with it, so the coupling is kept to
a minimum and written down here rather than left to be discovered.

**The contract under test is two things:**

| | |
|---|---|
| `haddock3 B.cfg --cache OLD-RUN-DIR` | repeatable; coverage is the union |
| `HADDOCK_CACHE_HARDLINK` | `1` force link and fail loudly, `0` force copy, unset best-effort |

Everything else is observed from outside: which files a run directory
contains, which of them share an inode with a file in a source run, what bytes
they hold, and how long the run took. The suite does not import the caching
implementation, call into it, monkeypatch it, or read the records it keeps in
order to decide any verdict.

**Three declared dependencies sit outside that line:**

1. **HADDOCK3 itself** — module names, step-folder naming, and the installed
   tree layout that Axis 1.3 and Axis 6.8/6.9 perturb. Also the public
   per-step `io.json`, which `runio.py` reads to learn *which input model* a
   refinement job consumed. That is the Axis 5 oracle and cannot be derived
   from outside; `io.json` is HADDOCK3's documented data-flow record and what
   `haddock3-traceback` consumes, not part of the cache.
2. **[`cachesuite/record_format.py`](cachesuite/record_format.py)** — the only
   implementation-coupled file. Six fixtures need it, because Axis 11.11–11.13
   and 10.6 are *about* malformed records and cannot be built without writing
   one. It is opt-in: if it goes stale, those six report themselves
   unavailable with a pointer to the file, and the other twenty-six fixtures
   and every other case are unaffected. **To port this suite, update that file
   and nothing else.**
3. **Two error messages** — `HADDOCK_CACHE_HARDLINK` must be named when an
   invalid value is rejected, and `--cache` must be named when it is refused
   outside local mode. Both assert the contract that a refusal says what was
   refused, not any particular wording.

---

## Layout

| | |
|---|---|
| `by-hand/` | **six scenarios to drive yourself**, with no framework at all |
| `by-hand/check_reuse.py` | which files were reused, for any two run directories |
| `cases/*.yaml` | **the test set** — one file per topic, written to be read |
| `list_cases.py` | render the whole matrix as prose, for reviewing what is claimed |
| `try_case.py` | explain, run, or export any single case |
| `build_corpus.py` | Phase 1: build the previous runs (see `corpus.md`) |
| `test_phase2_cases.py` | the runner that executes `cases/*.yaml` |
| `test_phase0_hardlink.py` | Phase 0: is a reused file detectable as one? |
| `test_phase0_fixtures.py` | Phase 0: does each fixture make the change its case claims? |
| `test_phase0_oracle.py` | Phase 0: does the expected-source search identify a job by content? |
| `test_axis0_determinism.py` | does CNS give the same answer twice? |
| `test_scope_boundaries.py` | what is deliberately *not* cached |
| `cachesuite/` | the machinery; you should not need to read it |
| `cachesuite/thresholds.py` | the timing assumptions, all in one file |
| `corpus/` | generated; not in git |

Why this suite exists, who reviews it and when, is
[`caching-publication-plan.md`](../../caching-publication-plan.md). What must
be reused and why, and how that is observed, is this directory: the axis files
in [`cases/`](cases/) each open by saying what their axis is, and every case
states its own reasoning.
