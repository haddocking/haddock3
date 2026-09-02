# Staged publication of the CNS caching feature

## Standing of this document

This note records **how** the CNS caching feature is being published and
reviewed, not what it does or how it works. It sits alongside
`caching-use-cases.md` (what must hit, what must miss, and why),
`caching-test-suite-plan.md` (how those verdicts are observed), and
`canonical-representation-consumption.md` (why the canonical form is not
executed in production).

## Why the work is staged

Caching is a **user-facing feature**, and the work behind it touches three
concerns that are best judged separately, by different people, against
different criteria.

**Usability.** Does the feature do what a HADDOCK user expects and wants — is
the reuse it offers the reuse they would ask for, and is it offered in a form
they can actually use? This is a domain judgement, made against behaviour, and
it requires no reading of the implementation. It is settled in Stage 2, where
the behaviour is written down as a readable contract, and confirmed in Stage 3,
where users exercise the feature as a black box.

**Code quality.** Is the implementation sound, maintainable, and a reasonable
thing to carry in HADDOCK3? This is a core developer judgement, made against
the diff. It is the content of Stage 3.

**Reproducibility.** Content-based caching rests on it entirely: reuse is only
meaningful if the same job, run again, produces the same result, and only safe
if a job's identity accounts for everything that can change that result.
Reproducibility is addressed twice — in Stage 1 independently of caching, as
fixes and a canonicalization library that stand on their own merits, and in
Stage 4 together with caching, where a job can be dumped and re-executed from
its declared dependencies alone.

Keeping these apart is what the staging buys. A user asked to review a diff
would be reviewing the wrong artifact; a core developer asked whether reusing a
refinement result across a changed cluster rank is scientifically acceptable
would be answering the wrong question. Reproducibility, meanwhile, is worth
having whether or not caching is ever merged, so it is not made to wait on a
review of caching.

The stages also order the reviews so that the promised behaviour is fixed
**before** the implementation that has to satisfy it exists.

## The stages

Each stage is the previous stage's branch with further work applied on top, so
the stages are cumulative branches rather than four independent pull requests.

| Stage | Adds | Concern | Reviewed by |
|---|---|---|---|
| 1 | Reproducibility fixes + canonicalization library | Reproducibility, independent of caching | Core devs |
| 2 | The caching test suite | Usability, as a written contract | HADDOCK users |
| 3 | The caching implementation | Usability in practice, then code quality | Users, then core devs |
| 4 | `seamless-run` job dumps | Reproducibility, together with caching | Core devs |

### Stage 1 — reproducibility and the canonicalization library

Branch: `bitwise-reproducibility`.

Two kinds of work, deliberately kept in separate commits:

- **General reproducibility fixes**, ten commits, each standing on its own:
  parallel scheduler results returned in submission order; deterministic
  coarse-grained topology generation and deterministic CNS seeding; a
  prefix-stable rigid-body sampling schedule; generated CNS input restricted to
  the parameters the recipes can actually read; run-volatile provenance
  normalized out of published PDB and PSF artifacts; complete output sets
  published atomically on the local, batch and grid paths; and three unrelated
  pre-existing defects corrected along the way (coarse-grain back-mapping input
  alignment, the unhonoured `sampling_factor` clamp, discarded grid worker
  exceptions). These contain no caching, and each is justifiable to a reviewer
  who never hears the word.
- **A canonicalization library**, two commits. It derives a stable canonical
  representation of a CNS transformation — the script with locators erased, plus
  the set of inputs bound to canonical pin names — from which a job identity can
  be computed, together with a committed golden canonical form for each of the
  nine CNS job shapes.

The library is **virtual** at this stage: nothing consumes its output. CNS
continues to execute the existing generated input through the existing
scheduler path, and the canonical form exists only to be checksummed. Making
CNS execute the canonical form directly was investigated and rejected; see
`canonical-representation-consumption.md` for the operational reasons.

Reviewed by core devs, on reproducibility grounds alone. The bar is that Stage 1
stands up **even if caching is never merged**: bitwise-reproducible CNS results
are worth having in their own right, and the canonicalization library is a
self-contained component. Stage 1 is not "caching, part one".

### Stage 2 — the behavioural contract

Stage 1 plus the caching test suite at `end-to-end_tests/caching/`.

The suite is nothing but ordinary invocations of
`haddock3 config.cfg --cache OLD-RUN-DIR`. Each case declares, per output file,
which earlier run's file it must be reused from, or that it must be recomputed.
There are no mocks, no test-only entry points, and no inspection of internals.

That constraint is what makes the user review possible. A case is legible to
someone who runs dockings and has never opened the source: it names a config
change, an earlier run, and the files that must or must not be recomputed as a
result.

Reviewed by **HADDOCK users, not core devs**. The questions put to them are
answerable from domain knowledge alone:

- Is anything declared a required hit that they would regard as a *different*
  calculation?
- Is anything declared a required miss that they would expect to be reused —
  and would therefore be surprised to pay for again?
- Is the feature offered in a shape they can use: opting in per run, naming
  earlier run directories directly, drawing on several of them at once?

The suite does not pass at Stage 2. There is nothing yet that could make it
pass. That is the point: freezing the expectations before the implementation
exists is what prevents the suite from being quietly written to describe
whatever the implementation happens to do.

Stage 2 exits when users agree the declared behaviour is the behaviour they
want.

### Stage 3 — the implementation

Stage 2 plus the caching implementation. Cache-key construction becomes the
canonicalizer's first production caller: the canonical mapping is computed
before cache lookup, a verified result is reused on a hit, and the existing CNS
path executes on a miss.

Stage 3 is reviewed twice, in order:

1. **The Stage 2 reviewers check it out and run it.** The test suite gives them
   a mechanical answer; their own configurations and molecular systems give
   them the answer that actually matters. This is the confirmation that the
   contract they signed off in Stage 2 has been met.
2. **On their green light, core devs review Stage 3 as code.**

The order is deliberate. Core developer review is the expensive, scarce
resource, and it is spent once the behaviour is settled rather than twice
across a moving specification.

### Stage 4 — dump to a workdir and run with `seamless-run`

The ability to dump a single CNS job to a working directory as a
self-contained `seamless-run` command.

Stage 4 returns to the reproducibility concern, this time bound to caching. A
dumped job is reproducible in the strongest available sense — re-executable
anywhere from its declared dependencies alone — and that same property is what
makes it a test of job identity.

Its primary purpose is therefore to test the Stage 1 canonicalization library
rigorously, which Stages 2 and 3 structurally cannot do. Their suite observes
the cache from outside: it can see whether a file was reused, never the key
that decided it. Two jobs that ought to differ but receive the same key are
invisible to it until some case happens to catch the resulting wrong hit.

A dump closes that gap, because it is self-contained by construction — it
carries precisely the declared dependencies and nothing else:

- running one in an isolated working directory **is** the strong test for an
  under-declared dependency set: an undeclared file is simply not there, so a
  silent key defect becomes a loud execution failure;
- the key and the pin table become directly inspectable, so canonical forms can
  be diffed and frozen against golden forms, with no CNS execution required;
- two dumps of one job that produce different results make a cache conflict
  *classifiable* — nondeterministic artifact versus incomplete key — rather
  than merely observable.

Stage 4 is sequenced last because it is an assurance instrument, not a
prerequisite: Stage 1 is reviewed on its own merits and already retains an
executable canonical-workspace probe for representative jobs, which Stage 4
generalizes into a suite. It does not change how HADDOCK executes CNS.
Executable canonical workspaces remain a test and audit instrument, never the
production execution path.

## What this staging is not

- **Not a schedule.** The stages are a review order, not a set of deadlines.
- **Not four independent pull requests.** Each stage is the previous branch plus
  work on top and is reviewed as a whole.
- **Not two reviews of the same thing.** The user review in Stages 2 and 3
  cannot substitute for core developer code review, and code review cannot
  substitute for it.

## Numbering note

`canonical-representation-consumption.md` was written under an earlier
numbering in which Stage 2 was the caching implementation. Under the scheme in
this document, that is Stage 3. Where that document says "the caching
implementation in Stage 2" or "Stage 2 must append a successful cache record
only after the complete declared output set has been published", read
**Stage 3**.
