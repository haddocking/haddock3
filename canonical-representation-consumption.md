# Consumption of the canonical CNS representation

## Standing of this document

This note records the investigation into whether CNS should execute the same
canonical representation that will be used to calculate CNS-job identities.
It belongs to the staged caching design discussion and captures a rejected
production architecture, the reasons for rejecting it, and a smaller related
change that remains desirable.

The conclusion is:

- the canonical representation should remain a **checksum-side identity
  representation**;
- normal CNS execution should continue to use the existing run and scheduler
  architecture;
- the caching implementation in Stage 2 should become the canonicalizer's
  production caller when it calculates a job key before cache lookup;
- executable canonical workspaces remain valuable as a test and audit
  instrument, but should not become the normal execution path; and
- incomplete CNS output publication should be solved independently using
  hidden partial output files and atomic rename.

## The question

The Stage 1 branch currently exposes canonicalization through
`CNSJob.canonical_mapping()`. Normal runs do not call it. `CNSJob.run()` sends
the original generated input to CNS, using the real input and output filenames,
and only normalizes selected output artifacts afterward.

This led to a reasonable question: instead of maintaining a canonical form only
for job identity, should CNS itself execute with canonical input and output
filenames?

Doing so would make the relationship between execution and identity especially
direct. The canonical script and canonical input bindings used to construct the
key would also be the script and bindings observed by CNS.

## Initially considered production architecture

The investigated design gave every CNS job a private execution directory. A job
would:

1. build its canonical mapping;
2. create a private job directory;
3. materialize job-specific dependencies under canonical names;
4. execute a canonical CNS input producing `canonical-output.pdb` and,
   optionally, `canonical-output.psf`;
5. validate and normalize the outputs; and
6. publish them under HADDOCK's user-facing filenames.

Only job-specific files need to be materialized. The CNS executable and the
installed module and toppar trees do not need to be copied for ordinary local
execution. Stable `MODULE:relative/path` and `TOPPAR:relative/path` references
can continue to resolve through the execution environment.

This distinction is important. A dependency has at least two relevant names:

- its canonical pin name, used in transformation identity; and
- the stable reference spelling CNS receives, which may be a `MODULE:` or
  `TOPPAR:` reference rather than a workspace path.

Conflating these names is one cause of the current `MODULE:module/...` rewrite
defect.

## Empirical feasibility

The basic execution model was tested with the real CNS executable and generated
HADDOCK inputs.

### Rigid-body job

A rigid-body job was run with only its two input PDB files and two input PSF
files materialized under canonical names. CNS, the module scripts, and toppar
remained in their installed locations.

The job completed successfully. After output normalization, its PDB was
byte-identical to the output of the original job.

### Topology job

A `topoaa` job was run with only its input PDB materialized under a canonical
name. Its canonical PDB and PSF outputs both completed successfully and, after
normalization, were byte-identical to the original outputs.

The probe also demonstrated that an executable canonical form would immediately
expose several current defects that a checksum-only form can hide:

- rewriting an input alias must not rewrite the CNS variable name on the
  assignment's left-hand side;
- `MODULE:relative` and `TOPPAR:relative` execution references must not acquire
  an extra `module/` or `toppar/` prefix;
- the current textual `canonical-count` replacement is not executable CNS; and
- output bindings must not be consumed by an input basename collision.

These findings establish that direct canonical execution is technically
possible after the known canonicalization defects are fixed.

## Why technical feasibility is not sufficient

HADDOCK commonly runs many jobs concurrently on HPC systems backed by a network
filesystem. A production design must account for filesystem metadata operations,
data movement, scheduler topology, local scratch lifetime, and many simultaneous
HADDOCK runs.

### Step-local workspaces

Putting one workspace per job in the step directory permits hardlinking inputs
and same-filesystem atomic rename of outputs. It also adds substantial metadata
traffic to the shared filesystem.

For 10,000 rigid-body jobs with two input molecules, a naive implementation can
add approximately:

- 10,000 workspace directories; and
- 40,000 input hardlinks for the two PDB/PSF pairs.

Creating and removing these entries is precisely the kind of workload that can
stress NFS, Lustre, GPFS, and similar metadata services, especially when several
HADDOCK runs do it concurrently.

This is a material regression from the current non-debug local path, where CNS
inputs may remain in memory and no per-job input tree is created.

### Node-local workspaces

Node-local scratch avoids shared-filesystem metadata pressure and makes CNS
working I/O local. However, it prevents hardlinking directly from network-hosted
inputs. A naive implementation would copy every input for every job.

That could be mitigated with a node-local content pool:

```text
shared filesystem -> one local copy per unique input per node
                              |
                              +-> hardlinks into private job directories
```

Checksumming, decompression, and local materialization would also need to be
fused so an input is not read from the network once for its checksum and again
for staging. Static module/toppar dependency scans and checksums would need
caching rather than being repeated for every job.

These optimizations are feasible, but together they constitute a new execution
and data-staging subsystem. They require decisions about:

- scheduler-provided scratch such as `SLURM_TMPDIR` and `TMPDIR`;
- scratch capacity and exhaustion;
- allocation- or node-scoped content pools;
- concurrent population and cleanup;
- cross-filesystem output publication;
- batch-job preparation on the submission host versus the compute node;
- MPI processes distributed across nodes;
- debug artifact retention; and
- grid payload construction, where invariant files really do need transport.

This is far beyond the intended responsibility of Stage 1.

### Scheduler integration

Local and MPI modes eventually call `CNSJob.run()`, so they could share a
central canonical executor. Batch mode bypasses that method and generates shell
commands that invoke CNS directly. It would require a preparation and
finalization protocol, including per-job status and output handling for
concatenated batch jobs.

Grid mode already creates remote workspaces, but it has its own input rewrite,
payload, and output-retrieval implementation. Converting it to the same model
would be another substantial refactor. Unlike ordinary local execution, a grid
payload must transport the executable and invariant dependencies because the
remote site cannot be assumed to contain the HADDOCK installation.

Thus direct consumption of the canonical form would not be a contained change
to `CNSJob.run()`. It would alter every CNS execution backend or make execution
semantics mode-dependent.

## Decision: do not execute the canonical representation in production

The proposed execution architecture is rejected for the staged caching work.
It is technically workable, but its operational cost and scope are not
proportionate to the identity problem it solves.

The intended separation is instead:

### Stage 1

- Produce a stable canonical representation of each CNS transformation.
- Normalize result artifacts so irrelevant output locators and volatile headers
  do not affect result identity.
- Test canonicalization against real generated CNS inputs and golden canonical
  forms.
- Keep an executable canonical-workspace synthesizer or equivalent companion
  test for representative jobs. This proves that the declared dependencies are
  complete and that the canonical form remains meaningful.
- Do not alter normal CNS input materialization or scheduler execution merely to
  make CNS observe the canonical pin names.

### Stage 2

- Make cache-key construction the canonicalizer's production caller.
- Calculate the canonical mapping before cache lookup.
- Reuse a verified result on a hit.
- Execute the existing CNS path on a miss.

The executed input and the checksum-side representation will therefore not be
byte-identical objects. Their equivalence is a contract maintained by focused
canonicalization tests, golden forms, adjacent MUST-MISS tests, and the
executable isolated-workspace companion suite.

This is the same deliberate separation already described in the caching design:
canonicalization erases irrelevant locators for identity without requiring a
change to the run-directory layout or normal CNS execution.

## Independent issue: incomplete CNS output publication

The workspace investigation highlighted a separate correctness problem that
does not require canonical input staging.

At present, CNS writes directly to the final HADDOCK output filename. If CNS or
its worker dies while writing, a truncated PDB can remain at that path.
`Persistent.is_present()` checks only whether the path exists, so later module
logic can treat that truncated file as a generated result. In addition,
`CNSJob.run()` currently does not use the subprocess return code as a success
condition.

This should be corrected independently.

### Proposed output-publication protocol

For each declared output, CNS writes to a hidden partial filename in the existing
step directory, retaining the logical suffix, for example:

```text
.rigidbody_1.partial.pdb
.molecule_haddock.partial.psf
```

Keeping `.pdb` or `.psf` as the final suffix matters because CNS recipes derive
auxiliary names from those suffixes.

After CNS terminates, the job must require all of the following before publishing
anything:

1. the subprocess exit status is acceptable;
2. no known CNS failure is present in its output;
3. every declared output exists; and
4. every declared output is non-empty.

The job then normalizes the partial PDB/PSF and publishes it with
same-filesystem `os.replace()` to the expected HADDOCK filename.

If CNS is interrupted while writing, only the hidden partial file is affected.
The public output remains absent. A subsequent attempt removes stale partial
files and treats the declared output set as a unit.

PDB and PSF cannot be made filesystem-atomic as a pair while retaining the flat
step-directory layout. They can nevertheless be validated as a pair before
publication, published only during job finalization, and cleared together before
a retry. Downstream modules do not start until the producing module finishes.
Stage 2 must append a successful cache record only after the complete declared
output set has been published.

### Filesystem impact

This targeted protocol does not create per-job directories, stage inputs, or
copy large data trees. Compared with current execution, it adds temporary output
names and one metadata rename per successfully published artifact. CNS still
writes the same amount of output data to the same filesystem.

Grid retrieval should use the analogous cross-filesystem-safe form: copy the
retrieved artifact to a hidden temporary file in the destination step directory,
verify and normalize it there, and then use `os.replace()`. An interrupted copy
cannot expose a truncated public output.

## Commit boundaries

The following concerns should remain separate:

1. Scheduler results returned in submission order. This is a pre-existing
   Scheduler correction whose importance was exposed by output normalization.
2. Canonicalization correctness: context-aware path rewriting, stable
   `MODULE:`/`TOPPAR:` references, output-binding validation, count handling,
   archive restraints, and indexed `cgtoaa` references.
3. Atomic CNS output publication through hidden partial outputs.
4. Coarse-grained random-number determinism.
5. CNS parameter inclusion and exclusion, including removal of `tolerance` and
   the recorded treatment of `log_level`.
6. Real generated-input golden canonical forms.
7. Stage 2 cache lookup and recording, which supplies the canonicalizer's normal
   production caller.

Direct production execution of the canonical representation should not be one
of these commits. It is intentionally not part of the staged caching feature.
