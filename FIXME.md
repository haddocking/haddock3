# FIXME — Code review findings

Identified by static code review of `src/haddock/`, then re-verified against the
current source. Each finding below records whether the original report held up,
a corrected severity, and whether it is a live bug or merely latent (no current
caller can trigger it). Style and test-coverage gaps are excluded.

**Verdict summary**

| # | Finding | Verified | Severity | Status |
|---|---------|----------|----------|--------|
| 1 | `randint` overshoots upper bound | yes | low | live, fix worthwhile |
| 3 | wrong dict key in string validation | yes | medium | live, fix worthwhile |
| 5 | malformed CNS block for list PSF | yes | low | latent |
| 2 | `open(None)` on missing `error_file` | yes | low | latent |
| 8 | `split_tasks([])` raises `ValueError` | yes | low | latent |
| 9 | double assignment to `cns_exec` | yes | trivial | cosmetic |
| 10 | vestigial `self.len` | yes | trivial | cosmetic |
| 11 | step-name underscore assumption | yes | trivial | fragile, safe today |
| 4 | unhelpful (not crashing) error msg | partial | trivial | mischaracterized |
| 7 | rigidbody sampling truncation | partial | trivial | mostly by-design |
| 6 | `lru_cache` on mutable dict | no | — | not a live bug |
| 12 | missing energy term → 0 | n/a | — | dropped, by design |

---

## 1. [OPEN — fix worthwhile] `randint` overshoots the upper bound

**File:** `src/haddock/libs/libmath.py:27`

```python
# current (wrong)
return int(self() * (upper_limit + 1)) + lower_limit
```

`self()` is called with no arguments, so it always returns `uniform(0, 1)` —
the `lower_limit`/`upper_limit` passed to `randint` never reach it. The result
therefore spans `[lower_limit, lower_limit + upper_limit]` instead of
`[lower_limit, upper_limit]`.

```python
# corrected
return lower_limit + int(self() * (upper_limit - lower_limit + 1))
```

**Callers:** `libcns.py:310` and `libcns.py:497`, both `RND.randint(100, 99999)`,
which currently yield `[100, 100099]`.

**Impact:** Low. These are CNS random seeds; CNS accepts any integer, so there
is no runtime failure — only a slightly shifted seed range. Worth fixing for
correctness, not urgency.

---

## 3. [OPEN — fix worthwhile] Wrong dict key in string-length validation

**File:** `src/haddock/gear/prepare_run.py:752`

```python
elif "minchars" in param.keys() and "maxchars" in param.keys():
    if (nbchars := len(val)) < param["minchars"] or nbchars > param["maxchars"]:
        _desc = "under" if nbchars < param["minitems"] else "exceeding"
        #                                  ^^^^^^^^^^^ should be param["minchars"]
```

This branch is reached only for string parameters that define `minchars`/
`maxchars` and **not** `minitems`. When such a parameter fails its length check,
`param["minitems"]` raises `KeyError`, masking the real validation message.

**Impact:** Medium — the only genuinely live correctness bug in this list,
though it fires solely on a validation failure path.

**Fix:** `param["minchars"]`.

---

## 5. [LATENT] Malformed CNS `structure...end` block when PSF is a list

**File:** `src/haddock/libs/libcns.py:287-297`

```python
if psf_input:
    # if isinstance(psf_input, str):           <- half-finished refactor
    input_str += f"structure{linesep}"
    input_str += f"  @@{psf_input}{linesep}"   # if a list, emits its repr
    input_str += f"end{linesep}"
    input_str += f"coor @@{pdb_input}{linesep}"
    if isinstance(psf_input, list):
        input_str += f"structure{linesep}"      # second, correct block
        for psf in psf_input:
            input_str += f"  @@{psf}{linesep}"
        input_str += f"end{linesep}"
```

If `psf_input` is a list, the first block emits the Python list repr
(`@@['a.psf', 'b.psf']`) — invalid CNS — followed by the correct block. The
commented-out `isinstance` check confirms this is an unfinished refactor.

**Status:** Latent. No current caller passes a list, so it never triggers.

**Fix:** Collapse into a single block that branches on `isinstance(psf_input, list)`.

---

## 2. [LATENT] `open(None, ...)` when `error_file` is omitted

**File:** `src/haddock/libs/libsubprocess.py:249`

```python
with open(self.error_file, "wb+") as errf:   # TypeError if error_file is None
```

`error_file` defaults to `None`. On a CNS error this would raise `TypeError`
instead of a useful diagnostic.

**Status:** Latent. Every `CNSJob` construction in the codebase passes
`err_fname` as the third positional argument, so `error_file` is never `None`
in practice. Severity downgraded from the original HIGH.

**Fix (defensive):** Guard with `if self.error_file is not None:` before opening.

---

## 8. [LATENT] `split_tasks([])` raises `ValueError`

**File:** `src/haddock/libs/libparallel.py:21`

```python
n = math.ceil(len(lst) / n)   # 0 when lst is empty → range(0, 0, 0)
```

An empty task list makes `n == 0`, and `range(0, 0, 0)` raises
`ValueError: range() arg 3 must not be zero`.

**Status:** Latent. Modules always submit at least one job.

**Fix:** Add `if not lst: return` at the top. Note `split_tasks` is a
**generator**, so use a bare `return`, not `return []`.

---

## 9. [COSMETIC] Double assignment to `self.params["cns_exec"]`

**File:** `src/haddock/modules/base_cns_module.py:112-114`

```python
self.params["cns_exec"] = shutil.copyfile(_cns_exec, new_cns)  # immediately overwritten
shutil.copystat(_cns_exec, new_cns)
self.params["cns_exec"] = Path("..", Path(_cns_exec).name)
```

The first assignment is dead. Harmless; clarity only.

---

## 10. [COSMETIC] Vestigial `self.len`

**File:** `src/haddock/libs/libontology.py:115`

```python
self.len = score   # score is NaN by default; .len is never read
```

Confirmed: `self.len` is the only occurrence of that attribute in the entire
tree — written on every `PDBFile`, never read. Dead code.

---

## 11. [COSMETIC] Step-folder parsing assumes no underscore in module names

**File:** `src/haddock/modules/analysis/caprieval/capri.py:69`

```python
sel_steps = [step for step in sel_steps if step.count("_") == 1]
```

Folder names are `<N>_<module_name>`; this silently drops any step whose module
name contains an underscore. Safe today (no such module), fragile for the future.

---

## 4. [MISCHARACTERIZED] Unhelpful `inject_in_modules` error message

**File:** `src/haddock/gear/prepare_run.py:1066-1068`

```python
raise ValueError(
    "key {key!r} already in {module!r} parameters. " "Can't inject."
)
```

The original report claimed this raises `NameError`. It does **not** — the
string is not an f-string, so `{key!r}`/`{module!r}` are printed literally and
`module` is never evaluated. The only real defect is an uninformative message
on a rarely-hit error path.

**Fix:** `raise ValueError(f"key {key!r} already in module parameters. Can't inject.")`

---

## 7. [MOSTLY BY DESIGN] Rigidbody sampling truncation

**File:** `src/haddock/modules/sampling/rigidbody/__init__.py:217`

```python
sampling_factor = int(self.params["sampling"] / len(models_to_dock))
if sampling_factor < 1:
    self.finish_with_error(...)
```

The original report said "no warning is emitted." In fact the `< 1` case is
already handled with `finish_with_error`. The only gap is that integer division
rounds the effective sample count down (e.g. 1000 / 3 → 999) with no notice —
inherent to sampling an integer number of copies per combination.

**Optional:** Log the effective total so the user isn't surprised by 999 vs 1000.

---

## 6. [NOT A LIVE BUG] `lru_cache` on `_read_defaults`

**File:** `src/haddock/gear/prepare_run.py:162`

The original report claimed the cached dict gets mutated, bleeding parameters
between steps. Re-checking disproves this:

- None of the four callers (`prepare_run.py:564, 646, 1251, 1292`) mutate the
  returned dict; they read keys/values only.
- The cache key is the raw per-step name (`caprieval.1` vs `caprieval.2` are
  distinct keys), so two `caprieval` steps would not share a cache entry anyway.

At most a theoretical footgun if a future caller mutates the result; not a
current defect. No action required beyond awareness.

---

## 12. [DROPPED] Missing energy term silently contributes 0 to HADDOCK score

Not an issue — intentional by design.
