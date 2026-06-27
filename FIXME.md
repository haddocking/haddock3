# FIXME — Code review findings

Identified by static code review of `src/haddock/`. All issues are bugs or
latent correctness problems; style, documentation, and test-coverage gaps are
excluded.

---

## 1. [OPEN] `randint` produces out-of-range values

**File:** `src/haddock/libs/libmath.py:27`

```python
# current (wrong)
return int(self() * (upper_limit + 1)) + lower_limit
```

For any nonzero `lower_limit` the result can reach `lower_limit + upper_limit`,
exceeding `upper_limit` by up to `lower_limit`. The correct formula is:

```python
return lower_limit + int(self() * (upper_limit - lower_limit + 1))
```

**Impact:** All CNS seed generation uses this function. Seeds can exceed the
intended upper bound, but CNS accepts any integer seed so there is no runtime
failure. The practical effect is a subtle bias in the sampled range rather than
a hard error.

---

## 2. [OPEN] `open(None, ...)` crash when `error_file` is omitted

**File:** `src/haddock/libs/libsubprocess.py:249`

```python
with open(self.error_file, "wb+") as errf:   # TypeError if error_file is None
```

`error_file` is `Optional[FilePath]` with default `None`. If a CNS error is
detected and the caller omitted `error_file`, this raises `TypeError` instead
of a useful diagnostic.

**Fix:** Guard with `if self.error_file is not None:` before opening.

---

## 3. [OPEN] Wrong key in string-length validation error message

**File:** `src/haddock/gear/prepare_run.py:752`

```python
_desc = "under" if nbchars < param["minitems"] else "exceeding"
#                                  ^^^^^^^^^^^ should be param["minchars"]
```

For string parameters with `minchars`/`maxchars` bounds, the direction label is
computed by comparing `nbchars` against `param["minitems"]` — a key that exists
only on list-type parameters. This raises `KeyError` on any string-length
validation failure, hiding the real error.

---

## 4. [OPEN] Missing `f`-string prefix in `inject_in_modules` error

**File:** `src/haddock/gear/prepare_run.py:1067`

```python
raise ValueError(
    "key {key!r} already in {module!r} parameters. Can't inject."
)
```

Not an f-string — `{key!r}` and `{module!r}` are emitted literally. In
addition, `module` is not defined in that scope (the loop variable is `params`),
so adding the `f` prefix would immediately raise `NameError`.

**Fix:**
```python
raise ValueError(
    f"key {key!r} already in parameters. Can't inject."
)
```

---

## 5. [OPEN] Malformed CNS `structure...end` block when `psf_input` is a list

**File:** `src/haddock/libs/libcns.py:290-296`

When `psf_input` is a list, the code emits a first `structure` block containing
the Python list repr (`['a.psf', 'b.psf']`), then a second correct block:

```python
input_str += f"  @@{psf_input}{linesep}"   # emits list repr — invalid CNS syntax
```

CNS would fail to parse the first block. Currently latent (no caller passes a
list), but the code path exists.

**Fix:** Merge both branches into one block:
```python
if psf_input:
    input_str += f"structure{linesep}"
    if isinstance(psf_input, list):
        for psf in psf_input:
            input_str += f"  @@{psf}{linesep}"
    else:
        input_str += f"  @@{psf_input}{linesep}"
    input_str += f"end{linesep}"
    input_str += f"coor @@{pdb_input}{linesep}"
```

---

## 6. [OPEN] `@lru_cache` on function returning a mutable dict

**File:** `src/haddock/gear/prepare_run.py:162`

```python
@lru_cache
def _read_defaults(module_name, default_only=True):
    ...  # returns a dict
```

The cached dict is returned by reference. Callers that mutate it corrupt the
cache for subsequent calls with the same key. Because the same module (e.g.
`caprieval`) can appear multiple times in a workflow, parameters from step N
can bleed into step N+k.

**Fix:** `return copy.deepcopy(result)` at the end of `_read_defaults`, or copy
on the caller side before mutating.

---

## 7. [OPEN] Silent sampling truncation in rigidbody module

**File:** `src/haddock/modules/sampling/rigidbody/__init__.py:217`

```python
sampling_factor = int(self.params["sampling"] / len(models_to_dock))
```

`int()` truncates: 1000 structures with 3 model combinations → factor 333 →
999 structures actually generated, not 1000. No warning is emitted.

**Fix:** Use `round()` or `math.ceil()`, or log the actual count that will be
produced.

---

## 8. [OPEN] `split_tasks([])` raises `ValueError`

**File:** `src/haddock/libs/libparallel.py:21`

```python
n = math.ceil(len(lst) / n)   # n=0 when lst is empty → range step of 0
```

An empty task list causes `range(0, 0, 0)` which raises `ValueError: range()
arg 3 must not be zero`. Latent crash for any module that produces no jobs
(e.g. no valid input models).

**Fix:** Add `if not lst: return []` at the top of `split_tasks`.

---

## 9. [OPEN] Double assignment to `self.params["cns_exec"]`

**File:** `src/haddock/modules/base_cns_module.py:112-114`

```python
self.params["cns_exec"] = shutil.copyfile(_cns_exec, new_cns)  # immediately overwritten
shutil.copystat(_cns_exec, new_cns)
self.params["cns_exec"] = Path("..", Path(_cns_exec).name)
```

The first assignment stores the `copyfile` return value (the destination path
as a string) then is immediately overwritten. Harmless today but masks the
intent and could cause confusion if the order is ever changed.

---

## 10. [OPEN] Vestigial `self.len = score` in `PDBFile.__init__`

**File:** `src/haddock/libs/libontology.py:115`

```python
self.len = score   # score is NaN by default; .len is never read
```

`.len` is set on every `PDBFile` object but never accessed anywhere in the
codebase. Likely a leftover from a refactor.

---

## 11. [OPEN] Step-folder name parsing assumes module names contain no underscores

**File:** `src/haddock/modules/analysis/caprieval/capri.py:69`

```python
sel_steps = [step for step in sel_steps if step.count("_") == 1]
```

Step folder names follow the pattern `<N>_<module_name>`. This filter silently
drops any step whose module name contains an underscore (e.g. a hypothetical
`rigid_body`). Currently safe, but fragile.

---

## 12. [DROPPED] Missing energy term silently contributes 0 to HADDOCK score

Not an issue — intentional by design.
