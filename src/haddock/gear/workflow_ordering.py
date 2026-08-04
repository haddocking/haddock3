"""
Validate the ordering of modules in a HADDOCK3 workflow.

Module- and parameter-level validation (see :mod:`haddock.gear.prepare_run`)
confirms that each requested module exists and that its parameters are valid,
but it does not enforce constraints on the *sequence* of modules. Some modules
only make sense when another module ran before them (e.g. ``clustrmsd`` needs
an RMSD matrix), some must not be chained to themselves (e.g. two consecutive
``seletop``), and a workflow must always start by building topologies.

Those constraints are declared as data in ``workflow_rules.yaml`` and applied
here against the ordered list of modules taken from the user configuration
file. Keeping the rules in a data file means new constraints can be added
without editing this logic.

Rule types (see ``workflow_rules.yaml`` for the exact syntax):

- ``disallowed_sequences``: module B must not directly follow module A.
- ``required_preceding``: module B must be immediately preceded by one of a set.
- ``required_prior``: module B must be preceded anywhere earlier by one of a set.
- ``required_first``: the first module must be one of a set.
"""
from pathlib import Path

from haddock import log
from haddock.core.exceptions import ConfigurationError
from haddock.libs.libio import read_from_yaml


DEFAULT_RULES = Path(Path(__file__).resolve().parent, "workflow_rules.yaml")


def read_workflow_rules(rules_file=DEFAULT_RULES):
    """
    Read the workflow ordering rules from a YAML file.

    Parameters
    ----------
    rules_file : str or pathlib.Path
        Path to the YAML file defining the ordering rules.
        Defaults to the bundled ``workflow_rules.yaml``.

    Returns
    -------
    dict
        The parsed rules, with every supported rule key present (empty
        containers for the ones not defined in the file).
    """
    rules = read_from_yaml(rules_file)
    # guarantee all keys exist so callers need not test for their presence
    rules.setdefault("disallowed_sequences", [])
    rules.setdefault("required_preceding", {})
    rules.setdefault("required_prior", {})
    rules.setdefault("required_first", [])
    return rules


def _check_required_first(modules, required_first):
    """Collect a violation if the first module is not an allowed one."""
    errors = []
    if required_first and modules and modules[0] not in required_first:
        errors.append(
            f"The first module of a workflow must be one of "
            f"{_join(required_first)}, but {modules[0]!r} was found."
        )
    return errors


def _check_disallowed_sequences(modules, disallowed_sequences):
    """Collect violations for modules directly following a disallowed one."""
    errors = []
    disallowed = {tuple(pair) for pair in disallowed_sequences}
    for position, (before, after) in enumerate(zip(modules, modules[1:]), start=2):
        if (before, after) in disallowed:
            errors.append(
                f"Module {after!r} (step {position}) cannot directly follow "
                f"module {before!r}."
            )
    return errors


def _check_required_preceding(modules, required_preceding):
    """Collect violations for modules lacking a required direct predecessor."""
    errors = []
    for position, module in enumerate(modules, start=1):
        allowed = required_preceding.get(module)
        if allowed is None:
            continue
        preceding = modules[position - 2] if position >= 2 else None
        if preceding not in allowed:
            found = f"{preceding!r}" if preceding is not None else "nothing"
            errors.append(
                f"Module {module!r} (step {position}) must be directly preceded "
                f"by one of {_join(allowed)}, but {found} was found."
            )
    return errors


def _check_required_prior(modules, required_prior):
    """Collect violations for modules lacking a required earlier module."""
    errors = []
    for position, module in enumerate(modules, start=1):
        allowed = required_prior.get(module)
        if allowed is None:
            continue
        earlier = modules[: position - 1]
        if not any(mod in earlier for mod in allowed):
            errors.append(
                f"Module {module!r} (step {position}) requires one of "
                f"{_join(allowed)} to run at some earlier step."
            )
    return errors


def _join(modules):
    """Render a collection of module names as a readable quoted list."""
    return ", ".join(repr(mod) for mod in modules)


def validate_workflow_order(modules, rules_file=DEFAULT_RULES):
    """
    Validate the order of the modules of a workflow against the rules.

    Parameters
    ----------
    modules : sequence of str
        The module names in workflow order (without their ``.N`` suffix).

    rules_file : str or pathlib.Path
        Path to the YAML file defining the ordering rules.
        Defaults to the bundled ``workflow_rules.yaml``.

    Raises
    ------
    haddock.core.exceptions.ConfigurationError
        If the workflow violates one or more ordering rules. All detected
        violations are reported together.
    """
    modules = list(modules)
    rules = read_workflow_rules(rules_file)

    errors = []
    errors.extend(_check_required_first(modules, rules["required_first"]))
    errors.extend(
        _check_disallowed_sequences(modules, rules["disallowed_sequences"])
    )
    errors.extend(
        _check_required_preceding(modules, rules["required_preceding"])
    )
    errors.extend(_check_required_prior(modules, rules["required_prior"]))

    if errors:
        msg = "Invalid workflow module order:" + "".join(
            f"\n  - {error}" for error in errors
        )
        raise ConfigurationError(msg)

    log.info("Workflow module order validated.")
