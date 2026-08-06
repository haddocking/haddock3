"""Test the workflow module ordering validation."""

import pytest

from haddock.core.exceptions import ConfigurationError
from haddock.gear.workflow_ordering import (
    read_workflow_rules,
    validate_workflow_order,
)


def test_read_workflow_rules_has_all_keys():
    """The bundled rules expose every supported rule key."""
    rules = read_workflow_rules()
    assert set(rules) >= {
        "disallowed_sequences",
        "required_preceding",
        "required_prior",
        "required_first",
    }


@pytest.mark.parametrize(
    "modules",
    [
        ["topoaa", "rigidbody", "flexref", "emref"],
        ["topoaa", "rigidbody", "rmsdmatrix", "clustrmsd", "seletopclusts"],
        ["topoaa", "rigidbody", "ilrmsdmatrix", "clustrmsd"],
        ["topoaa", "topocg", "rigidbody"],
        ["topoaa", "topocg", "rigidbody", "cgtoaa", "flexref"],
        ["topoaa", "rigidbody", "clustfcc", "seletopclusts", "clustfcc"],
    ],
)
def test_valid_workflows_pass(modules):
    """Well-ordered workflows validate without error."""
    validate_workflow_order(modules)


def test_required_first_violation():
    """A workflow not starting with topoaa is rejected."""
    with pytest.raises(ConfigurationError, match="first module"):
        validate_workflow_order(["rigidbody", "flexref"])


@pytest.mark.parametrize(
    "modules",
    [
        ["topoaa", "rmsdmatrix", "clustrmsd", "clustrmsd"],
        ["topoaa", "rigidbody", "clustfcc", "clustfcc"],
        ["topoaa", "rigidbody", "seletop", "seletop"],
        ["topoaa", "seletopclusts", "seletopclusts"],
    ],
)
def test_disallowed_sequence_violation(modules):
    """Repeating a non-repeatable module directly is rejected."""
    with pytest.raises(ConfigurationError, match="cannot directly follow"):
        validate_workflow_order(modules)


def test_required_preceding_violation_wrong_predecessor():
    """clustrmsd not directly preceded by an rmsd matrix is rejected."""
    with pytest.raises(ConfigurationError, match="directly preceded"):
        validate_workflow_order(["topoaa", "rigidbody", "clustrmsd"])


def test_required_preceding_violation_topocg():
    """topocg must directly follow topoaa."""
    with pytest.raises(ConfigurationError, match="directly preceded"):
        validate_workflow_order(["topoaa", "rigidbody", "topocg"])


def test_required_prior_violation():
    """cgtoaa requires topocg somewhere earlier."""
    # required_first also fires here, but the prior-module message must appear
    with pytest.raises(ConfigurationError, match="earlier step"):
        validate_workflow_order(["rigidbody", "cgtoaa"])


def test_multiple_violations_reported_together():
    """All ordering violations are collected in a single error."""
    with pytest.raises(ConfigurationError) as exc:
        validate_workflow_order(["rigidbody", "seletop", "seletop"])
    message = str(exc.value)
    assert "first module" in message
    assert "cannot directly follow" in message


def test_empty_workflow_does_not_raise():
    """An empty module list is not an ordering error on its own."""
    validate_workflow_order([])


def test_custom_rules_file(tmp_path):
    """Rules are read from the provided file."""
    rules_file = tmp_path / "rules.yaml"
    rules_file.write_text("required_first:\n  - rigidbody\n")
    # topoaa first now violates the custom rule
    with pytest.raises(ConfigurationError, match="first module"):
        validate_workflow_order(["topoaa"], rules_file=rules_file)
    # rigidbody first satisfies it
    validate_workflow_order(["rigidbody"], rules_file=rules_file)
