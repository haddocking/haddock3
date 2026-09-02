"""Test exceptions module."""
import pytest

from haddock.core.exceptions import (
    CNSRunningError,
    HaddockError,
    HaddockTaskExecutionError,
)


def test_cns_error():
    """Test CNS error."""
    with pytest.raises(CNSRunningError):
        raise CNSRunningError(b'something')


def test_cns_error_is_a_task_execution_error():
    """CNS failures can be tolerated by task schedulers."""
    assert issubclass(CNSRunningError, HaddockTaskExecutionError)
    assert issubclass(HaddockTaskExecutionError, HaddockError)
