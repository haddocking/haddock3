import pytest

from haddock.core.exceptions import CNSRunningError


def test_cns_error():
    with pytest.raises(CNSRunningError):
        raise CNSRunningError(b'something')
