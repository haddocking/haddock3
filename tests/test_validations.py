"""Test specific validations modules."""
from pathlib import Path

import pytest

from haddock.core.defaults import MODULE_DEFAULT_YAML
from haddock.core.exceptions import ConfigurationError
from haddock.gear.validations import v_rundir


@pytest.mark.parametrize(
    "rundir",
    [
        "some/folder",
        "some_folder/other",
        "some_folder_2/other",
        r"some\folder\on\windows",
        Path('some', 'folder', 'file.f'),
        ]
    )
def test_v_rundir_okay(rundir):
    """Test rundir validation okay."""
    v_rundir(rundir)


@pytest.mark.parametrize(
    "rundir",
    [
        "some/folder with spaces",
        "some_folder/ðßæ",
        "joão/folder",
        ]
    )
def test_v_rundir_invalid(rundir):
    """Test rundir validation okay."""
    with pytest.raises(ConfigurationError):
        v_rundir(rundir)
