"""Test general functions of haddock3 analysis modules."""

from haddock.modules.analysis import get_analysis_exec_mode


def test_get_analysis_exec_mode():
    """Test the get_analysis_exec_mode function."""
    assert get_analysis_exec_mode("local") == "local"
    assert get_analysis_exec_mode("batch") == "local"
    assert get_analysis_exec_mode("mpi") == "mpi"
