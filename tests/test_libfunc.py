"""Test libfunc."""

from haddock.libs.libfunc import chainf, chainfs


def test_chainf():
    result = chainf(2, *[str, int, float])
    assert isinstance(result, float)
    assert result == 2.0


def test_chainfs():
    do = chainfs(str, int, float)
    result = do(2)
    assert isinstance(result, float)
    assert result == 2.0
