import pytest

import haddock.clis.wrapper_haddock_restraints as wrapper


def test_wrapper_haddock_restraints(monkeypatch):

    monkeypatch.setattr("sys.argv", ["wrapper_haddock_restraints.py", "--help"])

    with pytest.raises(SystemExit) as e:
        wrapper.main()

    assert e.value.code == 0

    monkeypatch.setattr("sys.argv", ["wrapper_haddock_restraints.py", ""])

    with pytest.raises(SystemExit) as e:
        wrapper.main()

    assert e.value.code == 2
