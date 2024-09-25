import tempfile
from pathlib import Path

from tests import golden_data
from haddock.clis import cli_score

import io
from contextlib import redirect_stdout
import os


def test_cli_score():
    """Test the haddock3-score CLI."""
    pdb_f = Path(golden_data, "protprot_complex_1.pdb")
    # tempdir
    with tempfile.TemporaryDirectory(dir=".") as tmpdir:
        # parsing
        f = io.StringIO()
        with redirect_stdout(f):
            cli_score.main(pdb_f, tmpdir, full=True, keep_all=False)
        out = f.getvalue().split(os.linesep)
        assert out[-3].startswith("> HADDOCK-score (emscoring) = ")
        assert out[-2].startswith("> vdw")

        # now changing some weitght
        f = io.StringIO()
        with redirect_stdout(f):
            cli_score.main(pdb_f, tmpdir, full=True, keep_all=False, w_vdw=0.5)
        out = f.getvalue().split(os.linesep)
        assert out[4].startswith("> HADDOCK-score = (0.5 * vdw) + ")

        # now putting some non-standard emscoring parameters
        kwargs_dict = {
            "elecflag": "False",  # this is a boolean
            "epsilon": "2.0",  # this is a float
            "nemsteps": "100",  # this is an int
        }
        f = io.StringIO()
        with redirect_stdout(f):
            cli_score.main(pdb_f, tmpdir, full=True, keep_all=False, **kwargs_dict)
        out = f.getvalue().split(os.linesep)

        assert out[0].startswith(
            "* ATTENTION * Value (False) of parameter elecflag different from default"
        )
        assert out[1].startswith(
            "* ATTENTION * Value (2.0) of parameter epsilon different from default"
        )
        assert out[2].startswith(
            "* ATTENTION * Value (100) of parameter nemsteps different from default"
        )
        # check the used parameters
        assert out[4].startswith("used emscoring parameters: ")
        assert out[4].split("elecflag':")[1].split(",")[0].strip() == "False"

        assert out[-3].startswith("> HADDOCK-score (emscoring) = ")
