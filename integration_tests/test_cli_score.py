import tempfile
from pathlib import Path

from tests import golden_data
from haddock.clis import cli_score

import io
from contextlib import redirect_stdout
import os
from . import has_cns

@has_cns
def test_cli_score():
    """Test the haddock3-score CLI."""
    pdb_f = Path(golden_data, "protprot_complex_1.pdb")
    # tempdir
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
        assert out[3].startswith("> HADDOCK-score = (0.5 * vdw) + ")
