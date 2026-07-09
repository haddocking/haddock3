"""Make the repo root importable so tests here can reuse tests/conftest.py marks.

`end-to-end_tests` has a hyphen in its name, so it can never be a proper
Python package (unlike `tests`/`integration_tests`) and pytest's usual
rootdir-insertion via `__init__.py` doesn't apply here.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
