"""Make the repo root importable so tests here can reuse tests/conftest.py marks."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
