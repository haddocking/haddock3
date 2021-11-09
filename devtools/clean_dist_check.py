"""Clean dist files after tox -e check"""
from pathlib import Path

folder = Path(__file__).resolve().parents[1]

whls = folder.glob('dist/*.whl')

for f in whls:
    f.unlink()

targz = folder.glob('dist/*.tar.gz')

for f in targz:
    f.unlink()
