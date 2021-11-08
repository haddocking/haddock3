"""
Certify the developer has input all requirements for PR.

Situations tested:

* additions are reported in CHANGELOG.rst
"""
from pathlib import Path

folder = Path(__file__).resolve().parents[1]
changelog = Path(folder, 'CHANGELOG.rst')

with open(changelog, 'r') as fin:
    for line in fin:
        if line.startswith('v'):
            raise ValueError(
                'Please add a summary of your additions to CHANGELOG.rst. '
                'As described in: https://python-project-skeleton.readthedocs.io'
                '/en/latest/contributing.html#update-changelog.'
                )
        elif line.startswith('*'):
            break
