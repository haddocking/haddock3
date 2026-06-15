# HADDOCK3 Development Guide

This guide provides instructions for setting up a HADDOCK3 development environment.

## System Requirements

- **Python**: 3.9-3.14
- **MPI**: OpenMPI (optional, required for MPI support)

## Development Environment Setup

Clone the repository,setup Python environment and install in editable mode

> If you need help seting up your Python environment, look into [docs/PYTHON.md](docs/PYTHON.md)

```
git clone https://github.com/haddocking/haddock3.git
cd haddock3
python3.14 -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'
```

Please refer to `pyproject.toml` for additional dependencies.

## Running Tests

HADDOCK3 uses pytest for testing:

```bash
pytest tests/

pytest integration_tests/

pytest end-to-end_tests/
```

## Code Formatting

HADDOCK3 uses [ruff](https://docs.astral.sh/ruff/) to format Python code
(included in the `dev` extra). Before submitting a pull request, format any
files you changed:

```bash
ruff format <path/to/changed_file.py>
```

CI checks the formatting of changed files on every pull request and will
fail if `ruff format --check` reports differences.

## CNS Executable Troubleshooting

If you encounter CNS related errors please refer to [CNS.md](docs/CNS.md)
