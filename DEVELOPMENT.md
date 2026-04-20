# HADDOCK3 Development Guide

This guide provides instructions for setting up a HADDOCK3 development environment.

## System Requirements

- **Python**: 3.9-3.14
- **MPI**: OpenMPI (optional, required for MPI support)

## Development Environment Setup

### 1. Clone Repository

```bash
git clone https://github.com/haddocking/haddock3.git
cd haddock3
```

### 2. Python Environment

```bash
# Create virtual environment
python3.14 -m venv .venv
source .venv/bin/activate
```

### 3. Install Dependencies

```bash
# Basic development installation
pip install -e '.[dev,docs]'

# With MPI support (optional)
pip install -e '.[dev,docs,mpi]'
```

> **Note for macOS**: If mpi4py installation fails, run `brew install mpi4py` first.

## Running Tests

HADDOCK3 uses pytest for testing:

```bash
# Unit tests
pytest tests/

# Integration tests
pytest integration_tests/

# End-to-end tests
pytest end-to-end_tests/
```

## CNS Executable Troubleshooting

If you encounter CNS related errors please refer to [CNS.md](docs/CNS.md)
