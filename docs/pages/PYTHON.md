# Python Environment Setup for HADDOCK3

This guide provides instructions for setting up Python environments
for HADDOCK3 development and usage.

## Python Version Requirements

HADDOCK3 requires Python 3.10 or later:

- **Recommended**: Python 3.14+ for best compatibility
- **Supported versions**: 3.10-3.14
- **Minimum required**: Python 3.10

## Python Installation

HADDOCK3 requires Python 3.10 or later.

### Using uv (Recommended)

[uv](https://github.com/astral-sh/uv) is a modern Python package manager that
provides fast dependency resolution and virtual environment management:

```bash
# Install uv
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create and activate a virtual environment
uv venv .venv --python=3.14
source .venv/bin/activate
```

### System Python

Alternatively, use your system's package manager:

#### Linux (Ubuntu/Debian)

```bash
sudo apt-get update
sudo apt-get install python3.14 python3.14-dev python3.14-venv
```

#### macOS (Homebrew)

```bash
brew install python@3.14
```

---

```bash
python3.14 -m venv .venv
source .venv/bin/activate 
```

## Troubleshooting

### Common Issues

**Python version conflicts:**

```bash
# Check which Python is being used
which python
which python3
```

**Virtual environment activation issues:**

```bash
# Ensure proper shell integration
echo $VIRTUAL_ENV
```

**Permission errors:**

```bash
# Use --user flag or virtual environment
pip install --user package
```

### Debugging Tools

```bash
# Check installed packages
pip list

# Check package locations
pip show package_name

# Verify Python environment
python -c "import sys; print(sys.executable)"
```
