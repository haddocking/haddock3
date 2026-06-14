# Contributing to HADDOCK3

HADDOCK3 welcomes contributions to improve its
functionality, documentation, and code quality.
This guide provides the essential information for
contributing effectively.

## Contribution Workflow

### 1. Prerequisites

- Install HADDOCK3 following [INSTALL.md](docs/INSTALL.md)
- Python 3.10+ development environment
- Familiarity with Git and GitHub workflows

### 2. Development Process

1. **Fork the repository** and create a feature branch
2. **Implement changes** following coding standards
3. **Test thoroughly** using pytest framework
4. **Update documentation** as needed
5. **Submit Pull Request** for review

## Code Contribution Guidelines

### Project Structure

```text
src/haddock/
├── clis/          # Command-line interfaces
├── libs/          # General utility functions
├── gear/          # Plugin-like functionality modules
├── core/          # Physical constants and definitions
└── modules/       # HADDOCK3 simulation modules
```

### Coding Standards

- **Python Version**: Minimum 3.10 compatibility
- **Function Design**: Small, testable functions preferred over complex classes
- **Naming**: Use descriptive variable names
- **Comments**: Explain *why* not *how*
- **Documentation**: Update docstrings and markdown files

### Testing Requirements

- **Unit Tests**: Located in `tests/` directory
- **Integration Tests**: Located in `integration_tests/` directory
- **End-to-End**: Located in `end-to-end_tests/` directory
- **Coverage**: Aim for 100% test coverage for new code
- **Framework**: pytest

### Dependency Policy

HADDOCK3 maintains minimal dependencies:

1. Prefer Python standard library
2. NumPy allowed for numerical operations
3. Avoid adding new dependencies without discussion
4. Consider runtime dependencies for optional functionality

## Documentation Contributions

HADDOCK3 documentation uses Sphinx with Markdown and reStructuredText:

### Local Documentation Build

```bash
# Install documentation dependencies
pip install -e '.[docs]'

# Build documentation
sphinx-apidoc -f -e -o docs/ src/haddock -d 1
sphinx-build -b html docs haddock3-docs
```

### Documentation Structure

- Source files in `docs/` directory
- Follow existing organization and formatting
- Update both code docstrings and markdown files

## Pull Request Requirements

1. **Testing**: All tests must pass
2. **Documentation**: Updated for new features
3. **CHANGELOG**: Add entry for significant changes
4. **Version**: Update `pyproject.toml` if applicable
5. **Code Review**: Address all feedback

## Support and Communication

For questions or discussions:

- **Issues**: Report bugs or suggest features via GitHub Issues
- **Discussions**: Use GitHub Discussions for general questions
- **Development**: Contact maintainers for major changes
