# Contributing to HADDOCK3

HADDOCK3 welcomes contributions to improve its
functionality, documentation, and code quality.
This guide provides the essential information for
contributing effectively.

## Contribution Workflow

### 1. Prerequisites

- Install HADDOCK3 following [INSTALL.md](https://github.com/haddocking/haddock3/blob/main/docs/pages/INSTALL.md)
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
- **Reuse**: Existing code being functional doesn't mean it's the best
  implementation. If you ctrl-c ctrl-v code from elsewhere in the codebase evaluate it on its own merits and improve it if you can.
- **Documentation**: Update docstrings and markdown files
- **Formatting**: Code is formatted with [ruff](https://docs.astral.sh/ruff/)
  (`ruff format`); CI checks formatting of changed files on every pull request

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

- Update docstrings for any function, class, or module you add or change.
- Update the relevant pages under `docs/pages/` for user-facing changes.
- New modules need a title to be added manually in `docs/titles.yaml`.
- Confirm `docs/` still builds without introducing new warnings before submitting.

See [docs/README.md](https://github.com/haddocking/haddock3/blob/main/docs/README.md)
for how to build the docs locally and how docstrings should be formatted.

## Pull Request Requirements

1. **Testing**: All tests must pass
2. **Documentation**: Updated for new features
3. **CHANGELOG**: Add entry for significant changes
4. **Version**: Update `pyproject.toml` if applicable
5. **Code Review**: Address all feedback

## Create a new release

Here are the guidelines to create a new release:

1. On the `haddock3` repository:
 - Update the release version in the `pyproject.toml` file: `version = "YYYY.M.0"` (e.g.: version = "2026.8.0")
 - Create a pull request and merge it
 - Create a new release from [https://github.com/haddocking/haddock3/releases](https://github.com/haddocking/haddock3/releases)
 - The new version will automatically be published on PyPI thanks to the `.github/workflows/publish.yml`
2. On the `haddock3-user-manual` repository:
 - Create a new [release/tag](https://github.com/haddocking/haddock3-user-manual/releases) with the same version name, so both are matching 

## Use of AI tools

See [AI-POLICY.md](https://github.com/haddocking/haddock3/blob/main/AI-POLICY.md) for guidance on how AI coding assistants
may and may not be used when contributing.

## Support and Communication

For questions or discussions:

- **Issues**: Report bugs or suggest features via GitHub Issues
- **Discussions**: Use GitHub Discussions for general questions
- **Development**: Contact maintainers for major changes
