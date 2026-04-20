# Code Documentation Guide for HADDOCK3

This guide provides comprehensive information about HADDOCK3's code documentation standards, tools, and processes.

## Documentation Standards

### Docstring Format

HADDOCK3 uses [Google-style docstrings](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings) for all Python code:

```python
def calculate_distance(atom1, atom2):
    """Calculate Euclidean distance between two atoms.
    
    Args:
        atom1 (Atom): First atom object with x, y, z coordinates.
        atom2 (Atom): Second atom object with x, y, z coordinates.
        
    Returns:
        float: Euclidean distance between the atoms in Angstroms.
        
    Example:
        >>> distance = calculate_distance(atom_a, atom_b)
        >>> print(f"Distance: {distance:.2f} Å")
    """
    # Function implementation
    ...
```

### Key Documentation Elements

1. **Module Docstrings**: Every module should have a module-level docstring explaining its purpose
2. **Class Docstrings**: Include description, attributes, and methods
3. **Function Docstrings**: Include args, returns, raises, and examples where applicable
4. **Type Hints**: Use Python type hints for better code clarity

## Documentation Tools

### Sphinx Documentation

HADDOCK3 uses Sphinx to generate HTML documentation from:
- Python docstrings
- Markdown files (.md)
- reStructuredText files (.rst)

### Building Documentation Locally

```bash
# Install documentation dependencies
pip install -e '.[docs]'

# Generate API documentation
sphinx-apidoc -f -e -o docs/ src/haddock -d 1

# Build HTML documentation
sphinx-build -b html docs haddock3-docs

# View documentation
open haddock3-docs/index.html  # macOS
xdg-open haddock3-docs/index.html  # Linux
```

## Documentation Structure

### Source Code Documentation

```
src/haddock/
├── __init__.py          # Package-level documentation
├── modules/             # Module-specific documentation
│   ├── __init__.py      # Modules package docs
│   └── *.py             # Individual module docs
├── libs/                # Library function documentation
└── gear/                # Gear/component documentation
```

### Documentation Files

```
docs/
├── index.rst            # Main documentation entry point
├── user/                # User guide documentation
├── developer/           # Developer documentation
├── api/                 # Auto-generated API docs
└── *.md                  # Various markdown guides
```

## Writing Good Documentation

### Best Practices

1. **Be Concise but Complete**: Explain what the code does, not how it does it
2. **Use Examples**: Include usage examples where helpful
3. **Document Exceptions**: Use `Raises:` section for exceptions
4. **Update Regularly**: Keep documentation in sync with code changes
5. **Use Cross-References**: Link to related functions/classes

### Common Sections

```python
"""Function docstring template.

[Summary line]

[Extended description]

Args:
    param1 (type): Description
    param2 (type, optional): Description. Defaults to X.

Returns:
    type: Description of return value

Raises:
    ExceptionType: Description of when this exception is raised

Note:
    Additional notes or warnings

Example:
    >>> result = function_name(arg1, arg2)
    >>> print(result)
"""
```

## Documentation Workflow

### For New Features

1. Write code with proper docstrings
2. Add examples and usage notes
3. Update related documentation files
4. Test documentation builds locally
5. Submit with your pull request

### For Bug Fixes

1. Update affected docstrings if behavior changes
2. Add notes about fixed edge cases
3. Update examples if needed

## Documentation Testing

### Testing Your Documentation

```bash
# Test docstring examples (if any)
python -m doctest src/haddock/modules/your_module.py

# Check for documentation coverage
pydocstyle src/haddock/ --select=D100,D101,D102,D103,D104
```

### Common Documentation Issues

- **Missing docstrings**: All public functions/classes should have docstrings
- **Outdated examples**: Examples should work with current code
- **Inconsistent formatting**: Follow Google style guide
- **Broken references**: Check cross-references work

## Advanced Documentation

### Mathematical Notation

For complex algorithms, use LaTeX-style math notation:

```python
"""Calculate protein-ligand binding energy.

Uses the following formula:

.. math::
    \Delta G = \Delta H - T\Delta S

Where:
    - ΔG: Gibbs free energy change
    - ΔH: Enthalpy change  
    - T: Temperature in Kelvin
    - ΔS: Entropy change
"""
```

### Code References

```python
"""Perform molecular dynamics simulation.

See Also:
    :func:`~haddock.modules.simulation.setup_simulation`
    :class:`~haddock.core.simulation.SimulationEngine`

Note:
    This function uses the OpenMM engine for MD calculations.
    For alternative engines, see :mod:`haddock.modules.engines`.
"""
```

## Troubleshooting Documentation

### Common Problems

**Documentation won't build:**
- Check for syntax errors in .rst files
- Ensure all referenced modules can be imported
- Verify Sphinx extensions are installed

**Missing API documentation:**
- Run `sphinx-apidoc` to regenerate API docs
- Check module `__all__` variables

**Broken links:**
- Use `sphinx-build -b linkcheck docs output` to check links
- Fix any broken external references

