# Contributing to Genome Sequencing Analyzer

Thank you for your interest in contributing to this project! This document provides guidelines and instructions for contributing.

## Code of Conduct

Please be respectful and constructive in all interactions. We welcome contributors of all experience levels.

## How to Contribute

### Reporting Bugs

1. Check existing issues to avoid duplicates
2. Create a new issue with:
   - Clear, descriptive title
   - Steps to reproduce the problem
   - Expected vs. actual behavior
   - System information (Python version, OS)
   - Sample data if applicable (anonymized)

### Suggesting Features

1. Open an issue with the "enhancement" label
2. Describe the feature and its use case
3. Explain how it benefits the bioinformatics community

### Submitting Code

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/your-feature-name`
3. Make your changes following the style guidelines
4. Add tests for new functionality
5. Ensure all tests pass: `pytest tests/ -v`
6. Commit with clear messages: `git commit -m "Add feature: description"`
7. Push to your fork: `git push origin feature/your-feature-name`
8. Open a Pull Request

## Development Setup

```bash
# Clone your fork
git clone https://github.com/yourusername/genome-sequencing-analyzer.git
cd genome-sequencing-analyzer

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Run tests
pytest tests/ -v
```

## Code Style Guidelines

### Python

- Follow PEP 8 style guidelines
- Use type hints for function signatures
- Write docstrings for all public functions (Google style)
- Maximum line length: 100 characters

### Documentation

- Update docstrings when modifying functions
- Add examples for complex functionality
- Keep README.md current with new features

### Example Function

```python
def calculate_gc_content(sequence: str) -> float:
    """
    Calculate the GC content percentage of a DNA sequence.
    
    GC content is the proportion of guanine (G) and cytosine (C) bases
    in a nucleotide sequence.
    
    Args:
        sequence: DNA sequence string (case-insensitive)
    
    Returns:
        GC content as a percentage (0-100)
    
    Example:
        >>> calculate_gc_content("ATGCATGC")
        50.0
    """
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    total = len([n for n in sequence if n in 'ATCG'])
    return (gc_count / total * 100) if total > 0 else 0.0
```

## Testing

- Write tests for all new functionality
- Place tests in `tests/` directory
- Name test files as `test_<module>.py`
- Use descriptive test function names

```python
def test_gc_content_balanced():
    """Test GC content with balanced sequence."""
    assert calculate_gc_content("ATGC") == 50.0
```

## Project Structure

When adding new features:

- **Core analysis code**: Add to `app/core/`
- **UI components**: Add to `app/ui/`
- **Tests**: Add to `tests/`
- **Sample data**: Add to `sample_data/`

## Pull Request Process

1. Update documentation for any changed functionality
2. Add yourself to contributors in README if desired
3. Ensure CI checks pass
4. Request review from maintainers
5. Address feedback promptly

## Questions?

Open an issue with the "question" label or reach out to the maintainers.

Thank you for contributing!
