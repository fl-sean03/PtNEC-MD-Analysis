# Contributing to PtNEC MD Analysis

Thank you for your interest in contributing to the PtNEC MD Analysis project! This document provides guidelines and instructions for contributing.

## Development Environment Setup

1. Fork the repository on GitHub.
2. Clone your fork locally:
   ```bash
   git clone https://github.com/YOUR-USERNAME/MDAnalysis-PtNEC-Analysis.git
   cd MDAnalysis-PtNEC-Analysis
   ```

3. Create a virtual environment and install dependencies:
   ```bash
   python -m venv venv
   source venv/bin/activate  # Windows: venv\Scripts\activate
   pip install -r requirements.txt
   pip install -e .  # Install in development mode
   ```

4. Install development dependencies:
   ```bash
   pip install pytest pytest-cov flake8 black
   ```

## Code Style and Guidelines

- Use [Black](https://black.readthedocs.io/) for code formatting
- Follow [PEP 8](https://www.python.org/dev/peps/pep-0008/) style guidelines
- Use type hints where appropriate
- Write comprehensive docstrings following the NumPy/SciPy style

## Development Workflow

1. Create a new branch for your feature or bugfix:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes and test them:
   ```bash
   pytest tests/
   ```

3. Run the code formatter:
   ```bash
   black src/ tests/
   ```

4. Check for linting issues:
   ```bash
   flake8 src/ tests/
   ```

5. Commit your changes with a descriptive message:
   ```bash
   git commit -m "Add feature: your feature description"
   ```

6. Push your branch to GitHub:
   ```bash
   git push origin feature/your-feature-name
   ```

7. Create a pull request from your branch to the main repository.

## Pull Request Process

1. Update the README.md or documentation with details of changes if appropriate
2. Add or update tests for any new functionality
3. Make sure all tests are passing
4. Get your pull request reviewed and approved by maintainers
5. Your pull request will be merged once approved

## Adding New Features

When adding new features, please follow these guidelines:

1. Place new functionality in the appropriate module:
   - Fragment analysis related code in `src/fragment_analysis.py`
   - RDF analysis related code in `src/rdf_analysis.py`
   - Plotting functions in `src/plotting.py`
   - Utility functions in `src/utils.py`

2. Add unit tests for your feature in the `tests/` directory

3. Update the configuration file if needed

4. Document your feature in the code and in the documentation

## Reporting Issues

Please use the GitHub issue tracker to report bugs or suggest features.

When reporting bugs, please include:
- A clear description of the problem
- Steps to reproduce the issue
- Expected behavior
- Screenshots if applicable
- Your environment details (OS, Python version, etc.)

## License

By contributing to this project, you agree that your contributions will be licensed under the project's license.