# Development Guide for MD Analysis Repository

## Introduction

This repository is designed for analyzing MD simulations of Pt nanoparticles (cuboctahedra) interacting with various NEC species that differ in hydrogenation state. The codebase is organized to be modular, maintainable, and extensible.

## Repository Structure

```
/MDAnalysis-PtNEC-Analysis
├── README.md
├── docs/
│   └── DEVELOPMENT_GUIDE.md
├── data/ (references to external location)
│   ├── 0HPt/
│   ├── 12HPt/
│   ├── 4HPt/
│   └── 8HPt/
├── analysis/
│   ├── combined/
│   ├── 0HPt/
│   ├── 12HPt/
│   ├── 4HPt/
│   └── 8HPt/
├── src/
│   ├── __init__.py
│   ├── fragment_analysis.py
│   ├── rdf_analysis.py
│   ├── plotting.py
│   ├── utils.py
│   └── main.py
├── tests/
├── config/
└── requirements.txt
```

## Code Organization

### Fragment Analysis Module (`src/fragment_analysis.py`)

Functions for analyzing NEC fragments relative to Pt nanoparticles:
- Computing fragment metrics (min distance, COM distance, etc.)
- Identifying attached fragments
- Calculating contact fractions and orientation angles

### RDF Analysis Module (`src/rdf_analysis.py`)

Functions for computing radial distribution functions:
- COM-based RDF calculations
- Binning and normalization methods
- Functions to extract spatial distribution information

### Plotting Module (`src/plotting.py`)

Functions for visualization:
- Individual system plots
- Comparative plots across systems
- Standardized styling and formatting

### Utilities Module (`src/utils.py`)

Helper functions:
- Logging utilities
- File I/O operations
- Directory management

### Main Workflow (`src/main.py`)

Script that ties all modules together:
- Processing each system individually
- Generating combined analyses and comparisons
- Command-line interface

## Development Practices

- **Documentation:** Include docstrings for all functions and classes
- **Type Hints:** Use Python type hints for better code understanding
- **Testing:** Write unit tests for critical functionality
- **Modular Design:** Keep modules focused on specific functionality
- **Error Handling:** Use proper exception handling and logging
- **Configuration:** Externalize parameters in config files

## Data Handling

The repository references data stored at the following location:
```
/home/sf2/LabWork/Workspace/9-NAMDAnalysis/data/
```

Do not copy this data (12GB total) - always reference it in place.

## Extending the Codebase

When adding new features:
1. Place relevant functionality in the appropriate module
2. Update tests to cover new functionality
3. Document the feature in code and update this guide as needed
4. Ensure backward compatibility where possible