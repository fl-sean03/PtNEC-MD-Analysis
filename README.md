# MDAnalysis-PtNEC-Analysis

Analysis repository for MD simulations of Pt nanoparticles interacting with NEC species.

## Overview

This repository contains code for analyzing molecular dynamics simulations of platinum nanoparticles (cuboctahedra) interacting with various NEC species (differing in hydrogenation state). The analysis includes fragment metrics, COM-based RDF calculations, and visualization tools.

## Features

- **Fragment Metrics:** Quantifying contact versus overall positioning for attached NEC fragments
- **COM-Based RDF:** Calculating radial distribution function of NEC fragment centers of mass relative to Pt nanoparticles
- **Visualization:** Individual and comparative plots for key metrics

## Setup

```bash
# Clone the repository
git clone [repository-url]
cd MDAnalysis-PtNEC-Analysis

# Create and activate virtual environment
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt
```

## Usage

```bash
python src/main.py
```

See the documentation in `docs/` for detailed usage instructions.
