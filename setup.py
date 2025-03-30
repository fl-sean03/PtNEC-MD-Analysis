#!/usr/bin/env python
from setuptools import setup, find_packages

setup(
    name="ptnec_analysis",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "MDAnalysis>=2.0.0",
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "matplotlib>=3.4.0",
        "scipy>=1.7.0",
        "pyyaml>=5.1",
    ],
    entry_points={
        "console_scripts": [
            "ptnec-analysis=src.main:main",
        ],
    },
    python_requires=">=3.8",
    author="Lab Research Team",
    description="Analysis toolkit for MD simulations of Pt nanoparticles with NEC species",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)