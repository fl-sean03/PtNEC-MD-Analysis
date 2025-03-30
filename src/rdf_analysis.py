#!/usr/bin/env python
"""
RDF analysis module for PtNEC system analysis.
This module computes the center-of-mass based radial distribution function.
"""

import os
import numpy as np
import MDAnalysis as mda
from typing import Tuple, List, Dict, Any, Optional

def compute_com_rdf(run_dir: str, 
                    dcd_file: str, 
                    pdb_file: str, 
                    psf_file: str, 
                    rdf_range: Tuple[float, float] = (0, 20), 
                    nbins: int = 200, 
                    cutoff: float = 2.5) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute a COM-based radial distribution function g(r) for the entire system.
    Includes all NEC fragments (without Pt atoms).
    
    Parameters
    ----------
    run_dir : str
        Directory containing simulation files
    dcd_file : str
        Name of the DCD trajectory file
    pdb_file : str
        Name of the PDB structure file
    psf_file : str
        Name of the PSF topology file
    rdf_range : Tuple[float, float]
        Range for RDF calculation in Angstroms (default: (0, 20))
    nbins : int
        Number of bins for the histogram (default: 200)
    cutoff : float
        Distance cutoff in Angstroms (not used for filtering, kept for API compatibility)
    
    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        bin_centers: centers of the histogram bins
        rdf_values: normalized g(r) values
    """
    # Load the universe with the final frame
    dcd_path = os.path.join(run_dir, dcd_file)
    pdb_path = os.path.join(run_dir, pdb_file)
    psf_path = os.path.join(run_dir, psf_file)
    
    u = mda.Universe(psf_path, pdb_path)
    u.load_new(dcd_path)
    u.trajectory[-1]  # Go to last frame
    
    # Get Pt atoms and COM
    pt_atoms = u.select_atoms("name PT*")
    pt_com = pt_atoms.center_of_mass()
    
    # Get all fragments that don't contain Pt atoms
    fragments = [frag for frag in u.atoms.fragments if len(frag.select_atoms("name PT*")) == 0]
    
    # Calculate COM distances from Pt COM
    com_distances = []
    for frag in fragments:
        frag_com = frag.center_of_mass()
        dist = np.linalg.norm(frag_com - pt_com)
        com_distances.append(dist)
    com_distances = np.array(com_distances)
    
    # Create histogram
    counts, bin_edges = np.histogram(com_distances, bins=nbins, range=rdf_range)
    bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
    
    # Calculate shell volumes
    dr = bin_edges[1] - bin_edges[0]
    shell_volumes = 4 * np.pi * bin_centers**2 * dr
    
    # Normalize by density and shell volumes
    total_fragments = len(com_distances)
    r_min, r_max = rdf_range
    total_volume = (4/3) * np.pi * (r_max**3 - r_min**3)
    # Handle the case where r_min is 0
    if r_min == 0:
        total_volume = (4/3) * np.pi * r_max**3
    
    # Avoid division by zero
    density = total_fragments / total_volume if total_volume > 0 and total_fragments > 0 else 1
    rdf_values = counts / (density * shell_volumes)
    
    return bin_centers, rdf_values

def compute_rdf_statistics(rdf_values: np.ndarray, bin_centers: np.ndarray) -> Dict[str, float]:
    """
    Compute statistics from the RDF curve.
    
    Parameters
    ----------
    rdf_values : np.ndarray
        Array of g(r) values
    bin_centers : np.ndarray
        Centers of the histogram bins
        
    Returns
    -------
    Dict[str, float]
        Dictionary of computed statistics:
        - peak_height: Maximum value of g(r)
        - peak_position: r value at which g(r) is maximum
        - first_min_position: r value of the first minimum after the peak
        - coordination_number: Integral of g(r) up to the first minimum
    """
    # Find the peak
    peak_idx = np.argmax(rdf_values)
    peak_height = rdf_values[peak_idx]
    peak_position = bin_centers[peak_idx]
    
    # Find the first minimum after the peak
    # Look for the minimum of g(r) in the range after the peak
    if peak_idx < len(rdf_values) - 1:
        min_idx = peak_idx + 1
        try:
            while min_idx < len(rdf_values) - 1 and rdf_values[min_idx] >= rdf_values[min_idx + 1]:
                min_idx += 1
            while min_idx < len(rdf_values) - 1 and rdf_values[min_idx] <= rdf_values[min_idx + 1]:
                min_idx += 1
                
            first_min_position = bin_centers[min_idx]
        except:
            # If we can't find a minimum, use a default
            first_min_position = peak_position * 1.5
    else:
        # If the peak is at the end of the array, use a default
        first_min_position = peak_position * 1.5
    
    # Calculate coordination number (integral of g(r) up to first minimum)
    # For each bin, we integrate 4πr²ρg(r)dr
    # where ρ is the average number density
    r_indices = bin_centers <= first_min_position
    dr = bin_centers[1] - bin_centers[0]
    coordination_number = np.sum(rdf_values[r_indices] * 4 * np.pi * bin_centers[r_indices]**2 * dr)
    
    return {
        "peak_height": peak_height,
        "peak_position": peak_position,
        "first_min_position": first_min_position,
        "coordination_number": coordination_number
    }

def combine_rdfs(system_rdf_data: Dict[str, Tuple[np.ndarray, np.ndarray]]) -> Dict[str, np.ndarray]:
    """
    Combine RDFs from different systems for comparison.
    
    Parameters
    ----------
    system_rdf_data : Dict[str, Tuple[np.ndarray, np.ndarray]]
        Dictionary mapping system names to (bin_centers, rdf_values) tuples
        
    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary containing:
        - 'r': Common bin centers
        - For each system, system_name: RDF values
    """
    # Check if all systems have the same bin centers
    first_system = next(iter(system_rdf_data.values()))
    first_bins = first_system[0]
    
    # Create the combined data structure
    combined_data = {'r': first_bins}
    
    # Add each system's RDF values
    for system_name, (bins, rdf_values) in system_rdf_data.items():
        combined_data[system_name] = rdf_values
    
    return combined_data