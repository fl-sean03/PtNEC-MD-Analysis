#!/usr/bin/env python
"""
Fragment analysis module for PtNEC system analysis.
This module contains functions for analyzing NEC fragments relative to Pt nanoparticles.
"""

import os
import numpy as np
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from typing import Dict, List, Tuple, Any, Optional, Union

def compute_fragment_metrics(frag: mda.core.groups.AtomGroup, 
                             pt_atoms: mda.core.groups.AtomGroup,
                             pt_com: np.ndarray,
                             close_threshold: float = 2.5) -> Dict[str, Any]:
    """
    Compute distance statistics for a fragment relative to Pt atoms.
    
    Parameters
    ----------
    frag : mda.core.groups.AtomGroup
        Fragment atom group
    pt_atoms : mda.core.groups.AtomGroup
        Platinum atoms atom group
    pt_com : np.ndarray
        Center of mass of the platinum nanoparticle
    close_threshold : float
        Distance threshold for close contacts in Angstroms (default: 2.5)
        
    Returns
    -------
    Dict[str, Any]
        Dictionary with computed metrics:
        - min_dist: Minimum distance between any fragment atom and any Pt atom
        - avg_dist: Average distance between fragment atoms and Pt atoms
        - med_dist: Median distance between fragment atoms and Pt atoms
        - std_dist: Standard deviation of distances between fragment atoms and Pt atoms
        - com_dist: Distance between the fragment's COM and Pt COM
        - num_close_contacts: Count of atom pairs with distance < close_threshold
        - contact_fraction: Fraction of fragment atoms with a close contact to Pt
        - orientation_angle: Angle between fragment normal and COM-COM vector (if calculable)
    """
    # Calculate all pairwise distances
    dists = distances.distance_array(frag.positions, pt_atoms.positions)
    
    # Basic distance metrics
    min_dist = np.min(dists) if dists.size > 0 else np.inf
    avg_dist = np.mean(dists) if dists.size > 0 else np.nan
    med_dist = np.median(dists) if dists.size > 0 else np.nan
    std_dist = np.std(dists) if dists.size > 0 else np.nan
    
    # Fragment center of mass and distance to Pt COM
    frag_com = frag.center_of_mass()
    com_dist = np.linalg.norm(frag_com - pt_com)
    
    # Count close contacts and compute contact fraction
    num_close_contacts = np.sum(dists < close_threshold)
    
    # Calculate minimum distance per atom and count atoms in contact
    min_dists_per_atom = np.min(dists, axis=1) if dists.shape[0] > 0 else np.array([])
    atoms_in_contact = np.sum(min_dists_per_atom < close_threshold)
    contact_fraction = atoms_in_contact / len(frag) if len(frag) > 0 else 0.0
    
    # Calculate fragment orientation - may be None if can't compute
    try:
        orientation_angle = compute_orientation_angle(frag, pt_com)
    except Exception:
        orientation_angle = None
    
    return {
        "min_dist": min_dist,
        "avg_dist": avg_dist,
        "med_dist": med_dist,
        "std_dist": std_dist,
        "com_dist": com_dist,
        "num_close_contacts": int(num_close_contacts),
        "contact_fraction": contact_fraction,
        "orientation_angle": orientation_angle
    }

def compute_fragment_normal(frag: mda.core.groups.AtomGroup) -> np.ndarray:
    """
    Compute the normal vector of a fragment using SVD.
    
    Parameters
    ----------
    frag : mda.core.groups.AtomGroup
        Fragment atom group
        
    Returns
    -------
    np.ndarray
        Normal vector of the fragment
    
    Raises
    ------
    ValueError
        If normal vector cannot be computed (e.g., fragment has < 3 atoms)
    """
    if len(frag) < 3:
        raise ValueError("At least 3 atoms needed to compute normal vector")
    
    # Center the coordinates by subtracting the centroid
    coords = frag.positions
    centroid = np.mean(coords, axis=0)
    centered_coords = coords - centroid
    
    # Perform SVD
    u, s, vh = np.linalg.svd(centered_coords)
    
    # Normal vector is the right singular vector corresponding to the smallest singular value
    normal = vh[-1]
    
    # Normalize to unit length
    return normal / np.linalg.norm(normal)

def compute_orientation_angle(frag: mda.core.groups.AtomGroup, pt_com: np.ndarray) -> float:
    """
    Compute the orientation angle between fragment normal and COM-COM vector.
    
    Parameters
    ----------
    frag : mda.core.groups.AtomGroup
        Fragment atom group
    pt_com : np.ndarray
        Center of mass of the platinum nanoparticle
        
    Returns
    -------
    float
        Angle in degrees between fragment normal and COM-COM vector
        
    Raises
    ------
    ValueError
        If angle cannot be computed (e.g., insufficient atoms)
    """
    # Get fragment normal vector
    normal = compute_fragment_normal(frag)
    
    # Vector from fragment COM to Pt COM
    frag_com = frag.center_of_mass()
    com_vec = pt_com - frag_com
    
    # Normalize COM vector
    com_vec_norm = com_vec / np.linalg.norm(com_vec)
    
    # Calculate angle
    dot_product = np.dot(normal, com_vec_norm)
    # Clamp to [-1, 1] to handle numerical errors
    dot_product = np.clip(dot_product, -1.0, 1.0)
    angle_rad = np.arccos(dot_product)
    angle_deg = np.degrees(angle_rad)
    
    return angle_deg

def compute_contact_fraction(frag: mda.core.groups.AtomGroup, 
                             pt_atoms: mda.core.groups.AtomGroup, 
                             threshold: float = 2.5) -> float:
    """
    Compute the fraction of fragment atoms that are in contact with any Pt atom.
    
    Parameters
    ----------
    frag : mda.core.groups.AtomGroup
        Fragment atom group
    pt_atoms : mda.core.groups.AtomGroup
        Platinum atoms atom group
    threshold : float
        Distance threshold for close contacts in Angstroms (default: 2.5)
        
    Returns
    -------
    float
        Fraction of fragment atoms in contact with Pt (0.0 to 1.0)
    """
    dists = distances.distance_array(frag.positions, pt_atoms.positions)
    min_dists_per_atom = np.min(dists, axis=1) if dists.shape[0] > 0 else np.array([])
    atoms_in_contact = np.sum(min_dists_per_atom < threshold)
    return atoms_in_contact / len(frag) if len(frag) > 0 else 0.0

def analyze_run(run_dir: str, 
                dcd_file: str, 
                pdb_file: str, 
                psf_file: str, 
                outfile: Any, 
                cutoff: float = 2.5) -> Tuple[List[int], Dict[int, Dict[str, Any]], Dict[str, Any]]:
    """
    Analyze a simulation run and compute metrics for attached NEC fragments.
    
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
    outfile : Any
        File-like object for logging
    cutoff : float
        Distance threshold in Angstroms for identifying attached fragments (default: 2.5)
        
    Returns
    -------
    Tuple[List[int], Dict[int, Dict[str, Any]], Dict[str, Any]]
        attached_fragments: List of fragment indices meeting the criteria
        fragment_metrics: Dictionary mapping fragment index to computed metrics
        summary_stats: Dictionary of summary statistics
    """
    print(f"\nAnalyzing run in directory: {run_dir}")
    if hasattr(outfile, 'write'):
        outfile.write(f"\nAnalyzing run in directory: {run_dir}\n")
    
    # Load the universe
    dcd_path = os.path.join(run_dir, dcd_file)
    pdb_path = os.path.join(run_dir, pdb_file)
    psf_path = os.path.join(run_dir, psf_file)
    
    u = mda.Universe(psf_path, pdb_path)
    u.load_new(dcd_path)
    u.trajectory[-1]  # Go to last frame
    print(f"Using final frame: {u.trajectory.frame}")
    if hasattr(outfile, 'write'):
        outfile.write(f"Using final frame: {u.trajectory.frame}\n")
    
    # Get Pt atoms and COM
    pt_atoms = u.select_atoms("name PT*")
    print(f"Number of Pt atoms: {len(pt_atoms)}")
    if hasattr(outfile, 'write'):
        outfile.write(f"Number of Pt atoms: {len(pt_atoms)}\n")
    
    pt_com = pt_atoms.center_of_mass()
    print(f"Pt nanoparticle COM: {pt_com}")
    if hasattr(outfile, 'write'):
        outfile.write(f"Pt nanoparticle COM: {pt_com}\n")
    
    # Identify fragments
    fragments = list(u.atoms.fragments)
    print(f"Total number of fragments: {len(fragments)}")
    if hasattr(outfile, 'write'):
        outfile.write(f"Total number of fragments: {len(fragments)}\n")
    
    # Analyze attached fragments
    attached_fragments = []
    fragment_metrics = {}
    
    for frag_num, frag in enumerate(fragments):
        # Skip Pt-containing fragments
        if len(frag.select_atoms("name PT*")) > 0:
            continue
        
        # Check if any atom is within cutoff
        dists = distances.distance_array(frag.positions, pt_atoms.positions)
        min_dist = np.min(dists) if dists.size > 0 else np.inf
        
        if min_dist < cutoff:
            attached_fragments.append(frag_num)
            metrics = compute_fragment_metrics(frag, pt_atoms, pt_com, close_threshold=cutoff)
            fragment_metrics[frag_num] = metrics
            
            print(f"Fragment {frag_num}: min distance = {metrics['min_dist']:.2f} Å, " 
                  f"contact fraction = {metrics['contact_fraction']:.2f}")
            if hasattr(outfile, 'write'):
                outfile.write(f"Fragment {frag_num}: min distance = {metrics['min_dist']:.2f} Å, "
                              f"contact fraction = {metrics['contact_fraction']:.2f}\n")
    
    # Compute summary statistics
    avg_com_distance = np.mean([m["com_dist"] for m in fragment_metrics.values()]) if fragment_metrics else 0
    avg_min_distance = np.mean([m["min_dist"] for m in fragment_metrics.values()]) if fragment_metrics else 0
    avg_contact_fraction = np.mean([m["contact_fraction"] for m in fragment_metrics.values()]) if fragment_metrics else 0
    avg_orientation = np.mean([m["orientation_angle"] for m in fragment_metrics.values() 
                              if m["orientation_angle"] is not None]) if fragment_metrics else 0
    
    summary_stats = {
        "total_attached_fragments": len(attached_fragments),
        "avg_com_distance": avg_com_distance,
        "avg_min_distance": avg_min_distance,
        "avg_contact_fraction": avg_contact_fraction,
        "avg_orientation_angle": avg_orientation
    }
    
    print("\n--- Run Summary ---")
    print(f"Total attached fragments: {summary_stats['total_attached_fragments']}")
    print(f"Average fragment COM distance: {summary_stats['avg_com_distance']:.2f} Å")
    print(f"Average fragment min distance: {summary_stats['avg_min_distance']:.2f} Å")
    print(f"Average contact fraction: {summary_stats['avg_contact_fraction']:.2f}")
    print(f"Average orientation angle: {summary_stats['avg_orientation_angle']:.2f} degrees")
    
    if hasattr(outfile, 'write'):
        outfile.write("\n--- Run Summary ---\n")
        outfile.write(f"Total attached fragments: {summary_stats['total_attached_fragments']}\n")
        outfile.write(f"Average fragment COM distance: {summary_stats['avg_com_distance']:.2f} Å\n")
        outfile.write(f"Average fragment min distance: {summary_stats['avg_min_distance']:.2f} Å\n")
        outfile.write(f"Average contact fraction: {summary_stats['avg_contact_fraction']:.2f}\n")
        outfile.write(f"Average orientation angle: {summary_stats['avg_orientation_angle']:.2f} degrees\n")
    
    return attached_fragments, fragment_metrics, summary_stats