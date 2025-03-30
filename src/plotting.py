#!/usr/bin/env python
"""
Plotting module for PtNEC system analysis.
This module contains functions for visualizing analysis results.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Any, Optional

def set_plotting_style() -> None:
    """Set consistent plot style for all figures."""
    plt.style.use('seaborn-v0_8-whitegrid')
    plt.rcParams['figure.figsize'] = [10, 6]
    plt.rcParams['figure.dpi'] = 100
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['font.size'] = 12
    plt.rcParams['axes.titlesize'] = 14
    plt.rcParams['axes.labelsize'] = 12
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['figure.titlesize'] = 16

def create_plots(system_name: str, df: pd.DataFrame, system_folder: str) -> None:
    """
    Create individual plots for a system's fragment metrics.
    
    Parameters
    ----------
    system_name : str
        Name of the system
    df : pd.DataFrame
        DataFrame containing fragment metrics
    system_folder : str
        Folder where plots should be saved
    """
    plots_dir = os.path.join(system_folder, "plots")
    
    # Set consistent style
    set_plotting_style()
    
    # Plot 1: Histogram of minimum distances
    plt.figure()
    plt.hist(df['min_dist'], bins=20, edgecolor='black', alpha=0.7)
    plt.xlabel("Minimum Distance (Å)")
    plt.ylabel("Frequency")
    plt.title(f"{system_name}: Histogram of Minimum Distances")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = os.path.join(plots_dir, f"{system_name}_min_distance_hist.png")
    plt.savefig(filename)
    plt.close()
    
    # Plot 2: Histogram of COM distances
    plt.figure()
    plt.hist(df['com_dist'], bins=20, edgecolor='black', alpha=0.7)
    plt.xlabel("Fragment COM Distance to Pt COM (Å)")
    plt.ylabel("Frequency")
    plt.title(f"{system_name}: Histogram of Fragment COM Distances")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = os.path.join(plots_dir, f"{system_name}_com_distance_hist.png")
    plt.savefig(filename)
    plt.close()
    
    # Plot 3: Scatter plot of min vs COM distances
    plt.figure()
    plt.scatter(df['min_dist'], df['com_dist'], alpha=0.7, edgecolor='black')
    plt.xlabel("Minimum Distance (Å)")
    plt.ylabel("Fragment COM Distance (Å)")
    plt.title(f"{system_name}: Min vs COM Distance")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = os.path.join(plots_dir, f"{system_name}_scatter_min_vs_com.png")
    plt.savefig(filename)
    plt.close()
    
    # Plot 4: Histogram of contact fractions
    if 'contact_fraction' in df.columns:
        plt.figure()
        plt.hist(df['contact_fraction'], bins=20, edgecolor='black', alpha=0.7)
        plt.xlabel("Contact Fraction")
        plt.ylabel("Frequency")
        plt.title(f"{system_name}: Histogram of Contact Fractions")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        filename = os.path.join(plots_dir, f"{system_name}_contact_fraction_hist.png")
        plt.savefig(filename)
        plt.close()
    
    # Plot 5: Scatter plot of contact fraction vs COM distance
    if 'contact_fraction' in df.columns:
        plt.figure()
        plt.scatter(df['contact_fraction'], df['com_dist'], alpha=0.7, edgecolor='black')
        plt.xlabel("Contact Fraction")
        plt.ylabel("Fragment COM Distance (Å)")
        plt.title(f"{system_name}: Contact Fraction vs COM Distance")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        filename = os.path.join(plots_dir, f"{system_name}_scatter_contact_vs_com.png")
        plt.savefig(filename)
        plt.close()
    
    # Plot 6: Histogram of orientation angles (if available)
    if 'orientation_angle' in df.columns:
        # Filter out None values
        orientation_data = df['orientation_angle'].dropna()
        if not orientation_data.empty:
            plt.figure()
            plt.hist(orientation_data, bins=20, edgecolor='black', alpha=0.7)
            plt.xlabel("Orientation Angle (degrees)")
            plt.ylabel("Frequency")
            plt.title(f"{system_name}: Histogram of Orientation Angles")
            plt.grid(True, alpha=0.3)
            plt.tight_layout()
            filename = os.path.join(plots_dir, f"{system_name}_orientation_angle_hist.png")
            plt.savefig(filename)
            plt.close()
            
            # Plot 7: Scatter plot of orientation angle vs contact fraction
            if 'contact_fraction' in df.columns:
                # Filter out rows with None orientation angle
                filtered_df = df.dropna(subset=['orientation_angle'])
                if not filtered_df.empty:
                    plt.figure()
                    plt.scatter(filtered_df['orientation_angle'], filtered_df['contact_fraction'], 
                                alpha=0.7, edgecolor='black')
                    plt.xlabel("Orientation Angle (degrees)")
                    plt.ylabel("Contact Fraction")
                    plt.title(f"{system_name}: Orientation Angle vs Contact Fraction")
                    plt.grid(True, alpha=0.3)
                    plt.tight_layout()
                    filename = os.path.join(plots_dir, f"{system_name}_scatter_orientation_vs_contact.png")
                    plt.savefig(filename)
                    plt.close()

def create_individual_rdf_plot(system_name: str, 
                               bins: np.ndarray, 
                               rdf_values: np.ndarray, 
                               system_folder: str) -> None:
    """
    Create an individual RDF plot for a system.
    
    Parameters
    ----------
    system_name : str
        Name of the system
    bins : np.ndarray
        Bin centers for the RDF
    rdf_values : np.ndarray
        RDF values
    system_folder : str
        Folder where the plot should be saved
    """
    plots_dir = os.path.join(system_folder, "plots")
    
    # Set consistent style
    set_plotting_style()
    
    plt.figure()
    plt.plot(bins, rdf_values, marker='o', linestyle='-', markersize=4, alpha=0.7)
    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")
    plt.title(f"{system_name}: COM-Based RDF")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = os.path.join(plots_dir, f"{system_name}_rdf_plot.png")
    plt.savefig(filename)
    plt.close()

def create_comparative_plots(summary_list: List[Dict[str, Any]], 
                             system_dfs: Dict[str, pd.DataFrame], 
                             combined_folder: str) -> None:
    """
    Create comparative plots across systems.
    
    Parameters
    ----------
    summary_list : List[Dict[str, Any]]
        List of dictionaries containing summary stats for each system
    system_dfs : Dict[str, pd.DataFrame]
        Dictionary mapping system names to DataFrames with fragment metrics
    combined_folder : str
        Folder where the plots should be saved
    """
    # Set consistent style
    set_plotting_style()
    
    # Extract system names
    systems = [d["system"] for d in summary_list]
    
    # Plot 1: Bar chart of attached fragments
    attached_counts = [d["total_attached_fragments"] for d in summary_list]
    plt.figure()
    plt.bar(systems, attached_counts, color='skyblue', edgecolor='black', alpha=0.7)
    plt.xlabel("System")
    plt.ylabel("Total Attached Fragments")
    plt.title("Total Attached Fragments per System")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = os.path.join(combined_folder, "comparative_attached_fragments.png")
    plt.savefig(filename)
    plt.close()
    
    # Plot 2: Box plot of COM distances
    com_data = [system_dfs[system]["com_dist"].values for system in systems]
    plt.figure()
    plt.boxplot(com_data, labels=systems)
    plt.xlabel("System")
    plt.ylabel("Fragment COM Distance (Å)")
    plt.title("Distribution of Fragment COM Distances")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = os.path.join(combined_folder, "comparative_com_distance_boxplot.png")
    plt.savefig(filename)
    plt.close()
    
    # Plot 3: Box plot of minimum distances
    min_data = [system_dfs[system]["min_dist"].values for system in systems]
    plt.figure()
    plt.boxplot(min_data, labels=systems)
    plt.xlabel("System")
    plt.ylabel("Minimum Distance (Å)")
    plt.title("Distribution of Minimum Distances")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = os.path.join(combined_folder, "comparative_min_distance_boxplot.png")
    plt.savefig(filename)
    plt.close()
    
    # Plot 4: Combined scatter plot of min vs COM distances
    colors = ['red', 'green', 'blue', 'orange']
    plt.figure()
    for i, system in enumerate(systems):
        df = system_dfs[system]
        plt.scatter(df["min_dist"], df["com_dist"], color=colors[i % len(colors)], 
                    alpha=0.6, label=system, edgecolor='black')
    plt.xlabel("Minimum Distance (Å)")
    plt.ylabel("Fragment COM Distance (Å)")
    plt.title("Min vs COM Distance Across Systems")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = os.path.join(combined_folder, "comparative_min_vs_com_scatter.png")
    plt.savefig(filename)
    plt.close()
    
    # Plot 5: Box plot of contact fractions (if available)
    if all('contact_fraction' in system_dfs[system].columns for system in systems):
        contact_data = [system_dfs[system]["contact_fraction"].values for system in systems]
        plt.figure()
        plt.boxplot(contact_data, labels=systems)
        plt.xlabel("System")
        plt.ylabel("Contact Fraction")
        plt.title("Distribution of Contact Fractions")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        filename = os.path.join(combined_folder, "comparative_contact_fraction_boxplot.png")
        plt.savefig(filename)
        plt.close()
    
    # Plot 6: Bar chart of average contact fractions (if available)
    if all('avg_contact_fraction' in d for d in summary_list):
        avg_contact = [d.get('avg_contact_fraction', 0) for d in summary_list]
        plt.figure()
        plt.bar(systems, avg_contact, color='skyblue', edgecolor='black', alpha=0.7)
        plt.xlabel("System")
        plt.ylabel("Average Contact Fraction")
        plt.title("Average Contact Fraction per System")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        filename = os.path.join(combined_folder, "comparative_avg_contact_fraction.png")
        plt.savefig(filename)
        plt.close()
    
    # Plot 7: Bar chart of average orientation angles (if available)
    if all('avg_orientation_angle' in d for d in summary_list):
        avg_orientation = [d.get('avg_orientation_angle', 0) for d in summary_list]
        plt.figure()
        plt.bar(systems, avg_orientation, color='skyblue', edgecolor='black', alpha=0.7)
        plt.xlabel("System")
        plt.ylabel("Average Orientation Angle (degrees)")
        plt.title("Average Orientation Angle per System")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        filename = os.path.join(combined_folder, "comparative_avg_orientation_angle.png")
        plt.savefig(filename)
        plt.close()

def create_comparative_rdf_plot(system_rdf: Dict[str, Tuple[np.ndarray, np.ndarray]], 
                                combined_folder: str) -> None:
    """
    Create a comparative RDF plot across systems.
    
    Parameters
    ----------
    system_rdf : Dict[str, Tuple[np.ndarray, np.ndarray]]
        Dictionary mapping system names to (bins, rdf_values) tuples
    combined_folder : str
        Folder where the plot should be saved
    """
    # Set consistent style
    set_plotting_style()
    
    plt.figure()
    colors = ['red', 'green', 'blue', 'orange']
    
    for i, (system, (bins, rdf_values)) in enumerate(system_rdf.items()):
        plt.plot(bins, rdf_values, label=system, color=colors[i % len(colors)], 
                 alpha=0.7, linewidth=2)
    
    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")
    plt.title("COM-Based RDF Comparison")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    filename = os.path.join(combined_folder, "comparative_rdf_overlay.png")
    plt.savefig(filename)
    plt.close()

def save_fragment_metrics_csv(system_name: str, 
                              fragment_metrics: Dict[int, Dict[str, Any]], 
                              system_folder: str) -> pd.DataFrame:
    """
    Save fragment metrics to a CSV file.
    
    Parameters
    ----------
    system_name : str
        Name of the system
    fragment_metrics : Dict[int, Dict[str, Any]]
        Dictionary mapping fragment indices to metric dictionaries
    system_folder : str
        Folder where the CSV should be saved
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing the fragment metrics
    """
    # Create a list of dictionaries, one per fragment
    data = []
    for frag_num, metrics in fragment_metrics.items():
        row = {"fragment_num": frag_num}
        row.update(metrics)
        data.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(data)
    
    # Save to CSV
    csv_filename = os.path.join(system_folder, f"{system_name}_fragment_metrics.csv")
    df.to_csv(csv_filename, index=False)
    
    return df

def save_rdf_csv(system_name: str, 
                 bins: np.ndarray, 
                 rdf_values: np.ndarray, 
                 system_folder: str) -> pd.DataFrame:
    """
    Save RDF data to a CSV file.
    
    Parameters
    ----------
    system_name : str
        Name of the system
    bins : np.ndarray
        Bin centers for the RDF
    rdf_values : np.ndarray
        RDF values
    system_folder : str
        Folder where the CSV should be saved
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing the RDF data
    """
    # Create DataFrame
    df = pd.DataFrame({"r (Å)": bins, "g(r)": rdf_values})
    
    # Save to CSV
    csv_filename = os.path.join(system_folder, f"{system_name}_rdf.csv")
    df.to_csv(csv_filename, index=False)
    
    return df

def save_summary_csv(summary_list: List[Dict[str, Any]], combined_folder: str) -> pd.DataFrame:
    """
    Save summary statistics to a CSV file.
    
    Parameters
    ----------
    summary_list : List[Dict[str, Any]]
        List of dictionaries containing summary stats for each system
    combined_folder : str
        Folder where the CSV should be saved
        
    Returns
    -------
    pd.DataFrame
        DataFrame containing the summary stats
    """
    # Create DataFrame
    df = pd.DataFrame(summary_list)
    
    # Save to CSV
    csv_filename = os.path.join(combined_folder, "summary.csv")
    df.to_csv(csv_filename, index=False)
    
    return df