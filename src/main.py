#!/usr/bin/env python
"""
Main module for PtNEC MD analysis.
This module ties together the fragment analysis, RDF analysis, and plotting modules.
"""

import os
import sys
import argparse
import yaml
import logging
from datetime import datetime
import numpy as np
import pandas as pd

# Import local modules
from . import utils
from . import fragment_analysis
from . import rdf_analysis
from . import plotting

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='PtNEC MD Analysis')
    parser.add_argument('--config', type=str, default='config/analysis_config.yaml',
                        help='Path to configuration file')
    parser.add_argument('--recreate', action='store_true',
                        help='Recreate analysis folders if they exist')
    parser.add_argument('--systems', nargs='+', 
                        help='Specific systems to analyze (e.g., 0HPt 4HPt)')
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug logging')
    
    return parser.parse_args()

def main():
    """Main function to run the analysis workflow."""
    # Parse command line arguments
    args = parse_arguments()
    
    # Set up logging
    log_level = logging.DEBUG if args.debug else logging.INFO
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = f"analysis_{timestamp}.log"
    logger = utils.setup_logging(log_file, level=log_level)
    
    logger.info("Starting PtNEC MD analysis")
    
    # Load configuration
    config_path = args.config
    try:
        config = utils.load_config(config_path)
        logger.info(f"Loaded configuration from {config_path}")
    except Exception as e:
        logger.error(f"Error loading configuration: {e}")
        return 1
    
    # Get data and analysis directories
    base_dir = config.get('data_dir', "")
    base_analysis = config.get('analysis_dir', "analysis")
    
    # Validate data directory
    if not os.path.exists(base_dir):
        logger.error(f"Data directory {base_dir} does not exist")
        return 1
    
    # Create analysis folder structure
    utils.create_base_analysis_folder(base_analysis, recreate=args.recreate)
    logger.info(f"Created base analysis folder at {base_analysis}")
    
    # Get systems to analyze
    all_systems = config.get('systems', [])
    if args.systems:
        # Filter to only the specified systems
        systems = [s for s in all_systems if s['name'] in args.systems]
        if not systems:
            logger.warning(f"None of the specified systems {args.systems} found in config")
            return 1
    else:
        systems = all_systems
    
    # Analysis parameters
    contact_threshold = config.get('fragment_analysis', {}).get('contact_threshold', 2.5)
    rdf_range = tuple(config.get('rdf_analysis', {}).get('rdf_range', [0, 20]))
    nbins = config.get('rdf_analysis', {}).get('nbins', 200)
    
    # Process each system
    summary_list = []
    system_dfs = {}
    system_rdf = {}  # To store COM-based RDF data per system
    
    for system_info in systems:
        system_name = system_info['name']
        dcd_file = system_info['dcd_file']
        pdb_file = system_info['pdb_file']
        psf_file = system_info['psf_file']
        
        logger.info(f"Processing system: {system_name}")
        
        # Create system folder
        sys_dir = utils.create_system_folder(base_analysis, system_name)
        log_file_path = os.path.join(sys_dir, f"{system_name}_analysis.log")
        
        # Process the system
        with open(log_file_path, "w") as outfile:
            run_path = os.path.join(base_dir, system_name)
            outfile.write("=" * 50 + "\n")
            utils.log_message(f"Analyzing system: {system_name}", logger, outfile)
            
            try:
                # Fragment analysis
                utils.log_message("Running fragment analysis...", logger, outfile)
                attached_fragments, fragment_metrics, summary_stats = fragment_analysis.analyze_run(
                    run_path, dcd_file, pdb_file, psf_file, outfile, cutoff=contact_threshold)
                
                # Add system name to summary stats
                summary_stats["system"] = system_name
                summary_list.append(summary_stats)
                
                # Save fragment metrics and get DataFrame
                utils.log_message("Saving fragment metrics...", logger, outfile)
                df = plotting.save_fragment_metrics_csv(system_name, fragment_metrics, sys_dir)
                system_dfs[system_name] = df
                
                # Create individual plots
                utils.log_message("Creating individual plots...", logger, outfile)
                plotting.create_plots(system_name, df, sys_dir)
                
                # RDF analysis
                try:
                    utils.log_message("Running RDF analysis...", logger, outfile)
                    bins, rdf_values = rdf_analysis.compute_com_rdf(
                        run_path, dcd_file, pdb_file, psf_file,
                        rdf_range=rdf_range, nbins=nbins, cutoff=contact_threshold)
                    
                    # Save RDF data
                    plotting.save_rdf_csv(system_name, bins, rdf_values, sys_dir)
                    system_rdf[system_name] = (bins, rdf_values)
                    
                    # Create individual RDF plot
                    plotting.create_individual_rdf_plot(system_name, bins, rdf_values, sys_dir)
                    
                    # Compute RDF statistics
                    rdf_stats = rdf_analysis.compute_rdf_statistics(rdf_values, bins)
                    utils.log_message(f"RDF statistics: {rdf_stats}", logger, outfile)
                    
                except Exception as e:
                    utils.log_message(f"Error in RDF analysis: {e}", logger, outfile, level="error")
                
            except Exception as e:
                utils.log_message(f"Error processing system {system_name}: {e}", 
                                logger, outfile, level="error")
                continue
            
            utils.log_message(f"Completed analysis for system: {system_name}", logger, outfile)
            outfile.write("=" * 50 + "\n")
    
    # Create combined analysis
    combined_folder = os.path.join(base_analysis, "combined")
    try:
        logger.info("Creating comparative analysis...")
        
        # Save summary CSV
        plotting.save_summary_csv(summary_list, combined_folder)
        
        # Create comparative plots
        plotting.create_comparative_plots(summary_list, system_dfs, combined_folder)
        
        # Create comparative RDF plot
        if system_rdf:
            plotting.create_comparative_rdf_plot(system_rdf, combined_folder)
    
    except Exception as e:
        logger.error(f"Error in comparative analysis: {e}")
    
    logger.info(f"Analysis complete. Output saved to {base_analysis}")
    return 0

if __name__ == '__main__':
    sys.exit(main())