#!/usr/bin/env python
"""
Utility functions for the PtNEC MD analysis package.
"""

import os
import shutil
import logging
import yaml
from typing import Dict, Any, List, Optional, Tuple, Union, TextIO

def setup_logging(log_file: Optional[str] = None, level: int = logging.INFO) -> logging.Logger:
    """
    Set up logging for the analysis.
    
    Parameters
    ----------
    log_file : Optional[str]
        Path to the log file. If None, logs only to console.
    level : int
        Logging level (default: logging.INFO)
        
    Returns
    -------
    logging.Logger
        Configured logger instance
    """
    logger = logging.getLogger('pt_nec_analysis')
    logger.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    
    # Console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)
    
    return logger

def load_config(config_path: str) -> Dict[str, Any]:
    """
    Load configuration from a YAML file.
    
    Parameters
    ----------
    config_path : str
        Path to the configuration file.
        
    Returns
    -------
    Dict[str, Any]
        Configuration dictionary
    """
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    return config

def create_dir_structure(base_dir: str, subfolders: List[str]) -> None:
    """
    Create a directory structure with the specified base directory and subfolders.
    
    Parameters
    ----------
    base_dir : str
        Base directory path
    subfolders : List[str]
        List of subdirectories to create
    """
    os.makedirs(base_dir, exist_ok=True)
    for folder in subfolders:
        path = os.path.join(base_dir, folder)
        os.makedirs(path, exist_ok=True)

def create_base_analysis_folder(analysis_folder: str, recreate: bool = True) -> str:
    """
    Create the base analysis folder and combined subfolder.
    
    Parameters
    ----------
    analysis_folder : str
        Path to the analysis folder
    recreate : bool
        If True, delete the existing folder if it exists (default: True)
        
    Returns
    -------
    str
        Path to the created analysis folder
    """
    if recreate and os.path.exists(analysis_folder):
        shutil.rmtree(analysis_folder)
    
    os.makedirs(analysis_folder, exist_ok=True)
    combined_dir = os.path.join(analysis_folder, "combined")
    os.makedirs(combined_dir, exist_ok=True)
    
    return analysis_folder

def create_system_folder(analysis_folder: str, system_name: str) -> str:
    """
    Create a subdirectory for a system inside the analysis folder.
    Also create a 'plots' subdirectory for that system.
    
    Parameters
    ----------
    analysis_folder : str
        Path to the base analysis folder
    system_name : str
        Name of the system
        
    Returns
    -------
    str
        Path to the created system folder
    """
    sys_dir = os.path.join(analysis_folder, system_name)
    os.makedirs(sys_dir, exist_ok=True)
    plots_dir = os.path.join(sys_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    return sys_dir

def log_message(message: str, logger: Optional[logging.Logger] = None, 
                file_handle: Optional[TextIO] = None, level: str = "info") -> None:
    """
    Log a message to a logger and/or file handle.
    
    Parameters
    ----------
    message : str
        The message to log
    logger : Optional[logging.Logger]
        Logger instance to use
    file_handle : Optional[TextIO]
        File handle to write to
    level : str
        Log level (info, warning, error, critical, debug)
    """
    if logger:
        if level == "info":
            logger.info(message)
        elif level == "warning":
            logger.warning(message)
        elif level == "error":
            logger.error(message)
        elif level == "critical":
            logger.critical(message)
        elif level == "debug":
            logger.debug(message)
    
    # Also print to console and/or file if specified
    print(message)
    if file_handle:
        file_handle.write(message + "\n")