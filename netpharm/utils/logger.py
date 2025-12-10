"""
Logging configuration for the network pharmacology tool.
"""

import logging
import os
from datetime import datetime


def setup_logger(output_dir, compound_id):
    """
    Set up logging configuration.
    
    Args:
        output_dir: Directory to save log file
        compound_id: Compound identifier for log filename
    
    Returns:
        Logger instance
    """
    # Create logger
    logger = logging.getLogger('netpharm')
    logger.setLevel(logging.DEBUG)
    
    # Avoid duplicate handlers
    if logger.handlers:
        return logger
    
    # Console handler (INFO level)
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)
    
    # File handler (DEBUG level)
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(output_dir, f'{compound_id}_{timestamp}.log')
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_format = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler.setFormatter(file_format)
    
    # Add handlers
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)
    
    logger.info(f"Log file created: {log_file}")
    
    return logger
