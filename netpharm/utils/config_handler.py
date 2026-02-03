"""
Configuration handling - loading from YAML or prompting user.
"""

import yaml
import os


def load_config(config_path):
    """
    Load configuration from YAML file.
    
    Args:
        config_path: Path to YAML config file
    
    Returns:
        dict: Configuration dictionary
    
    Raises:
        FileNotFoundError: If config file doesn't exist
        yaml.YAMLError: If config file is invalid
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    return config


def prompt_user_config():
    """
    Interactively prompt user for configuration.
    
    Returns:
        dict: Configuration dictionary
    """
    print("\n" + "="*70)
    print("NETWORK PHARMACOLOGY ANALYSIS - CONFIGURATION")
    print("="*70 + "\n")
    
    config = {}
    
    # Compound input
    print("--- COMPOUND INPUT ---")
    input_type = input("Enter compound by CID or SMILES? [cid/smiles]: ").strip().lower()
    
    if input_type == 'cid':
        cid = input("Enter PubChem CID: ").strip()
        config['compound'] = {'cid': cid}
    elif input_type == 'smiles':
        smiles = input("Enter SMILES string: ").strip()
        config['compound'] = {'smiles': smiles}
    else:
        raise ValueError("Invalid input type. Please enter 'cid' or 'smiles'.")
    
    # Target prediction thresholds
    print("\n--- TARGET PREDICTION THRESHOLDS ---")
    print("📄 Paper defaults: SwissTargetPrediction=0.0, SuperPred=0.5")

    while True:
        swiss_thresh_input = input("SwissTargetPrediction threshold [0.0]: ").strip()
        if not swiss_thresh_input:
            swiss_thresh = 0.0
            break
        try:
            swiss_thresh = float(swiss_thresh_input)
            break
        except ValueError:
            print("❌ Invalid SwissTargetPrediction threshold. Please enter a number.")

    while True:
        superpred_thresh_input = input("SuperPred threshold [0.5]: ").strip()
        if not superpred_thresh_input:
            superpred_thresh = 0.5
            break
        try:
            superpred_thresh = float(superpred_thresh_input)
            break
        except ValueError:
            print("❌ Invalid SuperPred threshold. Please enter a number.")
    
    config['target_prediction'] = {
        'swiss_threshold': swiss_thresh,
        'superpred_threshold': superpred_thresh
    }
    
    # Pathway selection
    print("\n--- PATHWAY SELECTION ---")
    print("💡 Find pathways at: https://reactome.org/PathwayBrowser/")
    print("   You can enter keywords or specific pathway IDs (R-HSA-XXXXXX)")

    while True:
        pathway_input = input("\nEnter pathway keywords (comma-separated) or IDs: ").strip()
        if pathway_input:
            pathway_list = [p.strip() for p in pathway_input.split(',') if p.strip()]
            if pathway_list:
                break
        print("❌ Pathway input cannot be empty. Please provide at least one keyword or ID.")
    
    config['pathways'] = {'search_terms': pathway_list}
    
    # STRING network parameters
    print("\n--- STRING NETWORK PARAMETERS ---")
    while True:
        string_conf_input = input("STRING confidence threshold [0.700]: ").strip()
        if not string_conf_input:
            string_conf = 0.700
            break
        try:
            string_conf = float(string_conf_input)
            break
        except ValueError:
            print("❌ Invalid STRING confidence threshold. Please enter a number.")
    
    config['string'] = {'confidence': string_conf}
    
    # Enrichment analysis
    print("\n--- ENRICHMENT ANALYSIS ---")
    enrichment_method = input("Use DAVID (manual) or g:Profiler (API)? [david/gprofiler]: ").strip().lower()
    
    config['enrichment'] = {'method': enrichment_method if enrichment_method in ['david', 'gprofiler'] else 'david'}
    
    # Output directory
    print("\n--- OUTPUT SETTINGS ---")
    output_dir = input("Output directory [./outputs]: ").strip()
    output_dir = output_dir if output_dir else './outputs'
    
    config['output'] = {'directory': output_dir}
    
    print("\n" + "="*70)
    print("✓ Configuration complete!")
    print("="*70 + "\n")
    
    return config


def get_config(config_path=None):
    """
    Get configuration either from file or interactive prompt.
    
    Args:
        config_path: Optional path to config file
    
    Returns:
        dict: Configuration dictionary
    """
    if config_path:
        if os.path.exists(config_path):
            print(f"Loading configuration from: {config_path}")
            return load_config(config_path)
        raise FileNotFoundError(f"Config file not found: {config_path}")

    return prompt_user_config()
