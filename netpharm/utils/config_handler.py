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
        print("\n❌ Invalid input type. Please restart and enter 'cid' or 'smiles'.")
        exit(1)
    
    # Target prediction thresholds
    print("\n--- TARGET PREDICTION THRESHOLDS ---")
    print("📄 Paper defaults: SwissTargetPrediction=0.0, SuperPred=0.5")
    
    swiss_thresh = input("SwissTargetPrediction threshold [0.0]: ").strip()
    swiss_thresh = float(swiss_thresh) if swiss_thresh else 0.0
    
    superpred_thresh = input("SuperPred threshold [0.5]: ").strip()
    superpred_thresh = float(superpred_thresh) if superpred_thresh else 0.5
    
    config['target_prediction'] = {
        'swiss_threshold': swiss_thresh,
        'superpred_threshold': superpred_thresh
    }
    
    # Pathway selection
    print("\n--- PATHWAY SELECTION ---")
    print("💡 Find pathways at: https://reactome.org/PathwayBrowser/")
    print("   You can enter keywords or specific pathway IDs (R-HSA-XXXXXX)")
    
    pathway_input = input("\nEnter pathway keywords (comma-separated) or IDs: ").strip()
    pathway_list = [p.strip() for p in pathway_input.split(',')]
    
    config['pathways'] = {'search_terms': pathway_list}
    
    # STRING network parameters
    print("\n--- STRING NETWORK PARAMETERS ---")
    string_conf = input("STRING confidence threshold [0.700]: ").strip()
    string_conf = float(string_conf) if string_conf else 0.700
    
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
    if config_path and os.path.exists(config_path):
        print(f"Loading configuration from: {config_path}")
        return load_config(config_path)
    else:
        return prompt_user_config()
