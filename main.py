#!/usr/bin/env python3
"""
Network Pharmacology Tool - Main Entry Point

Usage:
    python main.py [--config path/to/config.yaml]
"""

import sys
import argparse
from netpharm import NetworkPharmacology
from netpharm.utils.config_handler import get_config


def main():
    """Main entry point for the network pharmacology tool."""
    
    parser = argparse.ArgumentParser(
        description='Network Pharmacology Analysis Tool v1.0.0',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Interactive mode (prompts for all inputs)
    python main.py
    
    # Using configuration file
    python main.py --config config.yaml
    
    # Specify output directory
    python main.py --config config.yaml --output ./my_results

For documentation, visit: [Your GitHub URL]
        """
    )
    
    parser.add_argument(
        '--config',
        type=str,
        help='Path to YAML configuration file (optional)'
    )
    
    parser.add_argument(
        '--output',
        type=str,
        default='./outputs',
        help='Base output directory (default: ./outputs)'
    )
    
    # DAVID Integration Arguments
    parser.add_argument(
        '--run_david',
        action='store_true',
        help='Run DAVID functional enrichment (requires registered email)'
    )
    
    parser.add_argument(
        '--david_email',
        type=str,
        help='Email address registered with DAVID Web Service'
    )
    
    args = parser.parse_args()
    
    try:
        # Get configuration (file or interactive)
        config = get_config(args.config)
        
        # Inject DAVID settings into config
        if args.run_david:
            if 'enrichment' not in config:
                config['enrichment'] = {}
            config['enrichment']['run_david'] = True
            if args.david_email:
                config['enrichment']['david_email'] = args.david_email
        
        # Extract compound info
        compound_config = config.get('compound', {})
        cid = compound_config.get('cid')
        smiles = compound_config.get('smiles')
        
        # Initialize pipeline
        npharma = NetworkPharmacology(
            cid=cid,
            smiles=smiles,
            config=config,
            output_base=args.output
        )
        
        # Run full pipeline
        npharma.run_full_pipeline()
        
        return 0
        
    except KeyboardInterrupt:
        print("\n\n⚠ Analysis interrupted by user")
        return 1
    except SystemExit as e:
        return e.code if isinstance(e.code, int) else 1
    except Exception as e:
        print(f"\n❌ Fatal error: {str(e)}")
        return 1


if __name__ == '__main__':
    sys.exit(main())
