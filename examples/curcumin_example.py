"""
Network Pharmacology Tool - Full Curcumin Example
Replicates the analysis from the paper.
"""

from netpharm import NetworkPharmacology

def main():
    """
    Complete curcumin analysis following the paper methodology.
    
    Reference:
    Limsakul et al. (2025). Immunomodulatory Effects of Curcumin 
    on CAR T-Cell Therapy. Antioxidants 14(4):454.
    """
    
    print("="*70)
    print("CURCUMIN NETWORK PHARMACOLOGY ANALYSIS")
    print("Replicating: Limsakul et al. (2025)")
    print("="*70)
    
    # Initialize with curcumin CID
    print("\nInitializing pipeline for curcumin (CID: 969516)...")
    npharma = NetworkPharmacology(
        cid=969516,
        output_base='./curcumin_results'
    )
    
    # Configuration matching the paper
    config = {
        'target_prediction': {
            'swiss_threshold': 0.0,
            'superpred_threshold': 0.5
        },
        'pathways': {
            'search_terms': [
                "cytokine signaling immune system",
                "adaptive immune system",
                "signaling by interleukins",
                "intracellular signaling second messengers"
            ]
        },
        'string': {
            'confidence': 0.700
        },
        'enrichment': {
            'method': 'gprofiler'
        }
    }
    
    npharma.config = config
    
    # Run full pipeline
    print("\n📊 Starting full analysis pipeline...")
    print("Note: Manual steps will require your interaction\n")
    
    npharma.run_full_pipeline()
    
    print("\n" + "="*70)
    print("✅ CURCUMIN ANALYSIS COMPLETE!")
    print("="*70)
    print(f"\nResults directory: {npharma.output_dir}")
    print("\nExpected findings (from paper):")
    print("  - ~149 predicted targets")
    print("  - Key overlapping proteins: CHUK, IKBKB, IKBKG (IKK complex)")
    print("  - Pathways: NF-κB signaling, MAPK/ERK, PI3K/AKT")
    print("  - Enriched in: Cytokine signaling, inflammatory response")
    print("\nCompare your results with the paper!")


if __name__ == '__main__':
    main()
