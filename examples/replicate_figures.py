"""
Network Pharmacology - Replicating Specific Paper Figures
Generates distinct networks for:
A. Cytokine Signaling
B. Adaptive Immune System
C. CAR Signaling
"""

import os
import pandas as pd
import sys

# Add parent directory to path to ensure we can import netpharm
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from netpharm import NetworkPharmacology

def main():
    print("="*70)
    print("REPLICATING FIGURES A, B, & C")
    print("======================================================================")
    
    # --- 1. SETUP PATHS (PORTABLE) ---
    # Determine where this script is running
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Define the relative path to the bundled data
    # Looks for: netpharmpy/data/curcumin
    data_dir = os.path.join(script_dir, '..', 'data', 'curcumin')
    
    # Define paths to the CSV files
    swiss_path = os.path.join(data_dir, "swiss_targets.csv")
    superpred_path = os.path.join(data_dir, "superpred_targets.csv")
    
    # Verify data exists
    if not os.path.exists(swiss_path) or not os.path.exists(superpred_path):
        print(f"❌ Critical Error: Data files not found in {data_dir}")
        print("   Please ensure 'swiss_targets.csv' and 'superpred_targets.csv'")
        print("   are in the 'netpharmpy/data/curcumin' directory.")
        return

    # --- 2. INITIALIZE PIPELINE ---
    # Create a new output directory for these specific figures
    npharma = NetworkPharmacology(
        cid=969516,
        output_base='./paper_replication_results' 
    )
    
    # Get Compound Info (Step 1)
    print("\n[1] Getting Compound Info...")
    npharma.get_compound_info()
    
    # --- 3. LOAD TARGETS (MANUAL INJECTION) ---
    print("\n[2] Loading bundled target data...")
    print(f"   Source: {data_dir}")
    
    # Manually load files into the target predictor
    try:
        npharma.target_predictor.swiss_targets = pd.read_csv(swiss_path)
        npharma.target_predictor.superpred_targets = pd.read_csv(superpred_path)
        
        # Apply the fix for parsing gene names (splitting 'IKBKG IKBKB CHUK')
        if 'Common name' in npharma.target_predictor.swiss_targets.columns:
            npharma.target_predictor.swiss_targets['gene_list'] = \
                npharma.target_predictor.swiss_targets['Common name'].astype(str).str.split()
            npharma.target_predictor.swiss_targets = \
                npharma.target_predictor.swiss_targets.explode('gene_list')
            npharma.target_predictor.swiss_targets['Target'] = \
                npharma.target_predictor.swiss_targets['gene_list']

        # Generate the full list of targets
        npharma.target_genes = npharma.target_predictor.get_all_target_genes()
        print(f"   ✅ Successfully loaded {len(npharma.target_genes)} total targets.")
        
    except Exception as e:
        print(f"❌ Error loading data: {e}")
        return

    # --- 4. DEFINE NETWORKS TO REPLICATE ---
    figures_config = {
        "Figure_A_Cytokine_Signaling": "R-HSA-1280215", # Cytokine Signaling in Immune system
        "Figure_B_Adaptive_Immune":    "R-HSA-1280218", # Adaptive Immune System
        "Figure_C_CAR_Signaling":      "R-HSA-9664407"  # Signaling by CAR-dependent...
    }
    
    # --- 5. RUN ANALYSIS FOR EACH FIGURE ---
    print("\n[3] Generating Individual Networks...")
    
    # Search pathways to find overlaps
    npharma.analyze_pathways(list(figures_config.values()))
    
    for figure_name, pathway_id in figures_config.items():
        print(f"\n--- Processing: {figure_name} ({pathway_id}) ---")
        
        # Filter: Which of our targets are in this specific pathway?
        subset_df = npharma.overlapping_targets[
            npharma.overlapping_targets['pathway_id'] == pathway_id
        ]
        
        specific_genes = subset_df['gene_name'].unique().tolist()
        
        if not specific_genes:
            print(f"   ❌ No overlapping genes found for {figure_name}")
            continue
            
        print(f"   Found {len(specific_genes)} genes.")
        
        # Create output directory
        fig_output_dir = os.path.join(npharma.output_dir, figure_name)
        os.makedirs(fig_output_dir, exist_ok=True)
        
        # Build Network
        print("   Building Network...")
        try:
            # Query STRING only for these specific genes
            npharma.network_analyzer.query_string_network(specific_genes, confidence=0.700)
            network = npharma.network_analyzer.build_network()
            
            # Save and Visualize
            npharma.network_analyzer.save_results(fig_output_dir)
            npharma.visualizer.create_all_visualizations(network, fig_output_dir)
            print(f"   ✅ Network saved to: {fig_output_dir}")
        except Exception as e:
            print(f"   ⚠ Network generation skipped: {e}")

        # Enrichment Analysis
        print("   Running Enrichment...")
        try:
            npharma.enrichment_analyzer.analyze_gprofiler(specific_genes)
            npharma.enrichment_analyzer.save_results(fig_output_dir)
            print(f"   ✅ Enrichment saved.")
        except Exception as e:
            print(f"   ⚠ Enrichment skipped: {e}")

    print("\n" + "="*70)
    print("✅ REPLICATION COMPLETE!")
    print(f"Results are in: {npharma.output_dir}")

if __name__ == '__main__':
    main()