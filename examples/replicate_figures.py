"""
Network Pharmacology - Replicating Specific Paper Figures
Generates:
1. Cytoscape-ready CSV files (string_interactions.csv) for Figures A, B, C.
2. Enrichment Bubble Plots (The 'Right Graph') for Figures A, B, C.
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Ensure we can import netpharm from the parent directory
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from netpharm import NetworkPharmacology

def plot_enrichment_bubble(enrichment_df, title, output_path):
    """
    Creates a publication-quality bubble plot for enrichment results.
    Matches the 'Right Graph' style in the paper.
    """
    if enrichment_df is None or len(enrichment_df) == 0:
        print(f"   ⚠ No enrichment results to plot for {title}")
        return

    # Filter for significant terms (p < 0.05) and take top 20
    df = enrichment_df[enrichment_df['P_value'] < 0.05].copy()
    if len(df) == 0:
        return
        
    top_terms = df.sort_values('P_value').head(20).copy()
    
    # Calculate -log10(p-value) for the color scale
    top_terms['-log10(P)'] = -1 * np.log10(top_terms['P_value'])
    
    plt.figure(figsize=(10, 8))
    sns.set_style("whitegrid")
    
    # Create bubble plot
    try:
        scatter = sns.scatterplot(
            data=top_terms,
            x="Intersection_Size",    # X-axis: Gene Count
            y="Term_Name",            # Y-axis: Pathway Name
            size="Intersection_Size", # Bubble Size
            hue="-log10(P)",          # Color
            palette="viridis",
            sizes=(50, 400),
            edgecolor="black",
            linewidth=0.5
        )
        
        plt.title(f"{title}", fontsize=14, fontweight='bold')
        plt.xlabel("Gene Count", fontsize=12)
        plt.ylabel("")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0)
        plt.tight_layout()
        
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"   📊 Scatter plot saved: {output_path}")
    except Exception as e:
        print(f"   ⚠ Could not create plot: {e}")

def main():
    print("="*70)
    print("REPLICATING FIGURES A, B, & C")
    print("======================================================================")
    
    # --- 1. SETUP PATHS (PORTABLE) ---
    # This finds the 'data' folder relative to where this script is located
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, '..', 'data', 'curcumin')
    
    swiss_path = os.path.join(data_dir, "swiss_targets.csv")
    superpred_path = os.path.join(data_dir, "superpred_targets.csv")
    
    # Verify data exists before starting
    if not os.path.exists(swiss_path):
        print(f"❌ Critical Error: Data files not found in: {data_dir}")
        print("   Did you run the 'cp' commands to copy your results there?")
        return

    # --- 2. INITIALIZE PIPELINE ---
    npharma = NetworkPharmacology(
        cid=969516,
        output_base='./paper_replication_results' 
    )
    
    print("\n[1] Getting Compound Info...")
    npharma.get_compound_info()
    
    # --- 3. LOAD TARGETS ---
    print("\n[2] Loading bundled target data...")
    try:
        # Load CSVs
        npharma.target_predictor.swiss_targets = pd.read_csv(swiss_path)
        npharma.target_predictor.superpred_targets = pd.read_csv(superpred_path)
        
        # Apply the fix for splitting gene names (e.g. "IKBKG IKBKB CHUK")
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

    # --- 4. DEFINE THE 3 NETWORKS (FIGURES) ---
    figures_config = {
        "Figure_A_Cytokine": "R-HSA-1280215", # Cytokine Signaling in Immune system
        "Figure_B_Adaptive": "R-HSA-1280218", # Adaptive Immune System
        "Figure_C_CAR":      "R-HSA-9664407"  # Signaling by CAR-dependent T cell activation
    }
    
    # --- 5. GENERATE NETWORKS & PLOTS ---
    print("\n[3] Generating Individual Figures...")
    
    # First, analyze pathways to identify overlaps
    npharma.analyze_pathways(list(figures_config.values()))
    
    for figure_name, pathway_id in figures_config.items():
        print(f"\n--- Processing: {figure_name} ({pathway_id}) ---")
        
        # Filter for targets in this specific pathway
        subset_df = npharma.overlapping_targets[
            npharma.overlapping_targets['pathway_id'] == pathway_id
        ]
        
        specific_genes = subset_df['gene_name'].unique().tolist()
        
        if not specific_genes:
            print(f"   ❌ No overlapping genes found. (Pathway ID might be too specific or no targets match)")
            continue
            
        print(f"   Found {len(specific_genes)} genes.")
        
        # Create output directory
        fig_output_dir = os.path.join(npharma.output_dir, figure_name)
        os.makedirs(fig_output_dir, exist_ok=True)
        
        # A. Build Network (For Cytoscape)
        print("   Building Network...")
        try:
            npharma.network_analyzer.query_string_network(specific_genes, confidence=0.700)
            network = npharma.network_analyzer.build_network()
            npharma.network_analyzer.save_results(fig_output_dir)
            # Create basic visualizations just in case
            npharma.visualizer.create_all_visualizations(network, fig_output_dir)
            print(f"   ✅ Cytoscape files saved to: {fig_output_dir}")
        except Exception as e:
            print(f"   ⚠ Network generation error: {e}")

        # B. Enrichment & Bubble Plot
        print("   Running Enrichment & Plotting...")
        try:
            results = npharma.enrichment_analyzer.analyze_gprofiler(specific_genes)
            npharma.enrichment_analyzer.save_results(fig_output_dir)
            
            # Generate the Bubble Plot
            plot_path = os.path.join(fig_output_dir, "enrichment_bubble_plot.png")
            plot_enrichment_bubble(results, figure_name.replace('_', ' '), plot_path)
            
        except Exception as e:
            print(f"   ⚠ Enrichment error: {e}")

    print("\n" + "="*70)
    print("✅ REPLICATION COMPLETE!")
    print(f"Results are in: {npharma.output_dir}")

if __name__ == '__main__':
    main()