"""
Network Pharmacology - Replicating Specific Paper Figures
Generates:
1. Cytoscape-ready CSV files for Figures A, B, C.
2. Enrichment Bubble Plots with SMART LABEL STACKING (No Overlap).
"""

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Optional import for adjustment (we will use manual stacking primarily)
try:
    from adjustText import adjust_text
except ImportError:
    adjust_text = None

# Ensure we can import netpharm from the parent directory
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from netpharm import NetworkPharmacology

def plot_enrichment_bubble_grouped(enrichment_df, title, output_path):
    """
    Creates a publication-quality grouped bubble plot with SMART LABEL STACKING.
    
    Specs:
    - X-Axis: Source Categories (GO, KEGG, REAC) with JITTER
    - Y-Axis: Significance (-log10 P-value)
    - Labels: Uses a 'Stacking' algorithm to force labels apart vertically
              if they share the same P-value.
    """
    if enrichment_df is None or len(enrichment_df) == 0:
        print(f"   ⚠ No enrichment results to plot for {title}")
        return

    # 1. Prepare Data
    df = enrichment_df.copy()
    
    # Map sources to Main Categories for X-axis
    source_map = {
        'GO:BP': 'GO', 'GO:CC': 'GO', 'GO:MF': 'GO',
        'KEGG': 'KEGG',
        'REAC': 'REAC'
    }
    df['Category'] = df['Source'].map(source_map)
    df = df[df['Category'].notna()].copy()
    
    # Calculate significance
    df['-log10(P)'] = -1 * np.log10(df['P_value'])
    
    # 2. Select Top 20 per Category
    final_df = pd.DataFrame()
    category_order = ['GO', 'KEGG', 'REAC']
    
    for cat in category_order:
        subset = df[df['Category'] == cat].sort_values('P_value').head(20)
        final_df = pd.concat([final_df, subset])
    
    if len(final_df) == 0:
        return

    # 3. Add Jitter to X-Axis (Safe spread 0.3 to prevent column overlap)
    cat_to_num = {cat: i for i, cat in enumerate(category_order)}
    final_df['x_base'] = final_df['Category'].map(cat_to_num)
    
    np.random.seed(42) 
    jitter_strength = 0.3
    final_df['x_jittered'] = final_df['x_base'] + np.random.uniform(-jitter_strength, jitter_strength, size=len(final_df))

    # 4. Plotting
    plt.figure(figsize=(12, 10)) # Increased height slightly for stacking space
    sns.set_style("whitegrid")
    
    palette = {'GO': '#4C72B0', 'KEGG': '#55A868', 'REAC': '#C44E52'}
    
    # Draw Scatter
    sns.scatterplot(
        data=final_df,
        x="x_jittered",
        y="-log10(P)",
        size="Intersection_Size",
        hue="Category",
        palette=palette,
        sizes=(50, 600),
        alpha=0.75,
        edgecolor="black",
        linewidth=1
    )
    
    plt.xticks(ticks=range(len(category_order)), labels=category_order, fontsize=12, fontweight='bold')
    plt.xlabel("Source Database", fontsize=12)
    plt.ylabel("-log10(P-value)", fontsize=12)
    plt.xlim(-0.5, 2.5)
    
    # Calculate safe vertical spacing (10% of total Y range)
    y_min, y_max = final_df['-log10(P)'].min(), final_df['-log10(P)'].max()
    y_range = y_max - y_min if y_max != y_min else 1.0
    safe_spacing = y_range * 0.10  # 10% spacing cushion
    
    # 5. Add Labels with SMART STACKING
    texts = []
    
    for cat in category_order:
        # Get top 3 most significant
        top3 = final_df[final_df['Category'] == cat].sort_values('P_value').head(3)
        
        # Prepare list of labels to place for this category
        labels_to_place = []
        for idx, row in top3.iterrows():
            # Label Text Formatting
            if cat == 'GO':
                label_text = f"[{row['Source']}] {row['Term_Name']}"
            else:
                label_text = row['Term_Name']
            
            if len(label_text) > 40:
                label_text = label_text[:37] + "..."
            
            labels_to_place.append({
                'text': label_text,
                'x': row['x_jittered'],
                'y_origin': row['-log10(P)'],
                'y_target': row['-log10(P)'] # Will be updated if collision
            })
        
        # Sort by Y position (ascending) so we stack from bottom to top
        labels_to_place.sort(key=lambda item: item['y_origin'])
        
        # STACKING LOGIC: Force separation
        for i in range(1, len(labels_to_place)):
            prev_label = labels_to_place[i-1]
            curr_label = labels_to_place[i]
            
            # If current is too close to previous (or below it due to sorting tie), push it up
            if curr_label['y_target'] < prev_label['y_target'] + safe_spacing:
                curr_label['y_target'] = prev_label['y_target'] + safe_spacing
        
        # Draw the stacked labels
        for label in labels_to_place:
            t = plt.annotate(
                text=label['text'],
                xy=(label['x'], label['y_origin']),        # Arrow points to Bubble
                xytext=(label['x'], label['y_target']),    # Text sits at stacked position
                fontsize=9,
                fontweight='bold',
                color='black',
                ha='center',
                arrowprops=dict(
                    arrowstyle='-', 
                    color='gray', 
                    lw=0.5,
                    shrinkA=0, 
                    shrinkB=5
                )
            )
            texts.append(t)

    # Use adjust_text only as a final polish if available
    if adjust_text:
        try:
            adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
        except Exception:
            pass

    plt.title(f"{title}", fontsize=16, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', title="Gene Count")
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"   📊 Grouped Bubble plot saved: {output_path}")

def main():
    print("="*70)
    print("REPLICATING FIGURES A, B, & C")
    print("======================================================================")
    
    # --- 1. SETUP PATHS ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(script_dir, '..', 'data', 'curcumin')
    
    swiss_path = os.path.join(data_dir, "swiss_targets.csv")
    superpred_path = os.path.join(data_dir, "superpred_targets.csv")
    
    if not os.path.exists(swiss_path):
        print(f"❌ Critical Error: Data files not found in: {data_dir}")
        print("   Please run: mkdir -p data/curcumin && cp curcumin_results/.../*.csv data/curcumin/")
        return

    # --- 2. INITIALIZE PIPELINE ---
    npharma = NetworkPharmacology(
        cid=969516,
        output_base='./paper_replication_results' 
    )
    npharma.get_compound_info()
    
    # --- 3. LOAD TARGETS ---
    print("\n[1] Loading bundled target data...")
    try:
        npharma.target_predictor.swiss_targets = pd.read_csv(swiss_path)
        npharma.target_predictor.superpred_targets = pd.read_csv(superpred_path)
        
        # Gene name parsing fix
        if 'Common name' in npharma.target_predictor.swiss_targets.columns:
            npharma.target_predictor.swiss_targets['gene_list'] = \
                npharma.target_predictor.swiss_targets['Common name'].astype(str).str.split()
            npharma.target_predictor.swiss_targets = \
                npharma.target_predictor.swiss_targets.explode('gene_list')
            npharma.target_predictor.swiss_targets['Target'] = \
                npharma.target_predictor.swiss_targets['gene_list']

        npharma.target_genes = npharma.target_predictor.get_all_target_genes()
        print(f"   ✅ Successfully loaded {len(npharma.target_genes)} total targets.")
    except Exception as e:
        print(f"❌ Error loading data: {e}")
        return

    # --- 4. DEFINE THE 3 NETWORKS ---
    car_manual_genes = [
        'AKT1', 'NFKB1', 'STAT3', 'RAF1', 'IKBKB', 
        'CHUK', 'IKBKG', 'BCL2', 'MAPK1', 'RELA'
    ]

    figures_config = {
        "Figure_A_Cytokine": {"type": "pathway", "id": "R-HSA-1280215"},
        "Figure_B_Adaptive": {"type": "pathway", "id": "R-HSA-1280218"},
        "Figure_C_CAR":      {"type": "manual",  "genes": car_manual_genes}
    }
    
    # --- 5. GENERATE FIGURES ---
    print("\n[2] Generating Individual Figures...")
    
    pathway_ids = [v['id'] for k, v in figures_config.items() if v['type'] == 'pathway']
    npharma.analyze_pathways(pathway_ids)
    
    for fig_name, config in figures_config.items():
        print(f"\n--- Processing: {fig_name} ---")
        
        # Determine gene list
        if config['type'] == 'pathway':
            subset_df = npharma.overlapping_targets[
                npharma.overlapping_targets['pathway_id'] == config['id']
            ]
            genes = subset_df['gene_name'].unique().tolist()
        else:
            genes = [g for g in config['genes'] if g in npharma.target_genes]
            print(f"   (Using manual gene list for CAR Signaling)")
        
        if not genes:
            print(f"   ❌ No overlapping genes found.")
            continue
            
        print(f"   Found {len(genes)} genes.")
        
        # Create output directory
        fig_output_dir = os.path.join(npharma.output_dir, fig_name)
        os.makedirs(fig_output_dir, exist_ok=True)
        
        # A. Build Network
        try:
            npharma.network_analyzer.query_string_network(genes, confidence=0.700)
            network = npharma.network_analyzer.build_network()
            npharma.network_analyzer.save_results(fig_output_dir)
            # npharma.visualizer.create_all_visualizations(network, fig_output_dir)
            print(f"   ✅ Cytoscape CSV saved.")
        except Exception as e:
            print(f"   ⚠ Network error: {e}")

        # B. Enrichment & Bubble Plot
        try:
            results = npharma.enrichment_analyzer.analyze_gprofiler(genes)
            npharma.enrichment_analyzer.save_results(fig_output_dir)
            
            # Generate the Custom Grouped JITTERED Bubble Plot
            plot_path = os.path.join(fig_output_dir, "enrichment_bubble_grouped.png")
            plot_enrichment_bubble_grouped(results, fig_name.replace('_', ' '), plot_path)
            
        except Exception as e:
            print(f"   ⚠ Enrichment error: {e}")

    print("\n" + "="*70)
    print("✅ REPLICATION COMPLETE!")
    print(f"Results are in: {npharma.output_dir}")

if __name__ == '__main__':
    main()