"""
docking/plots.py

Visualization utilities for molecular docking analysis.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def plot_delta_scores(
    df: pd.DataFrame,
    output_path: Path,
    reference_gene: str = "DYRK2",
    figsize: tuple = (8, 6)
):
    """
    Plot Delta Docking Scores as a horizontal bar chart.
    
    Visual Style:
    - High Affinity (Delta > 0): Teal/Green
    - Reference (Delta = 0): Black
    - Low Affinity (Delta < 0): Gray/Muted Red
    
    Args:
        df: DataFrame containing 'gene_name' and 'delta_score' columns.
            Should be sorted by delta_score descending.
        output_path: Path to save the figure.
        reference_gene: Name of the reference gene to highlight.
    """
    
    # Ensure sorted order for plotting (Top = Highest Affinity)
    # Reset index so plotting respects the order
    plot_df = df.copy()
    if not plot_df.index.is_monotonic_increasing:
        plot_df = plot_df.sort_values('delta_score', ascending=True) # Ascending for horiz bars (bottom-up)
    
    # Assign Colors
    colors = []
    for _, row in plot_df.iterrows():
        if row['gene_name'] == reference_gene:
            colors.append('#333333')  # Dark Gray/Black for Reference
        elif row['delta_score'] > 0:
            colors.append('#009E73')  # NPG Teal (Stronger)
        else:
            colors.append('#BDBDBD')  # Light Gray (Weaker)
            # Alternatively use a Muted Red like '#E64B35' if you want to emphasize "worse"
    
    plt.figure(figsize=figsize)
    
    # Plot Horizontal Bars
    bars = plt.barh(
        y=plot_df['gene_name'],
        width=plot_df['delta_score'],
        color=colors,
        edgecolor='none',
        height=0.7
    )
    
    # Add vertical reference line
    plt.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
    
    # Styling
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_color('#333333')
    ax.spines['bottom'].set_color('#333333')
    
    # Labels
    plt.xlabel(f'Δ Docking Score (vs {reference_gene})', fontsize=11, fontweight='bold', fontfamily='sans-serif')
    plt.title('Relative Molecular Docking Affinity', fontsize=14, fontweight='bold', fontfamily='sans-serif', pad=15)
    
    # Add annotations (optional, but helpful for papers)
    # "Stronger Binding" arrow or text
    max_val = plot_df['delta_score'].max()
    min_val = plot_df['delta_score'].min()
    
    # Text annotation for zones
    if max_val > 0:
        plt.text(max_val * 0.5, len(plot_df)-0.5, "Stronger Affinity", 
                 ha='center', va='center', color='#009E73', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved Figure 4A to: {output_path}")


def plot_enrichment_bubble(
    df: pd.DataFrame,
    output_path: Path,
    top_n: int = 15,
    figsize: tuple = (10, 8)
):
    """
    Plot functional enrichment as a Bubble Plot (Figure 4B).
    
    Dimensions:
    - X-axis: Gene Ratio (Recall)
    - Y-axis: Pathway Term
    - Size: Gene Count (Intersection Size)
    - Color: Significance (Adjusted P-value)
    
    Args:
        df: Enrichment results DataFrame (from g:Profiler).
        output_path: Path to save the figure.
        top_n: Number of top terms to display.
    """
    import numpy as np
    
    # 1. Prepare Data
    if df.empty:
        print("Warning: Enrichment DataFrame is empty. Skipping bubble plot.")
        return

    # Filter for relevant sources (GO:BP, KEGG, Reactome)
    # Adjust source names based on g:Profiler output format (usually 'GO:BP', 'KEGG', 'REAC')
    relevant_sources = ['GO:BP', 'KEGG', 'REAC']
    plot_df = df[df['Source'].isin(relevant_sources)].copy()
    
    if plot_df.empty:
        # Fallback if specific sources not found
        plot_df = df.copy()

    # Sort by significance and take top N
    plot_df = plot_df.sort_values('Adjusted_P_value', ascending=True).head(top_n)
    
    # Reverse for plotting (best at top)
    plot_df = plot_df.iloc[::-1]
    
    # 2. Plotting
    plt.figure(figsize=figsize)
    
    # Create scatter plot
    scatter = plt.scatter(
        x=plot_df['Recall'],          # X-axis: Gene Ratio
        y=plot_df['Term_Name'],       # Y-axis: Terms
        s=plot_df['Intersection_Size'] * 50,  # Size: Count (scaled)
        c=plot_df['Adjusted_P_value'],        # Color: P-value
        cmap='viridis_r',             # Reverse viridis (Dark = Low P = Good)
        alpha=0.8,
        edgecolors='white'
    )
    
    # 3. Styling
    plt.xlabel('Gene Ratio (Recall)', fontsize=12, fontweight='bold', fontfamily='sans-serif')
    plt.title(f'Functional Enrichment (Top {top_n})', fontsize=14, fontweight='bold', fontfamily='sans-serif', pad=15)
    
    # Add Colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label('Adjusted P-value', fontsize=10)
    
    # Add Size Legend (Dummy points)
    # Find min/max sizes for legend
    sizes = sorted(plot_df['Intersection_Size'].unique())
    if len(sizes) > 0:
        legend_sizes = [min(sizes), max(sizes)]
        if len(sizes) > 2:
            legend_sizes.insert(1, int(np.median(sizes)))
            
        legend_elements = []
        for s in legend_sizes:
            legend_elements.append(
                plt.scatter([], [], s=s*50, c='gray', alpha=0.6, label=str(s))
            )
        
        plt.legend(
            handles=legend_elements, 
            title="Gene Count",
            loc='lower right',
            labelspacing=1.0,
            frameon=True
        )

    # Grid and Spines
    ax = plt.gca()
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Saved Figure 4B to: {output_path}")
