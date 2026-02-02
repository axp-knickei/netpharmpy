"""
replicate_figure4.py

Script to reproduce the logic of Figure 4:
1. Load docking scores.
2. Normalize against DYRK2 reference.
3. Filter for high-affinity targets.
4. Run enrichment on filtered targets.
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))

from netpharm.docking.processor import DockingProcessor
from netpharm.docking.plots import plot_delta_scores, plot_enrichment_bubble
from netpharm.docking.adapters import parse_david_results
from netpharm.enrichment import EnrichmentAnalyzer

from netpharm.docking.processor import DockingProcessor
from netpharm.docking.plots import plot_delta_scores, plot_enrichment_bubble
from netpharm.docking.adapters import parse_david_results
from netpharm.enrichment import EnrichmentAnalyzer

# Simple logger for standalone script
class SimpleLogger:
    def info(self, msg): print(msg)
    def warning(self, msg): print(f"WARNING: {msg}")
    def error(self, msg): print(f"ERROR: {msg}")

def main():
    # Setup Output
    output_dir = Path("outputs/figure4_replication")
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Setup paths
    data_path = Path("data/reproduction/figure4_scores.csv")
    david_data_path = Path("data/reproduction/david_results") # Pre-calculated/Mock results
    
    print(f"--- Loading Docking Data from {data_path} ---")
    processor = DockingProcessor(data_path)
    
    # 2. Normalize and Rank
    print("\n--- Normalizing against DYRK2 ---")
    try:
        df_ranked = processor.normalize_scores(reference_gene="DYRK2")
        print(df_ranked[['gene_name', 'docking_score', 'delta_score']])
    except Exception as e:
        print(f"Error: {e}")
        return

    # 3. Filter High Affinity
    print("\n--- Filtering High Affinity Targets (Delta > 0) ---")
    targets = processor.get_high_affinity_targets(threshold=0.0)
    print(f"Found {len(targets)} targets stronger than DYRK2:")
    print(targets)
    
    # 4. Generate Figure 4A (Bar Plot)
    print("\n--- Generating Figure 4A (Biophysical Rank) ---")
    plot_delta_scores(
        df_ranked, 
        output_path=output_dir / "figure_4a_delta_docking.png",
        reference_gene="DYRK2"
    )
    
    # 5. Enrichment Analysis (DAVID)
    print("\n--- Running Enrichment Analysis (DAVID) ---")
    logger = SimpleLogger()
    enricher = EnrichmentAnalyzer(logger)
    
    if targets:
        # A. Generate the gene list for the user
        print("\n[Action] Generating gene list for DAVID...")
        # Use non_interactive=True to avoid blocking in automated script
        enricher.analyze_david_manual(
            targets, 
            output_dir=output_dir / "david_input",
            non_interactive=True
        )
        
        # B. Parse Results (Simulating the user having run DAVID and saved files)
        print(f"\n[Mock] Loading pre-calculated MOCK DAVID results from {david_data_path}...")
        try:
            enrichment_results = parse_david_results(david_data_path, prefix="mock_")
            
            if not enrichment_results.empty:
                print(f"Loaded {len(enrichment_results)} significant terms.")
                
                # C. Plot Bubble Chart
                print("\n--- Generating Figure 4B (Enrichment Bubble) ---")
                plot_enrichment_bubble(
                    enrichment_results,
                    output_path=output_dir / "figure_4b_enrichment_bubble.png",
                    top_n=15
                )
            else:
                print("No significant enrichment terms found in DAVID files.")
                
        except Exception as e:
            print(f"Error processing DAVID results: {e}")
            
    else:
        print("No high affinity targets found. Skipping enrichment.")

if __name__ == "__main__":
    main()

