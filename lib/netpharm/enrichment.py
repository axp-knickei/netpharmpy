# lib/netpharm/enrichment.py

from pathlib import Path
from gprofiler import GProfiler
import pandas as pd


def run_functional_enrichment(
    gene_list,
    output_dir,
    organism="hsapiens",
    fdr_threshold=0.05,
    label="gprofiler",
):
    """
    Run functional enrichment using g:Profiler.

    Parameters
    ----------
    gene_list : list[str]
        Gene symbols to analyze.
    output_dir : str or Path
        Directory where CSV results will be saved.
    organism : str
        g:Profiler organism code.
    fdr_threshold : float
        Adjusted p-value cutoff.
    label : str
        Prefix for output files.

    Returns
    -------
    pandas.DataFrame
        Full enrichment results table.
    """

    if not gene_list:
        raise ValueError("Gene list for enrichment is empty.")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nQuerying g:Profiler for {len(gene_list)} genes...")
    print("Databases: GO, KEGG, Reactome")
    print(f"Significance threshold: FDR < {fdr_threshold}")

    gp = GProfiler(return_dataframe=True)

    results = gp.profile(
        organism=organism,
        query=gene_list,
        sources=["GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"],
        user_threshold=fdr_threshold,
        significance_threshold_method="fdr",
    )

    if results.empty:
        print("⚠️ No significant enrichment found.")
        return results

    # Save full table
    full_path = output_dir / f"{label}_all_results.csv"
    results.to_csv(full_path, index=False)

    # Save per-database tables
    for source in results["source"].unique():
        df = results[results["source"] == source]
        out = output_dir / f"{label}_{source.lower().replace(':', '_')}.csv"
        df.to_csv(out, index=False)

    print(f"\nTotal significant terms: {len(results)}")

    return results
