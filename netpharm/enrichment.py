"""
Functional enrichment analysis using DAVID or g:Profiler.
"""

import pandas as pd
import os
from .utils.api_wrappers import query_gprofiler


class EnrichmentAnalyzer:
    """Handle functional enrichment analysis."""
    
    def __init__(self, logger):
        """
        Initialize EnrichmentAnalyzer.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger
        self.enrichment_results = None
    
    def analyze_david_manual(self, gene_list, output_dir="./data"):
        """
        Guide user through manual DAVID enrichment.
        
        Args:
            gene_list: List of gene names
            output_dir: Directory for downloaded files
        
        Returns:
            dict: Enrichment statistics
        """
        self.logger.info("\n" + "="*70)
        self.logger.info("[STEP 5] FUNCTIONAL ENRICHMENT ANALYSIS - DAVID (MANUAL)")
        self.logger.info("="*70)
        
        os.makedirs(output_dir, exist_ok=True)
        
        # Save gene list for user
        gene_list_file = os.path.join(output_dir, "gene_list_for_david.txt")
        with open(gene_list_file, 'w') as f:
            f.write('\n'.join(gene_list))
        
        self.logger.info(f"\n✓ Gene list saved to: {gene_list_file}")
        self.logger.info(f"   ({len(gene_list)} genes)")
        
        # DAVID instructions
        self.logger.info("\n--- DAVID Functional Annotation Tool ---")
        self.logger.info("Please follow these steps:")
        self.logger.info("\n1. Go to: https://david.ncifcrf.gov/")
        self.logger.info("2. Click 'Start Analysis'")
        self.logger.info("3. Paste gene list from the file above (or upload file)")
        self.logger.info("4. Select identifier: 'OFFICIAL_GENE_SYMBOL'")
        self.logger.info("5. Select list type: 'Gene List'")
        self.logger.info("6. Click 'Submit List'")
        self.logger.info("\n7. On results page, use 'Functional Annotation Tool'")
        self.logger.info("8. Download results:")
        self.logger.info("   - GO Biological Process")
        self.logger.info("   - KEGG Pathway")
        self.logger.info("   - Reactome Pathway")
        self.logger.info(f"\n9. Save downloaded files in: {output_dir}/")
        self.logger.info("   Suggested names:")
        self.logger.info("   - david_go_bp.csv")
        self.logger.info("   - david_kegg.csv")
        self.logger.info("   - david_reactome.csv")
        
        self.logger.info("\n📖 Default parameters used in paper:")
        self.logger.info("   - P-value threshold: 0.05")
        self.logger.info("   - Multiple testing correction: Default (Benjamini)")
        self.logger.info("   - Background: Whole genome")
        
        input("\n⏸ Press ENTER after downloading DAVID results...")
        
        # Try to load results
        stats = {'files_loaded': []}
        
        for filename in ['david_go_bp.csv', 'david_kegg.csv', 'david_reactome.csv']:
            filepath = os.path.join(output_dir, filename)
            if os.path.exists(filepath):
                try:
                    df = pd.read_csv(filepath)
                    stats['files_loaded'].append(filename)
                    self.logger.info(f"✓ Loaded: {filename} ({len(df)} terms)")
                except Exception as e:
                    self.logger.warning(f"⚠ Could not load {filename}: {str(e)}")
        
        if not stats['files_loaded']:
            self.logger.warning("\n⚠ No DAVID results loaded. Please check file names and formats.")
        else:
            self.logger.info(f"\n✓ Successfully loaded {len(stats['files_loaded'])} enrichment files")
        
        return stats
    
    def analyze_gprofiler(self, gene_list):
        """
        Perform enrichment analysis using g:Profiler API.
        
        Args:
            gene_list: List of gene names
        
        Returns:
            pd.DataFrame: Enrichment results
        """
        self.logger.info("\n" + "="*70)
        self.logger.info("[STEP 5] FUNCTIONAL ENRICHMENT ANALYSIS - g:Profiler (API)")
        self.logger.info("="*70)
        
        self.logger.info(f"\nQuerying g:Profiler for {len(gene_list)} genes...")
        self.logger.info("Databases: GO (BP/MF/CC), KEGG, Reactome")
        self.logger.info("Significance threshold: FDR < 0.05")
        
        try:
            response = query_gprofiler(gene_list, organism="hsapiens")
            
            if 'result' not in response or not response['result']:
                self.logger.error("\n❌ ERROR: No enrichment results found")
                self.logger.error("💡 This could mean:")
                self.logger.error("   - Gene list too small")
                self.logger.error("   - No significant enrichment")
                self.logger.error("   - Gene names not recognized")
                raise SystemExit(1)
            
            # Parse results
            results = response['result']
            
            enrichment_data = []
            for term in results:
                # FIX: Handle intersections robustly (ensure all items are strings)
                # The API might return lists of lists or non-string objects
                raw_intersections = term.get('intersections', [])
                cleaned_intersections = []
                
                for item in raw_intersections:
                    if isinstance(item, list):
                        # If nested list, flatten it or stringify elements
                        cleaned_intersections.extend([str(x) for x in item])
                    else:
                        cleaned_intersections.append(str(item))
                
                intersection_str = ','.join(cleaned_intersections)

                enrichment_data.append({
                    'Source': term['source'],
                    'Term_ID': term['native'],
                    'Term_Name': term['name'],
                    'P_value': term['p_value'],
                    'Adjusted_P_value': term.get('p_value', term['p_value']),  # FDR
                    'Term_Size': term['term_size'],
                    'Query_Size': term['query_size'],
                    'Intersection_Size': term['intersection_size'],
                    'Precision': term['precision'],
                    'Recall': term['recall'],
                    'Intersections': intersection_str
                })
            
            self.enrichment_results = pd.DataFrame(enrichment_data)
            self.enrichment_results = self.enrichment_results.sort_values('Adjusted_P_value')
            
            # Summary statistics
            source_counts = self.enrichment_results['Source'].value_counts()
            
            self.logger.info("\n" + "="*70)
            self.logger.info("ENRICHMENT ANALYSIS SUMMARY:")
            self.logger.info("="*70)
            self.logger.info(f"Total significant terms: {len(self.enrichment_results)}")
            self.logger.info("\nBy database:")
            for source, count in source_counts.items():
                self.logger.info(f"  {source}: {count} terms")
            
            # Top enriched terms
            self.logger.info("\nTop 10 enriched terms:")
            for idx, row in self.enrichment_results.head(10).iterrows():
                self.logger.info(f"  [{row['Source']}] {row['Term_Name']}")
                self.logger.info(f"     P-value: {row['P_value']:.2e}, Genes: {row['Intersection_Size']}/{row['Query_Size']}")
            
            return self.enrichment_results
            
        except Exception as e:
            self.logger.error(f"\n❌ ERROR: g:Profiler query failed")
            self.logger.error(f"   Reason: {str(e)}")
            raise SystemExit(1)
    
    def save_results(self, output_dir, method='gprofiler'):
        """
        Save enrichment results.
        
        Args:
            output_dir: Output directory
            method: 'gprofiler' or 'david'
        """
        os.makedirs(output_dir, exist_ok=True)
        
        if method == 'gprofiler' and self.enrichment_results is not None:
            # Save all results
            all_path = os.path.join(output_dir, "gprofiler_all_results.csv")
            self.enrichment_results.to_csv(all_path, index=False)
            self.logger.info(f"\n📄 All enrichment results saved to: {all_path}")
            
            # Save by source
            for source in self.enrichment_results['Source'].unique():
                source_df = self.enrichment_results[self.enrichment_results['Source'] == source]
                source_name = source.replace(':', '_').lower()
                source_path = os.path.join(output_dir, f"gprofiler_{source_name}.csv")
                source_df.to_csv(source_path, index=False)
                self.logger.info(f"📄 {source} results saved to: {source_path}")
            
            self.logger.info("\n   These files contain: Enriched pathways, p-values, gene overlaps")
        
        elif method == 'david':
            self.logger.info("\n📄 DAVID results should be in the data/ directory")
            self.logger.info("   Check for: david_go_bp.csv, david_kegg.csv, david_reactome.csv")