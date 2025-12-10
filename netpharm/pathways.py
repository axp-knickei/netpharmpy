"""
Reactome pathway analysis and target overlap identification.
"""

import pandas as pd
import os
from .utils.api_wrappers import query_reactome, get_pathway_proteins
from .utils.validators import validate_pathway_id


class PathwayAnalyzer:
    """Handle Reactome pathway analysis."""
    
    def __init__(self, logger):
        """
        Initialize PathwayAnalyzer.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger
        self.pathways = []
        self.pathway_proteins = None
        self.overlapping_targets = None
    
    def search_pathways(self, search_terms):
        """
        Search Reactome for pathways.
        
        Args:
            search_terms: List of search keywords or pathway IDs
        
        Returns:
            list: List of pathway dictionaries
        """
        self.logger.info("\n" + "="*70)
        self.logger.info("[STEP 3] REACTOME PATHWAY ANALYSIS")
        self.logger.info("="*70)
        
        all_pathways = []
        
        for term in search_terms:
            term = term.strip()
            self.logger.info(f"\nSearching for: '{term}'")
            
            # Check if it's a pathway ID or keyword
            try:
                if term.startswith('R-'):
                    validate_pathway_id(term)
                    # Direct pathway ID - get its info
                    self.logger.info(f"  Using pathway ID directly: {term}")
                    all_pathways.append({
                        'pathway_id': term,
                        'pathway_name': f'Pathway {term}',
                        'source': 'direct_id'
                    })
                else:
                    # Keyword search
                    results = query_reactome(term)
                    if results:
                        self.logger.info(f"  Found {len(results)} pathways")
                        for pathway in results[:5]:  # Top 5 results
                            all_pathways.append({
                                'pathway_id': pathway['stId'],
                                'pathway_name': pathway['name'],
                                'source': term
                            })
                    else:
                        self.logger.warning(f"  No pathways found for: '{term}'")
            except Exception as e:
                self.logger.error(f"  Error searching '{term}': {str(e)}")
        
        if not all_pathways:
            self.logger.error("\n❌ ERROR: No pathways found")
            self.logger.error("💡 Suggestions:")
            self.logger.error("   - Try broader keywords (e.g., 'immune system', 'signaling')")
            self.logger.error("   - Browse pathways at: https://reactome.org/PathwayBrowser/")
            self.logger.error("   - Use specific pathway IDs (e.g., R-HSA-1280218)")
            raise SystemExit(1)
        
        # Remove duplicates
        unique_pathways = []
        seen_ids = set()
        for p in all_pathways:
            if p['pathway_id'] not in seen_ids:
                unique_pathways.append(p)
                seen_ids.add(p['pathway_id'])
        
        self.pathways = unique_pathways
        
        self.logger.info(f"\n✓ Total unique pathways found: {len(self.pathways)}")
        for pathway in self.pathways:
            self.logger.info(f"  - {pathway['pathway_id']}: {pathway['pathway_name']}")
        
        return self.pathways
    
    def extract_pathway_proteins(self):
        """
        Extract proteins from identified pathways.
        
        Returns:
            pd.DataFrame: Pathway proteins
        """
        self.logger.info("\nExtracting proteins from pathways...")
        
        all_proteins = []
        
        for pathway in self.pathways:
            pathway_id = pathway['pathway_id']
            pathway_name = pathway['pathway_name']
            
            self.logger.info(f"  Processing: {pathway_id}")
            
            try:
                proteins = get_pathway_proteins(pathway_id)
                
                for protein in proteins:
                    protein['pathway_id'] = pathway_id
                    protein['pathway_name'] = pathway_name
                    all_proteins.append(protein)
                
                self.logger.info(f"    Found {len(proteins)} proteins")
                
            except Exception as e:
                self.logger.warning(f"    Error: {str(e)}")
        
        if not all_proteins:
            self.logger.error("\n❌ ERROR: No proteins extracted from pathways")
            raise SystemExit(1)
        
        self.pathway_proteins = pd.DataFrame(all_proteins)
        self.pathway_proteins = self.pathway_proteins.drop_duplicates(
            subset=['gene_name', 'pathway_id']
        )
        
        unique_genes = self.pathway_proteins['gene_name'].nunique()
        self.logger.info(f"\n✓ Extracted {unique_genes} unique proteins from pathways")
        
        return self.pathway_proteins
    
    def find_overlapping_targets(self, target_genes):
        """
        Find overlapping targets between predictions and pathways.
        
        Args:
            target_genes: List of predicted target gene names
        
        Returns:
            pd.DataFrame: Overlapping targets
        """
        if self.pathway_proteins is None:
            raise ValueError("Must extract pathway proteins first")
        
        self.logger.info("\nFinding overlapping targets...")
        
        # Normalize gene names
        target_genes_upper = set([g.upper() for g in target_genes])
        pathway_genes_upper = set(self.pathway_proteins['gene_name'].str.upper())
        
        overlapping_genes = target_genes_upper.intersection(pathway_genes_upper)
        
        if not overlapping_genes:
            self.logger.error("\n❌ ERROR: No overlapping targets found")
            self.logger.error("   Predicted targets don't match pathway proteins")
            self.logger.error("💡 This could mean:")
            self.logger.error("   - The compound doesn't target these pathways")
            self.logger.error("   - Try different pathways")
            self.logger.error("   - Check gene name formats")
            raise SystemExit(1)
        
        # Filter pathway proteins to overlapping ones
        self.overlapping_targets = self.pathway_proteins[
            self.pathway_proteins['gene_name'].str.upper().isin(overlapping_genes)
        ]
        
        self.logger.info("\n" + "="*70)
        self.logger.info("PATHWAY OVERLAP SUMMARY:")
        self.logger.info("="*70)
        self.logger.info(f"Predicted targets: {len(target_genes)}")
        self.logger.info(f"Pathway proteins: {len(pathway_genes_upper)}")
        self.logger.info(f"🎯 Overlapping targets: {len(overlapping_genes)}")
        self.logger.info("\nOverlapping genes:")
        for gene in sorted(overlapping_genes):
            self.logger.info(f"  - {gene}")
        
        return self.overlapping_targets
    
    def save_results(self, output_dir):
        """
        Save pathway analysis results.
        
        Args:
            output_dir: Output directory
        """
        os.makedirs(output_dir, exist_ok=True)
        
        # Save pathway list
        pathways_df = pd.DataFrame(self.pathways)
        pathways_path = os.path.join(output_dir, "pathways_searched.csv")
        pathways_df.to_csv(pathways_path, index=False)
        self.logger.info(f"\n📄 Pathway list saved to: {pathways_path}")
        self.logger.info("   Contains: Pathway IDs, names, and search sources")
        
        # Save pathway proteins
        if self.pathway_proteins is not None:
            proteins_path = os.path.join(output_dir, "pathway_proteins.csv")
            self.pathway_proteins.to_csv(proteins_path, index=False)
            self.logger.info(f"📄 Pathway proteins saved to: {proteins_path}")
            self.logger.info("   Contains: Gene names, UniProt IDs, pathway associations")
        
        # Save overlapping targets
        if self.overlapping_targets is not None:
            overlap_path = os.path.join(output_dir, "overlapping_targets.csv")
            self.overlapping_targets.to_csv(overlap_path, index=False)
            self.logger.info(f"📄 Overlapping targets saved to: {overlap_path}")
            self.logger.info("   Contains: Genes present in both predictions and pathways")
