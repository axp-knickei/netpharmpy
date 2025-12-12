"""
Reactome pathway analysis and target overlap identification.
"""

import pandas as pd
import os
import time
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
                    # Direct pathway ID
                    try:
                        validate_pathway_id(term)
                        self.logger.info(f"  Using pathway ID directly: {term}")
                        all_pathways.append({
                            'pathway_id': term,
                            'pathway_name': f'Pathway {term}',
                            'source': 'direct_id'
                        })
                    except ValueError as e:
                        self.logger.warning(f"  Invalid pathway ID format: {e}")
                        continue
                else:
                    # Keyword search with retry logic
                    max_retries = 3
                    retry_delay = 2
                    
                    for attempt in range(max_retries):
                        try:
                            results = query_reactome(term)
                            
                            if results:
                                self.logger.info(f"  Found {len(results)} pathways")
                                
                                for pathway in results[:5]:  # Top 5 results
                                    # Robust handling of API response
                                    pathway_id = pathway.get('stId') or pathway.get('dbId') or pathway.get('id')
                                    pathway_name = pathway.get('name') or pathway.get('displayName') or 'Unknown'
                                    
                                    if pathway_id:
                                        all_pathways.append({
                                            'pathway_id': str(pathway_id),
                                            'pathway_name': pathway_name,
                                            'source': term
                                        })
                                    else:
                                        self.logger.warning(f"    Skipping pathway without valid ID: {pathway}")
                                
                                break  # Success, exit retry loop
                            else:
                                self.logger.warning(f"  No pathways found for: '{term}'")
                                break  # No results, no need to retry
                                
                        except ConnectionError as e:
                            if attempt < max_retries - 1:
                                self.logger.warning(f"    Connection error (attempt {attempt + 1}/{max_retries}): {str(e)}")
                                self.logger.warning(f"    Retrying in {retry_delay} seconds...")
                                time.sleep(retry_delay)
                                retry_delay *= 2  # Exponential backoff
                            else:
                                self.logger.error(f"    Failed after {max_retries} attempts: {str(e)}")
                                
                        except KeyError as e:
                            self.logger.error(f"    API response missing expected field: {str(e)}")
                            self.logger.debug(f"    Response structure: {results[0] if results else 'empty'}")
                            break  # Don't retry on data format issues
                            
                        except Exception as e:
                            self.logger.error(f"    Unexpected error: {str(e)}")
                            break
                            
            except Exception as e:
                self.logger.error(f"  Error processing '{term}': {str(e)}")
                continue
        
        if not all_pathways:
            self.logger.error("\n❌ ERROR: No pathways found")
            self.logger.error("💡 Suggestions:")
            self.logger.error("   - Try simpler keywords (e.g., 'immune system', 'cytokine')")
            self.logger.error("   - Check your internet connection")
            self.logger.error("   - Browse pathways at: https://reactome.org/PathwayBrowser/")
            self.logger.error("   - Use specific pathway IDs (e.g., R-HSA-1280218)")
            self.logger.error("\n   📌 Example pathway IDs for your research:")
            self.logger.error("      R-HSA-1280218 (Adaptive Immune System)")
            self.logger.error("      R-HSA-1280215 (Cytokine Signaling in Immune system)")
            self.logger.error("      R-HSA-449147 (Signaling by Interleukins)")
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
            
            # Retry logic for protein extraction
            max_retries = 3
            retry_delay = 2
            
            for attempt in range(max_retries):
                try:
                    proteins = get_pathway_proteins(pathway_id)
                    
                    if proteins:
                        for protein in proteins:
                            protein['pathway_id'] = pathway_id
                            protein['pathway_name'] = pathway_name
                            all_proteins.append(protein)
                        
                        self.logger.info(f"    Found {len(proteins)} proteins")
                        break  # Success
                    else:
                        self.logger.warning(f"    No proteins found for this pathway")
                        break
                        
                except ConnectionError as e:
                    if attempt < max_retries - 1:
                        self.logger.warning(f"      Connection error (attempt {attempt + 1}/{max_retries})")
                        time.sleep(retry_delay)
                        retry_delay *= 2
                    else:
                        self.logger.error(f"      Failed after {max_retries} attempts: {str(e)}")
                        
                except Exception as e:
                    self.logger.warning(f"      Error: {str(e)}")
                    break
        
        if not all_proteins:
            self.logger.error("\n❌ ERROR: No proteins extracted from pathways")
            self.logger.error("💡 This could mean:")
            self.logger.error("   - Reactome API is temporarily unavailable")
            self.logger.error("   - Network connection issues")
            self.logger.error("   - Pathway IDs are invalid")
            self.logger.error("\n   Try running the pipeline again in a few minutes.")
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
        target_genes_upper = set([g.upper() for g in target_genes if g and isinstance(g, str)])
        pathway_genes_upper = set(self.pathway_proteins['gene_name'].str.upper())
        
        overlapping_genes = target_genes_upper.intersection(pathway_genes_upper)
        
        if not overlapping_genes:
            self.logger.error("\n❌ ERROR: No overlapping targets found")
            self.logger.error("   Predicted targets don't match pathway proteins")
            self.logger.error("💡 This could mean:")
            self.logger.error("   - The compound doesn't target these pathways")
            self.logger.error("   - Try different/broader pathways")
            self.logger.error("   - Check gene name formats in input files")
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
