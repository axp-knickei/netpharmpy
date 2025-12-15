"""
Target prediction handling (manual workflow).
"""

import pandas as pd
import os


class TargetPredictor:
    """Handle target prediction workflow."""
    
    def __init__(self, logger):
        """
        Initialize TargetPredictor.
        
        Args:
            logger: Logger instance
        """
        self.logger = logger
        self.swiss_targets = None
        self.superpred_targets = None
        self.all_targets = None
    
    def predict_targets_manual(self, smiles, swiss_threshold=0.0, superpred_threshold=0.5, 
                               output_dir="./data"):
        """
        Guide user through manual target prediction.
        
        Args:
            smiles: Canonical SMILES string
            swiss_threshold: SwissTargetPrediction probability threshold
            superpred_threshold: SuperPred probability threshold
            output_dir: Directory for downloaded files
        
        Returns:
            dict: Target prediction statistics
        """
        self.logger.info("\n" + "="*70)
        self.logger.info("[STEP 2] TARGET PREDICTION (MANUAL WORKFLOW)")
        self.logger.info("="*70)
        
        os.makedirs(output_dir, exist_ok=True)
        
        # SwissTargetPrediction instructions
        self.logger.info("\n--- SwissTargetPrediction ---")
        self.logger.info("Please follow these steps:")
        self.logger.info("1. Go to: http://www.swisstargetprediction.ch/")
        self.logger.info("2. Select organism: Homo sapiens (Human)")
        self.logger.info(f"3. Paste this SMILES: {smiles}")
        self.logger.info("4. Click 'Submit'")
        self.logger.info("5. Wait for results (~15-20 seconds)")
        self.logger.info("6. Click 'Download' button")
        self.logger.info(f"7. Save as: {output_dir}/swiss_results.csv")
        
        input("\n⏸ Press ENTER after downloading SwissTargetPrediction results...")
        
        # Load SwissTargetPrediction results
        swiss_file = os.path.join(output_dir, "swiss_results.csv")
        if not os.path.exists(swiss_file):
            self.logger.error(f"\n❌ ERROR: File not found: {swiss_file}")
            self.logger.error("   Please download the results and try again.")
            raise SystemExit(1)
        
        try:
            self.swiss_targets = pd.read_csv(swiss_file)
            
            # Handle column name variations
            prob_col = None
            for col in ['Probability', 'Probability*', 'Probability (%)']:
                if col in self.swiss_targets.columns:
                    prob_col = col
                    break
            
            if prob_col is None:
                self.logger.error(f"\n❌ ERROR: Could not find probability column")
                raise SystemExit(1)
            
            # Filter by threshold
            self.swiss_targets = self.swiss_targets[
                self.swiss_targets[prob_col] > swiss_threshold
            ].copy()
            
            # Standardize column name
            if prob_col != 'Probability':
                self.swiss_targets.rename(columns={prob_col: 'Probability'}, inplace=True)

            # FIX: Prioritize 'Common name' for gene symbols and split multiple genes
            if 'Common name' in self.swiss_targets.columns:
                self.logger.info("   Processing 'Common name' column for gene symbols...")
                # Split entries like "IKBKG IKBKB CHUK" into lists
                self.swiss_targets['gene_list'] = self.swiss_targets['Common name'].astype(str).str.split()
                # Explode the lists so each gene gets its own row
                self.swiss_targets = self.swiss_targets.explode('gene_list')
                # Use this new column as the primary 'Target'
                self.swiss_targets['Target'] = self.swiss_targets['gene_list']
            
            self.logger.info(f"\n✓ SwissTargetPrediction: {len(self.swiss_targets)} targets loaded")
            self.logger.info(f"   (filtered by probability > {swiss_threshold})")
            
        except Exception as e:
            self.logger.error(f"\n❌ ERROR: Failed to load SwissTargetPrediction results")
            self.logger.error(f"   Reason: {str(e)}")
            raise SystemExit(1)
        
        # SuperPred instructions
        self.logger.info("\n--- SuperPred ---")
        self.logger.info("Please follow these steps:")
        self.logger.info("1. Go to: https://prediction.charite.de/")
        self.logger.info("2. Click on 'Target Prediction' tab at the top")
        self.logger.info(f"3. Paste SMILES: {smiles}")
        self.logger.info("4. Click 'Start Calculation'")
        self.logger.info("5. Download 'Targets.csv' (Known) and 'Targets (1).csv' (Predicted)")
        self.logger.info(f"   Save to: {output_dir}/")
        
        input("\n⏸ Press ENTER after downloading BOTH SuperPred CSV files...")
        
        # Load and merge SuperPred results
        superpred_known_file = os.path.join(output_dir, "Targets.csv")
        superpred_pred_file = os.path.join(output_dir, "Targets (1).csv")
        
        all_superpred_targets = []
        
        # Helper to load superpred
        if os.path.exists(superpred_known_file):
            try:
                known_df = pd.read_csv(superpred_known_file)
                known_df['Probability'] = 1.0
                known_df['Source'] = 'Known'
                if 'Target Name' in known_df.columns:
                    known_df.rename(columns={'Target Name': 'Target'}, inplace=True)
                all_superpred_targets.append(known_df)
                self.logger.info(f"\n✓ Loaded {len(known_df)} known strong binders from SuperPred")
            except Exception as e:
                self.logger.warning(f"⚠ Could not load known binders: {e}")

        if os.path.exists(superpred_pred_file):
            try:
                pred_df = pd.read_csv(superpred_pred_file)
                if 'Probability' in pred_df.columns:
                    # Fix percentage strings
                    pred_df['Probability'] = pred_df['Probability'].astype(str).str.rstrip('%').astype('float') / 100
                
                pred_df = pred_df[pred_df['Probability'] > superpred_threshold]
                pred_df['Source'] = 'Predicted'
                if 'Target Name' in pred_df.columns:
                    pred_df.rename(columns={'Target Name': 'Target'}, inplace=True)
                
                all_superpred_targets.append(pred_df)
                self.logger.info(f"✓ Loaded {len(pred_df)} predicted targets from SuperPred")
            except Exception as e:
                self.logger.warning(f"⚠ Could not load predicted targets: {e}")

        if all_superpred_targets:
            self.superpred_targets = pd.concat(all_superpred_targets, ignore_index=True)
            # Remove duplicates
            self.superpred_targets = self.superpred_targets.sort_values('Probability', ascending=False)
            self.superpred_targets = self.superpred_targets.drop_duplicates(subset=['Target'], keep='first')
            self.logger.info(f"\n✓ SuperPred total: {len(self.superpred_targets)} unique targets")
        else:
            self.logger.warning("\n⚠ No SuperPred targets loaded")

        # Stats
        stats = {
            'swiss_count': len(self.swiss_targets) if self.swiss_targets is not None else 0,
            'superpred_count': len(self.superpred_targets) if self.superpred_targets is not None else 0,
            'total_unique': len(self.get_all_target_genes())
        }
        
        self.logger.info("\n" + "="*70)
        self.logger.info("TARGET PREDICTION SUMMARY:")
        self.logger.info("="*70)
        self.logger.info(f"SwissTargetPrediction: {stats['swiss_count']} targets")
        self.logger.info(f"SuperPred: {stats['superpred_count']} targets")
        self.logger.info(f"Total unique gene symbols: {stats['total_unique']}")
        
        return stats
    
    def save_targets(self, output_dir):
        """Save results to CSV."""
        os.makedirs(output_dir, exist_ok=True)
        if self.swiss_targets is not None:
            self.swiss_targets.to_csv(os.path.join(output_dir, "swiss_targets.csv"), index=False)
        if self.superpred_targets is not None:
            self.superpred_targets.to_csv(os.path.join(output_dir, "superpred_targets.csv"), index=False)
    
    def get_all_target_genes(self):
        """Get unique gene symbols from all sources."""
        all_genes = []
        
        if self.swiss_targets is not None and 'Target' in self.swiss_targets.columns:
            # Clean up names: upper case, strip spaces
            genes = self.swiss_targets['Target'].astype(str).str.upper().str.strip().tolist()
            all_genes.extend(genes)
        
        if self.superpred_targets is not None and 'Target' in self.superpred_targets.columns:
            # SuperPred targets are often descriptive names, but sometimes symbols.
            # We add them, but Swiss is the higher quality source for symbols here.
            genes = self.superpred_targets['Target'].astype(str).str.upper().str.strip().tolist()
            all_genes.extend(genes)
        
        # Remove duplicates and empty strings
        unique_genes = list(set([g for g in all_genes if g and g != 'NAN']))
        return unique_genes