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
            
            # Handle column name variations (Probability, Probability*, Probability (%))
            prob_col = None
            for col in ['Probability', 'Probability*', 'Probability (%)']:
                if col in self.swiss_targets.columns:
                    prob_col = col
                    break
            
            if prob_col is None:
                self.logger.error(f"\n❌ ERROR: Could not find probability column")
                self.logger.error(f"   Available columns: {self.swiss_targets.columns.tolist()}")
                raise SystemExit(1)
            
            # Filter by threshold
            self.swiss_targets = self.swiss_targets[
                self.swiss_targets[prob_col] > swiss_threshold
            ]
            
            # Standardize column name
            if prob_col != 'Probability':
                self.swiss_targets.rename(columns={prob_col: 'Probability'}, inplace=True)
            
            self.logger.info(f"\n✓ SwissTargetPrediction: {len(self.swiss_targets)} targets loaded")
            self.logger.info(f"   (filtered by probability > {swiss_threshold})")
        except Exception as e:
            self.logger.error(f"\n❌ ERROR: Failed to load SwissTargetPrediction results")
            self.logger.error(f"   Reason: {str(e)}")
            raise SystemExit(1)
        
        # SuperPred instructions (UPDATED WITH CLEAR NAVIGATION)
        self.logger.info("\n--- SuperPred ---")
        self.logger.info("Please follow these steps:")
        self.logger.info("1. Go to: https://prediction.charite.de/")
        self.logger.info("2. Click on 'Target Prediction' tab at the top")
        self.logger.info("   (or go directly to: https://prediction.charite.de/subpages/target_prediction.php)")
        self.logger.info(f"3. In the 'SMILES' field, paste this SMILES:")
        self.logger.info(f"   {smiles}")
        self.logger.info("4. Click 'Start Calculation' button")
        self.logger.info("5. Wait for results (~30-60 seconds)")
        self.logger.info("\n6. Download TWO CSV files:")
        self.logger.info("   a) Under 'Known strong binders' section:")
        self.logger.info("      → Click 'CSV' button → Save as 'Targets.csv'")
        self.logger.info("   b) Under 'Additionally predicted targets' section:")
        self.logger.info("      → Click 'CSV' button → Save as 'Targets (1).csv'")
        self.logger.info(f"\n   Save BOTH files to: {output_dir}/")
        self.logger.info("\n   📝 Note: You can skip 'Indications' section - not needed for pathway analysis")
        
        input("\n⏸ Press ENTER after downloading BOTH SuperPred CSV files...")
        
        # Load and merge SuperPred results
        superpred_known_file = os.path.join(output_dir, "Targets.csv")
        superpred_pred_file = os.path.join(output_dir, "Targets (1).csv")
        
        # Check if at least one file exists
        has_known = os.path.exists(superpred_known_file)
        has_predicted = os.path.exists(superpred_pred_file)
        
        if not has_known and not has_predicted:
            self.logger.error(f"\n❌ ERROR: No SuperPred files found")
            self.logger.error(f"   Expected files in: {output_dir}/")
            self.logger.error(f"   - Targets.csv (known strong binders)")
            self.logger.error(f"   - Targets (1).csv (predicted targets)")
            self.logger.error("\n   Please download at least one file and try again.")
            raise SystemExit(1)
        
        all_superpred_targets = []
        
        # Load known binders (if available)
        if has_known:
            try:
                known_df = pd.read_csv(superpred_known_file)
                
                # Standardize to probability format (known = 100% confidence)
                known_df['Probability'] = 1.0  # 100% - experimentally validated
                known_df['Source'] = 'Known'
                
                # Rename columns to match standard format
                if 'Target Name' in known_df.columns:
                    known_df.rename(columns={'Target Name': 'Target'}, inplace=True)
                
                # Keep relevant columns
                cols_to_keep = ['Target', 'Probability', 'Source']
                if 'UniProt ID' in known_df.columns:
                    cols_to_keep.append('UniProt ID')
                if 'ChEMBL-ID' in known_df.columns:
                    cols_to_keep.append('ChEMBL-ID')
                
                all_superpred_targets.append(known_df[cols_to_keep])
                self.logger.info(f"\n✓ Loaded {len(known_df)} known strong binders from SuperPred")
                
            except Exception as e:
                self.logger.warning(f"\n⚠ Could not load known binders: {str(e)}")
                self.logger.warning(f"   Continuing with predicted targets only...")
        
        # Load predicted targets (if available)
        if has_predicted:
            try:
                predicted_df = pd.read_csv(superpred_pred_file)
                
                # Handle probability column
                if 'Probability' in predicted_df.columns:
                    # Convert percentage string to float (e.g., "98.55%" -> 0.9855)
                    predicted_df['Probability'] = predicted_df['Probability'].str.rstrip('%').astype('float') / 100
                else:
                    self.logger.error("\n❌ ERROR: 'Probability' column not found in predicted targets")
                    self.logger.error(f"   Available columns: {predicted_df.columns.tolist()}")
                    if not has_known:
                        raise SystemExit(1)
                    else:
                        self.logger.warning("   Continuing with known targets only...")
                        has_predicted = False
                
                if has_predicted:
                    predicted_df['Source'] = 'Predicted'
                    
                    # Rename columns to match standard format
                    if 'Target Name' in predicted_df.columns:
                        predicted_df.rename(columns={'Target Name': 'Target'}, inplace=True)
                    
                    # Filter by threshold
                    predicted_df = predicted_df[predicted_df['Probability'] > superpred_threshold]
                    
                    # Keep relevant columns
                    cols_to_keep = ['Target', 'Probability', 'Source']
                    if 'UniProt ID' in predicted_df.columns:
                        cols_to_keep.append('UniProt ID')
                    if 'ChEMBL-ID' in predicted_df.columns:
                        cols_to_keep.append('ChEMBL-ID')
                    
                    all_superpred_targets.append(predicted_df[cols_to_keep])
                    self.logger.info(f"✓ Loaded {len(predicted_df)} predicted targets from SuperPred")
                    self.logger.info(f"   (filtered by probability > {superpred_threshold})")
                
            except Exception as e:
                self.logger.warning(f"\n⚠ Could not load predicted targets: {str(e)}")
                if not has_known:
                    self.logger.error("   No SuperPred data available")
                    raise SystemExit(1)
                else:
                    self.logger.warning("   Continuing with known targets only...")
        
        if not all_superpred_targets:
            self.logger.error("\n❌ ERROR: Failed to load any SuperPred results")
            raise SystemExit(1)
        
        # Combine known + predicted
        self.superpred_targets = pd.concat(all_superpred_targets, ignore_index=True)
        
        # Remove duplicates (prefer known over predicted based on higher probability)
        self.superpred_targets = self.superpred_targets.sort_values('Probability', ascending=False)
        self.superpred_targets = self.superpred_targets.drop_duplicates(subset=['Target'], keep='first')
        
        known_count = len(all_superpred_targets[0]) if has_known else 0
        predicted_count = len(all_superpred_targets[-1]) if has_predicted and len(all_superpred_targets) > 1 else 0
        
        self.logger.info(f"\n✓ SuperPred total: {len(self.superpred_targets)} unique targets")
        self.logger.info(f"   ({known_count} known + {predicted_count} predicted)")
        
        # Combine targets from both databases
        stats = {
            'swiss_count': len(self.swiss_targets),
            'superpred_count': len(self.superpred_targets),
            'superpred_known': known_count,
            'superpred_predicted': predicted_count,
            'total_unique': 0,
            'overlapping': 0
        }
        
        # Find overlapping targets
        if 'Target' in self.swiss_targets.columns and 'Target' in self.superpred_targets.columns:
            swiss_set = set(self.swiss_targets['Target'].str.upper())
            superpred_set = set(self.superpred_targets['Target'].str.upper())
            
            # Calculate statistics
            all_targets = swiss_set.union(superpred_set)
            overlap = swiss_set.intersection(superpred_set)
            
            stats['total_unique'] = len(all_targets)
            stats['overlapping'] = len(overlap)
        
        self.logger.info("\n" + "="*70)
        self.logger.info("TARGET PREDICTION SUMMARY:")
        self.logger.info("="*70)
        self.logger.info(f"SwissTargetPrediction: {stats['swiss_count']} targets")
        self.logger.info(f"SuperPred: {stats['superpred_count']} targets")
        self.logger.info(f"  └─ Known strong binders: {stats['superpred_known']}")
        self.logger.info(f"  └─ Predicted targets: {stats['superpred_predicted']}")
        self.logger.info(f"Overlapping: {stats['overlapping']} targets")
        self.logger.info(f"Total unique: {stats['total_unique']} targets")
        
        return stats
    
    def save_targets(self, output_dir):
        """
        Save target prediction results.
        
        Args:
            output_dir: Output directory
        """
        os.makedirs(output_dir, exist_ok=True)
        
        if self.swiss_targets is not None:
            swiss_path = os.path.join(output_dir, "swiss_targets.csv")
            self.swiss_targets.to_csv(swiss_path, index=False)
            self.logger.info(f"\n📄 SwissTargetPrediction results saved to: {swiss_path}")
            self.logger.info("   Columns: Target name, UniProt ID, probability, target class")
        
        if self.superpred_targets is not None:
            superpred_path = os.path.join(output_dir, "superpred_targets.csv")
            self.superpred_targets.to_csv(superpred_path, index=False)
            self.logger.info(f"📄 SuperPred results saved to: {superpred_path}")
            self.logger.info("   Columns: Target name, probability, source (Known/Predicted)")
            self.logger.info("   Contains: Known strong binders + predicted targets (merged)")
    
    def get_all_target_genes(self):
        """
        Get list of all unique target gene names.
        
        Returns:
            list: Unique gene names
        """
        all_genes = []
        
        if self.swiss_targets is not None and 'Target' in self.swiss_targets.columns:
            all_genes.extend(self.swiss_targets['Target'].tolist())
        
        if self.superpred_targets is not None and 'Target' in self.superpred_targets.columns:
            all_genes.extend(self.superpred_targets['Target'].tolist())
        
        # Remove duplicates and return
        return list(set(all_genes))
