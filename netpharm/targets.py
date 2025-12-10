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
            self.swiss_targets = self.swiss_targets[
                self.swiss_targets['Probability'] > swiss_threshold
            ]
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
        self.logger.info(f"2. Paste this SMILES: {smiles}")
        self.logger.info("3. Click 'Predict'")
        self.logger.info("4. Wait for results")
        self.logger.info("5. Copy/Download the results table")
        self.logger.info(f"6. Save as: {output_dir}/superpred_results.csv")
        self.logger.info("   (CSV format with columns: Target, Probability, etc.)")
        
        input("\n⏸ Press ENTER after downloading SuperPred results...")
        
        # Load SuperPred results
        superpred_file = os.path.join(output_dir, "superpred_results.csv")
        if not os.path.exists(superpred_file):
            self.logger.error(f"\n❌ ERROR: File not found: {superpred_file}")
            self.logger.error("   Please download the results and try again.")
            raise SystemExit(1)
        
        try:
            self.superpred_targets = pd.read_csv(superpred_file)
            # Assuming SuperPred has 'Probability' column (adjust if different)
            if 'Probability' in self.superpred_targets.columns:
                # Convert percentage to decimal if needed
                if self.superpred_targets['Probability'].max() > 1:
                    self.superpred_targets['Probability'] = self.superpred_targets['Probability'] / 100
                self.superpred_targets = self.superpred_targets[
                    self.superpred_targets['Probability'] > superpred_threshold
                ]
            self.logger.info(f"\n✓ SuperPred: {len(self.superpred_targets)} targets loaded")
            self.logger.info(f"   (filtered by probability > {superpred_threshold})")
        except Exception as e:
            self.logger.error(f"\n❌ ERROR: Failed to load SuperPred results")
            self.logger.error(f"   Reason: {str(e)}")
            raise SystemExit(1)
        
        # Combine targets
        stats = {
            'swiss_count': len(self.swiss_targets),
            'superpred_count': len(self.superpred_targets),
            'total_unique': len(set(list(self.swiss_targets.get('Target', [])) + 
                                   list(self.superpred_targets.get('Target', [])))),
            'overlapping': 0
        }
        
        # Find overlapping targets
        if 'Target' in self.swiss_targets.columns and 'Target' in self.superpred_targets.columns:
            swiss_set = set(self.swiss_targets['Target'])
            superpred_set = set(self.superpred_targets['Target'])
            overlap = swiss_set.intersection(superpred_set)
            stats['overlapping'] = len(overlap)
        
        self.logger.info("\n" + "="*70)
        self.logger.info("TARGET PREDICTION SUMMARY:")
        self.logger.info("="*70)
        self.logger.info(f"SwissTargetPrediction: {stats['swiss_count']} targets")
        self.logger.info(f"SuperPred: {stats['superpred_count']} targets")
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
            self.logger.info("   Columns: Target name, probability, model accuracy")
    
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
        
        return list(set(all_genes))
